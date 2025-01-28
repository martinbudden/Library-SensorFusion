#include "SensorFusionFilter.h"
#include <cstdint>


/*!
Reciprocal square root
Implementation of [fast inverse square root](http://en.wikipedia.org/wiki/Fast_inverse_square_root)
using [Pizer’s optimisation](https://pizer.wordpress.com/2008/10/12/fast-inverse-square-root/) and
using `union` rather than `reinterpret_cast` to avoid breaking strict-aliasing rules.

The Xtensa floating point coprocessor (used on the ESP32) has some hardware support for reciprocal square root: it has
an RSQRT0.S (single-precision reciprocal square root initial step) instruction.
However benchmarking shows that FAST_RECIPROCAL_SQUARE_ROOT is approximately 3.5 times faster than `1.0F / sqrtf()`
*/
inline float reciprocalSqrt(float x)
{
#if defined(USE_FAST_RECIPROCAL_SQUARE_ROOT) || defined(USE_FAST_RECIPROCAL_SQUARE_ROOT_TWO_ITERATIONS)
    union {
        float y;
        int32_t i;
    } u = {x};

// NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
    u.i = 0x5f1f1412 - (u.i >> 1); // Initial estimate for Newton–Raphson method
    // single iteration gives accuracy to 4.5 significant figures
    u.y *= 1.69000231F - 0.714158168F * x * u.y * u.y; // First iteration
#if defined(USE_FAST_RECIPROCAL_SQUARE_ROOT_TWO_ITERATIONS)
    // two iterations gives floating point accuracy to within 2 significant bits, and will pass platformio's Unity TEST_ASSERT_EQUAL_FLOAT
    u.y *= 1.5F - (0.5F * x * u.y * u.y); // Second iteration
#endif

    return u.y;
// NOLINTEND(cppcoreguidelines-pro-type-union-access)
#else
    return 1.0F / sqrtf(x);
#endif
}

/*!
Normalize a vector. Return the square of the magnitude.
*/
inline float normalize(xyz_t& v)
{
    const float magnitudeSquared = v.magnitude_squared();
    if (magnitudeSquared != 0.0F) { // [[likely]]
        v *= reciprocalSqrt(magnitudeSquared);
    }
    return magnitudeSquared;
}

/*!
Normalize a quaternion.
*/
inline void normalize(Quaternion& q)
{
    q *= reciprocalSqrt(q.magnitude_squared());
}

void SensorFusionFilterBase::reset()
{
    q0 = 0.0F;
    q1 = 0.0F;
    q2 = 0.0F;
    q3 = 0.0F;
}

/*!
Directly set the filter orientation quaternion. Used by unit test code.
*/
void SensorFusionFilterBase::_setAndNormalizeQ(float q0_, float q1_, float q2_, float q3_)
{
    q0 = q0_;
    q1 = q1_;
    q2 = q2_;
    q3 = q3_;
    const float magnitudeReciprocal = reciprocalSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    q0 *= magnitudeReciprocal;
    q1 *= magnitudeReciprocal;
    q2 *= magnitudeReciprocal;
    q3 *= magnitudeReciprocal;
}

/*!
Calculate quaternion derivative (qDot) from angular rate https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
Twice the actual value is returned to reduce the number of multiplications needed in subsequent code.
*/
Quaternion SensorFusionFilterBase::twoQdot(const xyz_t& gyroRPS) const
{
    return Quaternion(
        -q1*gyroRPS.x - q2*gyroRPS.y - q3*gyroRPS.z,
         q0*gyroRPS.x + q2*gyroRPS.z - q3*gyroRPS.y,
         q0*gyroRPS.y - q1*gyroRPS.z + q3*gyroRPS.x,
         q0*gyroRPS.z + q1*gyroRPS.y - q2*gyroRPS.x
    );
}


void ComplementaryFilter::setFreeParameters(float parameter0, float parameter1)
{
     _alpha = parameter0;
     (void)parameter1;
}

Quaternion ComplementaryFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    // Create q from q0, q1, q2, q3
    Quaternion q(q0, q1, q2, q3);

    // Calculate quaternion derivative (qDot) from angular rate https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
    // Twice the actual value is used to reduce the number of multiplications needed
    const Quaternion _2qDot = twoQdot(gyroRPS);

    // Update the attitude quaternion using simple Euler integration (qNew = qOld + qDot*deltaT).
    // Note: to reduce the number of multiplications, _2qDot and halfDeltaT are used, ie qNew = qOld +_2qDot*deltaT*0.5.
    q += _2qDot * (deltaT * 0.5F); // note brackets to ensure scalar multiplication is performed before vector multiplication

    // use the normalized accelerometer data to calculate an estimate of the attitude
    xyz_t acc = accelerometer;
    normalize(acc);
    const Quaternion a = Quaternion::fromEulerAnglesRadians(rollRadiansFromAccNormalized(acc), pitchRadiansFromAccNormalized(acc));

    // use a complementary filter to combine the gyro attitude estimate(q) with the accelerometer attitude estimate(a)
    q = (q - a)*_alpha + a; // optimized form of `_alpha * q + (1.0F - _alpha) * a` : uses fewer operations and can take advantage of multiply-add instruction

    // normalize the orientation quaternion
    normalize(q);

    // update the values of q0, q1, q2, and q3 for the next iteration
    q.getWXYZ(q0, q1, q2, q3);

    return q;
}

/*!
see https://ahrs.readthedocs.io/en/latest/filters/complementary.html
*/
Quaternion ComplementaryFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT)
{
    // Create q from q0, q1, q2, q3
    Quaternion q(q0, q1, q2, q3);

    // Calculate quaternion derivative (qDot) from angular rate https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
    // Twice the actual value is used to reduce the number of multiplications needed
    const Quaternion _2qDot = twoQdot(gyroRPS);

    // Update the attitude quaternion using simple Euler integration (qNew = qOld + qDot*deltaT).
    // Note: to reduce the number of multiplications, _2qDot and halfDeltaT are used, ie qNew = qOld +_2qDot*deltaT*0.5.
    q += _2qDot * (deltaT * 0.5F); // note brackets to ensure scalar multiplication is performed before vector multiplication

    xyz_t acc = accelerometer;
    (void)normalize(acc);

    // Calculate roll(phi) and pitch(theta) from normalized accelerometer values
    const float roll = rollRadiansFromAccNormalized(acc);
    const float pitch = pitchRadiansFromAccNormalized(acc);
    const float cosPhi = cos(roll);
    const float sinPhi = sin(roll);
    const float sinTheta = sin(pitch);
    const float cosTheta = cos(pitch);

    //  Calculate magnetic field vector, b. See https://ahrs.readthedocs.io/en/latest/filters/tilt.html#module-ahrs.filters.tilt
    xyz_t mag = magnetometer;
    (void)normalize(mag);
    const xyz_t b {
        .x =  mag.x*cosTheta + sinTheta*(mag.y*sinPhi + mag.z*cosPhi),
        .y =  mag.y*cosPhi - mag.z*sinPhi,
        // z not used, so set to zero to save calculation
        .z = 0.0F //-mag.x*sinTheta + cosTheta*(mag.y*sinPhi + mag.z*cosPhi);
    };

    // Calculate yaw from the magnetic field vector
    const float yaw = atan2(-b.y, b.x);

    // use the accelerometer and magnetometer data to calculate an estimate of the attitude, am
    const Quaternion am = Quaternion::fromEulerAnglesRadians(roll, pitch, yaw);

    // use a complementary filter to combine the gyro attitude estimate(q) with the accelerometer/magnetometer attitude estimate(am)
    q = (q - am)*_alpha + am; // optimized form of `_alpha * q + (1.0F - _alpha) * am` : uses fewer operations and can take advantage of multiply-add instruction

    // normalize the orientation quaternion
    normalize(q);

    // update the values of q0, q1, q2, and q3 for the next iteration
    q.getWXYZ(q0, q1, q2, q3);

    return q;
}


void MahonyFilter::setFreeParameters(float parameter0, float parameter1)
{
     _kp = parameter0;
     _ki = parameter1;
}

Quaternion MahonyFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    // Create q from q0, q1, q2, q3
    QuaternionG q(q0, q1, q2, q3);

    // Normalize acceleration
    xyz_t acc = accelerometer;
    (void)normalize(acc);

    // Calculate estimated direction of gravity in the sensor coordinate frame
    const xyz_t gravity = q.gravity();

    // Error is the cross product between direction measured by acceleration and estimated direction of gravity
    const xyz_t error = acc.cross_product(gravity);

    // Apply proportional feedback
    xyz_t gyro = gyroRPS;
    gyro += error * _kp;

    // Apply integral feedback if _ki set
    if (_ki > 0.0F) {
        _errorIntegral += error * (_ki * deltaT); // note brackets to ensure scalar multiplication is performed before vector multiplication
        gyro += _errorIntegral;
    }

    const Quaternion _2qDot = twoQdot(gyro);

    // Update the attitude quaternion using simple Euler integration (qNew = qOld + qDot*deltaT).
    // Note: to reduce the number of multiplications, _2qDot and deltaT*0.5 are used, ie qNew = qOld +_2qDot*deltaT*0.5F
    q += _2qDot * (deltaT * 0.5F); // note brackets to ensure scalar multiplication is performed before quaternion multiplication

    // Normalize the orientation quaternion
    normalize(q);
    // update the values of q0, q1, q2, and q3 for the next iteration
    q.getWXYZ(q0, q1, q2, q3);

    return q;
}

Quaternion MahonyFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT)
{
    (void)magnetometer;
    return update(gyroRPS, accelerometer, deltaT);
}


void MadgwickFilter::setFreeParameters(float parameter0, float parameter1)
{
     _beta = parameter0;
     (void)parameter1;
}

/*!
Madgwick AHRS algorithm, calculates orientation by fusing output from gyroscope and accelerometer.
(No magnetometer is used in this implementation.)

The orientation is calculated as the integration of the gyroscope measurements summed with the measurement from the accelerometer multiplied by a gain.
A low gain gives more weight to the gyroscope more and so is more susceptible to drift.
A high gain gives more weight to the accelerometer and so is more susceptible to accelerometer noise, lag, and other accelerometer errors.
A gain of zero means that orientation is determined by solely by the gyroscope.

See [Sebastian Madgwick's Phd thesis](https://x-io.co.uk/downloads/madgwick-phd-thesis.pdf)
and also x-io Technologies [sensor fusion library](https://github.com/xioTechnologies/Fusion)

For computation efficiency this code refactors the code used in many implementations (Arduino, Adafruit, M5Stack, Reefwing-AHRS),
[see MadgwickRefactoring](../../../documents/MadgwickRefactoring.md)
*/
Quaternion MadgwickFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    // Calculate quaternion derivative (qDot) from the angular rate
    // Twice the actual value is used to reduce the number of multiplications needed
    float _2qDot0 = -q1*gyroRPS.x - q2*gyroRPS.y - q3*gyroRPS.z;
    float _2qDot1 =  q0*gyroRPS.x + q2*gyroRPS.z - q3*gyroRPS.y;
    float _2qDot2 =  q0*gyroRPS.y - q1*gyroRPS.z + q3*gyroRPS.x;
    float _2qDot3 =  q0*gyroRPS.z + q1*gyroRPS.y - q2*gyroRPS.x;

    xyz_t a = accelerometer;
    // Normalize acceleration if it is non-zero
    const float accMagnitudeSquared = a.x*a.x + a.y*a.y + a.z*a.z;
    if (accMagnitudeSquared != 0.0F) { // [[likely]]
        const float accMagnitudeReciprocal = reciprocalSqrt(accMagnitudeSquared);
        a.x *= accMagnitudeReciprocal;
        a.y *= accMagnitudeReciprocal;
        a.z *= accMagnitudeReciprocal;
    }

    // Acceleration is an unreliable indicator of orientation when in high-g maneuvers,
    // so exclude it from the calculation in these cases
    if (accMagnitudeSquared <= _accMagnitudeSquaredMax) {
        // Auxiliary variables to avoid repeated arithmetic
        const float _2q1q1_plus_2q2q2 = 2.0F*(q1*q1 + q2*q2);
        const float common = 2.0F*(q0*q0 + q3*q3 - 1.0F + _2q1q1_plus_2q2q2 + a.z);

        // Gradient decent algorithm corrective step
        const float s0 = q0*(_2q1q1_plus_2q2q2) + q2*a.x - q1*a.y;
        const float s1 = q1*common              - q3*a.x - q0*a.y;
        const float s2 = q2*common              + q0*a.x - q3*a.y;
        const float s3 = q3*(_2q1q1_plus_2q2q2) - q1*a.x - q2*a.y;

        const float _2betaMagnitudeReciprocal = 2.0F * _beta * reciprocalSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3);

        // Add the corrective step to the quaternion derivative
        // Twice the actual value is used to reduce the number of multiplications needed
        _2qDot0 -= s0 * _2betaMagnitudeReciprocal;
        _2qDot1 -= s1 * _2betaMagnitudeReciprocal;
        _2qDot2 -= s2 * _2betaMagnitudeReciprocal;
        _2qDot3 -= s3 * _2betaMagnitudeReciprocal;
    }

    // Update the attitude quaternion using simple Euler integration (qNew = qOld + qDot*deltaT).
    // Note: to reduce the number of multiplications, _2qDot and halfDeltaT are used, ie qNew = qOld +_2qDot*halfDeltaT.
    const float halfDeltaT = deltaT * 0.5F;
    q0 += _2qDot0 * halfDeltaT;
    q1 += _2qDot1 * halfDeltaT;
    q2 += _2qDot2 * halfDeltaT;
    q3 += _2qDot3 * halfDeltaT;

    // Normalize the orientation quaternion
    const float magnitudeReciprocal = reciprocalSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    q0 *= magnitudeReciprocal;
    q1 *= magnitudeReciprocal;
    q2 *= magnitudeReciprocal;
    q3 *= magnitudeReciprocal;

    return Quaternion(q0, q1, q2, q3); // NOLINT(modernize-return-braced-init-list) false positive
}

/*!
For computation efficiency this code refactors the code used in many implementations (Arduino, Adafruit, M5Stack, Reefwing-AHRS),
[see MadgwickRefactoring](../../../documents/MadgwickRefactoring.md)
*/
Quaternion MadgwickFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT)
{
    xyz_t a = accelerometer;
    const float accMagnitudeSquared = normalize(a);
    // Acceleration is an unreliable indicator of orientation when in high-g maneuvers,
    // so exclude it from the calculation in these cases
    if (accMagnitudeSquared > _accMagnitudeSquaredMax) {
        a.x = 0.0F;
        a.y = 0.0F;
        a.z = 0.0F;
    }

    xyz_t m = magnetometer;
    (void)normalize(m);

    // Auxiliary variables to avoid repeated arithmetic
    const float q0q0 = q0*q0;
    const float q0q1 = q0*q1;
    const float q0q2 = q0*q2;
    const float q0q3 = q0*q3;
    const float q1q1 = q1*q1;
    const float q1q2 = q1*q2;
    const float q1q3 = q1*q3;
    const float q2q2 = q2*q2;
    const float q2q3 = q2*q3;
    const float q3q3 = q3*q3;

    const float q1q1_plus_q2q2 = q1q1 + q2q2;
    const float q2q2_plus_q3q3 = q2q2 + q3q3;

    // Reference direction of Earth's magnetic field
    const float hX = m.x*(q0q0 + q1q1 - q2q2_plus_q3q3) + 2.0F*(m.y*(q1q2 - q0q3) + m.z*(q0q2 + q1q3));
    const float hY = 2.0F*(m.x*(q0q3 + q1q2) + m.y*(q0q0 - q1q1 + q2q2 - q3q3) + m.z*(q2q3 - q0q1));

    const float bXbX = hX*hX + hY*hY;
    const float bX =   sqrt(bXbX);
    const float bZ =   2.0F*(m.x*(q1q3 - q0q2) + m.y*(q0q1 + q2q3)) + m.z*(q0q0 - q1q1_plus_q2q2 + q3q3);
    const float bZbZ = bZ*bZ;
    const float _4bXbZ = 4.0F * bX * bZ;

    const float mXbX = m.x*bX;
    const float mYbX = m.y*bX;
    const float mZbX = m.z*bX;
    const float mZbZ = m.z*bZ;

    const float aX_plus_mXbZ = a.x + m.x*bZ;
    const float aY_plus_mYbZ = a.y + m.y*bZ;

    const float sumSquaresMinusOne  = q0q0 + q1q1_plus_q2q2 + q3q3 - 1.0F;
    const float common              = sumSquaresMinusOne + q1q1_plus_q2q2 + a.z;

    // Gradient decent algorithm corrective step
    const float s0 =
        + q0 * 2.0F * (q1q1_plus_q2q2*(1.0F + bZbZ) + bXbX*q2q2_plus_q3q3)
        - q1 * aY_plus_mYbZ
        + q2 * (aX_plus_mXbZ - mZbX)
        + q3 * (mYbX  - _4bXbZ*q0q1);

    const float s1 =
        - q0 * aY_plus_mYbZ
        + q1 * 2.0F * (common + mZbZ + bXbX*q2q2_plus_q3q3 + bZbZ*(sumSquaresMinusOne + q1q1_plus_q2q2))
        - q2 * mYbX
        - q3 * (aX_plus_mXbZ + mZbX + _4bXbZ*(0.5F*sumSquaresMinusOne + q1q1));

    const float s2 =
        + q0 * (aX_plus_mXbZ - mZbX)
        - q1 * mYbX
        + q2 * 2.0F * (common + mZbZ + mXbX + bXbX*(sumSquaresMinusOne + q2q2_plus_q3q3) + bZbZ*(sumSquaresMinusOne + q1q1_plus_q2q2))
        - q3 * (aY_plus_mYbZ + _4bXbZ*q1q2);

    const float s3 =
        + q0 *  mYbX
        - q1 * (aX_plus_mXbZ + mZbX + _4bXbZ*(0.5F*sumSquaresMinusOne + q3q3))
        - q2 * aY_plus_mYbZ
        + q3 * 2.0F * (q1q1_plus_q2q2*(1.0F + bZbZ) + mXbX + bXbX*(sumSquaresMinusOne + q2q2_plus_q3q3));

    const float _2betaNormReciprocal = 2.0F * _beta * reciprocalSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3);

    // Calculate quaternion derivative (qDot) from the angular rate and the corrective step
    // Twice the actual value is used to reduce the number of multiplications needed
    const float _2qDot0 = -q1*gyroRPS.x - q2*gyroRPS.y - q3*gyroRPS.z - s0 * _2betaNormReciprocal;
    const float _2qDot1 =  q0*gyroRPS.x + q2*gyroRPS.z - q3*gyroRPS.y - s1 * _2betaNormReciprocal;
    const float _2qDot2 =  q0*gyroRPS.y - q1*gyroRPS.z + q3*gyroRPS.x - s2 * _2betaNormReciprocal;
    const float _2qDot3 =  q0*gyroRPS.z + q1*gyroRPS.y - q2*gyroRPS.x - s3 * _2betaNormReciprocal;

    // Update the attitude quaternion using simple Euler integration (qNew = qOld + qDot*deltaT).
    // Note: to reduce the number of multiplications, _2qDot and halfDeltaT are used, ie qNew = qOld +_2qDot*halfDeltaT.
    const float halfDeltaT = deltaT * 0.5F;
    q0 += _2qDot0 * halfDeltaT;
    q1 += _2qDot1 * halfDeltaT;
    q2 += _2qDot2 * halfDeltaT;
    q3 += _2qDot3 * halfDeltaT;

    // Normalize the orientation quaternion
    const float magnitudeReciprocal = reciprocalSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    q0 *= magnitudeReciprocal;
    q1 *= magnitudeReciprocal;
    q2 *= magnitudeReciprocal;
    q3 *= magnitudeReciprocal;

    return Quaternion(q0, q1, q2, q3); // NOLINT(modernize-return-braced-init-list) false positive
}