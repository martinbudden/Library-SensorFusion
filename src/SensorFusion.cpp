#include "SensorFusion.h"
#include <cstdint>

float normalize(xyz_t& v);
void normalize(Quaternion& q);

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
#if defined(LIBRARY_SENSOR_FUSION_USE_FAST_RECIPROCAL_SQUARE_ROOT) || defined(LIBRARY_SENSOR_FUSION_USE_FAST_RECIPROCAL_SQUARE_ROOT_TWO_ITERATIONS)
    union {
        float f;
        int32_t i;
    } u { .f = x };

// NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
    u.i = 0x5f1f1412 - (u.i >> 1); // Initial estimate for Newton–Raphson method
    // single iteration gives accuracy to 4.5 significant figures
    u.f *= 1.69000231F - 0.714158168F * x * u.f * u.f; // First iteration
#if defined(LIBRARY_SENSOR_FUSION_USE_FAST_RECIPROCAL_SQUARE_ROOT_TWO_ITERATIONS)
    // two iterations gives floating point accuracy to within 2 significant bits, and will pass platformio's Unity TEST_ASSERT_EQUAL_FLOAT
    u.f *= 1.5F - (0.5F * x * u.f * u.f); // Second iteration
#endif

    return u.f;
// NOLINTEND(cppcoreguidelines-pro-type-union-access)
#else
    return 1.0F / sqrtf(x);
#endif
}

/*!
Normalize a vector. Return the square of the magnitude.
*/
float normalize(xyz_t& v)
{
    const float magnitudeSquared = v.magnitudeSquared();
    if (magnitudeSquared != 0.0F) { // [[likely]]
        v *= reciprocalSqrt(magnitudeSquared);
    }
    return magnitudeSquared;
}

/*!
Normalize a quaternion.
*/
void normalize(Quaternion& q)
{
    q *= reciprocalSqrt(q.magnitudeSquared());
}

void SensorFusionFilterBase::reset()
{
    q0 = 1.0F;
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
    return Quaternion( // NOLINT(modernize-return-braced-init-list)
        -q1*gyroRPS.x - q2*gyroRPS.y - q3*gyroRPS.z,
         q0*gyroRPS.x + q2*gyroRPS.z - q3*gyroRPS.y,
         q0*gyroRPS.y - q1*gyroRPS.z + q3*gyroRPS.x,
         q0*gyroRPS.z + q1*gyroRPS.y - q2*gyroRPS.x
    );
}

bool SensorFusionFilterBase::requiresInitialization() const
{
    return false;
}

void SensorFusionFilterBase::setFreeParameters(float parameter0, float parameter1)
{
    (void)parameter0;
    (void)parameter1;
}

void ComplementaryFilter::setFreeParameters(float parameter0, float parameter1)
{
    _alpha = parameter0;
    (void)parameter1;
}

Quaternion ComplementaryFilter::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
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
Quaternion ComplementaryFilter::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT)
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
    const float cosPhi = cosf(roll);
    const float sinPhi = sinf(roll);
    const float sinTheta = sinf(pitch);
    const float cosTheta = cosf(pitch);

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
    const float yaw = atan2f(-b.y, b.x);

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

Quaternion MahonyFilter::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    // Create q from q0, q1, q2, q3
    QuaternionG q(q0, q1, q2, q3);

    // Normalize acceleration
    xyz_t acc = accelerometer;
    (void)normalize(acc);

    // Calculate estimated direction of gravity in the sensor coordinate frame
    const xyz_t gravity = q.gravity();

    // Error is the cross product between direction measured by acceleration and estimated direction of gravity
    const xyz_t error = acc.cross(gravity);

    // Quadratic Interpolation (From Attitude Representation and Kinematic Propagation for Low-Cost UAVs by Robert T. Casey, Equation 14)
    // See https://docs.rosflight.org/v1.3/algorithms/estimator/#modifications-to-original-passive-filter for a publicly available explanation
    xyz_t gyro; // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
    if (_useQuadraticInterpolation) {
        gyro =  (5.0F/12.0F)*gyroRPS + (8.0F/12.0F)*_gyroRPS_1 - (1.0F/12.0F)*_gyroRPS_2; 
        _gyroRPS_2 = _gyroRPS_1;
        _gyroRPS_1 = gyroRPS;
    } else {
        gyro = gyroRPS;
    }

    // Apply proportional feedback
    gyro += error * _kp;

    // Apply integral feedback if _ki set
    if (_ki > 0.0F) {
        _errorIntegral += error * (_ki * deltaT); // note brackets to ensure scalar multiplication is performed before vector multiplication
        gyro += _errorIntegral;
    }

    if (_useMatrixExponentialApproximation) {
        // Matrix Exponential Approximation (From Attitude Representation and Kinematic Propagation for Low-Cost UAVs by Robert T. Casey, Equation 12)
        const float gyroMagnitude = gyro.magnitude();
        const float theta = gyroMagnitude * 0.5F * deltaT;
#if defined(LIBRARY_VECTOR_SENSOR_FUSION_USE_FAST_TRIGONOMETRY) && defined(LIBRARY_VECTOR_QUATERNION_MATRIX_USE_FAST_TRIGONOMETRY)
        float sinTheta {};
        float cosCosTheta {};
        FastTrigonometry::sincos(theta, sinTheta, cosTheta);
        const float t1 = cosTheta;
        const float t2 = (1.0F / gyroMagnitude) * sinTheta;
#else
        const float t1 = cosf(theta);
        const float t2 = (1.0F / gyroMagnitude) * sinf(theta);
#endif

        const float qW = t1*q0 + t2*(-gyro.x * q1 - gyro.y * q2 - gyro.z * q3);
        const float qX = t1*q1 + t2*( gyro.x * q0 + gyro.z * q2 - gyro.y * q3);
        const float qY = t1*q2 + t2*( gyro.y * q0 - gyro.z * q1 + gyro.x * q3);
        const float qZ = t1*q3 + t2*( gyro.z * q0 + gyro.y * q1 - gyro.x * q2);
        q.set(qW, qX, qY, qZ);
    } else {
        const Quaternion _2qDot = twoQdot(gyro);
        // Update the attitude quaternion using simple Euler integration (qNew = qOld + qDot*deltaT).
        // Note: to reduce the number of multiplications, _2qDot and deltaT*0.5 are used, ie qNew = qOld +_2qDot*deltaT*0.5F
        q += _2qDot * (deltaT * 0.5F); // note brackets to ensure scalar multiplication is performed before quaternion multiplication
    }

    // Normalize the orientation quaternion
    normalize(q);
    // update the values of q0, q1, q2, and q3 for the next iteration
    q.getWXYZ(q0, q1, q2, q3);

    return q;
}

Quaternion MahonyFilter::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) // NOLINT(readability-convert-member-functions-to-static) false positive
{
    (void)magnetometer;
    return updateOrientation(gyroRPS, accelerometer, deltaT);
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
Quaternion MadgwickFilter::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    // Calculate quaternion derivative (qDot) from the angular rate
    // Twice the actual value is used to reduce the number of multiplications needed
    float _2qDot0 = -q1*gyroRPS.x - q2*gyroRPS.y - q3*gyroRPS.z;
    float _2qDot1 =  q0*gyroRPS.x + q2*gyroRPS.z - q3*gyroRPS.y;
    float _2qDot2 =  q0*gyroRPS.y - q1*gyroRPS.z + q3*gyroRPS.x;
    float _2qDot3 =  q0*gyroRPS.z + q1*gyroRPS.y - q2*gyroRPS.x;

    xyz_t a = accelerometer;
    // Normalize acceleration if it is non-zero
    const float accMagnitudeSquared = a.magnitudeSquared();
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

    return {q0, q1, q2, q3};
}

/*!
For computation efficiency this code refactors the code used in many implementations (Arduino, Adafruit, M5Stack, Reefwing-AHRS),
[see MadgwickRefactoring](../../../documents/MadgwickRefactoring.md)
*/
Quaternion MadgwickFilter::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT)
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
    const float bX =   sqrtf(bXbX);
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

    return {q0, q1, q2, q3};
}

bool MadgwickFilter::requiresInitialization() const
{
    return true;
}


FilterButterworthCompound::FilterButterworthCompound() :
    _initialized(false),
    _tau(0.0F),
    _deltaT(0.0F)
{}

FilterButterworthCompound::FilterButterworthCompound(float tau, float deltaT) :
    _initialized(false),
    _tau(tau),
    _deltaT(deltaT)
{
    setCoefficients(tau, deltaT);
}

void FilterButterworthCompound::setCoefficients(float tau, float deltaT)
{
    static constexpr float M_PI_F = 3.14159265358979323846F;
    static constexpr float M_SQRT2_F = 1.41421356237309504880F;

    _initialized = false;
    _tau = tau;
    _deltaT = deltaT;
    // second order Butterworth filter based on https://stackoverflow.com/a/52764064
    const float fc = (M_SQRT2_F / (2.0F*M_PI_F))/tau; // time constant of dampened, non-oscillating part of step response
    const float C = tanf(M_PI_F*fc*deltaT);
    const float C2 = C*C;
    const float D = C2 + M_SQRT2_F*C + 1.0F;
    const float b0 = C2/D;
    _coefficients.b0 = b0;
    _coefficients.b1 = 2*b0;
    _coefficients.b2 = b0;
    // a0 = 1.0F, implicitly
    _coefficients.a1 = 2.0F*(C2 - 1.0F)/D;
    _coefficients.a2 = (1.0F - M_SQRT2_F*C + C2)/D;
}

void FilterButterworthCompound::setState(state_t& state, float x0) const
{
    // initial state for steady state (equivalent to scipy.signal.lfilter_zi, obtained by setting y=x=x0 in the filter update equation)
    state.s0 = x0*(1.0F - _coefficients.b0); // a0 = 1.0F, implicitly
    state.s1 = x0*(_coefficients.b2 - _coefficients.a2);
}

float FilterButterworthCompound::filterStep(state_t& state, float x) const
{
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
    // difference equations based on scipy.signal.lfilter documentation
    // assumes that a0 == 1.0F
    const float y = _coefficients.b0*x + state.s0;
    state.s0 = _coefficients.b1*x - _coefficients.a1*y + state.s1;
    state.s1 = _coefficients.b2*x - _coefficients.a2*y;
    return y;
}

FilterButterworthXYZ::FilterButterworthXYZ(float tau, float deltaT)
    : FilterButterworthCompound(tau, deltaT)
{
    _state.fill({0.0F, 0.0F});
    setCoefficients(tau, deltaT);
}

xyz_t FilterButterworthXYZ::filter(const xyz_t& v)
{
    // to avoid depending on a single sample, average the first samples (for duration tau)
    // and then use this average to calculate the filter initial state
    if (!_initialized) {
        ++_state[0].s0; // ._state[0].s0 used as count
        _state[X].s1 += v.x; // ._state[index].s1 used as sum
        _state[Y].s1 += v.y;
        _state[Z].s1 += v.z;
        const float count = _state[0].s0;
        const xyz_t ret = {
            _state[X].s1 / count,
            _state[Y].s1 / count,
            _state[Z].s1 / count
        };
        if (count*_deltaT >= _tau) {
            _initialized = true;
            setState(_state[X], ret.x);
            setState(_state[Y], ret.y);
            setState(_state[Z], ret.z);
        }
        return ret;
    }

    // not in initialization phase, so just filter each component of v
    return xyz_t {
        filterStep(_state[X], v.x),
        filterStep(_state[Y], v.y),
        filterStep(_state[Z], v.z)
    };
}

/*!
Filter just X and Y values for computational efficiency.
*/
xyz_t FilterButterworthXYZ::filterXY(const xyz_t& v)
{
    // to avoid depending on a single sample, average the first samples (for duration tau)
    // and then use this average to calculate the filter initial state
    if (!_initialized) {
        ++_state[0].s0; // ._state[0].s0 used as count
        _state[X].s1 += v.x; // ._state[index].s1 used as sum
        _state[Y].s1 += v.y;
        const float count = _state[0].s0;
        const xyz_t ret = {
            _state[X].s1 / count,
            _state[Y].s1 / count,
            0.0F
        };
        if (count*_deltaT >= _tau) {
            _initialized = true;
            setState(_state[X], ret.x);
            setState(_state[Y], ret.y);
        }
        return ret;
    }

    // not in initialization phase, so just filter each component of v
    return xyz_t {
        filterStep(_state[X], v.x),
        filterStep(_state[Y], v.y),
        0.0F
    };
}

FilterButterworthMatrix3x3::FilterButterworthMatrix3x3(float tau, float deltaT)
    : FilterButterworthCompound(tau, deltaT)
{
    _state.fill({0.0F, 0.0F});
    setCoefficients(tau, deltaT);
}

Matrix3x3 FilterButterworthMatrix3x3::filter(const Matrix3x3& m)
{
    Matrix3x3 ret {};

    // to avoid depending on a single sample, average the first samples (for duration tau)
    // and then use this average to calculate the filter initial _state
    if (!_initialized) {
        ++_state[0].s0; // ._state[0].s0 used as count
        for (size_t ii = 0; ii < _state.size(); ++ii) {
            _state[ii].s1 += m[ii]; // ._state[ii].s1 used as sum
        }
        const float count = _state[0].s0;
        for (size_t ii = 0; ii < _state.size(); ++ii) {
            ret[ii] = _state[ii].s1 / count;
        }
        if (count*_deltaT >= _tau) {
            _initialized = true;
            for (size_t ii = 0; ii < _state.size(); ++ii) {
                setState(_state[ii], ret[ii]);
            }
        }
        return ret;
    }

    // not in initialization phase, so just filter each component of m
    // scope return value to allow eliding copy on return (return value optimization RVO)
    for (size_t ii = 0; ii < _state.size(); ++ii) {
        ret[ii] = filterStep(_state[ii], m[ii]);
    }
    return ret;
}

BasicVQF::BasicVQF(float gyroDeltaT, float accDeltaT, float magDeltaT) :
    _params {
        .tauAcc = 3.0F,
        .tauMag = 9.0F
    },
    _coeffs {
        .gyroDeltaT = gyroDeltaT,
        .accDeltaT = accDeltaT,
        .magDeltaT = magDeltaT,
        .kMag = gainFromTau(_params.tauMag, _coeffs.magDeltaT)
    }
{
    _state.accLPF.setCoefficients(_params.tauAcc, _coeffs.accDeltaT);
}

BasicVQF::BasicVQF(float deltaT) :
    BasicVQF(deltaT, deltaT, deltaT)
{
}

void BasicVQF::resetState()
{
    _state.gyroQuaternion.setToIdentity();
    _state.accQuaternion.setToIdentity();

    _state.accLPF.reset();

    _state.delta = 0.0F;
    _state.kMagInit = 1.0F;
}


float BasicVQF::gainFromTau(float tau, float deltaT)
{
    return (tau < 0.0F) ? 0.0F : (tau == 0.0F) ? 1.0F : 1.0F - expf(-deltaT/tau);  // fc = 1/(2*pi*tau)
}

void BasicVQF::updateGyro(const xyz_t& gyroRPS, float deltaT) // NOLINT(readability-convert-member-functions-to-static) false positive
{
    const float gyroMagnitude = gyroRPS.magnitude();
    if (gyroMagnitude > std::numeric_limits<float>::epsilon()) {
        // calculate integration step
        const float angle2 = gyroMagnitude * deltaT * 0.5F;
        const float c = cosf(angle2);
        const float s = sinf(angle2)/gyroMagnitude;
        const Quaternion integrationStep(c, s*gyroRPS.x, s*gyroRPS.y, s*gyroRPS.z);
        _state.gyroQuaternion *= integrationStep;
        normalize(_state.gyroQuaternion);
    }
}

Quaternion BasicVQF::updateAccelerometer(const xyz_t& accelerometer, [[maybe_unused]] float deltaT) // NOLINT(readability-convert-member-functions-to-static) false positive
{
    xyz_t acc = accelerometer * G_TO_METERS_PER_SECOND_PER_SECOND;
    normalize(acc);

    // rotate acc ito the auxiliary frame Ιi, and filter it
    xyz_t accAuxiliary = _state.gyroQuaternion.rotate(acc);
    accAuxiliary = _state.accLPF.filter(accAuxiliary);

    // rotate accAuxiliary to the almost-inertial frame ξi (called the Earth frame in Laidig and Seel's implementation) and normalize
    xyz_t accEarth = _state.accQuaternion.rotate(accAuxiliary);
    normalize(accEarth);

    // inclination correction
    const float qW = sqrtf((accEarth.z + 1.0F) * 0.5F);
    // calculate inclination correction, avoiding numeric issues when acc is close to [0 0 -1], ie the correction step is close (<= 0.00011°) to 180°
    const Quaternion correctionStep = (qW > 1e-6F) ? Quaternion(qW, accEarth.y/(2.0F*qW), -accEarth.x/(2.0F*qW), 0) : Quaternion(0.0F, 1.0F, 0.0F, 0.0F);

    _state.accQuaternion = correctionStep * _state.accQuaternion;
    normalize(_state.accQuaternion);
    _state.orientation6D = _state.accQuaternion * _state.gyroQuaternion;

    return _state.orientation6D;
}

Quaternion BasicVQF::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    updateGyro(gyroRPS, deltaT);
    return updateAccelerometer(accelerometer, deltaT);
}

Quaternion BasicVQF::updateMagnetometer(const xyz_t& magnetometer, float deltaT) // NOLINT(readability-convert-member-functions-to-static,readability-make-member-function-const) false positive
{
    // ignore [0 0 0] samples
    if (magnetometer.x == 0.0F && magnetometer.y == 0.0F && magnetometer.z == 0.0F) {
        return _state.orientation6D;
    }

    // bring magnetometer measurement into 6D earth frame
    const xyz_t magEarth = _state.orientation6D.rotate(magnetometer);

    // calculate disagreement angle based on current magnetometer measurement
    _state.magDisagreementAngle = atan2f(magEarth.x, magEarth.y);

    // update the disagreement angle and make sure it is in the range [-pi, pi]
    _state.magDisagreementAngle = wrapToPi(_state.magDisagreementAngle - _state.delta);

    float k = _coeffs.kMag;

    // ensure fast initial convergence
    if (_state.kMagInit != 0.0F) {
        // make sure that the gain k is at least 1/N, N=1,2,3,... in the first few samples
        if (k < _state.kMagInit) {
            k = _state.kMagInit;
        }

        // iterative expression to calculate 1/N
        _state.kMagInit /= _state.kMagInit + 1.0F;

        // disable if t > tauMag
        if (_state.kMagInit*_params.tauMag < deltaT) {
            _state.kMagInit = 0.0F;
        }
    }

    // first-order filter step
    _state.delta += k*_state.magDisagreementAngle;
    // make sure delta is in the range [-pi, pi]
    _state.delta = wrapToPi(_state.delta);

    _state.orientation9D = _state.orientation6D.rotateZ(_state.delta);
    return _state.orientation9D;
}

Quaternion BasicVQF::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT)
{
    updateGyro(gyroRPS, deltaT);
    updateAccelerometer(accelerometer, deltaT);

    return updateMagnetometer(magnetometer, deltaT);
}

VQF::VQF(float gyroDeltaT, float accDeltaT, float magDeltaT, bool restBiasEstimationEnabled, bool motionBiasEstimationEnabled, [[maybe_unused]] bool magDisturbanceRejectionEnabled) :
    _params {
        .tauAcc = 3.0F,
        .restBiasEstimationEnabled = restBiasEstimationEnabled,
        .motionBiasEstimationEnabled = motionBiasEstimationEnabled,
        .biasSigmaInit = 0.5F,
        .biasForgettingTime = 100.0F,
        .biasClipRPS = 2.0F * DEGREES_TO_RADIANS,
#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MOTION_BIAS_ESTIMATION)
        .biasSigmaMotion = 0.1F,
        .biasVerticalForgettingFactor = 0.0001F,
#endif
        .biasSigmaRest = 0.03F,
        .restMinT = 1.5F,
        .restFilterTau = 0.5F,
        .restThresholdGyroSquared = square(2.0F * DEGREES_TO_RADIANS),
        .restThresholdAccSquared = square(0.5F)
#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MAGNETOMETER)
        , .tauMag = 9.0F,
        .magDisturbanceRejectionEnabled = magDisturbanceRejectionEnabled,
        .magCurrentTau = 0.05F,
        .magRefTau = 20.0F,
        .magNormThreshold = 0.1F,
        .magDipThresholdRadians = 10.0F * DEGREES_TO_RADIANS,
        .magNewTime = 20.0F,
        .magNewFirstTime = 5.0F,
        .magNewMinGyroRPS = 20.0F * DEGREES_TO_RADIANS,
        .magMinUndisturbedTime = 0.5F,
        .magMaxRejectionTime = 60.0F,
        .magRejectionFactor = 2.0F
#endif
    },
    _coeffs {
        .gyroDeltaT = gyroDeltaT,
        .accDeltaT = accDeltaT,
#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MAGNETOMETER)
        .magDeltaT = magDeltaT,
        .kMag = gainFromTau(_params.tauMag, magDeltaT),
        .kMagRef = gainFromTau(_params.magRefTau, magDeltaT),
#endif
        .biasP0 = square(_params.biasSigmaInit*100.0F),
        // the system noise increases the variance from 0 to (0.1 °/s)^2 in biasForgettingTime seconds
        .biasV = square(0.1F*100.0F)*accDeltaT/_params.biasForgettingTime,
#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MOTION_BIAS_ESTIMATION)
        .biasMotionW = calculateBias(_params.biasSigmaMotion),
        .biasVerticalW = _coeffs.biasMotionW / std::max(_params.biasVerticalForgettingFactor, 1e-10F),
#endif
        .biasRestW = calculateBias(_params.biasSigmaRest)
    }
{
    _state.accLPF.setCoefficients(_params.tauAcc, _coeffs.accDeltaT);
    _state.restGyroLPF.setCoefficients(_params.restFilterTau, _coeffs.gyroDeltaT);
    _state.restAccLPF.setCoefficients(_params.restFilterTau, _coeffs.accDeltaT);

#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MOTION_BIAS_ESTIMATION)
    _state.motionBiasEstimateR_LPF.setCoefficients(_params.restFilterTau, _coeffs.accDeltaT);
    _state.motionBiasEstimateBiasLPF.setCoefficients(_params.restFilterTau, _coeffs.accDeltaT);
#endif

#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MAGNETOMETER)
    if (_params.magCurrentTau > 0) {
        _state.magNormDipLPF.setCoefficients(_params.magCurrentTau, _coeffs.magDeltaT);
    }
#else
    (void)magDeltaT;
#endif

    resetState();
}

VQF::VQF(float gyroDeltaT, float accDeltaT, float magDeltaT) :
    VQF(gyroDeltaT, accDeltaT, magDeltaT, true, true, true)
    {}

float VQF::calculateBias(float v) const
{
    const float p = square(v*100.0F);
    return square(p) / _coeffs.biasV + p;
}

float VQF::gainFromTau(float tau, float deltaT)
{
    return (tau < 0.0F) ? 0.0F : (tau == 0.0F) ? 1.0F : 1.0F - expf(-deltaT/tau);  // fc = 1/(2*pi*tau)
}

void VQF::resetState()
{
    _state.gyroQuaternion.setToIdentity();
    _state.accQuaternion.setToIdentity();

    _state.restDetected = false;

    _state.lastAccCorrectionAngularRate = 0.0F;

    _state.gyroBiasRPS = { 0.0F, 0.0F, 0.0F };
    _state.biasP.setToScaledIdentity(_coeffs.biasP0);

    _state.restT = 0.0F;
    _state.restLastGyro = { 0.0F, 0.0F, 0.0F };
    _state.restLastAcc = { 0.0F, 0.0F, 0.0F };

    // then initialize low pass filters
    _state.accLPF.reset();
    _state.restGyroLPF.reset();
    _state.restAccLPF.reset();

#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MOTION_BIAS_ESTIMATION)
    _state.motionBiasEstimateR_LPF.reset();
    _state.motionBiasEstimateBiasLPF.reset();
#endif

#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MAGNETOMETER)
    _state.delta = 0.0F;
    _state.kMagInit = 1.0F;
    _state.lastMagCorrectionAngularRate = 0.0F;
    _state.magDisturbanceDetected = false;
    _state.magDisagreementAngle = 0.0F;
    _state.magRefNorm = 0.0F;
    _state.magRefDip = 0.0F;
    _state.magUndisturbedT = 0.0F;
    _state.magRejectT = _params.magMaxRejectionTime;
    _state.magCandidateNorm = -1.0F;
    _state.magCandidateDip = 0.0F;
    _state.magCandidateT = 0.0F;
    _state.magNormDipLPF.reset();
#endif
}

void VQF::updateGyro(const xyz_t& gyroRPS, float deltaT) // NOLINT(readability-make-member-function-const) false positive
{
    // rest detection
    if (_params.restBiasEstimationEnabled){// || _params.magDisturbanceRejectionEnabled) {
        _state.restLastGyro = _state.restGyroLPF.filter(gyroRPS);

        _state.restLastGyroSquaredDeviation = (gyroRPS - _state.restLastGyro).magnitudeSquared();

        if (_state.restLastGyroSquaredDeviation >= _params.restThresholdGyroSquared
                || fabsf(_state.restLastGyro.x) > _params.biasClipRPS
                || fabsf(_state.restLastGyro.y) > _params.biasClipRPS
                || fabsf(_state.restLastGyro.z) > _params.biasClipRPS) {
            _state.restT = 0.0F;
            _state.restDetected = false;
        }
    }

    // procedure FilterUpdate, line 8, perform gyroscope strapdown integration
    // remove estimated gyro bias
    const xyz_t gyroNoBias = gyroRPS - _state.gyroBiasRPS;
    const float gyroNoBiasMagnitude = gyroNoBias.magnitude();
    if (gyroNoBiasMagnitude > std::numeric_limits<float>::epsilon()) {
        // calculate integration step
        const float angle2 = gyroNoBiasMagnitude * deltaT * 0.5F;
        const float c = cosf(angle2);
        const float s = sinf(angle2)/gyroNoBiasMagnitude;
        const Quaternion integrationStep(c, s*gyroNoBias.x, s*gyroNoBias.y, s*gyroNoBias.z);
        _state.gyroQuaternion *= integrationStep; // integrate
        normalize(_state.gyroQuaternion);
    }
}

Quaternion VQF::updateAccelerometer(const xyz_t& accelerometer, float deltaT) // NOLINT(readability-make-member-function-const)
{
    xyz_t acc = accelerometer * G_TO_METERS_PER_SECOND_PER_SECOND;
    const float accMagnitude = normalize(acc);
    if (accMagnitude == 0.0F) {
        // ignore [0 0 0] samples
        return _state.orientation6D;
    }

    if (_params.restBiasEstimationEnabled) {
        // rest detection
        _state.restLastAcc = _state.restAccLPF.filter(acc);
        _state.restLastAccSquaredDeviation = (acc - _state.restLastAcc).magnitudeSquared();
        if (_state.restLastAccSquaredDeviation >= _params.restThresholdAccSquared) {
            _state.restT = 0.0F;
            _state.restDetected = false;
        } else {
            _state.restT += deltaT;
            if (_state.restT >= _params.restMinT) {
                _state.restDetected = true;
            }
        }
    }

    // rotate acc into the auxiliary frame Ιi, and filter it
    xyz_t accAuxiliary = _state.gyroQuaternion.rotate(acc);
    accAuxiliary = _state.accLPF.filter(accAuxiliary);

    // rotate accAuxiliary to the almost-inertial frame ξi (called the Earth frame in Laidig and Seel's implementation) and normalize
    xyz_t accEarth = _state.accQuaternion.rotate(accAuxiliary);
    normalize(accEarth);

    // calculate correction angular rate to facilitate debugging
    _state.lastAccCorrectionAngularRate = acosf(accEarth.z)/deltaT;

    // calculate inclination correction, avoiding numeric issues when acc is close to [0 0 -1], ie the correction step is close (<= 0.00011°) to 180°
    const float qW = sqrtf((accEarth.z + 1.0F) * 0.5F);
    const Quaternion correctionStep = (qW > 1e-6F) ? Quaternion(qW, accEarth.y/(2.0F*qW), -accEarth.x/(2.0F*qW), 0) : Quaternion(0.0F, 1.0F, 0.0F, 0.0F);

    _state.accQuaternion = correctionStep * _state.accQuaternion;
    normalize(_state.accQuaternion);

    _state.orientation6D = _state.accQuaternion * _state.gyroQuaternion;

    // bias estimation
#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MOTION_BIAS_ESTIMATION)
    if (_params.restBiasEstimationEnabled || _params.motionBiasEstimationEnabled) {
        Matrix3x3 R(_state.orientation6D); // Create rotation matrix from orientation quaternion

        // set measurement error and covariance for the respective Kalman filter update
        xyz_t e {};
        xyz_t w {};
        if (_state.restDetected) {
            e = _state.restLastGyro - _state.gyroBiasRPS;
            w = { _coeffs.biasRestW, _coeffs.biasRestW, _coeffs.biasRestW };
        } else if (_params.motionBiasEstimationEnabled) {
            // calculate R*b_hat (only the x and y component, as z is not needed)
            xyz_t biasLp = {
                R[0]*_state.gyroBiasRPS.x + R[1]*_state.gyroBiasRPS.y + R[2]*_state.gyroBiasRPS.z,
                R[3]*_state.gyroBiasRPS.x + R[4]*_state.gyroBiasRPS.y + R[5]*_state.gyroBiasRPS.z,
                0.0F
            };
            // low-pass filter R and R*b_hat
            R = _state.motionBiasEstimateR_LPF.filter(R);
            biasLp = _state.motionBiasEstimateBiasLPF.filterXY(biasLp);
            e = R * (-_state.gyroBiasRPS); // negate vector before matrix multiplication for computational efficiency
            e.x += biasLp.x - accEarth.y/deltaT; // NOTE accEarth.y is not a typo
            e.y += biasLp.y + accEarth.x/deltaT; // NOTE accEarth.x is not a typo
            w = { _coeffs.biasMotionW, _coeffs.biasMotionW, _coeffs.biasVerticalW };
        }

        // Kalman filter update
        // STEP 1: P = P + V (also increase covariance if there is no measurement update!)
        if (_state.biasP[0] < _coeffs.biasP0) {
            _state.biasP[0] += _coeffs.biasV;
        }
        if (_state.biasP[4] < _coeffs.biasP0) {
            _state.biasP[4] += _coeffs.biasV;
        }
        if (_state.biasP[8] < _coeffs.biasP0) {
            _state.biasP[8] += _coeffs.biasV;
        }
        if (_state.restDetected) {
            // R is I, identity matrix
            // clip disagreement to -2..2 °/s
            // (this also effectively limits the harm done by the first inclination correction step)
            e.clampInPlace(-_params.biasClipRPS, _params.biasClipRPS);

            // STEP 2: K = P R^T inv(W + R P R^T), simplifies to K =  P inv(W + P)
            // K = W + P
            Matrix3x3 K = _state.biasP;
            // Add w to diagonal of K
            K[0] += w.x;
            K[4] += w.y;
            K[8] += w.z;
            K == _state.biasP * K.inverse();

            // STEP 3: bias = bias + K (y - R bias) = bias + K e
            _state.gyroBiasRPS += K * e;
            _state.gyroBiasRPS.clampInPlace(-_params.biasClipRPS, _params.biasClipRPS); // clip bias estimate to -2..2 °/s

            // STEP 4: P = P - K R P, simplifies to P = P - K P
            _state.biasP -= K * _state.biasP;
        } else if (_params.motionBiasEstimationEnabled) {
            // clip disagreement to -2..2 °/s
            // (this also effectively limits the harm done by the first inclination correction step)
            e.clampInPlace(-_params.biasClipRPS, _params.biasClipRPS);

            // STEP 2: K = P R^T inv(W + R P R^T)
            const Matrix3x3 RT = R.transpose();
            // K = W + R P R^T
            Matrix3x3 K = R * _state.biasP * RT;
            // Add w to diagonal of K
            K[0] += w.x;
            K[4] += w.y;
            K[8] += w.z;
            K == _state.biasP * (RT * K.inverse());

            // STEP 3: bias = bias + K (y - R bias) = bias + K e
            _state.gyroBiasRPS += K * e;
            _state.gyroBiasRPS.clampInPlace(-_params.biasClipRPS, _params.biasClipRPS); // clip bias estimate to -2..2 °/s

            // STEP 4: P = P - K R P
            _state.biasP -= K * R * _state.biasP;
        }
    }
#else
    // simplified implementation of bias estimation for the special case in which only resting bias estimation is enabled
    if (_params.restBiasEstimationEnabled) {
        if (_state.biasP[0] < _coeffs.biasP0) {
            _state.biasP[0] += _coeffs.biasV;
        }
        if (_state.restDetected) {
            xyz_t e = _state.restLastGyro - _state.gyroBiasRPS;
            e.clampInPlace(-_params.biasClipRPS, _params.biasClipRPS);

            // Kalman filter update, simplified scalar version for rest update
            // (this version only uses the first entry of P as P is diagonal and all diagonal elements are the same)
            // step 1: P = P + V (done above!)
            // step 2: K = P R^T inv(W + R P R^T)
            const float k = _state.biasP[0]/(_coeffs.biasRestW + _state.biasP[0]);
            // step 3: bias = bias + K (y - R bias) = bias + K e
            _state.gyroBiasRPS += k*e;
            _state.gyroBiasRPS.clampInPlace(-_params.biasClipRPS, _params.biasClipRPS);
            // step 4: P = P - K R P
            _state.biasP[0] -= k*_state.biasP[0];
        }
    }
#endif
    return _state.orientation6D;
}

#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MAGNETOMETER)
bool VQF::checkForMagneticDisturbance(const xyz_t& magEarth, float deltaT) // NOLINT(readability-make-member-function-const) false positive
{
    _state.magNormDip.x = magEarth.magnitude();
    _state.magNormDip.y = -asinf(magEarth.z/_state.magNormDip.x);

    if (_params.magCurrentTau > 0) {
        //filterVec(_state.magNormDip, 2, _params.magCurrentTau, _coeffs.magTs, _coeffs.magNormDipLpB, _coeffs.magNormDipLpA, _state.magNormDipLpState, _state.magNormDip);
        _state.magNormDip = _state.magNormDipLPF.filterXY(_state.magNormDip);
    }

    // magnetic disturbance detection
    if (fabsf(_state.magNormDip.x - _state.magRefNorm) < _params.magNormThreshold*_state.magRefNorm && fabsf(_state.magNormDip.y - _state.magRefDip) < _params.magDipThresholdRadians) {
        _state.magUndisturbedT += deltaT;
        if (_state.magUndisturbedT >= _params.magMinUndisturbedTime) {
            _state.magDisturbanceDetected = false;
            _state.magRefNorm += _coeffs.kMagRef*(_state.magNormDip.x - _state.magRefNorm);
            _state.magRefDip += _coeffs.kMagRef*(_state.magNormDip.y - _state.magRefDip);
        }
    } else {
        _state.magDisturbanceDetected = true;
        _state.magUndisturbedT = 0.0F;
    }

    // new magnetic field acceptance
    if (fabsf(_state.magNormDip.x - _state.magCandidateNorm) < _params.magNormThreshold*_state.magCandidateNorm && fabsf(_state.magNormDip.y - _state.magCandidateDip) < _params.magDipThresholdRadians) {
        if (sqrtf(_state.restLastGyro.magnitudeSquared()) >= _params.magNewMinGyroRPS) {
            _state.magCandidateT += deltaT;
        }
        _state.magCandidateNorm += _coeffs.kMagRef*(_state.magNormDip.x - _state.magCandidateNorm);
        _state.magCandidateDip += _coeffs.kMagRef*(_state.magNormDip.y - _state.magCandidateDip);

        if (_state.magDisturbanceDetected && (_state.magCandidateT >= _params.magNewTime || (_state.magRefNorm == 0.0F && _state.magCandidateT >= _params.magNewFirstTime))) {
            _state.magDisturbanceDetected = false;
            _state.magRefNorm = _state.magCandidateNorm;
            _state.magRefDip = _state.magCandidateDip;
            _state.magUndisturbedT = _params.magMinUndisturbedTime;
        }
    } else {
        _state.magCandidateT = 0.0F;
        _state.magCandidateNorm = _state.magNormDip.x;
        _state.magCandidateDip = _state.magNormDip.y;
    }
    return _state.magDisturbanceDetected;
}

Quaternion VQF::updateMagnetometer(const xyz_t& magnetometer, float deltaT) // NOLINT(readability-make-member-function-const) false positive
{
   // ignore [0 0 0] samples
    if (magnetometer.x == 0.0F && magnetometer.y == 0.0F && magnetometer.z == 0.0F) {
        return _state.orientation6D;
    }

    // bring magnetometer measurement into 6D earth frame
    const xyz_t magEarth = _state.orientation6D.rotate(magnetometer);

    if (_params.magDisturbanceRejectionEnabled) {
        checkForMagneticDisturbance(magEarth, deltaT);
    }

    // calculate disagreement angle based on current magnetometer measurement
    _state.magDisagreementAngle = atan2f(magEarth.x, magEarth.y);
    // update the disagreement angle and make sure it is in the range [-pi, pi]
    _state.magDisagreementAngle = wrapToPi(_state.magDisagreementAngle - _state.delta);


    // magnetic disturbance rejection (Algorithm 3 MagneticDisturbanceRejection : procedure MagDistRejection)
    float k = _coeffs.kMag;
    if (_state.magDisturbanceDetected) {
        if (_state.magRejectT <= _params.magMaxRejectionTime) {
            _state.magRejectT += deltaT;
            k = 0.0F; // do not perform magnetic correction
        } else {
            k /= _params.magRejectionFactor; // perform magnetic correction with reduced k
        }
    } else {
        // k unaltered ie perform magnetic correction
        _state.magRejectT = std::max(_state.magRejectT - _params.magRejectionFactor*deltaT, 0.0F);
    }

    // ensure fast initial convergence
    if (_state.kMagInit != 0.0F) {
        // make sure that the gain k is at least 1/N, N=1,2,3,... in the first few samples
        if (k < _state.kMagInit) {
            k = _state.kMagInit;
        }

        // iterative expression to calculate 1/N
        _state.kMagInit /= _state.kMagInit + 1.0F;

        // disable if t > tauMag
        if (_state.kMagInit*_params.tauMag < deltaT) {
            _state.kMagInit = 0.0F;
        }
    }

    // first-order filter step
    _state.delta += k*_state.magDisagreementAngle;
    // calculate correction angular rate to facilitate debugging
    _state.lastMagCorrectionAngularRate = k*_state.magDisagreementAngle/deltaT;

    // make sure delta is in the range [-pi, pi]
    _state.delta = wrapToPi(_state.delta);

    _state.orientation9D = _state.orientation6D.rotateZ(_state.delta);
    return _state.orientation9D;
}
#endif

Quaternion VQF::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) // NOLINT(readability-convert-member-functions-to-static)
{
    updateGyro(gyroRPS, deltaT);

#if defined(LIBRARY_SENSOR_FUSION_VQF_USE_MAGNETOMETER)
    updateAccelerometer(accelerometer, deltaT);
    return updateMagnetometer(magnetometer, deltaT);
#else
    (void)magnetometer;
    return updateAccelerometer(accelerometer, deltaT);
#endif
}

Quaternion VQF::updateOrientation(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)
{
    updateGyro(gyroRPS, deltaT);
    return updateAccelerometer(accelerometer, deltaT);
}
