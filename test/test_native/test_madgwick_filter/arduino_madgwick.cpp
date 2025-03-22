// https://github.com/arduino-libraries/MadgwickAHRS

#include "arduino_madgwick.h"
#include <cmath>
#include <cstdint>

inline static float invSqrt(float x)
{
#if defined(FAST_INVERSE_SQUARE_ROOT)
    //See: http://en.wikipedia.org/wiki/Fast_inverse_square_root
    const float halfX = 0.5F*x;
    float y = x;
    // cppcheck-suppress invalidPointerCast
    int32_t i = *reinterpret_cast<int32_t*>(&y); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    i = 0x5f3759dF - (i>>1);
    // cppcheck-suppress invalidPointerCast
    y = *reinterpret_cast<float*>(&i); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    y *= 1.5F - halfX*y*y;
    y *= 1.5F - halfX*y*y; // Second iteration
    return y;
#else
    return 1.0F / sqrtf(x);
#endif
}

// 88 multiplications, 27 additions, 17 subtractions, 3 comparisons, 13 auxiliary variables

void Madgwick::updateIMU(float gx, float gy, float gz, float ax, float ay, float az)
{
    float recipNorm; // cppcheck-suppress variableScope
    float s0, s1, s2, s3; // cppcheck-suppress variableScope NOLINT(readability-isolate-declaration)
    float qDot1, qDot2, qDot3, qDot4; // cppcheck-suppress variableScope NOLINT(readability-isolate-declaration)
    float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3; // cppcheck-suppress variableScope NOLINT(readability-isolate-declaration)

    // Convert gyroscope degrees/sec to radians/sec
    gx *= 0.0174533F;
    gy *= 0.0174533F;
    gz *= 0.0174533F;

    // Rate of change of quaternion from gyroscope
    qDot1 = 0.5F * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5F * (q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5F * (q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5F * (q0 * gz + q1 * gy - q2 * gx);

    // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalization)
    if(!((ax == 0.0F) && (ay == 0.0F) && (az == 0.0F))) { // NOLINT(readability-simplify-boolean-expr)

        // Normalize accelerometer measurement
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // Auxiliary variables to avoid repeated arithmetic
        _2q0 = 2.0F * q0;
        _2q1 = 2.0F * q1;
        _2q2 = 2.0F * q2;
        _2q3 = 2.0F * q3;
        _4q0 = 4.0F * q0;
        _4q1 = 4.0F * q1;
        _4q2 = 4.0F * q2;
        _8q1 = 8.0F * q1;
        _8q2 = 8.0F * q2;
        q0q0 = q0 * q0;
        q1q1 = q1 * q1;
        q2q2 = q2 * q2;
        q3q3 = q3 * q3;

        // Gradient decent algorithm corrective step
        s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
        s1 = _4q1 * q3q3 - _2q3 * ax + 4.0F * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
        s2 = 4.0F * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
        s3 = 4.0F * q1q1 * q3 - _2q1 * ax + 4.0F * q2q2 * q3 - _2q2 * ay;
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalize step magnitude
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;

        // Apply feedback step
        qDot1 -= beta * s0;
        qDot2 -= beta * s1;
        qDot3 -= beta * s2;
        qDot4 -= beta * s3;
    }

    // Integrate rate of change of quaternion to yield quaternion
    q0 += qDot1 * invSampleFreq;
    q1 += qDot2 * invSampleFreq;
    q2 += qDot3 * invSampleFreq;
    q3 += qDot4 * invSampleFreq;

    // Normalize quaternion
    recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;
    anglesComputed = 0;
}

void Madgwick::_setAndNormalizeQ(float q0_, float q1_, float q2_, float q3_)
{
    q0 = q0_;
    q1 = q1_;
    q2 = q2_;
    q3 = q3_;
    const float normReciprocal = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= normReciprocal;
    q1 *= normReciprocal;
    q2 *= normReciprocal;
    q3 *= normReciprocal;

    anglesComputed = 1;
}

