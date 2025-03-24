# Refactoring of Madgwick filter update

The Madgwick filter implementation in `MadgwickFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT)` has been refactored to
improve computational efficiency.

The steps taken in that refactoring are outlined here.

A similar refactoring was done for `MadgwickFilter::update(const xyz_t& gyroRPS, const xyz_t& accelerometer, xyz_t& magnetometer, float deltaT)`.

## Derivation

The original code as used in many implementations (Arduino, Adafruit, M5Stack, Reefwing-AHRS):

```cpp
s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
```

Reorder terms:

```cpp
s0 = _4q0 * q2q2 + _4q0 * q1q1                                                     + _2q2 * ax - _2q1 * ay;
s1 = _4q1 * q3q3 + 4.0f * q0q0 * q1 - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az - _2q3 * ax - _2q0 * ay;
s2 = 4.0f * q0q0 * q2 + _4q2 * q3q3 - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az + _2q0 * ax - _2q3 * ay;
s3 = 4.0f * q1q1 * q3  + 4.0f * q2q2 * q3                                          - _2q1 * ax - _2q2 * ay;
```

Factor out `4 * q<n>`:

```cpp
s0 = 4 * q0 * (q2q2 + q1q1)                            + _2q2 * ax - _2q1 * ay;
s1 = 4 * q1 * (q3q3 + q0q0 - 1 + 2*q1q1 + 2*q2q2 + az) - _2q3 * ax - _2q0 * ay;
s2 = 4 * q2 * (q0q0 + q3q3 - 1 + 2*q1q1 + 2*q2q2 + az) + _2q0 * ax - _2q3 * ay;
s3 = 4 * q3 * (q1q1 + q2q2)                            - _2q1 * ax - _2q2 * ay;
```

Remove common factor of 2:

```cpp
s0 = 2 * (q0 * 2 * (q2q2 + q1q1)                            + q2 * ax - q1 * ay);
s1 = 2 * (2 * q1 * (q3q3 + q0q0 - 1 + 2*q1q1 + 2*q2q2 + az) - q3 * ax - q0 * ay);
s2 = 2 * (2 * q2 * (q0q0 + q3q3 - 1 + 2*q1q1 + 2*q2q2 + az) + q0 * ax - q3 * ay);
s3 = 2 * (q3 * 2 * (q1q1 + q2q2)                            - q1 * ax - q2 * ay);
```

Substitute:
`_2q1q1_plus_2q2q2` = `2.0f * (q1*q1 + q2*q2)`:

```cpp
s0 = 2 * (q0 * _2q1q1_plus_2q2q2                              + q2 * ax - q1 * ay);
s1 = 2 * (q1 * 2 * (q3q3 + q0q0 - 1 + _2q1q1_plus_2q2q2 + az) - q3 * ax - q0 * ay);
s2 = 2 * (q2 * 2 * (q0q0 + q3q3 - 1 + _2q1q1_plus_2q2q2 + az) + q0 * ax - q3 * ay);
s3 = 2 * (q3 * _2q1q1_plus_2q2q2                              - q1 * ax - q2 * ay);
```

Substitute:
`common` = `2.0f * (q0*q0 + q3*q3 - 1.0f + _2q1q1_plus_2q2q2 + az)`:

```cpp
s0 = 2 * (q0 * _2q1q1_plus_2q2q2 + q2 * ax - q1 * ay);
s1 = 2 * (q1 * common            - q3 * ax - q0 * ay);
s2 = 2 * (q2 * common            + q0 * ax - q3 * ay);
s3 = 2 * (q3 * _2q1q1_plus_2q2q2 - q1 * ax - q2 * ay);
```

Instead, calculate half S-values:

```cpp
halfS0 = q0 * _2q1q1_plus_2q2q2 + q2 * ax - q1 * ay;
halfS1 = q1 * common            - q3 * ax - q0 * ay;
halfS2 = q2 * common            + q0 * ax - q3 * ay;
halfS3 = q3 * _2q1q1_plus_2q2q2 - q1 * ax - q2 * ay;
```

Observe that:
`s0 * invSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3)` == `halfS0 * invSqrt(halfS0*halfS0 + halfS1*halfS1 + halfS2*halfS2 + halfS3*halfS3)`

(or alternatively note that the norm of a quaternion is the same as the norm of half of that quaternion).
So, for the feedback step, we have:

```cpp
qDot0 -= beta * halfS0 * normReciprocal;
qDot1 -= beta * halfS1 * normReciprocal;
qDot2 -= beta * halfS2 * normReciprocal;
qDot3 -= beta * halfS3 * normReciprocal;
```

Rename `halfS<n>` back to `s<n>`, since `s` is being normalized and they are therefore equivalent:

```cpp
qDot0 -= beta * s0 * normReciprocal;
qDot1 -= beta * s1 * normReciprocal;
qDot2 -= beta * s2 * normReciprocal;
qDot3 -= beta * s3 * normReciprocal;
```

Substituting `betaNormReciprocal` = `beta * reciprocalSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3)`
we have:

```cpp
qDot0 -= s0 * betaNormReciprocal;
qDot1 -= s1 * betaNormReciprocal;
qDot2 -= s2 * betaNormReciprocal;
qDot3 -= s3 * betaNormReciprocal;
```

Substituting `_2betaNormReciprocal` = `2 * beta * reciprocalSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3)`
and doubling everything we have:

```cpp
_2qDot0 -= s0 * betaNormReciprocal;
_2qDot1 -= s1 * betaNormReciprocal;
_2qDot2 -= s2 * betaNormReciprocal;
_2qDot3 -= s3 * betaNormReciprocal;
```

Add in double the derivative from the gyro:

```cpp
_2qDot0 = -q1*gyroRPS.x - q2*gyroRPS.y - q3*gyroRPS.z - s0 * _2betaNormReciprocal;
_2qDot1 =  q0*gyroRPS.x + q2*gyroRPS.z - q3*gyroRPS.y - s1 * _2betaNormReciprocal;
_2qDot2 =  q0*gyroRPS.y - q1*gyroRPS.z + q3*gyroRPS.x - s2 * _2betaNormReciprocal;
_2qDot3 =  q0*gyroRPS.z + q1*gyroRPS.y - q2*gyroRPS.x - s3 * _2betaNormReciprocal;
```

Finally, we update the attitude quaternion using simple Euler integration `qNew = qOld + qDot*deltaT`.
To improve computation efficiency by avoiding multiplications, we use  `_2qDot` and `deltaT*0.5`,
ie `qNew = qOld + _2qDot*halfDeltaT`, that is `q += _2qDot*halfDeltaT`:

```cpp
halfDeltaT = deltaT * 0.5F;
q0 += _2qDot0 * halfDeltaT;
q1 += _2qDot1 * halfDeltaT;
q2 += _2qDot2 * halfDeltaT;
q3 += _2qDot3 * halfDeltaT;
```
