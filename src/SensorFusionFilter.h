#pragma once

#include <Quaternion.h>
#include <xyz_type.h>

/*!
Quaternion with added gravity function;
*/
class QuaternionG : public Quaternion {
public:
    explicit QuaternionG(const Quaternion& q) : Quaternion(q) {}
    QuaternionG(float w, float x, float y, float z) : Quaternion(w, x, y, z) {}
public:
    inline xyz_t halfGravity() {
        return xyz_t { .x = x*z - w*y, .y = w*x + y*z, .z = w*w - 0.5F + z*z };
    }
    inline xyz_t gravity() { return halfGravity()*2.0F; }
};

/*!
SensorFusionFilterBase

For Sensor fusion filters, Euler angles are defined in radians:
1. Roll, denoted by ϕ (phi), is rotation about the X axis
2. Pitch, denoted by θ (theta), is rotation about the Y axis
3. Yaw, denoted by ψ (psi), is rotation about the Z axis
*/
class SensorFusionFilterBase {
public:
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) = 0;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) = 0;
    virtual void setFreeParameters(float parameter0, float parameter1) = 0;
    void reset();
    inline Quaternion getOrientation() const { return Quaternion(q0, q1, q2, q3); }
    Quaternion twoQdot(const xyz_t& gyroRPS) const;
    // Attitude(tilt) from gravity, https://ahrs.readthedocs.io/en/latest/filters/tilt.html#module-ahrs.filters.tilt
     //! Calculate roll (theta) from the normalized accelerometer readings
    static inline float rollRadiansFromAccNormalized(const xyz_t& acc) { return atan2(acc.y, acc.z); }
    //! Calculate pitch (phi) from the normalized accelerometer readings
    static inline float pitchRadiansFromAccNormalized(const xyz_t& acc) { return atan2(-acc.x, sqrt(acc.y*acc.y + acc.z*acc.z)); }
public: // functions used for unit testing
    void _setAndNormalizeQ(float q0_, float q1_, float q2_, float q3_);
protected:
    // orientation quaternion
    float q0 { 1.0 };
    float q1 { 0.0 };
    float q2 { 0.0 };
    float q3 { 0.0 };
    float _accMagnitudeSquaredMax {4.0};
};

/*!
ComplementaryFilter
*/
class ComplementaryFilter : public SensorFusionFilterBase {
public:
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) override;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) override;
    virtual void setFreeParameters(float parameter0, float parameter1) override;
    inline void setAlpha(float alpha) { setFreeParameters(alpha, 0.0f); }
private:
    float _alpha { 0.96 };
};

/*!
MahonyFilter
*/
class MahonyFilter : public SensorFusionFilterBase {
public:
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) override;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) override;
    virtual void setFreeParameters(float parameter0, float parameter1) override;
    inline void setKpKi(float kp, float ki) { setFreeParameters(kp, ki); }
private:
    float _kp { 10.0 };
    float _ki { 0.0 };
    xyz_t _errorIntegral { 0.0, 0.0, 0.0 };
};

/*!
MadgwickFilter
*/
class MadgwickFilter : public SensorFusionFilterBase {
public:
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) override;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) override;
    virtual void setFreeParameters(float parameter0, float parameter1) override;
    inline void setBeta(float beta) { setFreeParameters(beta, 0.0f); }
    inline void setGyroMeasurementError(float gyroMeasurementError) { setBeta(gyroMeasurementError * sqrtf(3.0F / 4.0F)); }
private:
    float _beta { 1.0 }; // Initially gain is high, to give fast convergence
};
