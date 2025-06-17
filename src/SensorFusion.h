#pragma once

#include <Matrix3x3.h>


/*!
Quaternion with added gravity function;
*/
class QuaternionG : public Quaternion {
public:
    explicit QuaternionG(const Quaternion& q) : Quaternion(q) {}
    QuaternionG(float w_, float x_, float y_, float z_) : Quaternion(w_, x_, y_, z_) {}
public:
    inline xyz_t halfGravity() { return xyz_t {  x*z - w*y,  w*x + y*z, -0.5F + w*w + z*z }; } // ENU convention
    inline xyz_t gravity() { return halfGravity()*2.0F; }
};

/*!
SensorFusionFilterBase

The sensor fusion filters use the ENU (East North Up) coordinate frame.

For sensor fusion filters, Euler angles are defined in radians:
1. Roll, denoted by ϕ (phi), is rotation about the X axis
2. Pitch, denoted by θ (theta), is rotation about the Y axis
3. Yaw, denoted by ψ (psi), is rotation about the Z axis
*/
class SensorFusionFilterBase {
public:
    static constexpr float degreesToRadians = static_cast<float>(M_PI / 180.0);
    static constexpr float gToMetersPerSecondSquared = 9.80665F;
    static constexpr float M_PI_F = static_cast<float>(M_PI);
    static constexpr float M_SQRT2_F = static_cast<float>(M_SQRT2);
public:
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) = 0;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) = 0;
    virtual void setFreeParameters(float parameter0, float parameter1);
    virtual bool requiresInitialization() const;
    void reset();
    inline Quaternion getOrientation() const { return Quaternion(q0, q1, q2, q3); }
    Quaternion twoQdot(const xyz_t& gyroRPS) const;
    // Attitude(tilt) from gravity, https://ahrs.readthedocs.io/en/latest/filters/tilt.html#module-ahrs.filters.tilt
     //! Calculate roll (theta) from the normalized accelerometer readings
    static inline float rollRadiansFromAccNormalized(const xyz_t& acc) { return atan2f(acc.y, acc.z); }
    //! Calculate pitch (phi) from the normalized accelerometer readings
    static inline float pitchRadiansFromAccNormalized(const xyz_t& acc) { return atan2f(-acc.x, sqrtf(acc.y*acc.y + acc.z*acc.z)); }
    //! Wrap angle to [-pi, pi]
    static inline float wrapToPi(float theta) { return (theta > M_PI_F) ? theta - 2.0F*M_PI_F : (theta < -M_PI_F) ? theta + 2.0F*M_PI_F : theta; }
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
    float _alpha { 0.96F };
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
    virtual bool requiresInitialization() const override;
    virtual void setFreeParameters(float parameter0, float parameter1) override;
    inline void setBeta(float beta) { setFreeParameters(beta, 0.0f); }
    inline void setGyroMeasurementError(float gyroMeasurementError) { setBeta(gyroMeasurementError * sqrtf(3.0F / 4.0F)); }
private:
    float _beta { 1.0 }; // Initially gain is high, to give fast convergence
};

/*!
Compound Butterworth filters for use by VQF (Versatile Quaternion-based Filter)
*/
class FilterButterworthCompound {
public:
    struct state_t {
        float s0;
        float s1;
    };
    struct coefficients_t {
        float a1;
        float a2;
        float b0;
        float b1;
        float b2;
    };
public:
    FilterButterworthCompound();
    FilterButterworthCompound(float tau, float deltaT);
public:
    void setCoefficients(float tau, float deltaT);
    inline const coefficients_t& getCoefficients() const { return _coefficients; }
    void setState(state_t& state, float x0) const;
    float filterStep(state_t& state, float x) const;
    inline bool getInitialized() const { return _initialized; }
protected:
    coefficients_t _coefficients {};
    int _initialized;
    float _tau;
    float _deltaT;
};

class FilterButterworthXYZ : public FilterButterworthCompound {
public:
    enum {X=0, Y=1, Z=2};
public:
    FilterButterworthXYZ() = default;
    FilterButterworthXYZ(float tau, float deltaT);
public:
    xyz_t filter(const xyz_t& v);
    xyz_t filterXY(const xyz_t& v);
    inline void reset() { _initialized = false; _state.fill({0.0F, 0.0F}); }
    inline const std::array<state_t, 3>& getState() const { return _state; }
protected:
    std::array<state_t, 3> _state {};
};

class FilterButterworthMatrix3x3 : public FilterButterworthCompound {
public:
    FilterButterworthMatrix3x3() = default;
    FilterButterworthMatrix3x3(float tau, float deltaT);
public:
    Matrix3x3 filter(const Matrix3x3& m);
    inline void reset() { _initialized = false; _state.fill({0.0F, 0.0F}); }
    inline const std::array<state_t, 9>& getState() const { return _state; }
protected:
    std::array<state_t, 9> _state {};
};

/*!
Basic VQF (Versatile Quaternion-based Filter)
See [https://arxiv.org/pdf/2203.17024](VQF: Highly Accurate IMU Orientation Estimation with Bias Estimation and Magnetic Disturbance Rejection)

VQF is a Strapdown INS(Inertial Navigation System), that is the inertial sensors (gyroscope and accelerometer) are rigidly fixed to the body of the vehicle.
This contrasts to a Gimbaled INS where the inertial sensors are mounted on a stabilized platform that remains oriented in a fixed direction,
independent of the vehicle’s movements.

BasicVQF consists of gyroscope integration, accelerometer inclination correction, and optional magnetic heading correction.

The full version VQF additionally includes rest detection, gyroscope bias estimation, and magnetic disturbance rejection (which can be enabled or disabled independently).

Code adapted from https://github.com/dlaidig/vqf under MIT license.
*/
class BasicVQF : public SensorFusionFilterBase {
public:
    BasicVQF(float gyroDeltaT, float accDeltaT, float magDeltaT);
    inline explicit BasicVQF(float deltaT) : BasicVQF(deltaT, deltaT, deltaT) {}
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) override;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) override;
    void resetState();
protected:
    void updateGyro(const xyz_t& gyroRPS, float deltaT);
    Quaternion updateAccelerometer(const xyz_t& accelerometer, float deltaT);
    Quaternion updateMagnetometer(const xyz_t& magnetometer, float deltaT);
    static float gainFromTau(float tau, float deltaT);
public:
    struct params_t {
        float tauAcc;
        float tauMag;
    };
    struct coeffs_t {
        float gyroDeltaT;
        float accDeltaT;
        float magDeltaT;
        float kMag;
    };
    struct state_t {
        Quaternion gyroQuaternion {1.0F, 0.0F, 0.0F, 0.0F}; //<! gyro strapdown integration quaternion
        Quaternion accQuaternion {1.0F, 0.0F, 0.0F, 0.0F}; //<! accelerometer correction quaternion
        Quaternion orientation6D {1.0F, 0.0F, 0.0F, 0.0F}; //<! 6D orientation estimate
        Quaternion orientation9D {1.0F, 0.0F, 0.0F, 0.0F}; //<! 9D orientation estimate
        FilterButterworthXYZ accLPF {}; // accelerometer low pass filter
        float delta {0.0F}; // heading difference between the 9D inertial reference frame ξ and the 6D sensor specific almost-inertial reference frame ξi (magnetometer correction angle)
        float magDisagreementAngle {};
        float kMagInit {1.0F};
    };
    const params_t& getParams() { return _params; }
    const coeffs_t& getCoeffs() { return _coeffs; }
    const state_t& getState() { return _state; }
protected:
    const params_t _params;
    const coeffs_t _coeffs;
    state_t _state;
};

/*!
VQF (Versatile Quaternion-based Filter)
See [https://arxiv.org/pdf/2203.17024](VQF: Highly Accurate IMU Orientation Estimation with Bias Estimation and Magnetic Disturbance Rejection)
*/
class VQF : public SensorFusionFilterBase {
public:
    struct params_t {
        float tauAcc;
        bool restBiasEstimationEnabled;
        bool motionBiasEstimationEnabled;
        float biasSigmaInit;
        float biasForgettingTime;
        float biasClipRPS;
#if defined(VQF_MOTION_BIAS_ESTIMATION)
        float biasSigmaMotion;
        float biasVerticalForgettingFactor;
#endif
        float biasSigmaRest;
        float restMinT;
        float restFilterTau;
        float restThresholdGyroSquared;
        float restThresholdAccSquared;
#if defined(USE_MAGNETOMETER)
        float tauMag;
        bool magDisturbanceRejectionEnabled;
        float magCurrentTau;
        float magRefTau;
        float magNormThreshold;
        float magDipThresholdRadians;
        float magNewTime;
        float magNewFirstTime;
        float magNewMinGyroRPS;
        float magMinUndisturbedTime;
        float magMaxRejectionTime;
        float magRejectionFactor;
#endif
    };
    struct coeffs_t {
        float gyroDeltaT;
        float accDeltaT;
#if defined(USE_MAGNETOMETER)
        float magDeltaT;
        float kMag;
        float kMagRef;
#endif
        float biasP0;
        float biasV;
#if defined(VQF_MOTION_BIAS_ESTIMATION)
        float biasMotionW;
        float biasVerticalW;
#endif
        float biasRestW;
    };
    struct state_t {
        Quaternion gyroQuaternion {1.0F, 0.0F, 0.0F, 0.0F}; //<! gyro strapdown integration quaternion
        Quaternion accQuaternion {1.0F, 0.0F, 0.0F, 0.0F}; //<! accelerometer correction quaternion
        Quaternion orientation6D {1.0F, 0.0F, 0.0F, 0.0F}; //<! 6D orientation estimate
        Quaternion orientation9D {1.0F, 0.0F, 0.0F, 0.0F}; //<! 9D orientation estimate
        FilterButterworthXYZ accLPF {}; //<! accelerometer low pass filter
        float lastAccCorrectionAngularRate {}; // for debug
        xyz_t gyroBiasRPS { 0.0F, 0.0F, 0.0F }; // called bias in Laidig and Seel
        Matrix3x3 biasP {}; //<! Diagonal covariance matrix
#if defined(VQF_MOTION_BIAS_ESTIMATION)
        FilterButterworthMatrix3x3 motionBiasEstimateR_LPF {}; //<! Low-pass filter for rotation matrix
        FilterButterworthXYZ motionBiasEstimateBiasLPF {};
#endif
        bool restDetected {false};
        float restT {};
        xyz_t restLastGyro { 0.0F, 0.0F, 0.0F };
        float restLastGyroSquaredDeviation {};
        FilterButterworthXYZ restGyroLPF {};
        xyz_t restLastAcc { 0.0F, 0.0F, 0.0F };
        float restLastAccSquaredDeviation {};
        FilterButterworthXYZ restAccLPF {};
#if defined(USE_MAGNETOMETER)
        float delta {0.0F}; //<! δi, the heading difference between the 9D inertial reference frame ξ and the 6D sensor specific almost-inertial reference frame ξi (magnetometer correction angle)
        float kMagInit {1.0F};
        float lastMagCorrectionAngularRate {}; // for debug
        bool magDisturbanceDetected {false};
        float magDisagreementAngle {};
        float magRefNorm {};
        float magRefDip {};
        float magUndisturbedT {};
        float magRejectT {};
        float magCandidateNorm {};
        float magCandidateDip {};
        float magCandidateT {};
        xyz_t magNormDip { 0.0F, 0.0F, 0.0F };
        FilterButterworthXYZ magNormDipLPF {};
#endif
    };
    const params_t& getParams() { return _params; }
    const coeffs_t& getCoeffs() { return _coeffs; }
    const state_t& getState() { return _state; }
public:
    VQF(float gyroDeltaT, float accDeltaT, float magDeltaT, bool restBiasEstimationEnabled, bool motionBiasEstimationEnabled, bool magDisturbanceRejectionEnabled);
    VQF(float gyroDeltaT, float accDeltaT, float magDeltaT);
    inline explicit VQF(float deltaT) : VQF(deltaT, deltaT, deltaT) {}
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, float deltaT) override;
    virtual Quaternion update(const xyz_t& gyroRPS, const xyz_t& accelerometer, const xyz_t& magnetometer, float deltaT) override;
    void resetState();
    static float gainFromTau(float tau, float deltaT);
    float calculateBias(float v) const;
protected:
    static inline float square(float x) { return x*x; }
    void updateGyro(const xyz_t& gyroRPS, float deltaT);
    Quaternion updateAccelerometer(const xyz_t& accelerometer, float deltaT);
#if defined(USE_MAGNETOMETER)
    bool checkForMagneticDisturbance(const xyz_t& magEarth, float deltaT);
    Quaternion updateMagnetometer(const xyz_t& magnetometer, float deltaT);
#endif
private:
    const params_t _params;
    const coeffs_t _coeffs;
    state_t _state {};
};
