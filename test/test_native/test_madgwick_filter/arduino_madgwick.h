#pragma once

class Madgwick {
public:
    Madgwick() = default;
    void begin(float sampleFrequency, float beta_in) { invSampleFreq = 1.0f / sampleFrequency; beta = beta_in; }
    void _setAndNormalizeQ(float q0_in, float q1_in, float q2_in, float q3_in);
    void updateIMU(float gx, float gy, float gz, float ax, float ay, float az);
private:
    float invSampleFreq {};
    float beta {};
public:
    char anglesComputed {0};
    float q0 {1.0};
    float q1 {0.0};
    float q2 {0.0};
    float q3 {0.0};
};
