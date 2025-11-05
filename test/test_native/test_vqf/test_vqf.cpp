#include "SensorFusion.h"
#include <unity.h>


void setUp() {
}

void tearDown() {
}

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays,modernize-use-using,readability-non-const-parameter,cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
typedef float vqf_real_t;
#define EPS std::numeric_limits<vqf_real_t>::epsilon()

void matrix3SetToScaledIdentity(vqf_real_t scale, vqf_real_t out[9])
{
    out[0] = scale;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = scale;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = scale;
}

void matrix3Multiply(const vqf_real_t in1[9], const vqf_real_t in2[9], vqf_real_t out[9])
{
    vqf_real_t tmp[9];
    tmp[0] = in1[0]*in2[0] + in1[1]*in2[3] + in1[2]*in2[6];
    tmp[1] = in1[0]*in2[1] + in1[1]*in2[4] + in1[2]*in2[7];
    tmp[2] = in1[0]*in2[2] + in1[1]*in2[5] + in1[2]*in2[8];
    tmp[3] = in1[3]*in2[0] + in1[4]*in2[3] + in1[5]*in2[6];
    tmp[4] = in1[3]*in2[1] + in1[4]*in2[4] + in1[5]*in2[7];
    tmp[5] = in1[3]*in2[2] + in1[4]*in2[5] + in1[5]*in2[8];
    tmp[6] = in1[6]*in2[0] + in1[7]*in2[3] + in1[8]*in2[6];
    tmp[7] = in1[6]*in2[1] + in1[7]*in2[4] + in1[8]*in2[7];
    tmp[8] = in1[6]*in2[2] + in1[7]*in2[5] + in1[8]*in2[8];
    std::copy(tmp, tmp+9, out);
}

void matrix3MultiplyTpsFirst(const vqf_real_t in1[9], const vqf_real_t in2[9], vqf_real_t out[9])
{
    vqf_real_t tmp[9];
    tmp[0] = in1[0]*in2[0] + in1[3]*in2[3] + in1[6]*in2[6];
    tmp[1] = in1[0]*in2[1] + in1[3]*in2[4] + in1[6]*in2[7];
    tmp[2] = in1[0]*in2[2] + in1[3]*in2[5] + in1[6]*in2[8];
    tmp[3] = in1[1]*in2[0] + in1[4]*in2[3] + in1[7]*in2[6];
    tmp[4] = in1[1]*in2[1] + in1[4]*in2[4] + in1[7]*in2[7];
    tmp[5] = in1[1]*in2[2] + in1[4]*in2[5] + in1[7]*in2[8];
    tmp[6] = in1[2]*in2[0] + in1[5]*in2[3] + in1[8]*in2[6];
    tmp[7] = in1[2]*in2[1] + in1[5]*in2[4] + in1[8]*in2[7];
    tmp[8] = in1[2]*in2[2] + in1[5]*in2[5] + in1[8]*in2[8];
    std::copy(tmp, tmp+9, out);
}

void matrix3MultiplyTpsSecond(const vqf_real_t in1[9], const vqf_real_t in2[9], vqf_real_t out[9])
{
    vqf_real_t tmp[9];
    tmp[0] = in1[0]*in2[0] + in1[1]*in2[1] + in1[2]*in2[2];
    tmp[1] = in1[0]*in2[3] + in1[1]*in2[4] + in1[2]*in2[5];
    tmp[2] = in1[0]*in2[6] + in1[1]*in2[7] + in1[2]*in2[8];
    tmp[3] = in1[3]*in2[0] + in1[4]*in2[1] + in1[5]*in2[2];
    tmp[4] = in1[3]*in2[3] + in1[4]*in2[4] + in1[5]*in2[5];
    tmp[5] = in1[3]*in2[6] + in1[4]*in2[7] + in1[5]*in2[8];
    tmp[6] = in1[6]*in2[0] + in1[7]*in2[1] + in1[8]*in2[2];
    tmp[7] = in1[6]*in2[3] + in1[7]*in2[4] + in1[8]*in2[5];
    tmp[8] = in1[6]*in2[6] + in1[7]*in2[7] + in1[8]*in2[8];
    std::copy(tmp, tmp+9, out);
}

bool matrix3Inv(const vqf_real_t in[9], vqf_real_t out[9])
{
    // in = [a b c; d e f; g h i]
    const vqf_real_t A = in[4]*in[8] - in[5]*in[7]; // (e*i - f*h)
    const vqf_real_t D = in[2]*in[7] - in[1]*in[8]; // -(b*i - c*h)
    const vqf_real_t G = in[1]*in[5] - in[2]*in[4]; // (b*f - c*e)
    const vqf_real_t B = in[5]*in[6] - in[3]*in[8]; // -(d*i - f*g)
    const vqf_real_t E = in[0]*in[8] - in[2]*in[6]; // (a*i - c*g)
    const vqf_real_t H = in[2]*in[3] - in[0]*in[5]; // -(a*f - c*d)
    const vqf_real_t C = in[3]*in[7] - in[4]*in[6]; // (d*h - e*g)
    const vqf_real_t F = in[1]*in[6] - in[0]*in[7]; // -(a*h - b*g)
    const vqf_real_t I = in[0]*in[4] - in[1]*in[3]; // (a*e - b*d)

    const vqf_real_t det = in[0]*A + in[1]*B + in[2]*C; // a*A + b*B + c*C;

    if (det >= -EPS && det <= EPS) {
        std::fill(out, out+9, 0);
        return false;
    }

    // out = [A D G; B E H; C F I]/det
    out[0] = A/det;
    out[1] = D/det;
    out[2] = G/det;
    out[3] = B/det;
    out[4] = E/det;
    out[5] = H/det;
    out[6] = C/det;
    out[7] = F/det;
    out[8] = I/det;

    return true;
}

void test_kalman_acc()
{
    // step 2: K = P R^T inv(W + R P R^T)
    vqf_real_t eR[9] = {  2,  3,  5,  7, 11, 13, 17, 19, 23 };
    vqf_real_t eBiasP[9] = {29, 31, 37, 41, 43, 47, 53, 59, 61};
    vqf_real_t eK[9];
    matrix3MultiplyTpsSecond(eBiasP, eR, eK); // K = P R^T

    const Matrix3x3 R(eR);
    Matrix3x3 biasP(eBiasP);
    const Matrix3x3 RT = R.transpose();
    Matrix3x3 K = biasP * RT;
    TEST_ASSERT_TRUE(Matrix3x3(eK) == K);

    matrix3Multiply(eR, eK, eK); // K = R P R^T
    K = R * K;
    TEST_ASSERT_TRUE(Matrix3x3(eK) == K);
    TEST_ASSERT_TRUE(K == R * biasP * RT);

    const vqf_real_t ew[3] = { 2, 3, 5};
    eK[0] += ew[0];
    eK[4] += ew[1];
    eK[8] += ew[2]; // K = W + R P R^T
    const xyz_t w{2.0F, 3.0F, 5.0F};
    K[0] += w.x;
    K[4] += w.y;
    K[8] += w.z;
    TEST_ASSERT_TRUE(Matrix3x3(eK) == K);
    const Matrix3x3 KSave = K;

    matrix3Inv(eK, eK); // K = inv(W + R P R^T)
    K.invertInPlace();
    TEST_ASSERT_TRUE(Matrix3x3(eK) == K);

    matrix3MultiplyTpsFirst(eR, eK, eK); // K = R^T inv(W + R P R^T)
    K = RT * K;
    TEST_ASSERT_TRUE(Matrix3x3(eK) == K);
    TEST_ASSERT_TRUE(K == RT * KSave.inverse());

    matrix3Multiply(eBiasP, eK, eK); // K = P R^T inv(W + R P R^T)
    K = biasP * K;
    TEST_ASSERT_TRUE(Matrix3x3(eK) == K);
    TEST_ASSERT_TRUE(K == biasP * (RT * KSave.inverse()));

    // STEP 3: bias = bias + K (y - R bias) = bias + K e
    // STEP 4: P = P - K R P
    matrix3Multiply(eK, eR, eK); // K = K R
    matrix3Multiply(eK, eBiasP, eK); // K = K R P
    for(size_t i = 0; i < 9; i++) { // NOLINT(altera-unroll-loops)
        eBiasP[i] -= eK[i]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
    }
    biasP -= K * R * biasP;
    TEST_ASSERT_TRUE(Matrix3x3(eBiasP) == biasP);
}
// NOLINTEND(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays,modernize-use-using,readability-non-const-parameter,cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

void test_vqf()
{
    const xyz_t gyro0 { .x = 0.0F, .y = 0.0F, .z = 0.0F };
    const float deltaT = 0.1F;

    static VQF vqf(deltaT, deltaT, deltaT);
    const Quaternion q = vqf.updateOrientation(gyro0, xyz_t { .x = 0.0, .y = 0.0, .z = 1.0 }, deltaT);
    TEST_ASSERT_EQUAL_FLOAT(0, q.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(0, q.calculatePitchDegrees());

}

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
    UNITY_BEGIN();

    RUN_TEST(test_vqf);
    RUN_TEST(test_kalman_acc);

    UNITY_END();
}
