#include "SensorFusion.h"
#include <unity.h>


void setUp() {
}

void tearDown() {
}


void test_filter_butterworth_state()
{
    FilterButterworthXYZ f; // NOLINT(misc-const-correctness) false positive
    const float tau = 3.0F;
    const float deltaT = 0.01F;
    f.setCoefficients(tau, deltaT);
    const FilterButterworthCompound::coefficients_t coefficients = f.getCoefficients();

    FilterButterworthCompound::state_t state{}; // NOLINT(misc-const-correctness) false positive
    f.setState(state, 1.0F);
    TEST_ASSERT_EQUAL_FLOAT(1.0F - coefficients.b0, state.s0);
    TEST_ASSERT_EQUAL_FLOAT(coefficients.b2 - coefficients.a2, state.s1);
}

void test_filter_butterworth_xyz_1()
{
    FilterButterworthXYZ f; // NOLINT(misc-const-correctness) false positive
    const float tau = 2.0F;
    const float deltaT = 1.0F;
    f.setCoefficients(tau, deltaT);

    const xyz_t v0 = {2.0F, 3.0F, 5.0F};
    const xyz_t v0f = f.filter(v0);
    // first filter step just returns value unchanged
    TEST_ASSERT_FALSE(f.getInitialized());
    TEST_ASSERT_TRUE(v0 == v0f);
    TEST_ASSERT_EQUAL_FLOAT(v0.x, v0f.x);
    TEST_ASSERT_EQUAL_FLOAT(v0.y, v0f.y);
    TEST_ASSERT_EQUAL_FLOAT(v0.z, v0f.z);

    const xyz_t v1 = {4.0F, 5.0F, 7.0F};
    const xyz_t v1f = f.filter(v1);
    const xyz_t v1fe = (v0 + v1) / 2.0F;
    // second filter step is average of first two
    TEST_ASSERT_TRUE(f.getInitialized());
    TEST_ASSERT_TRUE(v1fe == v1f);
}

void test_filter_butterworth_xyz_2()
{
    FilterButterworthXYZ f; // NOLINT(misc-const-correctness) false positive
    const float tau = 2.0F;
    const float deltaT = 1.0F;
    f.setCoefficients(tau, deltaT);
    const xyz_t v0 = {2.0F, 3.0F, 5.0F};

    // first filter step just returns value unchanged
    const xyz_t v0f = f.filter(v0);
    TEST_ASSERT_FALSE(f.getInitialized());
    TEST_ASSERT_TRUE(v0 == v0f);

    // repeated filter of the same value should leave it unchanged
    const xyz_t v1f = f.filter(v0);
    TEST_ASSERT_TRUE(f.getInitialized());
    TEST_ASSERT_TRUE(v0 == v1f);

    // repeated filter of the same value should leave it unchanged
    const xyz_t v2f = f.filter(v0);
    TEST_ASSERT_TRUE(f.getInitialized());
    TEST_ASSERT_EQUAL_FLOAT(v0.x, v2f.x);
    TEST_ASSERT_EQUAL_FLOAT(v0.y, v2f.y);
    TEST_ASSERT_EQUAL_FLOAT(v0.z, v2f.z);

    const xyz_t v3 = {4.0F, 5.0F, 7.0F};
    std::array<FilterButterworthCompound::state_t, 3> state = f.getState();
    const xyz_t v3f = f.filter(v3);
    TEST_ASSERT_TRUE(f.getInitialized());

    const xyz_t v3fe = {
        f.filterStep(state[FilterButterworthXYZ::X], v3.x),
        f.filterStep(state[FilterButterworthXYZ::Y], v3.y),
        f.filterStep(state[FilterButterworthXYZ::Z], v3.z)
    };
    TEST_ASSERT_EQUAL_FLOAT(v3fe.x, v3f.x);
    TEST_ASSERT_EQUAL_FLOAT(v3fe.y, v3f.y);
    TEST_ASSERT_EQUAL_FLOAT(v3fe.z, v3f.z);
}

void test_filter_butterworth_matrix3x3_1()
{
    FilterButterworthMatrix3x3 f;
    const float tau = 2.0F;
    const float deltaT = 1.0F;
    f.setCoefficients(tau, deltaT);

    const Matrix3x3 m0( 2,  3,  5,  7, 11, 13, 17, 19, 23);

    // first filter step just returns value unchanged
    const Matrix3x3 m0f = f.filter(m0);
    TEST_ASSERT_FALSE(f.getInitialized());
    TEST_ASSERT_TRUE(m0 == m0f);

    // second filter step is average of first two
    const Matrix3x3 m1(29, 31, 37, 41, 43, 47, 53, 59, 61);
    const Matrix3x3 m1f = f.filter(m1);
    const Matrix3x3 m1fe = (m0 + m1) / 2.0F;
    TEST_ASSERT_TRUE(f.getInitialized());
    TEST_ASSERT_TRUE(m1fe == m1f);
}

void test_filter_butterworth_matrix3x3_2()
{
    FilterButterworthMatrix3x3 f;
    const float tau = 2.0F;
    const float deltaT = 1.0F;
    f.setCoefficients(tau, deltaT);

    const Matrix3x3 m0( 2,  3,  5,  7, 11, 13, 17, 19, 23);

    // first filter step just returns value unchanged
    const Matrix3x3 m0f = f.filter(m0);
    TEST_ASSERT_FALSE(f.getInitialized());
    TEST_ASSERT_TRUE(m0 == m0f);

    // repeated filter of the same value should leave it unchanged
    const Matrix3x3 m1f = f.filter(m0);
    TEST_ASSERT_TRUE(f.getInitialized());
    TEST_ASSERT_TRUE(m0 == m1f);

    // repeated filter of the same value should leave it unchanged
    const Matrix3x3 m2f = f.filter(m0);
    TEST_ASSERT_TRUE(f.getInitialized());
    TEST_ASSERT_EQUAL_FLOAT(m0[0], m2f[0]);
    TEST_ASSERT_EQUAL_FLOAT(m0[1], m2f[1]);
    TEST_ASSERT_EQUAL_FLOAT(m0[2], m2f[2]);
    TEST_ASSERT_EQUAL_FLOAT(m0[3], m2f[3]);
    TEST_ASSERT_EQUAL_FLOAT(m0[4], m2f[4]);
    TEST_ASSERT_EQUAL_FLOAT(m0[5], m2f[5]);
    TEST_ASSERT_EQUAL_FLOAT(m0[6], m2f[6]);
    TEST_ASSERT_EQUAL_FLOAT(m0[7], m2f[7]);
    TEST_ASSERT_EQUAL_FLOAT(m0[8], m2f[8]);

    const Matrix3x3 m3(29, 31, 37, 41, 43, 47, 53, 59, 61);
    std::array<FilterButterworthCompound::state_t, 9> state = f.getState();
    const Matrix3x3 m3f = f.filter(m3);
    TEST_ASSERT_TRUE(f.getInitialized());

    const Matrix3x3 m3fe = {
        f.filterStep(state[0], m3[0]),
        f.filterStep(state[1], m3[1]),
        f.filterStep(state[2], m3[2]),
        f.filterStep(state[3], m3[3]),
        f.filterStep(state[4], m3[4]),
        f.filterStep(state[5], m3[5]),
        f.filterStep(state[6], m3[6]),
        f.filterStep(state[7], m3[7]),
        f.filterStep(state[8], m3[8]),
    };
    TEST_ASSERT_EQUAL_FLOAT(m3fe[0], m3f[0]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[1], m3f[1]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[2], m3f[2]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[3], m3f[3]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[4], m3f[4]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[5], m3f[5]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[6], m3f[6]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[7], m3f[7]);
    TEST_ASSERT_EQUAL_FLOAT(m3fe[8], m3f[8]);
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
    UNITY_BEGIN();

    RUN_TEST(test_filter_butterworth_state);
    RUN_TEST(test_filter_butterworth_xyz_1);
    RUN_TEST(test_filter_butterworth_xyz_2);
    RUN_TEST(test_filter_butterworth_matrix3x3_1);
    RUN_TEST(test_filter_butterworth_matrix3x3_2);

    UNITY_END();
}
