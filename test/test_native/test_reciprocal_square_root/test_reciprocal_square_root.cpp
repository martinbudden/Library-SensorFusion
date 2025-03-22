#include <cmath>
#include <unity.h>

void setUp() {
}

void tearDown() {
}

inline static float reciprocalSqrt1(float x)
{
    const float halfX = 0.5F * x;
    union {
        float y;
        long i;
    } u = {x};

    u.i = 0x5f3759dF - (u.i >> 1); // Initial estimate for Newton–Raphson method
    u.y *= 1.5F - (halfX * u.y * u.y); // Single interation

    return u.y;
}

inline static float reciprocalSqrt2(float x)
{
    const float halfX = 0.5F * x;
    union {
        float y;
        long i;
    } u = {x};

    u.i = 0x5f3759dF - (u.i >> 1); // Initial estimate for Newton–Raphson method
    u.y *= 1.5F - (halfX * u.y * u.y); // First iteration
    u.y *= 1.5F - (halfX * u.y * u.y); // Second iteration

    return u.y;
}

// [Pizer’s optimisation](https://pizer.wordpress.com/2008/10/12/fast-inverse-square-root/)
float reciprocalSqrtPizer1(float x) {
    union {
        float y;
        long i;
    } u = {x};

    u.i = 0x5f1f1412 - (u.i >> 1); // Initial estimate for Newton–Raphson method
    u.y *= 1.69000231F - 0.714158168F * x * u.y * u.y; // Single iteration

    return u.y;
}

float reciprocalSqrtPizer2(float x) {
    union {
        float y;
        long i;
    } u = {x};

    u.i = 0x5f1f1412 - (u.i >> 1); // Initial estimate for Newton–Raphson method
    u.y *= 1.69000231F - 0.714158168F * x * u.y * u.y; // First iteration
    u.y *= 1.5F - (0.5F * x * u.y * u.y); // Second iteration

    return u.y;
}

void test_square_root_reciprocal()
{

    TEST_ASSERT_FLOAT_WITHIN(0.0009F,    1.41421356F, reciprocalSqrt1(0.5F));
    TEST_ASSERT_FLOAT_WITHIN(0.0007F,    1.41421356F, reciprocalSqrtPizer1(0.5F));
    TEST_ASSERT_FLOAT_WITHIN(0.0000003F, 1.41421356F, reciprocalSqrt2(0.5F));
    TEST_ASSERT_EQUAL_FLOAT(1.41421356F, reciprocalSqrtPizer2(0.5F));

    TEST_ASSERT_FLOAT_WITHIN(0.0009F,    0.5F, reciprocalSqrt1(4));
    TEST_ASSERT_FLOAT_WITHIN(0.00006F,   0.5F, reciprocalSqrtPizer1(4));
    TEST_ASSERT_FLOAT_WITHIN(0.000003F,  0.5F, reciprocalSqrt2(4));
    TEST_ASSERT_EQUAL_FLOAT(0.5F, reciprocalSqrtPizer2(4));

    TEST_ASSERT_FLOAT_WITHIN(0.0003F,    0.125F, reciprocalSqrt1(64));
    TEST_ASSERT_FLOAT_WITHIN(0.000015F,  0.125F, reciprocalSqrtPizer1(64));
    TEST_ASSERT_EQUAL_FLOAT(0.125F, reciprocalSqrtPizer2(64));

    TEST_ASSERT_FLOAT_WITHIN(0.00016F,   0.1F, reciprocalSqrt1(100));
    TEST_ASSERT_FLOAT_WITHIN(0.00006F,   0.1F, reciprocalSqrtPizer1(100));
    TEST_ASSERT_FLOAT_WITHIN(0.00000037F,0.1F, reciprocalSqrt2(100));
    TEST_ASSERT_EQUAL_FLOAT(0.1F, reciprocalSqrtPizer2(100));

    TEST_ASSERT_FLOAT_WITHIN(0.00002F,   0.01F, reciprocalSqrt2(10000));
    TEST_ASSERT_EQUAL_FLOAT(0.01F, reciprocalSqrtPizer2(10000));
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
    UNITY_BEGIN();

    RUN_TEST(test_square_root_reciprocal);

    UNITY_END();
}
