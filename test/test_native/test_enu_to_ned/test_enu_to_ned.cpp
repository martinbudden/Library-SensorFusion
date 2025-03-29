#include "Quaternion.h"
#include <unity.h>


void setUp() {
}

void tearDown() {
}

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays,modernize-use-using,readability-non-const-parameter,cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)

void test_enu_to_ned()
{
    // _E prefix denotes expected value

    const xyz_t vENU   { .x = 3.0F, .y = 5.0F, .z =  7.0F };
    const xyz_t vNED_E { .x = 5.0F, .y = 3.0F, .z = -7.0F };

    const Quaternion qENUtoNED(0.0F, sqrtf(2.0F)/2.0F, sqrtf(2.0F)/2.0F, 0.0F);

    const xyz_t vNED = qENUtoNED.rotate(vENU);
    TEST_ASSERT_EQUAL_FLOAT(vNED_E.x, vNED.x);
    TEST_ASSERT_EQUAL_FLOAT(vNED_E.y, vNED.y);
    TEST_ASSERT_EQUAL_FLOAT(vNED_E.z, vNED.z);

    const float rollNED = 19.0F;
    const float pitchNED = 43.0F;
    const float yawNED = 67.0F;
    const Quaternion qENU = Quaternion::fromEulerAnglesDegrees(rollNED , pitchNED, yawNED);
    TEST_ASSERT_EQUAL_FLOAT(rollNED, qENU.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(pitchNED, qENU.calculatePitchDegrees());
    TEST_ASSERT_EQUAL_FLOAT(yawNED, qENU.calculateYawDegrees());

    const Quaternion qNED = qENUtoNED * qENU;
    TEST_ASSERT_EQUAL_FLOAT(rollNED - 180.0F, qNED.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(-pitchNED, qNED.calculatePitchDegrees());
    TEST_ASSERT_EQUAL_FLOAT(90.0F - yawNED, qNED.calculateYawDegrees());

    const Quaternion qNED2 = qENU * qENUtoNED;
    TEST_ASSERT_EQUAL_FLOAT(-135.3967, qNED2.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(-13.77475, qNED2.calculatePitchDegrees());
    TEST_ASSERT_EQUAL_FLOAT(143.7846, qNED2.calculateYawDegrees());

}
// NOLINTEND(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays,modernize-use-using,readability-non-const-parameter,cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)


int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
    UNITY_BEGIN();

    RUN_TEST(test_enu_to_ned);

    UNITY_END();
}
