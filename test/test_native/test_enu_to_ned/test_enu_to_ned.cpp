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
    const Quaternion qNEDtoENU(0.0F, -sqrtf(2.0F)/2.0F, -sqrtf(2.0F)/2.0F, 0.0F);

    const xyz_t vNED = qENUtoNED.rotate(vENU);
    TEST_ASSERT_EQUAL_FLOAT(vNED_E.x, vNED.x);
    TEST_ASSERT_EQUAL_FLOAT(vNED_E.y, vNED.y);
    TEST_ASSERT_EQUAL_FLOAT(vNED_E.z, vNED.z);

    const float rollENU = 19.0F;
    const float pitchENU = 43.0F;
    const float yawENU = 67.0F;

    const Quaternion qENU = Quaternion::fromEulerAnglesDegrees(rollENU , pitchENU, yawENU);
    TEST_ASSERT_EQUAL_FLOAT(rollENU, qENU.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(pitchENU, qENU.calculatePitchDegrees());
    TEST_ASSERT_EQUAL_FLOAT(yawENU, qENU.calculateYawDegrees());

    const Quaternion qNED = qENUtoNED * qENU;
    TEST_ASSERT_EQUAL_FLOAT(rollENU - 180.0F, qNED.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(180.0F + pitchENU, qNED.calculatePitchDegrees());
    TEST_ASSERT_EQUAL_FLOAT(90.0F - yawENU, qNED.calculateYawDegrees());

    const Quaternion qENU2 = qNEDtoENU * qNED;
    TEST_ASSERT_EQUAL_FLOAT(rollENU, qENU2.calculateRollDegrees());
    TEST_ASSERT_EQUAL_FLOAT(pitchENU, qENU2.calculatePitchDegrees());
    TEST_ASSERT_EQUAL_FLOAT(yawENU, qENU2.calculateYawDegrees());
}
// NOLINTEND(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays,modernize-use-using,readability-non-const-parameter,cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)


int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    UNITY_BEGIN();

    RUN_TEST(test_enu_to_ned);

    UNITY_END();
}
