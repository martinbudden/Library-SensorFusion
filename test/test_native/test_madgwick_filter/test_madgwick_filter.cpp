#include "arduino_madgwick.h"
#include "SensorFusion.h"
#include <cmath>
#include <unity.h>


void setUp() {
}

void tearDown() {
}

void test_arduino_madgwick() {
    static Madgwick arduinoMadgwick;
    static MadgwickFilter madgwickFilter;

    const float gx = 0.33F;
    const float gy = 0.55F;
    const float gz = -0.19F;
    xyz_t acc = { .x =0.66F, .y = 0.77F, .z = 0.88F };
    acc *= 1.0F / sqrtf(acc.x*acc.x + acc.y*acc.y +acc.z*acc.z);

    const float q0 = 3.0F;
    const float q1 = 5.0F;
    const float q2 = 7.0F;
    const float q3 = 9.0F;

    const xyz_t gyro = { .x = gx*0.0174533F, .y = gy*0.0174533F, .z = gz*0.0174533F };
    float beta { 1.0 };
    float deltaT = 0.1F;

    arduinoMadgwick.begin(1.0F / deltaT, beta);
    madgwickFilter.setBeta(beta);

    arduinoMadgwick._setAndNormalizeQ(q0, q1, q2, q3);
    arduinoMadgwick.updateIMU(gx, gy, gz, acc.x, acc.y, acc.z);
    TEST_ASSERT_EQUAL(0, arduinoMadgwick.anglesComputed);
    float qE0 = arduinoMadgwick.q0;
    float qE1 = arduinoMadgwick.q1;
    float qE2 = arduinoMadgwick.q2;
    float qE3 = arduinoMadgwick.q3;

    madgwickFilter._setAndNormalizeQ(q0, q1, q2, q3);
    Quaternion q = madgwickFilter.updateOrientation(gyro, acc, deltaT);
    TEST_ASSERT_EQUAL_FLOAT(qE0, q.getW());
    TEST_ASSERT_EQUAL_FLOAT(qE1, q.getX());
    TEST_ASSERT_EQUAL_FLOAT(qE2, q.getY());
    TEST_ASSERT_EQUAL_FLOAT(qE3, q.getZ());


    beta = 0.1F;
    deltaT = 0.01F;
    arduinoMadgwick.begin(1.0F / deltaT, beta);
    madgwickFilter.setBeta(beta);
    acc = { .x = -0.66F, .y = 1.77F, .z = 3.88F };
    acc *= 1.0F / sqrtf(acc.x*acc.x + acc.y*acc.y +acc.z*acc.z);

    arduinoMadgwick._setAndNormalizeQ(q0, q1, q2, q3);
    arduinoMadgwick.updateIMU(gx, gy, gz, acc.x, acc.y, acc.z);
    qE0 = arduinoMadgwick.q0;
    qE1 = arduinoMadgwick.q1;
    qE2 = arduinoMadgwick.q2;
    qE3 = arduinoMadgwick.q3;

    madgwickFilter._setAndNormalizeQ(q0, q1, q2, q3);
    q = madgwickFilter.updateOrientation(gyro, acc, deltaT);
    TEST_ASSERT_EQUAL_FLOAT(qE0, q.getW());
    TEST_ASSERT_EQUAL_FLOAT(qE1, q.getX());
    TEST_ASSERT_EQUAL_FLOAT(qE2, q.getY());
    TEST_ASSERT_EQUAL_FLOAT(qE3, q.getZ());
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
    UNITY_BEGIN();

    RUN_TEST(test_arduino_madgwick);

    UNITY_END();
}
