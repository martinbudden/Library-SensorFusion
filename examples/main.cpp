#include <SensorFusion.h>
#include <cstdio>


int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
    static MadgwickFilter sensorFusionFilter;
    sensorFusionFilter.setBeta(0.1F);
    const float deltaT = 0.01F; // deltaT should be set to the interval between IMU readings

    // this loop should repeat each time new IMU data is available
    while (true) {
        const xyz_t gyroRPS = {0.0, 0.0, 0.0}; // replace this with actual gyroscope data in radians/s
        const xyz_t accelerometer = {0.0, 0.0, 1.0}; // replace this with actual accelerometer data in g

        const Quaternion orientation = sensorFusionFilter.update(gyroRPS, accelerometer, deltaT);

        const float rollDegrees = orientation.calculateRollDegrees();
        const float pitchDegrees = orientation.calculatePitchDegrees();
        const float yawDegrees = orientation.calculateYawDegrees();

        printf("Roll %0.1f, Pitch %0.1f, Yaw %0.1f\n", rollDegrees, pitchDegrees, yawDegrees);
    }
}
