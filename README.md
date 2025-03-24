# Sensor Fusion

This library contains [sensor fusion](https://en.wikipedia.org/wiki/Sensor_fusion) algorithms to combine
output from a gyroscope, accelerometer, and optionally a magnetometer to give output that has less uncertainty
than the output of the individual sensors.

Three sensor fusion algorithms are provided:

1. Complementary Filter
2. Mahony Filter
3. Madgwick Filter

The Madgwick filter has been refactored to be more computationally efficient (and so faster) than
the standard version used in many implementations (Arduino, Adafruit, M5Stack, Reefwing-AHRS),
[see MadgwickRefactoring](documents/MadgwickRefactoring.md) for details.
