# Terminology and Notation used by VQF Sensor fusion

See [VQF: Highly Accurate IMU Orientation Estimation with Bias Estimation and Magnetic Disturbance Rejection](https://arxiv.org/pdf/2203.17024)

Daniel Laidig's implementation [here](https://github.com/dlaidig/vqf)

| Term | meaning |
| ---- | ------- |
| ENU | East North Up - x east, y north, z down |
| IOE | Inertial Orientation Estimation |
| 6D IOE | IOE using only gyroscopes and accelerometers |
| 9D IOE | IEO additionally using magnetometers |

| Notation | meaning |
| -------- | ------- |
| **ω** | gyroscope reading (omega) |
| **a** | accelerometer reading |
| **m** | magnetometer reading |
| [**a**]<sub>ε</sub> | accelerometer reading transformed to ε frame |
| (t<sub>k</sub>) | reading at time t<sub>k</sub> eg **ω**(t<sub>k</sub>) |
| ς<sub>i</sub> | moving sensor frame. **Mnemonic:** ς **s**igma **s**ensor|
| ε<sub>i</sub> | sensor-specific almost-inertial reference frame used in 6D IOE. **Mnemonic:** ε **e**psilon **e**arth |
| ε | sensor-specific almost-inertial reference frame used in 9D IOE |
| Ι<sub>i</sub> | auxiliary frame representing orientation obtained by pure gyroscopic integration (slowly drifting frame)<br>Ι<sub>i</sub>(t<sub>0</sub>) = ς<sub>i</sub>(t<sub>0</sub>). **Mnemonic:** Ι **i**ota **i**ntegration|
| **q** | quaternion |
| ⊗ | quaternion multiplication |
| <sup>a</sup><sub>b</sub>**q** | quaternion that rotates from frame a to frame b|
| <sup>a</sup><sub>b</sub>**q** ⊗ **v** ⊗ <sup>a</sup><sub>b</sub>**q**<sup>-1</sup> | rotate vector **v** from frame a to frame b using quaternion <sup>a</sup><sub>b</sub>**q** |
| ( α @ **v** ) | quaternion that represents a rotation of angle α around axis **v**.<br>Equates to the quaternion [cos(α/2) **v**<sup>T</sup>sin(α/2)]<sup>T</sup> |

| Value    | variable | description |
| -------- | -------- | ----------- |
| T<sub>s</sub> | deltaT | time step |
| <sup>ςi</sup><sub>Ιi</sub>**q** | _state.gyroQuaternion |gyroscopic strapdown integration quaternion |
| <sup>Ii</sup><sub>εi</sub>**q** | _state.accQuaternion |inclination correction quaternion |
| <sup>ςi</sup><sub>εi</sub>**q** | _state.orientation6D |6D orientation estimate  = <sup>Ii</sup><sub>εi</sub>**q** ⊗ <sup>ςi</sup><sub>Ιi</sub>**q** |
| <sup>εi</sup><sub>ε</sub>**q** | n/a | heading correction rotation. Not used: δ<sub>i</sub> used directly instead |
| δ<sub>i</sub> | _state.delta |scalar heading offset, which represents the heading correction rotation <sup>εi</sup><sub>ε</sub>**q** |
| [**a**]<sub>Ιi</sub> | accAuxiliary | acceleration transformed into Ιi (auxiliary) frame |
| [**a**<sub>LP</sub>]<sub>Ιi</sub> | accAuxiliary (re-used) | filtered acceleration in Ιi frame |
| [**a**<sub>LP</sub>]<sub>εi</sub> | accEarth | filtered acceleration transformed into εi (earth) frame|
| **z** | n/a | unit vector in z-direction = [ 0 0 1 ]<sup>T</sup> |
| <sup>ςi</sup><sub>ε</sub>**q** | _state.orientation9D| 9D orientation estimate = <sup>εi</sup><sub>ε</sub>**q** ⊗ <sup>ςi</sup><sub>εi</sub>**q** = (δ<sub>i</sub> @ **z**) ⊗ <sup>ςi</sup><sub>εi</sub>**q** |

## Basic VQF Algorithm 1: **InitializeFilter** and **FilterUpdate**

### In BasicVQF::resetState()

Steps 1 to 6

<sup>ςi</sup><sub>Ιi</sub>**q** = [1 0 0 0]<sup>T</sup><br>
<sup>Ii</sup><sub>εi</sub>**q** = [1 0 0 0]<sup>T</sup><br>
δ<sub>i</sub> = 0<br>
initialize low-pass filter state

```cpp
_state.gyroQuaternion.setToIdentity();
_state.accQuaternion.setToIdentity();

_state.accLPF.reset();

_state.delta = 0.0;
_state.kMagInit = 1.0;
```

### In BasicVQF::updateGyro()

8\. Perform gyroscope strapdown integration

<sup>ςi</sup><sub>Ιi</sub>**q** = <sup>ςi</sup><sub>Ιi</sub>**q** ⊗ ( T<sub>s</sub> ||**ω**|| @ **v** )

```cpp
const float angle2 = gyroMagnitude * deltaT * 0.5F;
const float c = cosf(angle2);
const float s = sinf(angle2)/gyroMagnitude;
const Quaternion integrationStep(c, s*gyroRPS.x, s*gyroRPS.y, s*gyroRPS.z);
_state.gyroQuaternion *= integrationStep
```

### In BasicVQF::updateAccelerometer()

9\. Transform acceleration to Ι<sub>i</sub> frame

[**a**]<sub>Ιi</sub> = <sup>ςi</sup><sub>Ιi</sub>**q** ⊗ **a** ⊗ <sup>ςi</sup><sub>Ιi</sub>**q**<sup>-1</sup>

```cpp
xyz_t accAuxiliary = _state.gyroQuaternion.rotate(acc);
```

10\. Apply low-pas filter

[**a**<sub>LP</sub>]<sub>Ιi</sub> = lpfStep([**a**]<sub>Ιi</sub>, f<sub>c</sub> = f<sub>c,acc</sub>)

```cpp
accAuxiliary = _state.accLPF.filter(accAuxiliary);
```

11\. Transform to ε<sub>i</sub> frame

[**a**<sub>LP</sub>]<sub>εi</sub> = <sup>Ιi</sup><sub>εi</sub>**q** ⊗ [**a**<sub>LP</sub>]<sub>Ιi</sub> ⊗ <sup>Ιi</sup><sub>εi</sub>**q**<sup>-1</sup>

```cpp
xyz_t accEarth = _state.accQuaternion.rotate(accAuxiliary);
```

12\. Normalize

[ a<sub>x</sub> a<sub>y</sub> a<sub>z</sub> ]<sup>T</sup> = [**a**<sub>LP</sub>]<sub>εi</sub> / ||[**a**<sub>LP</sub>]<sub>εi</sub>||

```cpp
normalize(accEarth);
```

13\.

q<sub>w</sub> = sqrt((a<sub>z</sub> + 1)/2)

```cpp
const float qW = sqrtf((accEarth.z + 1.0F) * 0.5F);
```

14\. Update correction quaternion

<sup>Ii</sup><sub>εi</sub>**q** = [q<sub>w</sub> a<sub>y</sub>/(2q<sub>w</sub>) -a<sub>x</sub>/(2q<sub>w</sub>) 0]<sup>T</sup> ⊗ <sup>Ii</sup><sub>εi</sub>**q**

```cpp
const Quaternion correctionStep = (qW > 1e-6F) ? Quaternion(qW, accEarth.y/(2.0F*qW), -accEarth.x/(2.0F*qW), 0) : Quaternion(0.0F, 1.0F, 0.0F, 0.0F);
state.accQuaternion = correctionStep * _state.accQuaternion;
normalize(_state.accQuaternion);
```

15\. Calculate 6D orientation estimate

<sup>ςi</sup><sub>εi</sub>**q** = <sup>Ii</sup><sub>εi</sub>**q** ⊗ <sup>ςi</sup><sub>Ii</sub>**q**

```cpp
_state.orientation6D = _state.accQuaternion * _state.gyroQuaternion;
```

### In BasicVQF::updateMagnetometer()

16\.

**if** **m** is given **then**

```cpp
if (magnetometer.x == 0.0F && magnetometer.y == 0.0F && magnetometer.z == 0.0F) {
    return _state.orientation6D;
}
```

17\. Transform magnetometer sample to ε<sub>i</sub> frame

[ m<sub>x</sub> m<sub>y</sub> m<sub>z</sub> ]<sup>T</sup> = <sup>ςi</sup><sub>εi</sub>**q** ⊗ **m** ⊗ <sup>ςi</sup><sub>εi</sub>**q**<sup>-1</sup>

```cpp
const xyz_t magEarth = _state.orientation6D.rotate(magnetometer);
```

18\.Calculate heading offset from magnetometer sample

δ<sub>mag</sub> = atan2(m<sub>x</sub>, m<sub>y</sub>)

```cpp
_state.magDisagreementAngle = atan2f(magEarth.x, magEarth.y);
```

19\. Update correction angle

δ<sub>i</sub> = δ<sub>i</sub> + k<sub>mag</sub>wrapToPi(δ<sub>mag</sub> - δ<sub>i</sub>)

```cpp
_state.magDisagreementAngle = wrapToPi(_state.magDisagreementAngle - _state.delta);
_state.delta += k*_state.magDisagreementAngle;
```

20\.
**end if**


21\. Calculate 9D orientation estimate

<sup>ςi</sup><sub>ε</sub>**q** = [cos(δ<sub>i</sub>/2) 0 0 sin(δ<sub>i</sub>/2)]<sup>T</sup> ⊗ <sup>ςi</sup><sub>εi</sub>**q**

```cpp
_state.orientation9D = _state.orientation6D.rotateZ(_state.delta);
```
