# Multidimensional Robust Controller

Cpp source code for the robust controller (based on URDF convention) for modular reconfigurable robots. It includes a small cpp script to parse some controller parameters in order to quickly test the code.

The developed code is based on controller presented in *Nainer, C. and Giusti, A., 2022. Automatically deployable robust control of modular reconfigurable robot manipulators. IEEE Robotics and Automation Letters, 7(2), pp.5286-5293.*

## Build and Run
Download required dependencies:
```bash
sudo apt-get install nlohmann-json3-dev
```
In order to build the code
```
cmake .
make
```
Then the executable can be simply run as
```
./robust_controller_test
```

## Notes
This controller now works for any degrees of freedom since it uses ```coder::array``` for the parameters and states, instead of the fixed-size C-array.

In this version, gravity can be set in any direction: the parameter ```g``` is an array.