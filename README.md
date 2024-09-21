# fhi_manipulator_control_stack

Author: [Marco Magri](marco.magri@fraunhofer.it) Fraunhofer Italia 2023-2024

This repository contains a framework for developing and implementing controllers for modular and reconfigurable robots. It allows automatic controller generation and tuning for different mechanical configurations of robots. The framework is designed to provide high stability and performance in various robot assemblies.

### Overview
This project tackles the challenges of controlling reconfigurable robots by offering a modular software architecture that adapts to different hardware configurations. The framework decouples the low-level controller logic from hardware-specific implementations, such as ROS 2 Control or XBot2, providing flexibility and scalability for a wide range of robotic platforms.
Key Components:
- Robot Dynamics Library: Handles kinematic and dynamic calculations, including forward/inverse kinematics, Jacobians, and dynamic equations.
- Manipulator Control Library: Standardizes low-level controller logic across different hardware abstraction frameworks.
- ControllerCore Class: The central interface for all controller implementations.

### Features
- Modular Control Framework: A software architecture that supports seamless integration with different hardware abstraction layers.
- Automatic Controller Generation: Controllers are automatically generated and tuned based on robot modules.
- Library for Robot Dynamics: A custom library for performing efficient kinematic and dynamic calculations.
- Manipulator Control Library: Provides a standardized interface for low-level controller implementations.
- Framework Support: The project supports multiple frameworks, including ROS 2 Control and XBot2, ensuring consistent control logic implementation across platforms.

### Installation
Create catkin workspace and clone package:
```bash
mkdir -vp ros1_ws/src
cd ros1_ws/src
git clone https://github.com/FraunhoferItalia/fhi_manipulator_control_stack.git --recurse-submodules
```
Build the ROS1 workspace:
```bash
catkin_make -DXBOT2_ENABLE_XENO=ON
source devel/setup.bash
```

### Usage
This repository provides configurable controllers that can be deployed on different robots. Here’s how you can use the framework:
- Configure your robot's modules and their kinematic/dynamic properties using the Robot Dynamics Library.
- Deploy the low-level controller using the Manipulator Control Library.
- Use the framework-specific controller implementation to interface with your robot's hardware abstraction layer (e.g., ROS 2 Control or XBot2).

### Licence
fhi_manipulator_control_stack is licensed under the terms of the Apache License 2.0. The project has received financial support by the Horizon 2020 EU Project [CONCERT](https://concertproject.eu/).
