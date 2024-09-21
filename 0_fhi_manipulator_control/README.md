# Fhi Manipulator Control (FMC)

Author: [Marco Magri](marco.magri@fraunhofer.it) Fraunhofer Italia 2023-2024

The Manipulator Control Library offers a standardized interface for low-level controller implementation through the `ControllerCore` class. This library facilitates the seamless integration of control logic across different hardware abstraction frameworks, such as XBot, ROS 2 control, and others. By abstracting these frameworks, it ensures consistent control logic implementation regardless of the underlying hardware.

## Features
- **Standardized Interface:** Provides a consistent interface for low-level controller implementation.
- **Framework Abstraction:** Compatible with multiple hardware abstraction frameworks including XBot and ROS 2 control.
- **Reusable Controllers:** Allows reuse of the same low-level controller across different frameworks without redeveloping the control logic.
- **Customizable Methods:** Includes methods that can be overridden for specific controller implementations.

## ControllerCore Class

The `ControllerCore` class is at the heart of the Manipulator Control Library. It serves as an interface for the implementation of low-level controllers and acts as the foundation for all controller implementations. This class allows the reuse of the same low-level controller implementation across different frameworks, such as XBot and ROS 2 control, without needing to redevelop the low-level control logic.

### Key Methods

The `ControllerCore` class provides several methods intended to be called by the framework-specific implementation:

- **`set_reference`**: Sets the value of an interface (e.g., position, velocity, torque) in the desired chain attribute.
- **`set_state`**: Sets the value of an interface (e.g., position, velocity, torque) in the current chain attribute.
- **`update`**: Executes the update loop of the low-level controller. This method ensures that all interfaces (both reference and state) required by the controller have been correctly provided. It updates the current chain and desired chain attributes based on the values provided by the `set_reference` and `set_state` calls.
- **`get_output`**: Retrieves the result of the update computations, which serves as the control input for the robot.

### Virtual Methods for Custom Implementations

The `ControllerCore` class also defines a set of virtual methods that must be overridden for specific controller implementations:

- **`update_impl`**: To be overridden with controller-specific logic. The controller developer can access the state of the robot and the references via the current chain and desired chain attributes.
- **`on_activate`**: Allows custom operations upon controller activation.
- **`on_deactivate`**: Allows custom operations upon controller deactivation.
- **`set_output`**: Used to set the output after computations.
- **`export_parameters`** [Optional]: Allows the specification of controller parameters that can be changed at runtime.
- **`export_register_variables`** [Optional]: For debugging purposes, this allows specifying a set of controller variables to be logged while a controller is running.

## Installation
This library can be installed a standard ROS1 package using catkin.

## Usage

To use the Manipulator Control Library, follow these steps:

1. Implement the `ControllerCore` class and override the necessary virtual methods.
2. Integrate your implementation with the desired hardware abstraction framework (e.g., XBot, ROS 2 control).
3. Use the provided methods to set references, update states, and retrieve outputs during the control loop.

### Licence
Fhi Manipulator Control is licensed under the terms of the Apache License 2.0. The project has received financial support by the Horizon 2020 EU Project [CONCERT](https://concertproject.eu/).