# FhI Robot Dynamic Library (frdl)

Author: [Marco Magri](marco.magri@fraunhofer.it) Fraunhofer Italia 2023-2024

This repository contains a custom robot dynamics library designed to perform kinematic and dynamic calculations for a robot model. Unlike other libraries such as [KDL](https://www.orocos.org/wiki/Kinematic_and_Dynamic_Solvers.html) and [Pinocchio](https://stack-of-tasks.github.io/pinocchio/), this library accepts joint module descriptions as input for the robot model, providing a more tailored solution for our low-level controllers.

## Overview

The library is a core component of the control framework, enabling efficient and accurate calculations of robot dynamics and kinematics. It includes two primary classes:

- **ChainModules**: Manages the dynamic and kinematic information for each robot module.
- **Chain**: Wraps kinematic and dynamic methods for a specific module composition. This class ensures computational efficiency by reusing results when no state updates are detected between function calls. It also provides user-friendly methods to set and get both Cartesian and joint states of the kinematic chain.

## Features

### Kinematic Functions

The library implements several essential kinematic functions:

- **`get_forward_kinematic(chain_modules, q, link_index)`**: Computes the Cartesian pose of the specified link's reference frame given the current joint positions `q` relative to the kinematic chain's base link.
  
- **`get_jacobian(chain_modules, q)`**: Computes the Jacobian matrix of the kinematic chain given the joint positions `q`.

- **`get_jacobian_dot(chain_modules, q, dq)`**: Calculates the time derivative of the Jacobian matrix given the joint positions `q` and joint velocities `dq`.

- **`get_inverse_kinematics(chain_modules, X)`**: Computes the joint state values `q` required to place the last link's reference frame in the Cartesian pose defined by the homogeneous transformation `X`.

### Dynamics Functions

The library also provides key robot dynamics functions:

- **`get_rne(chain_modules, q, dq, ddq)`**: Uses the Recursive Newton-Euler (RNE) algorithm to compute the joint torque `τ` based on the joint positions `q`, velocities `dq`, and accelerations `ddq`.

- **`get_M_matrix(chain_modules, q)`**: Computes the inertia matrix `M` for the module assembly at the specified joint positions `q`.

- **`get_C_matrix(chain_modules, q, dq)`**: Calculates the Coriolis and centrifugal contribution matrix `C` for the module assembly given the joint positions `q` and velocities `dq`.

- **`get_gravity(chain_modules, q)`**: Computes the joint torque needed to compensate for the gravity effect on the kinematic chain at the given joint positions `q`.


## Installation
This library can be installed a standard ROS1 package using catkin.

## Usage
1. Create a `chain_config.yaml` file starting from your urdf by running:
    ``` bash
    rosrun fhi_robot_dynamic_library generate_chain_configuration --urdf_file_path <your_urdf_file_path> 
    ```
2. Provide the file content to the `ChainModules` class:
   ```cpp
    std::ifstream fin("chain_config.yaml");
    std::stringstream buffer;
    buffer << fin.rdbuf();
    std::string yamlContent = buffer.str();
    frdl::ChainModules chain_modules = frdl::ChainModules(yamlContent);
   ```
3. Use the `ChainModules` instance to instate a `Chain` object or use it directly with the functions:
    ```cpp
    auto chain = frdl::Chain(chain_modules);
    // or 
    Eigen::VectorXd q;
    Eigen::Affine3d x;
    frdl::Kinematics::get_forward_kinematic(modules, q, x);
    ```

### Licence
FhI Robot Dynamic Library is licensed under the terms of the Apache License 2.0. The project has received financial support by the Horizon 2020 EU Project [CONCERT](https://concertproject.eu/).
