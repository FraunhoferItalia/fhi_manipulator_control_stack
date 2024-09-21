/*
 * Copyright 2022-2024 Fraunhofer Italia Research
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef CONTROLLER_CORE_H
#define CONTROLLER_CORE_H

#pragma once
#include <atomic>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <fhi_manipulator_control/parameters.hpp>
#include <fhi_robot_dynamic_library/chain.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

namespace fmc
{

bool vector_to_eigen(const std::vector<double> & vector, Eigen::VectorXd & eigen);
bool eigen_to_vector(const Eigen::VectorXd & eigen, std::vector<double> & vector);

enum ReturnType { OK, ERROR };

enum InterfaceType {
  NONE = 0,
  POSITION = 1 << 0,
  VELOCITY = 1 << 1,
  ACCELERATION = 1 << 2,
  EFFORT = 1 << 3,
  CURRENT = 1 << 4,
  RESIDUAL_TORQUES = 1 << 5,
  TCP_WRENCH = 1 << 6,
  STIFFNESS = 1 << 7,
  DAMPING = 1 << 8,
};

InterfaceType operator|(InterfaceType lhs, InterfaceType rhs);
InterfaceType operator&(InterfaceType lhs, InterfaceType rhs);

using Interfaces = std::map<InterfaceType, std::vector<double>>;

using ControllerParameterDescriptions = std::map<unsigned int, parameters::Description>;
using ControllerParameters = std::map<unsigned int, parameters::Parameter::SharedPtr>;

class ControllerCore
{
public:
  using SharedPtr = std::shared_ptr<ControllerCore>;
  using UniquePtr = std::unique_ptr<ControllerCore>;

  ControllerCore(
    const InterfaceType & required_state_interfaces,
    const InterfaceType & required_reference_interfaces,
    const InterfaceType & required_command_interfaces);
  ~ControllerCore() = default;

  ReturnType init_chains(const std::string & chain_config_file_path);

  virtual ReturnType init_params(const std::string & parameters_file_path);
  virtual fmc::ControllerParameterDescriptions export_parameters();
  virtual std::map<std::string, double *> export_register_variables();

  ReturnType set_reference(
    const InterfaceType & interface_type, const std::vector<double> & reference);
  ReturnType set_reference(const InterfaceType & interface_type, const Eigen::VectorXd & reference);

  ReturnType set_state(const InterfaceType & interface_type, const std::vector<double> & state);
  ReturnType set_state(const InterfaceType & interface_type, const Eigen::VectorXd & state);

  ReturnType set_output(const InterfaceType & interface_type, const std::vector<double> & output);
  ReturnType set_output(const InterfaceType & interface_type, const Eigen::VectorXd & output);

  bool interfaces_valid(const InterfaceType & required, const InterfaceType & provided);

  ReturnType update(const std::chrono::nanoseconds & dt);

  ReturnType get_output(const InterfaceType & interface_type, std::vector<double> & out_vector);
  ReturnType get_output(const InterfaceType & interface_type, Eigen::VectorXd & out);

  template <
    typename T, typename = std::enable_if_t<std::is_same_v<T, bool> || std::is_same_v<T, double>>>
  T get_parameter(const unsigned int & param_enum)
  {
    return parameters::get_parameter<T>(parameters_.at(param_enum));
  }

  const InterfaceType type;
  std::string error_msg;
  ControllerParameters parameters_;
  const InterfaceType required_state_interfaces_;
  const InterfaceType required_reference_interfaces_;
  const InterfaceType required_command_interfaces_;

protected:
  frdl::ChainModules::SharedPtr modules_;
  frdl::Chain::UniquePtr current_chain_;
  frdl::Chain::UniquePtr desired_chain_;
  frdl::JointState::SharedPtr current_chain_state_;
  frdl::JointState::SharedPtr desired_chain_state_;
  virtual ReturnType update_impl(const std::chrono::nanoseconds & dt) = 0;
  virtual ReturnType on_activate() = 0;
  std::vector<double> vector_temp_;

private:
  std::string chain_config_file;
  Interfaces output_;

  InterfaceType provided_state_interfaces_ = InterfaceType::NONE;
  InterfaceType provided_reference_interfaces_ = InterfaceType::NONE;
  InterfaceType provided_command_interfaces_ = InterfaceType::NONE;

  ReturnType controller_impl_result_;
};

}  // namespace fmc

#endif  // CONTROLLER_CORE_H