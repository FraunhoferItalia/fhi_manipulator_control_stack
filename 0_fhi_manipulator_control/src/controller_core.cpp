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
#include <fhi_manipulator_control/controller_core.hpp>

namespace fmc
{

bool vector_to_eigen(const std::vector<double> & vector, Eigen::VectorXd & eigen)
{
  if (static_cast<long int>(vector.size()) != eigen.size()) return false;
  for (size_t i = 0; i < vector.size(); ++i) {
    eigen[i] = vector[i];
  }
  return true;
}

bool eigen_to_vector(const Eigen::VectorXd & eigen, std::vector<double> & vector)
{
  if (static_cast<long int>(vector.size()) != eigen.size()) return false;
  for (size_t i = 0; i < vector.size(); ++i) {
    vector[i] = eigen[i];
  }
  return true;
}

InterfaceType operator|(InterfaceType lhs, InterfaceType rhs)
{
  return static_cast<InterfaceType>(static_cast<int>(lhs) | static_cast<int>(rhs));
}
InterfaceType operator&(InterfaceType lhs, InterfaceType rhs)
{
  return static_cast<InterfaceType>(static_cast<int>(lhs) & static_cast<int>(rhs));
}

ControllerCore::ControllerCore(
  const InterfaceType & required_state_interfaces,
  const InterfaceType & required_reference_interfaces,
  const InterfaceType & required_command_interfaces)
: type{required_command_interfaces},
  required_state_interfaces_{required_state_interfaces},
  required_reference_interfaces_{required_reference_interfaces},
  required_command_interfaces_{required_command_interfaces}
{
}

ReturnType ControllerCore::init_chains(const std::string & chain_config_file_path)
{
  if (chain_config_file_path.empty()) {
    return ReturnType::ERROR;
  }
  std::ifstream fin(chain_config_file_path);
  if (!fin) {
    std::cerr << "Failed to open file chain file: " << chain_config_file_path << std::endl;
    return ReturnType::ERROR;
  }

  std::stringstream buffer;
  buffer << fin.rdbuf();
  std::string yamlContent = buffer.str();
  modules_ = std::make_shared<frdl::ChainModules>(yamlContent);
  current_chain_ = std::make_unique<frdl::Chain>(modules_);
  desired_chain_ = std::make_unique<frdl::Chain>(modules_);
  current_chain_state_ = current_chain_->get_state();
  desired_chain_state_ = desired_chain_->get_state();
  vector_temp_.resize(current_chain_->get_size());
  return ReturnType::OK;
}

ReturnType ControllerCore::init_params(const std::string & /* parameters_file_path */)
{
  for (auto & [param_key, param_description] : export_parameters()) {
    parameters_.insert({param_key, std::make_shared<parameters::Parameter>(param_description)});
  }
  return ReturnType::OK;
}

fmc::ControllerParameterDescriptions ControllerCore::export_parameters() { return {}; }
std::map<std::string, double *> ControllerCore::export_register_variables() { return {}; }

ReturnType ControllerCore::set_reference(
  const InterfaceType & interface_type, const std::vector<double> & reference)
{
  provided_reference_interfaces_ = provided_reference_interfaces_ | interface_type;
  switch (interface_type) {
    case fmc::InterfaceType::POSITION:
      vector_to_eigen(reference, desired_chain_state_->q);
      break;
    case fmc::InterfaceType::VELOCITY:
      vector_to_eigen(reference, desired_chain_state_->dq);
      break;
    case fmc::InterfaceType::ACCELERATION:
      vector_to_eigen(reference, desired_chain_state_->ddq);
      break;
    case fmc::InterfaceType::EFFORT:
      vector_to_eigen(reference, desired_chain_state_->tau);
      current_chain_state_->is_tau_residual = false;
      break;
    case fmc::InterfaceType::RESIDUAL_TORQUES:
      vector_to_eigen(reference, desired_chain_state_->tau);
      current_chain_state_->is_tau_residual = true;
      break;
    case fmc::InterfaceType::TCP_WRENCH:
      vector_to_eigen(reference, desired_chain_state_->wrench);
      break;
    default:
      break;
  }
  return ReturnType::OK;
}

ReturnType ControllerCore::set_reference(
  const InterfaceType & interface_type, const Eigen::VectorXd & reference)
{
  provided_reference_interfaces_ = provided_reference_interfaces_ | interface_type;
  switch (interface_type) {
    case fmc::InterfaceType::POSITION:
      desired_chain_state_->q = reference;
      break;
    case fmc::InterfaceType::VELOCITY:
      desired_chain_state_->dq = reference;
      break;
    case fmc::InterfaceType::ACCELERATION:
      desired_chain_state_->ddq = reference;
      break;
    case fmc::InterfaceType::EFFORT:
      desired_chain_state_->tau = reference;
      current_chain_state_->is_tau_residual = false;
      break;
    case fmc::InterfaceType::RESIDUAL_TORQUES:
      desired_chain_state_->tau = reference;
      current_chain_state_->is_tau_residual = true;
      break;
    case fmc::InterfaceType::TCP_WRENCH:
      desired_chain_state_->wrench = reference;
      break;
    default:
      break;
  }
  return ReturnType::OK;
}

ReturnType ControllerCore::set_state(
  const InterfaceType & interface_type, const std::vector<double> & state)
{
  provided_state_interfaces_ = provided_state_interfaces_ | interface_type;
  switch (interface_type) {
    case fmc::InterfaceType::POSITION:
      vector_to_eigen(state, current_chain_state_->q);
      break;
    case fmc::InterfaceType::VELOCITY:
      vector_to_eigen(state, current_chain_state_->dq);
      break;
    case fmc::InterfaceType::ACCELERATION:
      vector_to_eigen(state, current_chain_state_->ddq);
      break;
    case fmc::InterfaceType::EFFORT:
      vector_to_eigen(state, current_chain_state_->tau);
      current_chain_state_->is_tau_residual = false;
      break;
    case fmc::InterfaceType::RESIDUAL_TORQUES:
      vector_to_eigen(state, current_chain_state_->tau);
      current_chain_state_->is_tau_residual = true;
      break;
    case fmc::InterfaceType::TCP_WRENCH:
      vector_to_eigen(state, current_chain_state_->wrench);
      break;
    default:
      break;
  }
  return ReturnType::OK;
}

ReturnType ControllerCore::set_state(
  const InterfaceType & interface_type, const Eigen::VectorXd & state)
{
  provided_state_interfaces_ = provided_state_interfaces_ | interface_type;
  switch (interface_type) {
    case fmc::InterfaceType::POSITION:
      current_chain_state_->q = state;
      break;
    case fmc::InterfaceType::VELOCITY:
      current_chain_state_->dq = state;
      break;
    case fmc::InterfaceType::ACCELERATION:
      current_chain_state_->ddq = state;
      break;
    case fmc::InterfaceType::EFFORT:
      current_chain_state_->tau = state;
      current_chain_state_->is_tau_residual = false;
      break;
    case fmc::InterfaceType::RESIDUAL_TORQUES:
      current_chain_state_->tau = state;
      current_chain_state_->is_tau_residual = true;
      break;
    case fmc::InterfaceType::TCP_WRENCH:
      current_chain_state_->wrench = state;
      break;
    default:
      break;
  }
  return ReturnType::OK;
}

ReturnType ControllerCore::set_output(
  const InterfaceType & interface_type, const std::vector<double> & output)
{
  provided_command_interfaces_ = provided_command_interfaces_ | interface_type;
  output_[interface_type] = output;
  return ReturnType::OK;
}

ReturnType ControllerCore::set_output(
  const InterfaceType & interface_type, const Eigen::VectorXd & output)
{
  eigen_to_vector(output, vector_temp_);
  return set_output(interface_type, vector_temp_);
}

bool ControllerCore::interfaces_valid(
  const InterfaceType & required, const InterfaceType & provided)
{
  if ((required & provided) != required) {
    return false;
  }
  error_msg = "Error during interfaces validation. Missing interfaces: ";
  // for (auto & i : missing_interfaces) {
  //   error_msg += "\n" + i;
  // };
  return true;
}

ReturnType ControllerCore::update(const std::chrono::nanoseconds & dt)
{
  if (!interfaces_valid(required_state_interfaces_, provided_state_interfaces_)) {
    return ReturnType::ERROR;
  }
  if (!interfaces_valid(required_reference_interfaces_, provided_reference_interfaces_)) {
    return ReturnType::ERROR;
  }

  current_chain_->set_state(current_chain_state_);
  desired_chain_->set_state(desired_chain_state_);
  controller_impl_result_ = update_impl(dt);

  if (!interfaces_valid(required_command_interfaces_, provided_command_interfaces_)) {
    return ReturnType::ERROR;
  }
  provided_state_interfaces_ = InterfaceType::NONE;
  provided_reference_interfaces_ = InterfaceType::NONE;
  provided_command_interfaces_ = InterfaceType::NONE;
  return controller_impl_result_;
}

ReturnType ControllerCore::get_output(
  const InterfaceType & interface_type, std::vector<double> & out_vector)
{
  out_vector.resize(output_.at(interface_type).size());
  std::copy(
    output_.at(interface_type).begin(), output_.at(interface_type).end(), out_vector.begin());
  return ReturnType::OK;
}

ReturnType ControllerCore::get_output(const InterfaceType & interface_type, Eigen::VectorXd & out)
{
  return (vector_to_eigen(output_.at(interface_type), out)) ? ReturnType::OK : ReturnType::ERROR;
}

}  // namespace fmc