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
#pragma once
#include "fhi_manipulator_control/controller_core.hpp"
#include "robust_controller.hpp"

namespace robust_controller
{

class RobustControllerCore : public fmc::ControllerCore
{
public:
  RobustControllerCore()
  : fmc::ControllerCore(
      fmc::InterfaceType::POSITION | fmc::InterfaceType::VELOCITY,
      fmc::InterfaceType::POSITION | fmc::InterfaceType::VELOCITY |
        fmc::InterfaceType::ACCELERATION,
      fmc::InterfaceType::EFFORT){};

  fmc::ReturnType init_params(const std::string & parameters_file_path) override
  {
    robust_controller_.init(parameters_file_path);
    return fmc::ReturnType::OK;
  }

  fmc::ReturnType on_activate() override
  {
    q_.resize(modules_->size());
    dq_.resize(modules_->size());
    q_ref_.resize(modules_->size());
    dq_ref_.resize(modules_->size());
    ddq_ref_.resize(modules_->size());
    std::cout << "Robust controller: on_activate" << std::endl;
    return fmc::ReturnType::OK;
  }

  fmc::ReturnType update_impl(const std::chrono::nanoseconds & dt) override
  {
    dt_sec = dt.count() / 1e9;
    fmc::eigen_to_vector(current_chain_state_->q, q_);
    fmc::eigen_to_vector(current_chain_state_->dq, dq_);
    fmc::eigen_to_vector(desired_chain_state_->q, q_ref_);
    fmc::eigen_to_vector(desired_chain_state_->dq, dq_ref_);
    fmc::eigen_to_vector(desired_chain_state_->ddq, ddq_ref_);

    robust_controller_.set_states(q_, dq_);
    robust_controller_.set_reference(q_ref_, dq_ref_, ddq_ref_);
    robust_controller_.update();
    set_output(fmc::InterfaceType::EFFORT, robust_controller_.get_torque());
    return fmc::ReturnType::OK;
  };

  std::map<std::string, double *> export_register_variables() override
  {
    std::map<std::string, double *> registered_variables;
    for (size_t i = 0; i < modules_->size(); i++) {
      registered_variables.insert({"tau_" + std::to_string(i), &robust_controller_.tau_[i]});
      registered_variables.insert(
        {"tau_robust_" + std::to_string(i), &robust_controller_.tauRobust_[i]});
    }
    registered_variables.insert({"gain", &robust_controller_.gainTuneContr_});
    return registered_variables;
  }

  bool gravity_mode{false};
  double dt_sec;

private:
  RobustController robust_controller_ = RobustController();
  std::vector<double> q_, dq_;
  std::vector<double> q_ref_, dq_ref_, ddq_ref_;
};

}  // namespace robust_controller