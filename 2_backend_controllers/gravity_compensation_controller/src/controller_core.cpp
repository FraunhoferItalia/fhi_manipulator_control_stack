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
#include "gravity_compensation_controller/controller_core.hpp"

namespace gravity_compensation_controller
{
GraivityCompensationController::GraivityCompensationController()
: fmc::ControllerCore(
    fmc::InterfaceType::POSITION, fmc::InterfaceType::NONE, fmc::InterfaceType::EFFORT)
{
}

fmc::ControllerParameterDescriptions GraivityCompensationController::export_parameters()
{
  return {
    {Parameters::FRICTION_COMPENSATION_RATIO,
     fmc::parameters::Description("friction_compensation_ratio", 0.5, 0.0, 1.0)},
    {Parameters::COMPENSATE_CORIOLIS, fmc::parameters::Description("compensate_coriolis", false)},
    {Parameters::PERFORMANCE_CONSTRAINTS,
     fmc::parameters::Description("performance_constraints", false)}};
}

fmc::ReturnType GraivityCompensationController::on_activate()
{
  tau_.resize(modules_->size());
  if (!limits_handler_) {
    limits_handler_ =
      std::make_shared<fmc::LimitsHandler>(current_chain_->get_modules(), 5.0 * 3.14 / 180);
  }
  if (!performance_constrainer_) {
    performance_constrainer_ =
      std::make_shared<fmc::PerformanceConstrainer>(current_chain_->get_modules());
  }
  return fmc::ReturnType::OK;
}

fmc::ReturnType GraivityCompensationController::update_impl(
  const std::chrono::nanoseconds & /* dt */)
{
  tau_ = (get_parameter<bool>(Parameters::COMPENSATE_CORIOLIS))
           ? current_chain_->get_rne()
           : current_chain_->get_torque_gravity();
  tau_ += get_parameter<double>(Parameters::FRICTION_COMPENSATION_RATIO) *
          current_chain_->get_torque_friction();

  if (get_parameter<bool>(Parameters::PERFORMANCE_CONSTRAINTS)) {
    tau_ += performance_constrainer_->get_performance_contraint_torque(
      current_chain_->get_state()->q, tau_);
  }

  tau_ += limits_handler_->get_tau_limits(
    current_chain_->get_state()->q, current_chain_->get_state()->dq, tau_);
  set_output(fmc::InterfaceType::EFFORT, tau_);
  return fmc::ReturnType::OK;
}
}  // namespace gravity_compensation_controller