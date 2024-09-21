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
#include "fhi_manipulator_control/utils/limits_handler.hpp"
#include "fhi_manipulator_control/utils/performance_constrainer.hpp"

namespace gravity_compensation_controller
{

enum Parameters { FRICTION_COMPENSATION_RATIO, COMPENSATE_CORIOLIS, PERFORMANCE_CONSTRAINTS };

class GraivityCompensationController : public fmc::ControllerCore
{
public:
  GraivityCompensationController();

  fmc::ControllerParameterDescriptions export_parameters() override;

  fmc::ReturnType on_activate() override;

  fmc::ReturnType update_impl(const std::chrono::nanoseconds & /* dt */) override;

private:
  Eigen::VectorXd tau_;
  fmc::LimitsHandler::SharedPtr limits_handler_;
  fmc::PerformanceConstrainer::SharedPtr performance_constrainer_;
};
}  // namespace gravity_compensation_controller
