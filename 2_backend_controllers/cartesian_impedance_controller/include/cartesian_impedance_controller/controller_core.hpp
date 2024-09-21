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
#include "fhi_manipulator_control/utils/performance_constrainer.hpp"
#include "fhi_manipulator_control/utils/limits_handler.hpp"

namespace cartesian_impedance_controller
{

static double DAMPING_EPS = 1e-06;

Eigen::Matrix3d skew(const Eigen::Vector3d & v)
{
  Eigen::Matrix3d skewMatrix;
  skewMatrix << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

  return skewMatrix;
}

enum Parameters {
  COMPENSATE_GRAVITY,
  FRICTION_COMPENSATION_RATIO,
  USE_DESIRED_FRAME,
  COMPENSATE_FEXT,
  PERFORMANCE_CONSTRAINTS,
  USE_BASE_FRAME,
  M_x,
  M_y,
  M_z,
  M_rx,
  M_ry,
  M_rz,
  K_x,
  K_y,
  K_z,
  K_rx,
  K_ry,
  K_rz,
  V_x,
  V_y,
  V_z,
  V_rx,
  V_ry,
  V_rz
};

class CartesianImpedanceControllerCore : public fmc::ControllerCore
{
public:
  CartesianImpedanceControllerCore();

  fmc::ControllerParameterDescriptions export_parameters() override;

  fmc::ReturnType on_activate() override;

  fmc::ReturnType update_impl(const std::chrono::nanoseconds & dt) override;

  void update_compliance_matrices();

private:
  Eigen::VectorXd tau_out, tau_ns, zero_vec;
  Eigen::MatrixXd identity;
  Eigen::MatrixXd J, J_sharp, J_sharp_T, M, delta_lambda, lambda;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_;

  Eigen::EigenSolver<Eigen::MatrixXd> solver;
  Eigen::MatrixXd U_, V_, Mu_, Ku_, Du_, Dmat_;

  frdl::CartesianState::SharedPtr cartesian_error;
  frdl::CartesianState::SharedPtr current_cartesian_state;
  frdl::CartesianState::SharedPtr desired_cartesian_state;

  Eigen::MatrixXd Md_ = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd Dd_ = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd Kd_ = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd Kns_ = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd Dns_ = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd H_t = Eigen::MatrixXd::Identity(6, 6);
  Eigen::MatrixXd H_tp = Eigen::MatrixXd::Identity(6, 6);

  fmc::LimitsHandler::SharedPtr limits_handler_;
  fmc::PerformanceConstrainer::SharedPtr performance_constrainer_;
};

}  // namespace cartesian_impedance_controller
