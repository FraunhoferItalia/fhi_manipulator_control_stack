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
#include "cartesian_impedance_controller/controller_core.hpp"
namespace cartesian_impedance_controller
{

CartesianImpedanceControllerCore::CartesianImpedanceControllerCore()
: fmc::ControllerCore(
    fmc::InterfaceType::POSITION | fmc::InterfaceType::VELOCITY | fmc::InterfaceType::TCP_WRENCH,
    fmc::InterfaceType::POSITION | fmc::InterfaceType::VELOCITY, fmc::InterfaceType::EFFORT)
{
  cartesian_error = std::make_shared<frdl::CartesianState>();
}

fmc::ControllerParameterDescriptions CartesianImpedanceControllerCore::export_parameters()
{
  return {
    {Parameters::COMPENSATE_GRAVITY,
     fmc::parameters::Description("options.compensate_gravity", true)},
    {Parameters::USE_DESIRED_FRAME,
     fmc::parameters::Description("options.use_desired_frame", true)},
    {Parameters::COMPENSATE_FEXT, fmc::parameters::Description("options.compensate_Fext", false)},
    {Parameters::PERFORMANCE_CONSTRAINTS,
     fmc::parameters::Description("options.performance_constraints", false)},
    {Parameters::USE_BASE_FRAME, fmc::parameters::Description("options.use_base_frame", false)},
    {Parameters::FRICTION_COMPENSATION_RATIO,
     fmc::parameters::Description("friction_compensation_ratio", 0.5, 0.0, 1.0)},
    {Parameters::M_x, fmc::parameters::Description("M.x", 1.0, 0.0, 10000.0)},
    {Parameters::M_y, fmc::parameters::Description("M.y", 1.0, 0.0, 10000.0)},
    {Parameters::M_z, fmc::parameters::Description("M.z", 1.0, 0.0, 10000.0)},  // 1
    {Parameters::M_rx, fmc::parameters::Description("M.rx", 1.0, 0.0, 10000.0)},
    {Parameters::M_ry, fmc::parameters::Description("M.ry", 1.0, 0.0, 10000.0)},
    {Parameters::M_rz, fmc::parameters::Description("M.rz", 1.0, 0.0, 10000.0)},
    {Parameters::K_x, fmc::parameters::Description("K.x", 1200.0, 0.0, 10000.0)},
    {Parameters::K_y, fmc::parameters::Description("K.y", 1200.0, 0.0, 10000.0)},
    {Parameters::K_z, fmc::parameters::Description("K.z", 1200.0, 0.0, 10000.0)},  //100
    {Parameters::K_rx, fmc::parameters::Description("K.rx", 600.0, 0.0, 10000.0)},
    {Parameters::K_ry, fmc::parameters::Description("K.ry", 600.0, 0.0, 10000.0)},
    {Parameters::K_rz, fmc::parameters::Description("K.rz", 100.0, 0.0, 10000.0)},
    {Parameters::V_x, fmc::parameters::Description("gravity_mode.D.x", 0.001, 0.0, 10000.0)},
    {Parameters::V_y, fmc::parameters::Description("gravity_mode.D.y", 0.001, 0.0, 10000.0)},
    {Parameters::V_z, fmc::parameters::Description("gravity_mode.D.z", 0.001, 0.0, 10000.0)},
    {Parameters::V_rx, fmc::parameters::Description("gravity_mode.D.rx", 0.001, 0.0, 10000.0)},
    {Parameters::V_ry, fmc::parameters::Description("gravity_mode.D.ry", 0.001, 0.0, 10000.0)},
    {Parameters::V_rz, fmc::parameters::Description("gravity_mode.D.rz", 0.001, 0.0, 10000.0)},
  };
}

fmc::ReturnType CartesianImpedanceControllerCore::on_activate()
{
  update_compliance_matrices();
  identity = Eigen::MatrixXd::Identity(6, 6);
  zero_vec = Eigen::VectorXd::Zero(current_chain_->get_size());
  tau_out = Eigen::VectorXd::Zero(current_chain_->get_size());

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

fmc::ReturnType CartesianImpedanceControllerCore::update_impl(
  const std::chrono::nanoseconds & /* dt */)
{
  update_compliance_matrices();
  current_cartesian_state = current_chain_->get_cartesian_state();
  desired_cartesian_state = desired_chain_->get_cartesian_state();
  cartesian_state_difference(desired_cartesian_state, current_cartesian_state, cartesian_error);
  tau_out = current_chain_->get_rne();

  J = current_chain_->get_jacobian();
  M = current_chain_->get_M();

  if (get_parameter<bool>(Parameters::USE_BASE_FRAME)) {
    H.setIdentity();
    H_t.setIdentity();
    H_tp.setIdentity();
    H_tp.bottomRightCorner(3, 3) = current_cartesian_state->x.rotation();
  } else if (get_parameter<bool>(Parameters::USE_DESIRED_FRAME)) {
    H.topLeftCorner(3, 3) = desired_cartesian_state->x.rotation();
    H.bottomRightCorner(3, 3) = desired_cartesian_state->x.rotation();
    H_t.topLeftCorner(3, 3) = desired_cartesian_state->x.rotation().transpose();
    H_t.bottomRightCorner(3, 3) = desired_cartesian_state->x.rotation().transpose();
    H_tp.topLeftCorner(3, 3) = desired_cartesian_state->x.rotation().transpose();
    H_tp.bottomRightCorner(3, 3).setIdentity();
  } else {
    H.topLeftCorner(3, 3) = current_cartesian_state->x.rotation();
    H.bottomRightCorner(3, 3) = current_cartesian_state->x.rotation();
    H_t.topLeftCorner(3, 3) = current_cartesian_state->x.rotation().transpose();
    H_t.bottomRightCorner(3, 3) = current_cartesian_state->x.rotation().transpose();
    H_tp.topLeftCorner(3, 3) = current_cartesian_state->x.rotation().transpose();
    H_tp.bottomRightCorner(3, 3).setIdentity();
  }

  // Tau impedance contribution
  if (get_parameter<bool>(Parameters::COMPENSATE_FEXT)) {
    tau_out += M * J.transpose() *
               (J * J.transpose()).completeOrthogonalDecomposition().pseudoInverse() * H *
               Md_.inverse() *
               (Kd_ * H_tp * cartesian_error->vector() + Dd_ * H_t * cartesian_error->dx);
    // Tau ext contribution
    tau_out +=
      J.transpose() * (lambda * Md_.inverse() - identity) * H * current_chain_->get_state()->wrench;
  } else {
    tau_out += J.transpose() * H *
               (Kd_ * H_tp * cartesian_error->vector() + Dd_ * H_t * cartesian_error->dx);
  }

  tau_out += get_parameter<double>(Parameters::FRICTION_COMPENSATION_RATIO) *
             current_chain_->get_torque_friction();

  if (get_parameter<bool>(Parameters::PERFORMANCE_CONSTRAINTS)) {
    tau_out += performance_constrainer_->get_performance_contraint_torque(
      current_chain_->get_state()->q, tau_out);
  }
  tau_out += limits_handler_->get_tau_limits(
    current_chain_->get_state()->q, current_chain_->get_state()->dq, tau_out);
  set_output(fmc::InterfaceType::EFFORT, tau_out);
  return fmc::ReturnType::OK;
}

void CartesianImpedanceControllerCore::update_compliance_matrices()
{
  current_chain_->rne_options.a0[2] =
    (get_parameter<bool>(Parameters::COMPENSATE_GRAVITY)) ? 9.81 : 0.0;
  Md_.setIdentity();
  Kd_.setIdentity();
  Dd_.setIdentity();
  Md_.diagonal() << get_parameter<double>(Parameters::M_x), get_parameter<double>(Parameters::M_y),
    get_parameter<double>(Parameters::M_z), get_parameter<double>(Parameters::M_rx),
    get_parameter<double>(Parameters::M_ry), get_parameter<double>(Parameters::M_rz);
  Kd_.diagonal() << get_parameter<double>(Parameters::K_x), get_parameter<double>(Parameters::K_y),
    get_parameter<double>(Parameters::K_z), get_parameter<double>(Parameters::K_rx),
    get_parameter<double>(Parameters::K_ry), get_parameter<double>(Parameters::K_rz);
  Dd_ = 2 * (Md_ * Kd_).cwiseSqrt();
  Dd_.diagonal()[0] =
    (Dd_.diagonal()[0] < DAMPING_EPS) ? get_parameter<double>(Parameters::V_x) : Dd_.diagonal()[0];
  Dd_.diagonal()[1] =
    (Dd_.diagonal()[1] < DAMPING_EPS) ? get_parameter<double>(Parameters::V_y) : Dd_.diagonal()[1];
  Dd_.diagonal()[2] =
    (Dd_.diagonal()[2] < DAMPING_EPS) ? get_parameter<double>(Parameters::V_z) : Dd_.diagonal()[2];
  Dd_.diagonal()[3] =
    (Dd_.diagonal()[3] < DAMPING_EPS) ? get_parameter<double>(Parameters::V_rx) : Dd_.diagonal()[3];
  Dd_.diagonal()[4] =
    (Dd_.diagonal()[4] < DAMPING_EPS) ? get_parameter<double>(Parameters::V_rx) : Dd_.diagonal()[4];
  Dd_.diagonal()[5] =
    (Dd_.diagonal()[5] < DAMPING_EPS) ? get_parameter<double>(Parameters::V_rx) : Dd_.diagonal()[5];
}

}  // namespace cartesian_impedance_controller
