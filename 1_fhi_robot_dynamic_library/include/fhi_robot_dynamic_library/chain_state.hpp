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
#include <memory>

namespace frdl
{

class JointState
{
public:
  using SharedPtr = std::shared_ptr<JointState>;
  JointState(
    const Eigen::VectorXd & position, const Eigen::VectorXd & velocity = Eigen::VectorXd::Zero(6),
    const Eigen::VectorXd & acceleration = Eigen::VectorXd::Zero(6),
    const Eigen::VectorXd & torques = Eigen::VectorXd::Zero(6),
    const Eigen::VectorXd & external_wrench = Eigen::VectorXd::Zero(6),
    const bool & use_residual_torques = false)
  : q{position},
    dq{velocity},
    ddq{acceleration},
    tau{torques},
    wrench{external_wrench},
    is_tau_residual{use_residual_torques}
  {
    dq = Eigen::VectorXd::Zero(q.size());
    ddq = Eigen::VectorXd::Zero(q.size());
    tau = Eigen::VectorXd::Zero(q.size());
  };
  Eigen::VectorXd q;
  Eigen::VectorXd dq;
  Eigen::VectorXd ddq;
  Eigen::VectorXd tau;
  Eigen::VectorXd wrench;
  bool is_tau_residual = false;
};

class CartesianState
{
public:
  using SharedPtr = std::shared_ptr<CartesianState>;
  CartesianState(
    const Eigen::Affine3d & position = Eigen::Affine3d::Identity(),
    const Eigen::VectorXd & velocity = Eigen::VectorXd::Zero(6),
    const Eigen::VectorXd & acceleration = Eigen::VectorXd::Zero(6),
    const Eigen::VectorXd & forces_torques = Eigen::VectorXd::Zero(6))
  : x{position}, dx{velocity}, ddx{acceleration}, wrench{forces_torques} {};
  Eigen::VectorXd vector()
  {
    Eigen::Quaterniond q(x.rotation());
    vector_state.head(3) << x.translation();
    vector_state.tail(3) << q.coeffs().head(3);
    return vector_state;
  }
  Eigen::Affine3d x = Eigen::Affine3d::Identity();
  Eigen::VectorXd dx = Eigen::VectorXd::Zero(6);
  Eigen::VectorXd ddx = Eigen::VectorXd::Zero(6);
  Eigen::VectorXd wrench = Eigen::VectorXd::Zero(6);
  Eigen::VectorXd vector_state = Eigen::VectorXd::Zero(6);
};

static Eigen::Quaterniond lsupport_quaternion;
static Eigen::Quaterniond rsupport_quaternion;

void cartesian_state_difference(
  const CartesianState & lcartesian_state, const CartesianState & rcartesian_state,
  CartesianState & result)
{
  lsupport_quaternion = lcartesian_state.x.rotation();
  rsupport_quaternion = rcartesian_state.x.rotation();

  result.x.translation() = lcartesian_state.x.translation() - rcartesian_state.x.translation();
  // To interpolate the "short way"
  // https://math.stackexchange.com/questions/3572459/how-to-compute-the-orientation-error-between-two-3d-coordinate-frames
  if (lsupport_quaternion.coeffs().dot(rsupport_quaternion.coeffs()) < 0.0)
    rsupport_quaternion.coeffs() << -rsupport_quaternion.coeffs();

  result.x.linear() = (rsupport_quaternion.inverse() * lsupport_quaternion).toRotationMatrix();
  result.dx = lcartesian_state.dx - rcartesian_state.dx;
  result.ddx = lcartesian_state.ddx - rcartesian_state.ddx;
  result.wrench = lcartesian_state.wrench - rcartesian_state.wrench;
}

void cartesian_state_difference(
  const CartesianState::SharedPtr & lcartesian_state,
  const CartesianState::SharedPtr & rcartesian_state, const CartesianState::SharedPtr & result)
{
  cartesian_state_difference(*lcartesian_state.get(), *rcartesian_state.get(), *result.get());
}

void cartesian_state_sum(
  const CartesianState & lcartesian_state, const CartesianState & rcartesian_state,
  CartesianState & result)
{
  lsupport_quaternion = lcartesian_state.x.rotation();
  rsupport_quaternion = rcartesian_state.x.rotation();

  result.x.translation() = lcartesian_state.x.translation() + rcartesian_state.x.translation();
  result.x.linear() = (rsupport_quaternion * lsupport_quaternion).toRotationMatrix();
  result.dx = lcartesian_state.dx + rcartesian_state.dx;
  result.ddx = lcartesian_state.ddx + rcartesian_state.ddx;
  result.wrench = lcartesian_state.wrench + rcartesian_state.wrench;
}

void cartesian_state_sum(
  const CartesianState::SharedPtr & lcartesian_state,
  const CartesianState::SharedPtr & rcartesian_state, const CartesianState::SharedPtr & result)
{
  cartesian_state_sum(*lcartesian_state.get(), *rcartesian_state.get(), *result.get());
}

}  // namespace frdl
