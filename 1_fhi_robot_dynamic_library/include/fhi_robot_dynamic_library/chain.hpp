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
#include <eigen3/Eigen/Dense>
#include <fhi_robot_dynamic_library/chain_state.hpp>
#include <fhi_robot_dynamic_library/module.hpp>

namespace frdl
{
namespace Kinematics
{

bool get_forward_kinematic(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, int link_index,
  Eigen::Affine3d & result, bool include_fixed_joints = true)
{
  if (
    link_index < 1 || static_cast<size_t>(link_index) > chain_modules.size() ||
    chain_modules.size() != static_cast<size_t>(q.size())) {
  }
  result.setIdentity();
  int i = 0;
  auto modules = (include_fixed_joints) ? chain_modules.modules : chain_modules.active_modules;
  for (auto module : modules) {
    switch (module.type) {
      case ModuleType::REVOLUTE:
        result = result * module.origin *
                 Eigen::Affine3d(Eigen::AngleAxisd(q[i], module.axis.normalized()));
        break;
      case ModuleType::PRISMATIC:
        break;
      case ModuleType::FIXED:
        if (include_fixed_joints) {
          result = result * module.origin;
        }
        break;
      default:
        std::cout << "Unrecognized module type " << module.type << std::endl;
        break;
    }
    if (i >= link_index) return true;
    i++;
  }
  return true;
}

bool get_forward_kinematic(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, Eigen::Affine3d & result)
{
  return get_forward_kinematic(chain_modules, q, q.size(), result);
}

bool get_jacobian(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, Eigen::MatrixXd & J)
{
  Eigen::Affine3d T_a;
  Eigen::Vector3d p_tilde;
  std::vector<Eigen::Vector3d> p;
  std::vector<Eigen::Vector3d> z;
  Eigen::Vector3d z0 = Eigen::Vector3d::Zero();
  z0[2] = 1.0;
  p.push_back(chain_modules.modules[0].origin.translation());
  z.push_back(chain_modules.modules[0].origin.rotation().matrix() * z0);
  for (size_t i = 0; i < chain_modules.size(); i++) {
    get_forward_kinematic(chain_modules, q, i + 1, T_a);
    p.push_back(T_a.translation());
    z.push_back(T_a.rotation() * z0);
  }
  J.setZero();
  for (size_t i = 0; i < chain_modules.size(); i++) {
    if (chain_modules.active_modules[i].type == ModuleType::REVOLUTE) {
      J.block<3, 1>(0, i) = z[i].cross(p.back() - p[i]);
      J.block<3, 1>(3, i) = z[i];
    } else if (chain_modules.active_modules[i].type == ModuleType::PRISMATIC) {
      J.block<3, 1>(0, i) = z[i];
      J.block<3, 1>(3, i) = Eigen::Vector3d::Zero();
    } else {
      std::cout << "Unrecognized module type " << static_cast<int>(chain_modules.modules[i].type)
                << std::endl;
    }
  }

  return true;
}

struct IKParameters
{
  using SharedPtr = std::shared_ptr<IKParameters>;
  IKParameters(const Eigen::VectorXd seed) : q0{seed} {}
  int max_iterations{10000};
  double max_error_pos{1e-4};
  double max_error_ori{1e-4};
  double K{0.8};
  double k_null_rev{100.};
  double k_null_pris{7.};
  bool use_transpose{false};
  bool position_only{false};
  double dls_threshold{0.0};
  double max_distance_from_seed{std::numeric_limits<double>::infinity()};
  Eigen::VectorXd q0;
};

bool get_inverse_kinematic(
  const ChainModules & chain_modules, const Eigen::Affine3d & x, const IKParameters & ik_parameters,
  Eigen::VectorXd & solution)
{
  Eigen::VectorXd q = ik_parameters.q0;
  Eigen::Affine3d x_sol;
  Eigen::MatrixXd J_sol;
  Eigen::MatrixXd J_sol_pinv;
  Eigen::VectorXd error_sol(6);
  J_sol.resize(6, chain_modules.size());
  J_sol.setZero();
  error_sol.setZero();
  Eigen::VectorXd dq_null(chain_modules.size());
  dq_null.setZero();
  Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(chain_modules.size(), chain_modules.size());
  int count{0};

  Eigen::Quaterniond x_quat(x.rotation());
  Eigen::Quaterniond x_sol_quat;
  Eigen::Matrix3d reference_skew;
  reference_skew << 0, -x_quat.z(), x_quat.y(), x_quat.z(), 0, -x_quat.x(), -x_quat.y(), x_quat.x(),
    0;

  double lower = -3.14;
  double upper = 3.14;
  while (true) {
    get_forward_kinematic(chain_modules, q, x_sol);
    get_jacobian(chain_modules, q, J_sol);
    error_sol.head(3) = x.translation() - x_sol.translation();
    x_sol_quat = x_sol.rotation();
    error_sol.tail(3) = x_sol_quat.w() * x_quat.coeffs().head(3) -
                        x_quat.w() * x_sol_quat.coeffs().head(3) -
                        reference_skew * x_sol_quat.coeffs().head(3);
    double jlm;
    double k_null;
    bool violated{false};
    std::vector<int> violated_indexes{};
    for (size_t i = 0; i < chain_modules.size(); i++) {
      jlm = (chain_modules.active_modules[i].position_limits[0] +
             chain_modules.active_modules[i].position_limits[1]) /
            2;
      k_null = (chain_modules.active_modules[i].type == ModuleType::REVOLUTE)
                 ? ik_parameters.k_null_rev
                 : ik_parameters.k_null_pris;
      dq_null[i] = k_null * (-1.0 / chain_modules.active_modules.size()) * (q[i] - jlm) /
                   std::pow(
                     chain_modules.active_modules[i].position_limits[0] -
                       chain_modules.active_modules[i].position_limits[1],
                     2);
      if (
        q[i] >= chain_modules.active_modules[i].position_limits[1] ||
        q[i] <= chain_modules.active_modules[i].position_limits[0]) {
        violated = true;
        violated_indexes.push_back(i);
      }
    }
    if (ik_parameters.position_only) {
      J_sol.bottomRows(3).setZero();
    }
    if (ik_parameters.use_transpose || ik_parameters.position_only) {
      q = q + J_sol.transpose() * ik_parameters.K * error_sol;
    } else {
      J_sol_pinv = J_sol.completeOrthogonalDecomposition().pseudoInverse();
      q = q + J_sol_pinv * ik_parameters.K * error_sol + (eye - J_sol_pinv * J_sol) * dq_null;
    }

    // Wraps solution to a given limit (i.e [-pi,pi])
    // (AndreaG: why not joint limits?)
    for (int i = 0; i < q.size(); ++i) {
      q(i) = lower + std::fmod(q(i) - lower, upper - lower);
    }

    if (
      error_sol.head(3).norm() < ik_parameters.max_error_pos &&
      (error_sol.tail(3).norm() < ik_parameters.max_error_ori || ik_parameters.position_only) &&
      !violated && (q - ik_parameters.q0).norm() < ik_parameters.max_distance_from_seed) {
      solution = q;
      return true;
    } else if (count > ik_parameters.max_iterations) {
      std::cout << "pos error: " << error_sol.head(3).norm() << std::endl;
      std::cout << "ori error: " << error_sol.tail(3).norm() << std::endl;
      std::cout << "seed distance: " << (q - ik_parameters.q0).norm() << std::endl;
      std::cout << "violated: " << violated << std::endl;
      return false;
    }
    count++;
  }
  return true;
}

bool get_jacobian_dot(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, const Eigen::VectorXd & dq,
  Eigen::MatrixXd & J, Eigen::MatrixXd & J_dot)
{
  Eigen::Affine3d T_a;
  Eigen::Vector3d p_tilde;
  std::vector<Eigen::Vector3d> p;
  std::vector<Eigen::Vector3d> z;
  std::vector<Eigen::Vector3d> omega;
  std::vector<Eigen::Vector3d> dp;
  Eigen::Vector3d z0 = Eigen::Vector3d::Zero();
  z0[2] = 1.0;
  p.push_back(chain_modules.modules[0].origin.translation());
  z.push_back(chain_modules.modules[0].origin.rotation().matrix() * z0);
  omega.push_back(Eigen::Vector3d::Zero());
  dp.push_back(Eigen::Vector3d::Zero());
  for (size_t i = 0; i < chain_modules.size(); i++) {
    get_forward_kinematic(chain_modules, q, i + 1, T_a, false);
    p.push_back(T_a.translation());
    z.push_back(T_a.rotation() * z0);
  }
  if (chain_modules.modules.back().type == ModuleType::FIXED) {
    p.push_back(
      T_a.rotation() * chain_modules.modules.back().origin.translation() + T_a.translation());
    z.push_back(T_a.rotation() * z0);
  }
  J.setZero();
  for (size_t i = 0; i < chain_modules.size(); i++) {
    if (chain_modules.active_modules[i].type == ModuleType::REVOLUTE) {
      J.block<3, 1>(0, i) = z[i].cross(p.back() - p[i]);
      J.block<3, 1>(3, i) = z[i];
      omega.push_back(z[i] * dq[i] + omega[i]);
      dp.push_back(omega.back().cross(p[i + 1] - p[i]) + dp[i]);
    } else if (chain_modules.active_modules[i].type == ModuleType::PRISMATIC) {
      J.block<3, 1>(0, i) = z[i];
      J.block<3, 1>(3, i) = Eigen::Vector3d::Zero();
      omega.push_back(omega[i]);
      dp.push_back(omega.back().cross(p[i + 1] - p[i]) + dp[i] + z[i] * dq[i]);
    } else {
      std::cout << "Unrecognized module type "
                << static_cast<int>(chain_modules.active_modules[i].type) << std::endl;
    }
  }
  J_dot.setZero();
  Eigen::Vector3d dp_e = J.middleRows(0, 3) * dq;
  for (size_t i = 0; i < chain_modules.size(); i++) {
    if (chain_modules.active_modules[i].type == ModuleType::REVOLUTE) {
      J_dot.block<3, 1>(0, i) =
        (omega[i].cross(z[i])).cross(p.back() - p[i]) + z[i].cross(dp_e - dp[i]);
      J_dot.block<3, 1>(3, i) = omega[i].cross(z[i]);
    } else if (chain_modules.active_modules[i].type == ModuleType::PRISMATIC) {
      J_dot.block<3, 1>(0, i) = omega[i].cross(z[i]);
      J_dot.block<3, 1>(3, i) = Eigen::Vector3d::Zero();
    } else {
      std::cout << "Unrecognized module type "
                << static_cast<int>(chain_modules.active_modules[i].type) << std::endl;
    }
  }
  return true;
}
}  // namespace Kinematics
namespace Dynamics
{

struct RNEOptions
{
  Eigen::Vector3d w0 = Eigen::Vector3d::Zero();
  Eigen::Vector3d w_aux0 = Eigen::Vector3d::Zero();
  Eigen::Vector3d dw0 = Eigen::Vector3d::Zero();
  Eigen::Vector3d a0 = Eigen::Vector3d(0.0, 0.0, 9.81);
  Eigen::VectorXd ext_wrench = Eigen::VectorXd::Zero(6);
  Eigen::Affine3d ext_forces_rotation = Eigen::Affine3d::Identity();

  static RNEOptions zero_gravity()
  {
    RNEOptions iframe_info;
    iframe_info.a0.setZero();
    return iframe_info;
  }
};

// https://bitbucket.org/MatthiasAlthoff/improv/src/master/software/dyn_fcns/NE_mod.m
bool get_rne(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, const Eigen::VectorXd & dq,
  const Eigen::VectorXd & dq_aux, const Eigen::VectorXd & ddq, Eigen::VectorXd & torque,
  const RNEOptions & rne_options = RNEOptions())
{
  // Calculate rotation matrices for frame i to i-i and i+1 to i
  if (static_cast<size_t>(torque.size()) != chain_modules.size()) {
    return false;
  }
  size_t k = chain_modules.modules.size();  // number of joint arms + 1
  std::vector<Eigen::Vector3d> omega, omegaDot, omegaA, a_j, a_c, r_link;
  std::vector<Eigen::Matrix3d> R, R_T;
  omega.reserve(k);
  omegaDot.reserve(k);
  omegaA.reserve(k);
  a_j.reserve(k);
  a_c.reserve(k);
  r_link.reserve(k);
  R.reserve(k + 1);
  R_T.reserve(k);

  Eigen::Vector3d z0 = Eigen::Vector3d(0.0, 0.0, 1.0);
  Eigen::Matrix3d identity_matrix3d = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d rDof = Eigen::Matrix3d::Identity();
  Eigen::Affine3d H;
  for (size_t i = 0; i < k; i++) {
    H = chain_modules.modules[i].origin;
    if (chain_modules.modules[i].type == ModuleType::REVOLUTE) {
      rDof = Eigen::AngleAxisd(q[i], chain_modules.modules[i].axis.normalized());
    } else if (chain_modules.modules[i].type == ModuleType::FIXED) {
      rDof = identity_matrix3d;
    }
    R[i] = H.rotation() * rDof;
    R_T[i] = R[i].transpose();
    r_link[i] = H.translation();
  }

  // Forward recursion
  for (size_t i = 0; i < k; i++) {
    if (i == 0) {
      a_j[i] = rne_options.a0;
      omega[i] = rne_options.w0;
      omegaA[i] = rne_options.w_aux0;
      omegaDot[i] = rne_options.dw0;
    } else {
      if (chain_modules.modules[i - 1].type == ModuleType::REVOLUTE) {
        // omega(:,i) = R_T(:,:,i-1)*(omega(:,i-1)) + dq(i-1)*z0; % dq(i-1) because it starts from i = 1
        omega[i] = R_T[i - 1] * omega[i - 1] + dq[i - 1] * z0;
        omegaA[i] = R_T[i - 1] * omegaA[i - 1] + dq_aux[i - 1] * z0;
        // omegaDot(:,i) = R_T(:,:,i-1)*(omegaDot(:,i-1)) + Xcross(R_T(:,:,i-1)*omega(:,i-1))*dq(i-1)*z0 + ddq(i-1)*z0;
        omegaDot[i] = R_T[i - 1] * omegaDot[i - 1] + ddq[i - 1] * z0 +
                      (R_T[i - 1] * omega[i - 1]).cross(dq_aux[i - 1] * z0);
        // a_j(:,i) = R_T(:,:,i-1)*(a_j(:,i-1) + Xcross(omegaDot(:,i-1))*r_link(:,i-1) + Xcross(omega(:,i-1))*Xcross(omega(:,i-1))*r_link(:,i-1));
        a_j[i] = R_T[i - 1] * (a_j[i - 1] + omegaDot[i - 1].cross(r_link[i - 1]) +
                               omegaA[i - 1].cross(omega[i - 1].cross(r_link[i - 1])));
      } else if (chain_modules.modules[i - 1].type == ModuleType::PRISMATIC) {
        std::cout << "Unsupported module" << std::endl;
        /* omega[i] = R_T[i - 1] * omega[i - 1];
        omegaA[i] = R_T[i - 1] * omegaA[i - 1];
        omegaDot[i] = R_T[i - 1] * omegaDot[i - 1];
        a_j[i] = R_T[i - 1] * a_j[i - 1] + omegaDot[i - 1].cross(r_link[i - 1]) +
                 omegaA[i - 1].cross(omega[i - 1].cross(r_link[i - 1])) +
                 omegaA[i - 1].cross(dq[i - 1] * z0) + omega[i - 1].cross(dq_aux[i - 1] * z0) +
                 ddq[i - 1] * z0; */
      } else if (chain_modules.modules[i - 1].type == ModuleType::FIXED) {
        omega[i] = R_T[i - 1] * omega[i - 1];
        omegaA[i] = R_T[i - 1] * omegaA[i - 1];
        omegaDot[i] = R_T[i - 1] * omegaDot[i - 1];
        a_j[i] = R_T[i - 1] * (a_j[i - 1] + omegaDot[i - 1].cross(r_link[i - 1]) +
                               omegaA[i - 1].cross(omega[i - 1].cross(r_link[i - 1])));
      } else {
        std::cout << "Unrecognized module type " << static_cast<int>(chain_modules.modules[i].type)
                  << std::endl;
      }
    }
    a_c[i] = a_j[i] + omegaDot[i].cross(chain_modules.modules[i].CoM) +
             omegaA[i].cross(omega[i].cross(chain_modules.modules[i].CoM));
  }

  // Backward recursion
  Eigen::Vector3d mu, f_next, f;
  for (size_t i = k; i >= 1; i--) {
    if (i == k) {
      f_next = rne_options.ext_wrench.head(3);
      mu = rne_options.ext_wrench.head(3);
      R[i] = rne_options.ext_forces_rotation.rotation();
    }
    f = R[i] * f_next + chain_modules.modules[i].mass * a_c[i];
    mu = R[i] * mu + (R[i] * f_next).cross(chain_modules.modules[i].CoM - r_link[i]) -
         f.cross(chain_modules.modules[i].CoM) +
         chain_modules.modules[i].inertia_matrix * omegaDot[i] +
         omega[i].cross(Eigen::Vector3d(chain_modules.modules[i].inertia_matrix * omegaA[i]));
    f_next = f;
    if (chain_modules.modules[i - 1].type == ModuleType::REVOLUTE) {
      torque[i - 1] = mu.dot(z0);
    } else if (chain_modules.modules[i - 1].type == ModuleType::PRISMATIC) {
      torque[i - 1] = f.dot(z0);
    } else if (chain_modules.modules[i - 1].type == ModuleType::FIXED) {
    } else {
      std::cout << "Module " << i - 1 << "Unrecognized module type "
                << static_cast<int>(chain_modules.modules[i - 1].type) << std::endl;
    }
  }

  return true;
}

bool get_torque_gravity(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, Eigen::VectorXd & torque_gravity)
{
  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(chain_modules.size());
  return get_rne(chain_modules, q, zero_vec, zero_vec, zero_vec, torque_gravity);
}

bool get_torque_friction(
  const ChainModules & chain_modules, const Eigen::VectorXd & dq, Eigen::VectorXd & torque_friction,
  bool const & include_static)
{
  if (
    (chain_modules.size() != static_cast<size_t>(dq.size())) ||
    chain_modules.size() != static_cast<size_t>(torque_friction.size())) {
    return false;
  }
  for (size_t i = 0; i < chain_modules.size(); i++) {
    torque_friction[i] = chain_modules.modules[i].dynamic_friction * dq[i];
    if (std::abs(dq[i]) > 0.0 && include_static)
      torque_friction[i] += chain_modules.modules[i].static_friction * dq[i] / std::abs(dq[i]);
  }
  return true;
}

bool get_M(const ChainModules & chain_modules, const Eigen::VectorXd & q, Eigen::MatrixXd & M)
{
  // Check that the provided size of M is correct or do something
  Eigen::VectorXd ddq_aux = Eigen::VectorXd::Zero(M.cols());
  Eigen::VectorXd M_col = Eigen::VectorXd::Zero(M.cols());
  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(M.cols());
  for (long int i = 0; i < M.cols(); i++) {
    ddq_aux[i] = 1.0;
    get_rne(chain_modules, q, zero_vec, zero_vec, ddq_aux, M_col, RNEOptions::zero_gravity());
    M.col(i) = M_col;
    ddq_aux[i] = 0.0;
  }

  return true;
}

bool get_C(
  const ChainModules & chain_modules, const Eigen::VectorXd & q, const Eigen::VectorXd & dq,
  Eigen::MatrixXd & C)
{
  Eigen::VectorXd dq_aux(C.cols());
  Eigen::VectorXd C_col(C.cols());
  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(C.cols());
  dq_aux.setZero();
  C_col.setZero();
  for (long int i = 0; i < C.cols(); i++) {
    dq_aux[i] = 1;
    get_rne(chain_modules, q, dq, dq_aux, zero_vec, C_col, RNEOptions::zero_gravity());
    C.col(i) = C_col;
    dq_aux[i] = 0;
  }

  return true;
}

}  // namespace Dynamics

struct MethodUpToDate
{
  void reset()
  {
    J = false;
    M = false;
    C = false;
    g = false;
    rne = false;
    ik = false;
    friction = false;
  }
  bool J{false};  // Always J_dot is also computed
  bool M{false};
  bool C{false};
  bool g{false};
  bool rne{false};
  bool ik{false};
  bool friction{false};
};

class Chain
{
public:
  using SharedPtr = std::shared_ptr<Chain>;
  using UniquePtr = std::unique_ptr<Chain>;

  Chain(ChainModules chain_modules) : modules_{chain_modules}
  {
    zero_vec_ = Eigen::VectorXd::Zero(chain_modules.size());
    ik_parameters = std::make_shared<frdl::Kinematics::IKParameters>(zero_vec_);
    state_ = std::make_shared<JointState>(Eigen::VectorXd::Zero(chain_modules.size()));
    cartesian_state_ = std::make_shared<CartesianState>();

    J = Eigen::MatrixXd::Zero(6, chain_modules.size());
    J_pinv = Eigen::MatrixXd::Zero(chain_modules.size(), 6);
    J_dot = Eigen::MatrixXd::Zero(6, chain_modules.size());
    M = Eigen::MatrixXd::Zero(chain_modules.size(), chain_modules.size());
    C = Eigen::MatrixXd::Zero(chain_modules.size(), chain_modules.size());
    g = Eigen::VectorXd::Zero(chain_modules.size());
    friction_torque = Eigen::VectorXd::Zero(chain_modules.size());
    rne_torque = Eigen::VectorXd::Zero(chain_modules.size());
    temp_ = Eigen::VectorXd::Zero(chain_modules.size());
  }

  Chain(ChainModules::SharedPtr const & chain_modules) : Chain(*chain_modules.get()) {}

  /* Chain(std::vector<Module> chain_modules, JointState::SharedPtr state) : modules_{chain_modules}, state_{state}
  {
    zero_vec_ = Eigen::VectorXd::Zero(chain_modules.size());
    ik_parameters = std::make_shared<frdl::Kinematics::IKParameters>(zero_vec_);
  };
  Chain(std::vector<Module> chain_modules, CartesianState::SharedPtr state) : modules_{chain_modules}
  {
    zero_vec_ = Eigen::VectorXd::Zero(chain_modules.size());
    ik_parameters = std::make_shared<frdl::Kinematics::IKParameters>(zero_vec_);
    set_state(state);
  }; */

  bool set_state(CartesianState::SharedPtr const & state)
  {
    method_up_to_date.reset();
    if (!frdl::Kinematics::get_inverse_kinematic(
          modules_, state->x, *ik_parameters.get(), state_->q)) {
      std::cout << "Ik failed" << std::endl;
      return false;
    }
    frdl::Kinematics::get_jacobian(modules_, state_->q, J);
    J_svd.compute(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
    J_pinv =
      J_svd.matrixV() * J_svd.singularValues().asDiagonal().inverse() * J_svd.matrixU().transpose();

    state_->dq = J_pinv * state->dx;

    frdl::Kinematics::get_jacobian_dot(modules_, state_->q, state_->dq, J, J_dot);
    state_->ddq = J_pinv * (state->ddx - J_dot * state_->dq);
    state_->tau = J.transpose() * state->wrench;
    state_->is_tau_residual = true;

    method_up_to_date.J = true;
    method_up_to_date.ik = true;
    return true;
  }

  bool set_state(CartesianState::SharedPtr const & state, const Eigen::VectorXd & seed)
  {
    ik_parameters->q0 = seed;
    return set_state(state);
  }

  bool set_state(JointState::SharedPtr const & state)
  {
    method_up_to_date.reset();
    state_ = state;
    return true;
  }

  JointState::SharedPtr & get_state() { return state_; }

  CartesianState::SharedPtr & get_cartesian_state()
  {
    frdl::Kinematics::get_forward_kinematic(modules_, state_->q, cartesian_state_->x);
    if (!method_up_to_date.J) {
      frdl::Kinematics::get_jacobian_dot(modules_, state_->q, state_->dq, J, J_dot);
      method_up_to_date.J = true;
    }
    cartesian_state_->dx = J * state_->dq;
    cartesian_state_->ddx = J_dot * state_->dq + J * state_->ddq;
    if (!state_->is_tau_residual) {
      // temp_ =
    } else {
      temp_ = state_->tau;
    }
    J_svd.compute(J.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    J_pinv =
      J_svd.matrixV() * J_svd.singularValues().asDiagonal().inverse() * J_svd.matrixU().transpose();
    cartesian_state_->wrench = J_pinv * temp_;
    return cartesian_state_;
  }

  Eigen::MatrixXd & get_jacobian()
  {
    if (method_up_to_date.J) return J;
    frdl::Kinematics::get_jacobian_dot(modules_, state_->q, state_->dq, J, J_dot);
    method_up_to_date.J = true;
    return J;
  }

  Eigen::MatrixXd & get_jacobian_dot()
  {
    if (method_up_to_date.J) return J_dot;
    frdl::Kinematics::get_jacobian_dot(modules_, state_->q, state_->dq, J, J_dot);
    method_up_to_date.J = true;
    return J_dot;
  }

  Eigen::MatrixXd & get_M()
  {
    if (method_up_to_date.M) return M;
    frdl::Dynamics::get_M(modules_, state_->q, M);
    method_up_to_date.M = true;
    return M;
  }

  Eigen::MatrixXd & get_C()
  {
    if (method_up_to_date.C) return C;
    frdl::Dynamics::get_C(modules_, state_->q, state_->dq, C);
    method_up_to_date.C = true;
    return C;
  }

  Eigen::VectorXd & get_torque_gravity()
  {
    if (method_up_to_date.g) return g;
    frdl::Dynamics::get_torque_gravity(modules_, state_->q, g);
    method_up_to_date.g = true;
    return g;
  }

  Eigen::VectorXd & get_torque_friction(bool const & include_static = true)
  {
    if (method_up_to_date.friction) return friction_torque;
    frdl::Dynamics::get_torque_friction(modules_, state_->dq, friction_torque, include_static);
    method_up_to_date.friction = true;
    return friction_torque;
  }

  Eigen::VectorXd & get_rne()
  {
    if (method_up_to_date.rne) return rne_torque;
    frdl::Dynamics::get_rne(
      modules_, state_->q, state_->dq, zero_vec_, state_->ddq, rne_torque, rne_options);
    method_up_to_date.rne = true;
    return rne_torque;
  }

  int get_size() { return modules_.size(); }

  const ChainModules get_modules() { return modules_; }

  frdl::Kinematics::IKParameters::SharedPtr ik_parameters;
  Dynamics::RNEOptions rne_options;

protected:
  ChainModules modules_;
  JointState::SharedPtr state_;
  CartesianState::SharedPtr cartesian_state_;

  Eigen::MatrixXd J;
  Eigen::MatrixXd J_pinv;
  Eigen::MatrixXd J_dot;
  Eigen::MatrixXd M;
  Eigen::MatrixXd C;
  Eigen::VectorXd g;
  Eigen::VectorXd rne_torque;
  Eigen::VectorXd friction_torque;
  Eigen::VectorXd temp_;
  Eigen::VectorXd zero_vec_;
  Eigen::JacobiSVD<Eigen::MatrixXd> J_svd;

private:
  MethodUpToDate method_up_to_date;
};
}  // namespace frdl