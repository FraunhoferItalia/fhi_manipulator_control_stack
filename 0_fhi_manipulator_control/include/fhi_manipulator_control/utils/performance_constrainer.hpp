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
#include <atomic>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fhi_robot_dynamic_library/chain.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

namespace fmc
{

class PerformanceConstrainer
{
public:
  using SharedPtr = std::shared_ptr<PerformanceConstrainer>;
  using UniquePtr = std::unique_ptr<PerformanceConstrainer>;

  struct Parameters
  {
    double lambda = 100;
    double w_cr = 0.01;
    double w_th = 0.025;
    double dt = 0.001;
    double k_max = 1000000;
  };

  PerformanceConstrainer(frdl::ChainModules const & modules) : modules_{modules}
  {
    A_.setZero(6);
    V_.setZero(6);
    F_.setZero(6);
    tau_.setZero(modules.size());
    tau_limits_.setZero(modules.size());
    for (size_t i = 0; i < modules.size(); i++) {
      tau_limits_[i] = modules_.active_modules[i].effort_limit;
    }
    J_ = Eigen::MatrixXd::Zero(6, modules.size());
    J_current_ = Eigen::MatrixXd::Zero(6, modules.size());
    J_pinv_ = Eigen::MatrixXd::Zero(modules.size(), 6);
  }

  Eigen::VectorXd & get_performance_contraint_torque(
    Eigen::VectorXd const & q, Eigen::VectorXd const & actual_tau)
  {
    tau_ = J_current_.transpose() * get_performance_contraint_force(q);
    scaling_factor_ =
      ((tau_limits_ - actual_tau.cwiseAbs()).cwiseQuotient(tau_.cwiseAbs())).minCoeff();
    scaling_factor_ = (scaling_factor_ > 0.0) ? scaling_factor_ : 0.0;
    tau_ = (scaling_factor_ > 1.0) ? tau_ : scaling_factor_ * tau_;
    return tau_;
  }

  Eigen::VectorXd & get_performance_contraint_force(Eigen::VectorXd const & q)
  {
    w = evaluate_measure(q);
    frdl::Kinematics::get_jacobian(modules_, q, J_current_);
    F_.setZero();
    if (w > parameters_.w_th) {
      return F_;
    }
    k = (w <= parameters_.w_cr) ? parameters_.k_max
                                : parameters_.lambda * (1 / (w - parameters_.w_cr) -
                                                        1 / (parameters_.w_th - parameters_.w_cr));
    k = (k > parameters_.k_max) ? parameters_.k_max : k;

    J_pinv_ = J_current_.completeOrthogonalDecomposition().pseudoInverse();
    for (int i = 0; i < 6; i++) {
      V_.setZero();
      V_[i] = 1.0;
      q_plus = q + J_pinv_ * parameters_.dt * V_;
      q_minus = q - J_pinv_ * parameters_.dt * V_;
      w_plus = evaluate_measure(q_plus);
      w_minus = evaluate_measure(q_minus);
      if (w_plus <= w && w_minus <= w) {
        A_[i] = 0;
      } else {
        A_[i] = (w_minus > w_plus) ? -(w_minus - w) : (w_plus - w);
      }
    }
    // std::cout << "w: " << w <<std::endl;
    // std::cout << "k: " << k <<std::endl;
    // std::cout << "A_: " << A_.transpose() <<std::endl;

    F_ = k * A_;
    return F_;
  }

  double evaluate_measure(Eigen::VectorXd const & q)
  {
    frdl::Kinematics::get_jacobian(modules_, q, J_);
    // solver.compute(J_.topLeftCorner(3,6)*J_.topLeftCorner(3,6).transpose());
    // solver.compute(J_*J_.transpose());
    J_T_ = J_.topLeftCorner(3, 6);
    J_R_ = J_.bottomRightCorner(3, 6);
    solver.compute(
      J_T_ * (identity6 - (J_R_.transpose() * (J_R_ * J_R_.transpose()).inverse()) * J_R_) *
      J_T_.transpose());
    return solver.eigenvalues().real().minCoeff();
  }

private:
  frdl::ChainModules modules_;
  Eigen::MatrixXd identity6 = Eigen::MatrixXd::Identity(6, 6);
  double k, w, w_plus, w_minus, scaling_factor_;
  Eigen::VectorXd A_, V_, F_, tau_, tau_limits_;
  Eigen::MatrixXd J_, J_pinv_, J_current_, J_R_, J_T_;
  Eigen::VectorXd q_plus, q_minus;
  Eigen::EigenSolver<Eigen::MatrixXd> solver;
  Parameters parameters_ = Parameters();
};
}  // namespace fmc