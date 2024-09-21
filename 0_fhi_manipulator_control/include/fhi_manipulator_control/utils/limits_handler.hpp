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
class LimitsHandler
{
public:
  using SharedPtr = std::shared_ptr<LimitsHandler>;
  using UniquePtr = std::unique_ptr<LimitsHandler>;

  LimitsHandler(frdl::ChainModules const & modules, double const & threshold)
  : modules_{modules}, threshold_{threshold}
  {
    tau_joint_limits_.resize(modules_.size());
    K_.resize(modules_.size());
    D_ = Eigen::VectorXd::Constant(modules_.size(), 1.0);
    // std::cout << "K: " << K_.transpose() << std::endl;
    // std::cout << "D: " << D_.transpose() << std::endl;
  }

  Eigen::VectorXd & get_tau_limits(
    Eigen::VectorXd const & q, Eigen::VectorXd const & dq, Eigen::VectorXd const & actual_tau)
  {
    for (int i = 0; i < tau_joint_limits_.size(); i++) {
      K_[i] = (std::abs(actual_tau[i]) < modules_.active_modules[i].effort_limit)
                ? (modules_.active_modules[i].effort_limit - std::abs(actual_tau[i])) /
                    (1.01 * threshold_)
                : 0.0;
    }
    tau_joint_limits_.setZero();
    for (int i = 0; i < tau_joint_limits_.size(); i++) {
      delta_min_ = std::abs(modules_.active_modules[i].position_limits.x() - q[i]);
      delta_max_ = std::abs(modules_.active_modules[i].position_limits.y() - q[i]);
      if (delta_min_ < threshold_) {
        tau_joint_limits_[i] = K_[i] * (threshold_ - delta_min_) - D_[i] * dq[i];
      } else if (delta_max_ < threshold_) {
        tau_joint_limits_[i] = -K_[i] * (threshold_ - delta_max_) - D_[i] * dq[i];
      }
    }
    // std::cout << tau_joint_limits_.transpose() << std::endl;
    return tau_joint_limits_;
  }

protected:
  frdl::ChainModules modules_;
  double threshold_;
  double delta_min_, delta_max_;
  Eigen::VectorXd K_;
  Eigen::VectorXd D_;
  Eigen::VectorXd tau_joint_limits_;
};
}  // namespace fmc
