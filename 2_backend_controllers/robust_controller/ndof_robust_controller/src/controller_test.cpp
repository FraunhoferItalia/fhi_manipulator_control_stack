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
#include "robust_controller.hpp"

int main(int argc, char **argv)
{

    RobustController my_controller;
    // my_controller.init("/home/user/ros1_ws/src/fhi_xbot_plugins/robust_controller_plugin/ndof_robust_controller/data/frankaParamData.json");
    my_controller.init("/home/user/ros1_ws/src/fhi_xbot_plugins/robust_controller_plugin/ndof_robust_controller/data/ModularBot_4DOF_PD.json");
    // my_controller.init("/home/user/MyGitRepos/ndof_robust_controller/data/ModularBot_4DOF_PD.json");

    std::vector<double> q{1, 0.2, 2, 0.6, -1, 0.3, 0.1};
    std::vector<double> dq{0.1, 0.2, 0.3, -0.2, 0.3, 0.5, 0.4};

    std::vector<double> qref{1.05, 0.22, 2.4, 0.58, -1.1, 0.25, 0.11};
    std::vector<double> dqref{0.11, 0.25, 0.35, -0.28, 0.31, 0.55, 0.35};
    std::vector<double> ddqref{0.4, 0.3, -0.2, 0.1, -0.6, 0.2, 0.3};

    my_controller.set_states(q, dq);
    my_controller.set_reference(qref, dqref, ddqref);

    my_controller.update();

    // 21.08 17.2326 20.7409 -25.0036 3.65878 5.20032 3.97796

    std::vector<double> generated_torque;
    generated_torque = my_controller.get_torque();

    std::vector<double> gravity_torque;
    gravity_torque = my_controller.get_gravity();

    std::cout << "gen torque " << generated_torque[0] << " " << generated_torque[1] << " " << generated_torque[2] << " " << generated_torque[3] << " " << generated_torque[4] << " " << generated_torque[5] << " " << generated_torque[6] << " "
              << "\n";

    std::cout << "gravity torque " << gravity_torque[0] << " " << gravity_torque[1] << " " << gravity_torque[2] << " " << gravity_torque[3] << " " << gravity_torque[4] << " " << gravity_torque[5] << " " << gravity_torque[6] << " "
              << "\n";

    return 0;
}

// int main() {
//     return EXIT_SUCCESS;
// }