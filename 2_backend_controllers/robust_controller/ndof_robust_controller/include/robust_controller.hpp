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
#include <iostream>
#include <string>
#include <fstream>
#include <thread>
#include <vector>

#include <chrono> // to measure execution time

#include <algorithm>

#include "utils.hpp"

// Include Files for parameter tuning and robust controller
#include "matlab_interface.h"
#include "coder_bounded_array.h"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

class RobustController
{
    // Initialize variables

public:
    void init(std::string const &data_file);

    void update();

    void tuning(coder::array<double, 1U> dqState, double normRho,
                coder::array<double, 1U> tauRobust, double sampleT,
                double dqThreshold, double winSize, double winPerc,
                double PSDMax, double PSDMin, double kpUnique,
                double gainMax, double freqLim, double lambdaGain,
                coder::array<double, 1U> avgPSDdata);

    std::vector<double> get_torque();

    std::vector<double> get_torque_robust();

    std::vector<double> get_gravity();

    double get_controller_gain();

    void set_states(std::vector<double> const & q, std::vector<double> const & dq);

    void set_reference(std::vector<double> const & q_ref, std::vector<double> const & dq_ref, std::vector<double> const & ddq_ref);

    double PSDMax;
    double PSDMin;
    coder::array<double, 1U> avgPSDdata;
    double gainTuneContr_{0.0};
    std::vector<double> tau_;
    std::vector<double> tauRobust_;

private:
    void read_param(std::string const &data_file);

    coder::array<double, 1U> qState;  //{0, 0, 0, 0, 0, 0, 0}};
    coder::array<double, 1U> dqState; //{0, 0, 0, 0, 0, 0, 0}};
    coder::array<double, 1U> qRef;    //{0, 0, 0, 0, 0, 0, 0}};
    coder::array<double, 1U> dqRef;   //{0, 0, 0, 0, 0, 0, 0}};
    coder::array<double, 1U> ddqRef;  //{0, 0, 0, 0, 0, 0, 0}};
    double sampleT{0.001};

    int n_joints;
    double g[3];
    coder::array<double, 2U> B;
    coder::array<double, 2U> P;
    double Kd, Kp;
    coder::array<double, 1U> tau;
    coder::array<double, 1U> tauRobust;
    coder::array<double, 1U> tauGravity;
    double normRho;
    double dqThreshold;
    double winSize;
    double winPerc;
    double kpUnique;
    double gainMax;
    double freqLim;
    double lambdaGain;
    struct1_T KinParExt;
    struct0_T DynPar;
    struct2_T DynPar_inf;
    struct2_T DynPar_sup;


    // for parallel computing
    std::atomic<double> atomic_gainTuneContr{1.0}; // initialized to 1 for safety
    std::atomic<bool> is_thread_running{false};
};
