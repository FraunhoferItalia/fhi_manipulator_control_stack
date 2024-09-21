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

// using Duration = std::chrono::duration<double, std::milli>;

void RobustController::read_param(std::string const &data_file)
{
    // read paramters json
    std::ifstream df(data_file.c_str());
    json data = json::parse(df);

    // ------ read number of degrees of freedom robot ------

    read_scalar(n_joints, data, "nDof");

    // ------ set size for outputs ------

    tau.set_size(n_joints);
    tauRobust.set_size(n_joints);
    tauGravity.set_size(n_joints);
    avgPSDdata.set_size(n_joints);

    qState.set_size(n_joints);
    dqState.set_size(n_joints);
    qRef.set_size(n_joints);
    dqRef.set_size(n_joints);
    ddqRef.set_size(n_joints);

    // ------ read double type parameters ------

    read_scalar(Kd, data, "Kd");
    read_scalar(Kp, data, "Kp");
    read_scalar(PSDMax, data, "PSDmax");
    read_scalar(PSDMin, data, "PSDmin");
    read_scalar(freqLim, data, "freqLim");
    read_scalar(gainMax, data, "gainMax");
    read_scalar(kpUnique, data, "kpUnique");
    read_scalar(lambdaGain, data, "lambdaK");
    read_scalar(sampleT, data, "sampleT");
    read_scalar(winSize, data, "winSize");
    read_scalar(winPerc, data, "winPerc");
    read_scalar(dqThreshold, data, "dqThreshold");

    read_double_array(g, data, "g", 3);

    // ------ read matrix ------
    coder::array<double, 2U> jType;
    read_vect(jType, data, "jType", 1, n_joints);
    read_vect(B, data, "B", 2 * n_joints, n_joints);
    read_vect(P, data, "P", 2 * n_joints, 2 * n_joints);

    // ------ read structure ------
    read_struct(KinParExt.xyz, data, "KinPar", "xyz", 3, n_joints + 1);
    read_struct(KinParExt.rpy, data, "KinPar", "rpy", 3, n_joints + 1);

    read_struct(DynPar.Fc, data, "DynPar", "Fc", 1, n_joints);
    read_struct(DynPar.Fv, data, "DynPar", "Fv", 1, n_joints);
    read_struct(DynPar.b_I, data, "DynPar", "I", 3, 3, n_joints + 1);
    read_struct(DynPar.Im, data, "DynPar", "Im", 1, n_joints);
    DynPar.Jtype = jType;
    // read_struct(DynPar.Jtype, data, "DynPar", "Jtype", 1, n_joints);
    read_struct(DynPar.kr, data, "DynPar", "kr", 1, n_joints);
    read_struct(DynPar.m, data, "DynPar", "m", 1, n_joints + 1);
    read_struct(DynPar.r, data, "DynPar", "r", 3, n_joints + 1);

    read_struct(DynPar_inf.Fc, data, "DynPar_inf", "Fc", 1, n_joints);
    read_struct(DynPar_inf.Fv, data, "DynPar_inf", "Fv", 1, n_joints);
    read_struct(DynPar_inf.b_I, data, "DynPar_inf", "I", 3, 3, n_joints + 1);
    read_struct(DynPar_inf.Im, data, "DynPar_inf", "Im", 1, n_joints);
    DynPar_inf.Jtype = jType;
    // read_struct(DynPar_inf.Jtype, data, "DynPar_inf", "Jtype", 1, n_joints);
    read_struct(DynPar_inf.kr, data, "DynPar_inf", "kr", 1, n_joints);
    read_struct(DynPar_inf.m, data, "DynPar_inf", "m", 1, n_joints + 1);
    read_struct(DynPar_inf.r, data, "DynPar_inf", "r", 3, n_joints + 1);

    read_struct(DynPar_sup.Fc, data, "DynPar_sup", "Fc", 1, n_joints);
    read_struct(DynPar_sup.Fv, data, "DynPar_sup", "Fv", 1, n_joints);
    read_struct(DynPar_sup.b_I, data, "DynPar_sup", "I", 3, 3, n_joints + 1);
    read_struct(DynPar_sup.Im, data, "DynPar_sup", "Im", 1, n_joints);
    DynPar_sup.Jtype = jType;
    // read_struct(DynPar_sup.Jtype, data, "DynPar_sup", "Jtype", 1, n_joints);
    read_struct(DynPar_sup.kr, data, "DynPar_sup", "kr", 1, n_joints);
    read_struct(DynPar_sup.m, data, "DynPar_sup", "m", 1, n_joints + 1);
    read_struct(DynPar_sup.r, data, "DynPar_sup", "r", 3, n_joints + 1);
}

std::vector<double> RobustController::get_torque()
{
    return tau;
}

std::vector<double> RobustController::get_torque_robust()
{
    return tauRobust;
}

double RobustController::get_controller_gain()
{
    return atomic_gainTuneContr.load();
}

std::vector<double> RobustController::get_gravity()
{
    return tauGravity;
}

void RobustController::set_states(std::vector<double> const & q, std::vector<double> const & dq)
{
    qState = q;
    dqState = dq;
}

void RobustController::set_reference(std::vector<double> const & q_ref, std::vector<double> const & dq_ref, std::vector<double> const & ddq_ref)
{
    qRef = q_ref;
    dqRef = dq_ref;
    ddqRef = ddq_ref;
}

void RobustController::init(std::string const &data_file)
{
    read_param(data_file);
    tau_.resize(n_joints);
    tauRobust_.resize(n_joints);
}

void RobustController::tuning(coder::array<double, 1U> dqState, double normRho,
                              coder::array<double, 1U> tauRobust, double sampleT,
                              double dqThreshold, double winSize, double winPerc,
                              double PSDMax, double PSDMin, double kpUnique,
                              double gainMax, double freqLim, double lambdaGain,
                              coder::array<double, 1U> avgPSDdata)

{

    is_thread_running.store(true);
    double gainTuneContr = atomic_gainTuneContr.load();
    gainTuneContr_ = gainTuneContr;
    parameterTuning(dqState, normRho, tauRobust, sampleT, dqThreshold,
                    winSize, winPerc, PSDMax, PSDMin,
                    kpUnique, gainMax, freqLim, lambdaGain,
                    &gainTuneContr, avgPSDdata);

    atomic_gainTuneContr.store(gainTuneContr);
    is_thread_running.store(false);
}

void RobustController::update()
{

    // auto const begin_controller = std::chrono::high_resolution_clock::now();
    robustControllerURDF(qRef, dqRef, ddqRef, qState, dqState,
                         g, &DynPar, &KinParExt, Kd, Kp, B,
                         P, &DynPar_inf, &DynPar_sup, sampleT, atomic_gainTuneContr.load(),
                         tau, &normRho, tauRobust, tauGravity);
    for (size_t i=0; i < tau_.size(); i++)
    {
        tau_[i] = tau[i];
        tauRobust_[i] = tauRobust[i];
    }
    // auto const end_controller = std::chrono::high_resolution_clock::now();
    // auto const elapsed_controller = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end_controller - begin_controller);
    //  std::cout << "execution time controller" << elapsed_controller.count() << "ms \n";

    // auto const begin_tuning = std::chrono::high_resolution_clock::now();

    if (is_thread_running.load() == false)
    {

        std::thread t{&RobustController::tuning, this, dqState, normRho, tauRobust, sampleT,
                      dqThreshold, winSize, winPerc, PSDMax, PSDMin, kpUnique,
                      gainMax, freqLim, lambdaGain, avgPSDdata};
        t.detach();
    }

    // auto const end_tuning = std::chrono::high_resolution_clock::now();
    // auto const elapsed_tuning = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end_tuning - begin_tuning);
    // std::cout << "execution time tuning" << elapsed_tuning.count() << "ms \n";
}