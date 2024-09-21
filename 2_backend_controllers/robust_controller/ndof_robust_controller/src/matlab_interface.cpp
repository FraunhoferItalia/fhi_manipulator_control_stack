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
#include "matlab_interface.h"
#include "parameterTuning.h"
#include "parameterTuning_terminate.h"
#include "robustControllerURDF.h"
#include "robustControllerURDF_terminate.h"
#include "robustControllerURDF_types.h"
#include "coder_bounded_array.h"

// Function Definitions
// ----------- for both packages -----------

static double argInit_real_T()
{
    return 0.0;
}

static coder::array<double, 1U> argInit_Unboundedx1_real_T()
{
    coder::array<double, 1U> result;
    // Set the size of the array.
    // Change this size to the value that the application requires.
    result.set_size(2);
    // Loop over the array to initialize each element.
    for (int idx0{0}; idx0 < result.size(0); idx0++)
    {
        // Set the value of the array element.
        // Change this value to the value that the application requires.
        result[idx0] = argInit_real_T();
    }
    return result;
}

// ----------- for parameter tuning -----------

static void main_parameterTuning()
{
    coder::array<double, 1U> avgPSDdata;
    coder::array<double, 1U> dq_tmp;
    double gainTuneContr;
    double normRho_tmp;
    // Initialize function 'parameterTuning' input arguments.
    // Initialize function input argument 'dq'.
    dq_tmp = argInit_Unboundedx1_real_T();
    normRho_tmp = argInit_real_T();
    // Initialize function input argument 'robTorque'.
    // Call the entry-point 'parameterTuning'.
    parameterTuning(dq_tmp, normRho_tmp, dq_tmp, normRho_tmp, normRho_tmp,
                    normRho_tmp, normRho_tmp, normRho_tmp, normRho_tmp,
                    normRho_tmp, normRho_tmp, normRho_tmp, normRho_tmp,
                    &gainTuneContr, avgPSDdata);
}

// ----------- for robust controller -----------

static void argInit_1x3_real_T(double result[3])
{
    // Loop over the array to initialize each element.
    for (int idx1{0}; idx1 < 3; idx1++)
    {
        // Set the value of the array element.
        // Change this value to the value that the application requires.
        result[idx1] = argInit_real_T();
    }
}

static coder::array<double, 2U> argInit_1xUnbounded_real_T()
{
    coder::array<double, 2U> result;
    // Set the size of the array.
    // Change this size to the value that the application requires.
    result.set_size(1, 2);
    // Loop over the array to initialize each element.
    for (int idx0{0}; idx0 < 1; idx0++)
    {
        for (int idx1{0}; idx1 < result.size(1); idx1++)
        {
            // Set the value of the array element.
            // Change this value to the value that the application requires.
            result[idx1] = argInit_real_T();
        }
    }
    return result;
}

static coder::array<double, 3U> argInit_3x3xUnbounded_real_T()
{
    coder::array<double, 3U> result;
    // Set the size of the array.
    // Change this size to the value that the application requires.
    result.set_size(3, 3, 2);
    // Loop over the array to initialize each element.
    for (int idx0{0}; idx0 < 3; idx0++)
    {
        for (int idx1{0}; idx1 < 3; idx1++)
        {
            for (int idx2{0}; idx2 < result.size(2); idx2++)
            {
                // Set the value of the array element.
                // Change this value to the value that the application requires.
                result[(idx0 + 3 * idx1) + 9 * idx2] = argInit_real_T();
            }
        }
    }
    return result;
}

static coder::array<double, 2U> argInit_3xUnbounded_real_T()
{
    coder::array<double, 2U> result;
    // Set the size of the array.
    // Change this size to the value that the application requires.
    result.set_size(3, 2);
    // Loop over the array to initialize each element.
    for (int idx0{0}; idx0 < 3; idx0++)
    {
        for (int idx1{0}; idx1 < result.size(1); idx1++)
        {
            // Set the value of the array element.
            // Change this value to the value that the application requires.
            result[idx0 + 3 * idx1] = argInit_real_T();
        }
    }
    return result;
}

static coder::array<double, 2U> argInit_UnboundedxUnbounded_real_T()
{
    coder::array<double, 2U> result;
    // Set the size of the array.
    // Change this size to the value that the application requires.
    result.set_size(2, 2);
    // Loop over the array to initialize each element.
    for (int idx0{0}; idx0 < result.size(0); idx0++)
    {
        for (int idx1{0}; idx1 < result.size(1); idx1++)
        {
            // Set the value of the array element.
            // Change this value to the value that the application requires.
            result[idx0 + result.size(0) * idx1] = argInit_real_T();
        }
    }
    return result;
}

static void argInit_struct0_T(struct0_T *result)
{
    coder::array<double, 2U> result_tmp;
    // Set the value of each structure field.
    // Change this value to the value that the application requires.
    result_tmp = argInit_1xUnbounded_real_T();
    result->m = result_tmp;
    result->r = argInit_3xUnbounded_real_T();
    result->kr = result_tmp;
    result->Im = result_tmp;
    result->Fv = result_tmp;
    result->Fc = result_tmp;
    result->Jtype = result_tmp;
    result->b_I = argInit_3x3xUnbounded_real_T();
}

static struct1_T argInit_struct1_T()
{
    coder::array<double, 2U> result_tmp;
    struct1_T result;
    // Set the value of each structure field.
    // Change this value to the value that the application requires.
    result_tmp = argInit_3xUnbounded_real_T();
    result.xyz = result_tmp;
    result.rpy = result_tmp;
    return result;
}

static void argInit_struct2_T(struct2_T *result)
{
    coder::array<double, 2U> result_tmp;
    // Set the value of each structure field.
    // Change this value to the value that the application requires.
    result_tmp = argInit_1xUnbounded_real_T();
    result->Jtype = result_tmp;
    result->kr = result_tmp;
    result->Im = result_tmp;
    result->m = result_tmp;
    result->b_I = argInit_3x3xUnbounded_real_T();
    result->Fv = result_tmp;
    result->Fc = result_tmp;
    result->r = argInit_3xUnbounded_real_T();
}

static void main_robustControllerURDF()
{
    coder::array<double, 2U> B_tmp;
    coder::array<double, 1U> qRef_tmp;
    coder::array<double, 1U> tau;
    coder::array<double, 1U> tauGravity;
    coder::array<double, 1U> tauRobust;
    struct0_T DynPar;
    struct1_T KinPar;
    struct2_T DynPar_inf_tmp;
    struct2_T DynPar_sup;
    double dv[3];
    double Kd_tmp;
    double normRho;
    // Initialize function 'robustControllerURDF' input arguments.
    // Initialize function input argument 'qRef'.
    qRef_tmp = argInit_Unboundedx1_real_T();
    // Initialize function input argument 'dqRef'.
    // Initialize function input argument 'ddqRef'.
    // Initialize function input argument 'q'.
    // Initialize function input argument 'dq'.
    // Initialize function input argument 'g'.
    // Initialize function input argument 'DynPar'.
    argInit_struct0_T(&DynPar);
    // Initialize function input argument 'KinPar'.
    KinPar = argInit_struct1_T();
    Kd_tmp = argInit_real_T();
    // Initialize function input argument 'B'.
    B_tmp = argInit_UnboundedxUnbounded_real_T();
    // Initialize function input argument 'P'.
    // Initialize function input argument 'DynPar_inf'.
    argInit_struct2_T(&DynPar_inf_tmp);
    // Initialize function input argument 'DynPar_sup'.
    DynPar_sup = DynPar_inf_tmp;
    // Call the entry-point 'robustControllerURDF'.
    argInit_1x3_real_T(dv);
    robustControllerURDF(qRef_tmp, qRef_tmp, qRef_tmp, qRef_tmp, qRef_tmp, dv,
                         &DynPar, &KinPar, Kd_tmp, Kd_tmp, B_tmp, B_tmp,
                         &DynPar_inf_tmp, &DynPar_sup, Kd_tmp, Kd_tmp, tau,
                         &normRho, tauRobust, tauGravity);
}
