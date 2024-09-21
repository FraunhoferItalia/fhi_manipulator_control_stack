//
// File: robustControllerURDF.h
//
// MATLAB Coder version            : 5.4
// C/C++ source code generated on  : 20-Feb-2023 17:14:04
//

#ifndef ROBUSTCONTROLLERURDF_H
#define ROBUSTCONTROLLERURDF_H

// Include Files
#include "robustControllerURDF_types.h"
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void robustControllerURDF(
    const coder::array<double, 1U> &qRef, const coder::array<double, 1U> &dqRef,
    const coder::array<double, 1U> &ddqRef, const coder::array<double, 1U> &q,
    const coder::array<double, 1U> &dq, const double g[3],
    const struct0_T *DynPar, const struct1_T *KinPar, double Kd, double Kp,
    const coder::array<double, 2U> &B, const coder::array<double, 2U> &P,
    const struct2_T *DynPar_inf, struct2_T *DynPar_sup, double sampleT,
    double gainTuneContr, coder::array<double, 1U> &tau, double *normRho,
    coder::array<double, 1U> &tauRobust, coder::array<double, 1U> &tauGravity);

#endif
//
// File trailer for robustControllerURDF.h
//
// [EOF]
//
