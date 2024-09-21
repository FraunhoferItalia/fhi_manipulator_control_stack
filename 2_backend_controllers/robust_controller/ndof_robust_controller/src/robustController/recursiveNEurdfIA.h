//
// File: recursiveNEurdfIA.h
//
// MATLAB Coder version            : 5.4
// C/C++ source code generated on  : 20-Feb-2023 17:14:04
//

#ifndef RECURSIVENEURDFIA_H
#define RECURSIVENEURDFIA_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct struct2_T;

// Function Declarations
void recursiveNEurdfIA(const coder::array<double, 2U> &q,
                       const coder::array<double, 2U> &dq,
                       const coder::array<double, 2U> &ddq,
                       const coder::array<double, 2U> &KinPar_xyz,
                       const coder::array<double, 2U> &KinPar_rpy,
                       const coder::array<double, 2U> &DynParInf_Jtype,
                       const coder::array<double, 2U> &DynParInf_kr,
                       const coder::array<double, 2U> &DynParInf_Im,
                       const coder::array<double, 2U> &DynParInf_m,
                       const coder::array<double, 3U> &DynParInf_I,
                       const coder::array<double, 2U> &DynParInf_Fv,
                       const coder::array<double, 2U> &DynParInf_r,
                       const struct2_T *DynParSup, const double g[3],
                       coder::array<double, 2U> &u_out_inf,
                       coder::array<double, 2U> &u_out_sup);

#endif
//
// File trailer for recursiveNEurdfIA.h
//
// [EOF]
//
