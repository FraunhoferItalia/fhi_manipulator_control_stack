//
// File: recursiveNEurdf.h
//
// MATLAB Coder version            : 5.4
// C/C++ source code generated on  : 20-Feb-2023 17:14:04
//

#ifndef RECURSIVENEURDF_H
#define RECURSIVENEURDF_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
void recursiveNEurdf(const coder::array<double, 2U> &q,
                     const coder::array<double, 2U> &dq,
                     const coder::array<double, 2U> &ddq,
                     const coder::array<double, 2U> &KinPar_xyz,
                     const coder::array<double, 2U> &KinPar_rpy,
                     const coder::array<double, 2U> &DynPar_m,
                     const coder::array<double, 2U> &DynPar_r,
                     const coder::array<double, 2U> &DynPar_kr,
                     const coder::array<double, 2U> &DynPar_Im,
                     const coder::array<double, 2U> &DynPar_Fv,
                     const coder::array<double, 2U> &DynPar_Jtype,
                     const coder::array<double, 3U> &DynPar_I,
                     const double g[3], coder::array<double, 2U> &u_out);

#endif
//
// File trailer for recursiveNEurdf.h
//
// [EOF]
//
