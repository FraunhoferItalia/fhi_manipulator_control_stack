//
// File: robustControllerURDF_types.h
//
// MATLAB Coder version            : 5.4
// C/C++ source code generated on  : 20-Feb-2023 17:14:04
//

#ifndef ROBUSTCONTROLLERURDF_TYPES_H
#define ROBUSTCONTROLLERURDF_TYPES_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"

// Type Definitions
struct struct1_T {
  coder::array<double, 2U> xyz;
  coder::array<double, 2U> rpy;
};

struct struct2_T {
  coder::array<double, 2U> Jtype;
  coder::array<double, 2U> kr;
  coder::array<double, 2U> Im;
  coder::array<double, 2U> m;
  coder::array<double, 3U> b_I;
  coder::array<double, 2U> Fv;
  coder::array<double, 2U> Fc;
  coder::array<double, 2U> r;
};

struct struct0_T {
  coder::array<double, 2U> m;
  coder::array<double, 2U> r;
  coder::array<double, 2U> kr;
  coder::array<double, 2U> Im;
  coder::array<double, 2U> Fv;
  coder::array<double, 2U> Fc;
  coder::array<double, 2U> Jtype;
  coder::array<double, 3U> b_I;
};

#endif
//
// File trailer for robustControllerURDF_types.h
//
// [EOF]
//
