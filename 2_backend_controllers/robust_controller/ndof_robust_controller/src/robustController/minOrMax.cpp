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

// Include Files
#include "minOrMax.h"

// Function Definitions
//
// Arguments    : const double x[2]
// Return Type  : double
//
namespace control_coder {
namespace coder {
namespace internal {
double b_maximum(const double x[2])
{
  double ex;
  if (x[0] < x[1]) {
    ex = x[1];
  } else {
    ex = x[0];
  }
  return ex;
}

//
// Arguments    : const double x[2]
// Return Type  : double
//
double b_minimum(const double x[2])
{
  double ex;
  if (x[0] > x[1]) {
    ex = x[1];
  } else {
    ex = x[0];
  }
  return ex;
}

//
// Arguments    : const double x[4]
// Return Type  : double
//
double maximum(const double x[4])
{
  double ex;
  ex = x[0];
  if (x[0] < x[1]) {
    ex = x[1];
  }
  if (ex < x[2]) {
    ex = x[2];
  }
  if (ex < x[3]) {
    ex = x[3];
  }
  return ex;
}

//
// Arguments    : const double x[4]
// Return Type  : double
//
double minimum(const double x[4])
{
  double ex;
  ex = x[0];
  if (x[0] > x[1]) {
    ex = x[1];
  }
  if (ex > x[2]) {
    ex = x[2];
  }
  if (ex > x[3]) {
    ex = x[3];
  }
  return ex;
}

} // namespace internal
} // namespace coder
} // namespace control_coder

//
// File trailer for minOrMax.cpp
//
// [EOF]
//
