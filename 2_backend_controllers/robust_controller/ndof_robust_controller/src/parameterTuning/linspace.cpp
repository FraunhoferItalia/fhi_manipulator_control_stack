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
#include "linspace.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// Arguments    : double d2
//                short n
//                ::coder::array<double, 2U> &y
// Return Type  : void
//
namespace coderTuner {
void linspace(double d2, short n, ::coder::array<double, 2U> &y)
{
  if (n < 0) {
    y.set_size(1, 0);
  } else {
    y.set_size(1, static_cast<int>(n));
    if (n >= 1) {
      int y_tmp;
      y_tmp = n - 1;
      y[n - 1] = d2;
      if (y.size(1) >= 2) {
        y[0] = 0.0;
        if (y.size(1) >= 3) {
          if ((-d2 == 0.0) && (n > 2)) {
            double delta1;
            delta1 = d2 / (static_cast<double>(n) - 1.0);
            for (int k{2}; k <= y_tmp; k++) {
              y[k - 1] = (static_cast<double>((k << 1) - n) - 1.0) * delta1;
            }
            if ((n & 1) == 1) {
              y[n >> 1] = 0.0;
            }
          } else if ((d2 < 0.0) && (std::abs(d2) > 8.9884656743115785E+307)) {
            double delta1;
            delta1 = d2 / (static_cast<double>(y.size(1)) - 1.0);
            y_tmp = y.size(1);
            for (int k{0}; k <= y_tmp - 3; k++) {
              y[k + 1] = delta1 * (static_cast<double>(k) + 1.0);
            }
          } else {
            double delta1;
            delta1 = d2 / (static_cast<double>(y.size(1)) - 1.0);
            y_tmp = y.size(1);
            for (int k{0}; k <= y_tmp - 3; k++) {
              y[k + 1] = (static_cast<double>(k) + 1.0) * delta1;
            }
          }
        }
      }
    }
  }
}

} // namespace coderTuner

//
// File trailer for linspace.cpp
//
// [EOF]
//
