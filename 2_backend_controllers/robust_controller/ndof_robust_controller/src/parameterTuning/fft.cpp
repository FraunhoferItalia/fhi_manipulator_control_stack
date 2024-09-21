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
#include "fft.h"
#include "FFTImplementationCallback.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const ::coder::array<double, 2U> &x
//                ::coder::array<creal_T, 2U> &y
// Return Type  : void
//
namespace coderTuner {
void fft(const ::coder::array<double, 2U> &x, ::coder::array<creal_T, 2U> &y)
{
  ::coder::array<double, 2U> costab;
  ::coder::array<double, 2U> costab1q;
  ::coder::array<double, 2U> sintab;
  ::coder::array<double, 2U> sintabinv;
  if ((x.size(0) == 0) || (x.size(1) == 0)) {
    int pmax;
    y.set_size(x.size(0), x.size(1));
    pmax = x.size(0) * x.size(1);
    for (int pow2p{0}; pow2p < pmax; pow2p++) {
      y[pow2p].re = 0.0;
      y[pow2p].im = 0.0;
    }
  } else {
    double e;
    int k;
    int n;
    int pmax;
    int pmin;
    int pow2p;
    boolean_T useRadix2;
    useRadix2 = ((x.size(0) & (x.size(0) - 1)) == 0);
    pmin = 1;
    if (useRadix2) {
      pmax = x.size(0);
    } else {
      n = (x.size(0) + x.size(0)) - 1;
      pmax = 31;
      if (n <= 1) {
        pmax = 0;
      } else {
        boolean_T exitg1;
        pmin = 0;
        exitg1 = false;
        while ((!exitg1) && (pmax - pmin > 1)) {
          k = (pmin + pmax) >> 1;
          pow2p = 1 << k;
          if (pow2p == n) {
            pmax = k;
            exitg1 = true;
          } else if (pow2p > n) {
            pmax = k;
          } else {
            pmin = k;
          }
        }
      }
      pmin = 1 << pmax;
      pmax = pmin;
    }
    e = 6.2831853071795862 / static_cast<double>(pmax);
    n = pmax / 2 / 2;
    costab1q.set_size(1, n + 1);
    costab1q[0] = 1.0;
    pmax = n / 2 - 1;
    for (k = 0; k <= pmax; k++) {
      costab1q[k + 1] = std::cos(e * (static_cast<double>(k) + 1.0));
    }
    pow2p = pmax + 2;
    pmax = n - 1;
    for (k = pow2p; k <= pmax; k++) {
      costab1q[k] = std::sin(e * static_cast<double>(n - k));
    }
    costab1q[n] = 0.0;
    if (!useRadix2) {
      n = costab1q.size(1) - 1;
      pmax = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, pmax + 1);
      sintab.set_size(1, pmax + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, pmax + 1);
      for (k = 0; k < n; k++) {
        sintabinv[k + 1] = costab1q[(n - k) - 1];
      }
      pow2p = costab1q.size(1);
      for (k = pow2p; k <= pmax; k++) {
        sintabinv[k] = costab1q[k - n];
      }
      for (k = 0; k < n; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(n - k) - 1];
      }
      pow2p = costab1q.size(1);
      for (k = pow2p; k <= pmax; k++) {
        costab[k] = -costab1q[pmax - k];
        sintab[k] = -costab1q[k - n];
      }
    } else {
      n = costab1q.size(1) - 1;
      pmax = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, pmax + 1);
      sintab.set_size(1, pmax + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      for (k = 0; k < n; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(n - k) - 1];
      }
      pow2p = costab1q.size(1);
      for (k = pow2p; k <= pmax; k++) {
        costab[k] = -costab1q[pmax - k];
        sintab[k] = -costab1q[k - n];
      }
      sintabinv.set_size(1, 0);
    }
    if (useRadix2) {
      coder::internal::fft::FFTImplementationCallback::r2br_r2dit_trig(
          x, x.size(0), costab, sintab, y);
    } else {
      coder::internal::fft::FFTImplementationCallback::dobluesteinfft(
          x, pmin, x.size(0), costab, sintab, sintabinv, y);
    }
  }
}

} // namespace coderTuner

//
// File trailer for fft.cpp
//
// [EOF]
//
