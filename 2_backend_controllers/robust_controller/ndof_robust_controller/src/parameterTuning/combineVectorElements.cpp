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
#include "combineVectorElements.h"
#include "coder_array.h"

// Function Definitions
//
// Arguments    : const ::coder::array<double, 2U> &x
//                ::coder::array<double, 2U> &y
// Return Type  : void
//
namespace coderTuner {
void combineVectorElements(const ::coder::array<double, 2U> &x,
                           ::coder::array<double, 2U> &y)
{
  if ((x.size(0) == 0) || (x.size(1) == 0)) {
    int npages;
    y.set_size(1, x.size(1));
    npages = x.size(1);
    for (int firstBlockLength{0}; firstBlockLength < npages;
         firstBlockLength++) {
      y[firstBlockLength] = 0.0;
    }
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    int npages;
    npages = x.size(1);
    y.set_size(1, x.size(1));
    if (x.size(0) <= 1024) {
      firstBlockLength = x.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x.size(0) / 1024;
      lastBlockLength = x.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    for (int xi{0}; xi < npages; xi++) {
      int xpageoffset;
      xpageoffset = xi * x.size(0);
      y[xi] = x[xpageoffset];
      for (int k{2}; k <= firstBlockLength; k++) {
        y[xi] = y[xi] + x[(xpageoffset + k) - 1];
      }
      for (int ib{2}; ib <= nblocks; ib++) {
        double bsum;
        int hi;
        int xblockoffset;
        xblockoffset = xpageoffset + ((ib - 1) << 10);
        bsum = x[xblockoffset];
        if (ib == nblocks) {
          hi = lastBlockLength;
        } else {
          hi = 1024;
        }
        for (int k{2}; k <= hi; k++) {
          bsum += x[(xblockoffset + k) - 1];
        }
        y[xi] = y[xi] + bsum;
      }
    }
  }
}

} // namespace coderTuner

//
// File trailer for combineVectorElements.cpp
//
// [EOF]
//
