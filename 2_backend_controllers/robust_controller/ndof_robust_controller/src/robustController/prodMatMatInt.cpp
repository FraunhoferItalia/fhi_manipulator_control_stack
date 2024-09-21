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
#include "prodMatMatInt.h"

// Function Definitions
//
// function [ CInf, CSup ] = prodMatMatInt( AInf, ASup, BInf, BSup )
//
// UNTITLED3 Summary of this function goes here
//    Detailed explanation goes here
//
// Arguments    : const double AInf[9]
//                const double ASup[9]
//                const double BInf[3]
//                const double BSup[3]
//                double CInf[3]
//                double CSup[3]
// Return Type  : void
//
void prodMatMatInt(const double AInf[9], const double ASup[9],
                   const double BInf[3], const double BSup[3], double CInf[3],
                   double CSup[3])
{
  double temp1Inf;
  double temp1Sup;
  // 'prodMatMatInt:5' CInf = zeros(size(AInf,1),size(BInf,2));
  // 'prodMatMatInt:6' CSup = zeros(size(AInf,1),size(BInf,2));
  // 'prodMatMatInt:8' temp1Inf = 0;
  temp1Inf = 0.0;
  //  inizialization for simulink
  // 'prodMatMatInt:9' temp1Sup = 0;
  temp1Sup = 0.0;
  //  inizialization for simulink
  // 'prodMatMatInt:11' for i = 1:size(AInf,1)
  for (int i{0}; i < 3; i++) {
    // 'prodMatMatInt:12' for k = 1:size(BInf,2)
    // 'prodMatMatInt:13' for j = 1:size(AInf,2)
    for (int j{0}; j < 3; j++) {
      // 'prodMatMatInt:14' if j == 1
      if (j + 1 == 1) {
        double b_varargin_1_tmp_tmp;
        double varargin_1_tmp_idx_0;
        double varargin_1_tmp_idx_1;
        double varargin_1_tmp_idx_2;
        // 'prodMatMatInt:15' [temp1Inf, temp1Sup] =
        // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3 Summary
        // of this function goes here
        //    Detailed explanation goes here
        //  if bInf > bSup || aInf > aSup
        //      display('inf value greater than sup')
        //  end
        // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf, aSup*bSup]);
        varargin_1_tmp_idx_0 = AInf[i] * BInf[0];
        varargin_1_tmp_idx_1 = AInf[i] * BSup[0];
        varargin_1_tmp_idx_2 = ASup[i] * BInf[0];
        b_varargin_1_tmp_tmp = ASup[i] * BSup[0];
        temp1Inf = varargin_1_tmp_idx_0;
        // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
        // aSup*bSup]);
        if (varargin_1_tmp_idx_0 > varargin_1_tmp_idx_1) {
          temp1Inf = varargin_1_tmp_idx_1;
        }
        if (varargin_1_tmp_idx_0 < varargin_1_tmp_idx_1) {
          varargin_1_tmp_idx_0 = varargin_1_tmp_idx_1;
        }
        if (temp1Inf > varargin_1_tmp_idx_2) {
          temp1Inf = varargin_1_tmp_idx_2;
        }
        if (varargin_1_tmp_idx_0 < varargin_1_tmp_idx_2) {
          varargin_1_tmp_idx_0 = varargin_1_tmp_idx_2;
        }
        if (temp1Inf > b_varargin_1_tmp_tmp) {
          temp1Inf = b_varargin_1_tmp_tmp;
        }
        if (varargin_1_tmp_idx_0 < b_varargin_1_tmp_tmp) {
          varargin_1_tmp_idx_0 = b_varargin_1_tmp_tmp;
        }
        temp1Sup = varargin_1_tmp_idx_0;
        //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
        //                     display('inf value greater than sup')
        //                 end
      } else {
        double b_varargin_1_tmp_tmp;
        double ex;
        double varargin_1_tmp_idx_0;
        double varargin_1_tmp_idx_1;
        double varargin_1_tmp_idx_2;
        int varargin_1_tmp_tmp;
        // 'prodMatMatInt:19' else
        // 'prodMatMatInt:20' [temp2Inf, temp2Sup] =
        // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3 Summary
        // of this function goes here
        //    Detailed explanation goes here
        //  if bInf > bSup || aInf > aSup
        //      display('inf value greater than sup')
        //  end
        // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf, aSup*bSup]);
        varargin_1_tmp_tmp = i + 3 * j;
        b_varargin_1_tmp_tmp = AInf[varargin_1_tmp_tmp];
        varargin_1_tmp_idx_0 = b_varargin_1_tmp_tmp * BInf[j];
        varargin_1_tmp_idx_1 = b_varargin_1_tmp_tmp * BSup[j];
        b_varargin_1_tmp_tmp = ASup[varargin_1_tmp_tmp];
        varargin_1_tmp_idx_2 = b_varargin_1_tmp_tmp * BInf[j];
        b_varargin_1_tmp_tmp *= BSup[j];
        ex = varargin_1_tmp_idx_0;
        // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
        // aSup*bSup]);
        if (varargin_1_tmp_idx_0 > varargin_1_tmp_idx_1) {
          ex = varargin_1_tmp_idx_1;
        }
        if (varargin_1_tmp_idx_0 < varargin_1_tmp_idx_1) {
          varargin_1_tmp_idx_0 = varargin_1_tmp_idx_1;
        }
        if (ex > varargin_1_tmp_idx_2) {
          ex = varargin_1_tmp_idx_2;
        }
        if (varargin_1_tmp_idx_0 < varargin_1_tmp_idx_2) {
          varargin_1_tmp_idx_0 = varargin_1_tmp_idx_2;
        }
        if (ex > b_varargin_1_tmp_tmp) {
          ex = b_varargin_1_tmp_tmp;
        }
        if (varargin_1_tmp_idx_0 < b_varargin_1_tmp_tmp) {
          varargin_1_tmp_idx_0 = b_varargin_1_tmp_tmp;
        }
        // 'prodMatMatInt:21' [temp1Inf, temp1Sup] =
        // sumInt(temp1Inf,temp1Sup,temp2Inf,temp2Sup); UNTITLED3 Summary of
        // this function goes here
        //    Detailed explanation goes here
        //  if bInf > bSup || aInf > aSup
        //      display('inf value greater than sup')
        //  end
        // 'sumInt:9' cInf = aInf + bInf;
        temp1Inf += ex;
        // 'sumInt:10' cSup = aSup + bSup;
        temp1Sup += varargin_1_tmp_idx_0;
        //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
        //                     display('inf value greater than sup')
        //                 end
      }
    }
    // 'prodMatMatInt:27' CInf(i,k) = temp1Inf;
    CInf[i] = temp1Inf;
    // 'prodMatMatInt:28' CSup(i,k) = temp1Sup;
    CSup[i] = temp1Sup;
  }
}

//
// File trailer for prodMatMatInt.cpp
//
// [EOF]
//
