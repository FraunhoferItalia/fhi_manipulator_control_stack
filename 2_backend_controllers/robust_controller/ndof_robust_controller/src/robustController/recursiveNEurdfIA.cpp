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
#include "recursiveNEurdfIA.h"
#include "minOrMax.h"
#include "prodMatMatInt.h"
#include "robustControllerURDF_data.h"
#include "robustControllerURDF_types.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// function [ u_out_inf, u_out_sup ] = recursiveNEurdfIA(
// q,dq,ddq,KinPar,DynParInf, DynParSup,g )
//
// UNTITLED18 Summary of this function goes here
//    Detailed explanation goes here
//
// Arguments    : const coder::array<double, 2U> &q
//                const coder::array<double, 2U> &dq
//                const coder::array<double, 2U> &ddq
//                const coder::array<double, 2U> &KinPar_xyz
//                const coder::array<double, 2U> &KinPar_rpy
//                const coder::array<double, 2U> &DynParInf_Jtype
//                const coder::array<double, 2U> &DynParInf_kr
//                const coder::array<double, 2U> &DynParInf_Im
//                const coder::array<double, 2U> &DynParInf_m
//                const coder::array<double, 3U> &DynParInf_I
//                const coder::array<double, 2U> &DynParInf_Fv
//                const coder::array<double, 2U> &DynParInf_r
//                const struct2_T *DynParSup
//                const double g[3]
//                coder::array<double, 2U> &u_out_inf
//                coder::array<double, 2U> &u_out_sup
// Return Type  : void
//
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
                       coder::array<double, 2U> &u_out_sup)
{
  static const signed char b[3]{0, 0, 1};
  static const signed char b_iv[3]{1, 0, 0};
  coder::array<double, 3U> R;
  coder::array<double, 3U> R_T;
  coder::array<double, 2U> a_c_inf;
  coder::array<double, 2U> a_c_sup;
  coder::array<double, 2U> a_j;
  coder::array<double, 2U> f_inf;
  coder::array<double, 2U> f_sup;
  coder::array<double, 2U> mu_inf;
  coder::array<double, 2U> mu_sup;
  coder::array<double, 2U> omega;
  coder::array<double, 2U> omegaDot;
  coder::array<double, 2U> r_link;
  coder::array<double, 2U> u_inf;
  coder::array<double, 2U> u_sup;
  double c_Ry_tmp[9];
  double c_Rz_tmp[9];
  double dv[9];
  double e_Rz_tmp[9];
  double xw_I_sup[9];
  double Fi_inf[3];
  double Fi_sup[3];
  double Rf_inf[3];
  double Rf_sup[3];
  double Rx_tmp;
  double Ry_tmp;
  double Rz_tmp;
  double b_Ry_tmp;
  double b_Rz_tmp;
  double d;
  double d1;
  double d_Rz_tmp;
  double temp1Inf;
  double temp1Sup;
  int b_i;
  int b_k;
  int d_i;
  int exitg2;
  int i;
  unsigned int k;
  int loop_ub;
  //  recursive Newton-Euler algorithm for URDF convention with interval
  //  arithmetic
  // 'recursiveNEurdfIA:8' jointType = DynParInf.Jtype;
  // 'recursiveNEurdfIA:9' m_inf = DynParInf.m;
  // 'recursiveNEurdfIA:10' I_inf = DynParInf.I;
  // 'recursiveNEurdfIA:11' r_inf = DynParInf.r;
  // 'recursiveNEurdfIA:13' m_sup = DynParSup.m;
  // 'recursiveNEurdfIA:14' I_sup = DynParSup.I;
  // 'recursiveNEurdfIA:15' r_sup = DynParSup.r;
  // 'recursiveNEurdfIA:17' k = length(jointType) + 1;
  k = static_cast<unsigned int>(DynParInf_Jtype.size(1));
  //  number of arms = number of joints + 1
  // 'recursiveNEurdfIA:19' xyz = KinPar.xyz;
  //  transposed
  // 'recursiveNEurdfIA:20' rpy = KinPar.rpy;
  //  transposed
  // 'recursiveNEurdfIA:22' kr = DynParInf.kr;
  // 'recursiveNEurdfIA:23' Im = DynParInf.Im;
  // 'recursiveNEurdfIA:24' Fv_inf = DynParInf.Fv;
  // 'recursiveNEurdfIA:25' Fc_inf = DynParInf.Fc;
  // 'recursiveNEurdfIA:27' Fv_sup = DynParSup.Fv;
  // 'recursiveNEurdfIA:28' Fc_sup = DynParSup.Fc;
  //  initialize parameters
  // 'recursiveNEurdfIA:32' omega       = zeros(3,k);
  omega.set_size(3, DynParInf_Jtype.size(1) + 1);
  loop_ub = 3 * (DynParInf_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    omega[i] = 0.0;
  }
  // 'recursiveNEurdfIA:33' omegaDot    = zeros(3,k);
  omegaDot.set_size(3, DynParInf_Jtype.size(1) + 1);
  loop_ub = 3 * (DynParInf_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    omegaDot[i] = 0.0;
  }
  // 'recursiveNEurdfIA:34' a_j         = zeros(3,k);
  a_j.set_size(3, DynParInf_Jtype.size(1) + 1);
  loop_ub = 3 * (DynParInf_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    a_j[i] = 0.0;
  }
  // 'recursiveNEurdfIA:35' a_c_inf     = zeros(3,k);
  a_c_inf.set_size(3, DynParInf_Jtype.size(1) + 1);
  loop_ub = 3 * (DynParInf_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    a_c_inf[i] = 0.0;
  }
  // 'recursiveNEurdfIA:36' a_c_sup     = zeros(3,k);
  a_c_sup.set_size(3, DynParInf_Jtype.size(1) + 1);
  loop_ub = 3 * (DynParInf_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    a_c_sup[i] = 0.0;
  }
  // 'recursiveNEurdfIA:37' R           = zeros(3,3,k+1);
  R.set_size(3, 3,
             static_cast<int>(
                 (static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) + 1.0));
  loop_ub = 9 * static_cast<int>(
                    (static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    R[i] = 0.0;
  }
  //  from i to i-1
  // 'recursiveNEurdfIA:38' R_T         = zeros(3,3,k);
  R_T.set_size(3, 3, DynParInf_Jtype.size(1) + 1);
  loop_ub = 9 * (DynParInf_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    R_T[i] = 0.0;
  }
  //  from i-1 to i
  // 'recursiveNEurdfIA:39' z           = zeros(3,k);
  // 'recursiveNEurdfIA:40' u_inf       = zeros(1,k);
  u_inf.set_size(1, DynParInf_Jtype.size(1) + 1);
  loop_ub = DynParInf_Jtype.size(1) + 1;
  for (i = 0; i < loop_ub; i++) {
    u_inf[i] = 0.0;
  }
  // 'recursiveNEurdfIA:41' u_sup       = zeros(1,k);
  u_sup.set_size(1, DynParInf_Jtype.size(1) + 1);
  loop_ub = DynParInf_Jtype.size(1) + 1;
  for (i = 0; i < loop_ub; i++) {
    u_sup[i] = 0.0;
  }
  // 'recursiveNEurdfIA:42' z0          = [0;0;1];
  // 'recursiveNEurdfIA:43' f_inf       = zeros(3,k+1);
  f_inf.set_size(
      3, static_cast<int>((static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) +
                          1.0));
  loop_ub = 3 * static_cast<int>(
                    (static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    f_inf[i] = 0.0;
  }
  //  forces (k+1 for end effector forces)
  // 'recursiveNEurdfIA:44' f_sup       = zeros(3,k+1);
  f_sup.set_size(
      3, static_cast<int>((static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) +
                          1.0));
  loop_ub = 3 * static_cast<int>(
                    (static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    f_sup[i] = 0.0;
  }
  //  forces (k+1 for end effector forces)
  // 'recursiveNEurdfIA:45' mu_inf      = zeros(3,k+1);
  mu_inf.set_size(
      3, static_cast<int>((static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) +
                          1.0));
  loop_ub = 3 * static_cast<int>(
                    (static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    mu_inf[i] = 0.0;
  }
  //  torques (k+1 for end effector torques)
  // 'recursiveNEurdfIA:46' mu_sup      = zeros(3,k+1);
  mu_sup.set_size(
      3, static_cast<int>((static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) +
                          1.0));
  loop_ub = 3 * static_cast<int>(
                    (static_cast<double>(DynParInf_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    mu_sup[i] = 0.0;
  }
  //  torques (k+1 for end effector torques)
  // 'recursiveNEurdfIA:48' r_link = zeros(3,k);
  r_link.set_size(3, DynParInf_Jtype.size(1) + 1);
  // 'recursiveNEurdfIA:50' for i = 1:k
  i = DynParInf_Jtype.size(1);
  for (b_i = 0; b_i <= i; b_i++) {
    // 'recursiveNEurdfIA:51' if i < k
    if (static_cast<unsigned int>(b_i) < k) {
      double f_Rz_tmp[16];
      //  all but last link
      // 'recursiveNEurdfIA:53' RzCurrent = RzMat(q(i));
      // UNTITLED3 Summary of this function goes here
      //    Detailed explanation goes here
      // 'RzMat:5' Rz = [cos(x), -sin(x), 0 ,0;
      // 'RzMat:6'       sin(x), cos(x), 0, 0;
      // 'RzMat:7'       0, 0, 1, 0;
      // 'RzMat:8'       0, 0, 0, 1];
      temp1Sup = std::sin(q[b_i]);
      Rz_tmp = std::cos(q[b_i]);
      // 'recursiveNEurdfIA:54' R(:,:,i) =
      // rotMatURDF(rpy(:,i))*RzCurrent(1:3,1:3); UNTITLED2 Summary of this
      // function goes here
      //    Detailed explanation goes here
      // 'rotMatURDF:5' r = angles(1);
      //  roll
      // 'rotMatURDF:6' p = angles(2);
      //  pitch
      // 'rotMatURDF:7' y = angles(3);
      //  yaw
      // 'rotMatURDF:9' Rx = [1, 0, 0; 0, cos(r), -sin(r); 0, sin(r), cos(r)];
      d = KinPar_rpy[3 * b_i];
      Rx_tmp = std::sin(d);
      temp1Inf = std::cos(d);
      // 'rotMatURDF:10' Ry = [cos(p), 0, -sin(p); 0, 1, 0; sin(p), 0, cos(p)];
      d = KinPar_rpy[3 * b_i + 1];
      Ry_tmp = std::sin(d);
      b_Ry_tmp = std::cos(d);
      // 'rotMatURDF:11' Rz = [cos(y), -sin(y), 0; sin(y), cos(y), 0; 0, 0, 1];
      d = KinPar_rpy[3 * b_i + 2];
      d_Rz_tmp = std::sin(d);
      b_Rz_tmp = std::cos(d);
      // 'rotMatURDF:13' R = Rz * Ry * Rx;
      c_Rz_tmp[0] = b_Rz_tmp;
      c_Rz_tmp[3] = -d_Rz_tmp;
      c_Rz_tmp[6] = 0.0;
      c_Rz_tmp[1] = d_Rz_tmp;
      c_Rz_tmp[4] = b_Rz_tmp;
      c_Rz_tmp[7] = 0.0;
      c_Ry_tmp[0] = b_Ry_tmp;
      c_Ry_tmp[3] = 0.0;
      c_Ry_tmp[6] = -Ry_tmp;
      c_Rz_tmp[2] = 0.0;
      c_Ry_tmp[1] = 0.0;
      c_Rz_tmp[5] = 0.0;
      c_Ry_tmp[4] = 1.0;
      c_Rz_tmp[8] = 1.0;
      c_Ry_tmp[7] = 0.0;
      c_Ry_tmp[2] = Ry_tmp;
      c_Ry_tmp[5] = 0.0;
      c_Ry_tmp[8] = b_Ry_tmp;
      for (b_k = 0; b_k < 3; b_k++) {
        d = c_Rz_tmp[b_k];
        d1 = c_Rz_tmp[b_k + 3];
        loop_ub = static_cast<int>(c_Rz_tmp[b_k + 6]);
        for (d_i = 0; d_i < 3; d_i++) {
          e_Rz_tmp[b_k + 3 * d_i] =
              (d * c_Ry_tmp[3 * d_i] + d1 * c_Ry_tmp[3 * d_i + 1]) +
              static_cast<double>(loop_ub) * c_Ry_tmp[3 * d_i + 2];
        }
        dv[3 * b_k] = b_iv[b_k];
      }
      dv[1] = 0.0;
      dv[4] = temp1Inf;
      dv[7] = -Rx_tmp;
      dv[2] = 0.0;
      dv[5] = Rx_tmp;
      dv[8] = temp1Inf;
      f_Rz_tmp[0] = Rz_tmp;
      f_Rz_tmp[4] = -temp1Sup;
      f_Rz_tmp[8] = 0.0;
      f_Rz_tmp[12] = 0.0;
      f_Rz_tmp[1] = temp1Sup;
      f_Rz_tmp[5] = Rz_tmp;
      f_Rz_tmp[9] = 0.0;
      f_Rz_tmp[13] = 0.0;
      f_Rz_tmp[2] = 0.0;
      f_Rz_tmp[3] = 0.0;
      f_Rz_tmp[6] = 0.0;
      f_Rz_tmp[7] = 0.0;
      f_Rz_tmp[10] = 1.0;
      f_Rz_tmp[11] = 0.0;
      f_Rz_tmp[14] = 0.0;
      f_Rz_tmp[15] = 1.0;
      // 'recursiveNEurdfIA:55' R_T(:,:,i) = R(:,:,i)';
      // 'recursiveNEurdfIA:57' z(:,i) = R_T(:,:,i)*z0;
      // 'recursiveNEurdfIA:59' r_link(:,i) = xyz(:,i);
      for (b_k = 0; b_k < 3; b_k++) {
        d = e_Rz_tmp[b_k];
        d1 = e_Rz_tmp[b_k + 3];
        b_Rz_tmp = e_Rz_tmp[b_k + 6];
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          c_Rz_tmp[b_k + 3 * loop_ub] =
              (d * dv[3 * loop_ub] + d1 * dv[3 * loop_ub + 1]) +
              b_Rz_tmp * dv[3 * loop_ub + 2];
        }
        d = c_Rz_tmp[b_k];
        d1 = c_Rz_tmp[b_k + 3];
        b_Rz_tmp = c_Rz_tmp[b_k + 6];
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          d_i = loop_ub << 2;
          Ry_tmp = (d * f_Rz_tmp[d_i] + d1 * f_Rz_tmp[d_i + 1]) +
                   b_Rz_tmp * f_Rz_tmp[d_i + 2];
          R[(b_k + 3 * loop_ub) + 9 * b_i] = Ry_tmp;
          R_T[(loop_ub + 3 * b_k) + 9 * b_i] = Ry_tmp;
        }
        r_link[b_k + 3 * b_i] = KinPar_xyz[b_k + 3 * b_i];
      }
    } else {
      // 'recursiveNEurdfIA:61' else
      //  last link
      // 'recursiveNEurdfIA:62' R(:,:,i) = rotMatURDF(rpy(:,i));
      // UNTITLED2 Summary of this function goes here
      //    Detailed explanation goes here
      // 'rotMatURDF:5' r = angles(1);
      //  roll
      // 'rotMatURDF:6' p = angles(2);
      //  pitch
      // 'rotMatURDF:7' y = angles(3);
      //  yaw
      // 'rotMatURDF:9' Rx = [1, 0, 0; 0, cos(r), -sin(r); 0, sin(r), cos(r)];
      d = KinPar_rpy[3 * b_i];
      Rx_tmp = std::sin(d);
      temp1Inf = std::cos(d);
      // 'rotMatURDF:10' Ry = [cos(p), 0, -sin(p); 0, 1, 0; sin(p), 0, cos(p)];
      d = KinPar_rpy[3 * b_i + 1];
      Ry_tmp = std::sin(d);
      b_Ry_tmp = std::cos(d);
      // 'rotMatURDF:11' Rz = [cos(y), -sin(y), 0; sin(y), cos(y), 0; 0, 0, 1];
      d = KinPar_rpy[3 * b_i + 2];
      temp1Sup = std::sin(d);
      Rz_tmp = std::cos(d);
      // 'rotMatURDF:13' R = Rz * Ry * Rx;
      c_Rz_tmp[0] = Rz_tmp;
      c_Rz_tmp[3] = -temp1Sup;
      c_Rz_tmp[6] = 0.0;
      c_Rz_tmp[1] = temp1Sup;
      c_Rz_tmp[4] = Rz_tmp;
      c_Rz_tmp[7] = 0.0;
      c_Ry_tmp[0] = b_Ry_tmp;
      c_Ry_tmp[3] = 0.0;
      c_Ry_tmp[6] = -Ry_tmp;
      c_Rz_tmp[2] = 0.0;
      c_Ry_tmp[1] = 0.0;
      c_Rz_tmp[5] = 0.0;
      c_Ry_tmp[4] = 1.0;
      c_Rz_tmp[8] = 1.0;
      c_Ry_tmp[7] = 0.0;
      c_Ry_tmp[2] = Ry_tmp;
      c_Ry_tmp[5] = 0.0;
      c_Ry_tmp[8] = b_Ry_tmp;
      for (b_k = 0; b_k < 3; b_k++) {
        d = c_Rz_tmp[b_k];
        d1 = c_Rz_tmp[b_k + 3];
        loop_ub = static_cast<int>(c_Rz_tmp[b_k + 6]);
        for (d_i = 0; d_i < 3; d_i++) {
          e_Rz_tmp[b_k + 3 * d_i] =
              (d * c_Ry_tmp[3 * d_i] + d1 * c_Ry_tmp[3 * d_i + 1]) +
              static_cast<double>(loop_ub) * c_Ry_tmp[3 * d_i + 2];
        }
        dv[3 * b_k] = b_iv[b_k];
      }
      dv[1] = 0.0;
      dv[4] = temp1Inf;
      dv[7] = -Rx_tmp;
      dv[2] = 0.0;
      dv[5] = Rx_tmp;
      dv[8] = temp1Inf;
      // 'recursiveNEurdfIA:63' R_T(:,:,i) = R(:,:,i)';
      // 'recursiveNEurdfIA:65' z(:,i) = R_T(:,:,i)*z0;
      // 'recursiveNEurdfIA:67' r_link(:,i) = xyz(:,i);
      for (b_k = 0; b_k < 3; b_k++) {
        d = e_Rz_tmp[b_k];
        d1 = e_Rz_tmp[b_k + 3];
        b_Rz_tmp = e_Rz_tmp[b_k + 6];
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          Ry_tmp = (d * dv[3 * loop_ub] + d1 * dv[3 * loop_ub + 1]) +
                   b_Rz_tmp * dv[3 * loop_ub + 2];
          R[(b_k + 3 * loop_ub) + 9 * b_i] = Ry_tmp;
          R_T[(loop_ub + 3 * b_k) + 9 * b_i] = Ry_tmp;
        }
        r_link[b_k + 3 * b_i] = KinPar_xyz[b_k + 3 * b_i];
      }
    }
    // 'recursiveNEurdfIA:70' z(:,i) = R_T(:,:,i)*z0;
  }
  // 'recursiveNEurdfIA:74' B = eye(3);
  //  forward recursion
  // 'recursiveNEurdfIA:78' for i = 1:k
  b_i = 0;
  do {
    exitg2 = 0;
    if (b_i <= static_cast<int>(k)) {
      boolean_T guard1{false};
      //  initialization
      // 'recursiveNEurdfIA:80' if i == 1
      guard1 = false;
      if (static_cast<double>(b_i) + 1.0 == 1.0) {
        // 'recursiveNEurdfIA:81' z(:,i) = z0;
        // 'recursiveNEurdfIA:82' a_j(:,i) = -B'*[g(1);g(2);g(3)];
        d = g[0];
        d1 = g[1];
        b_Rz_tmp = g[2];
        for (i = 0; i < 3; i++) {
          a_j[i] = (static_cast<double>(iv[i]) * d +
                    static_cast<double>(iv[i + 3]) * d1) +
                   static_cast<double>(iv[i + 6]) * b_Rz_tmp;
        }
        guard1 = true;
      } else {
        // 'recursiveNEurdfIA:84' else
        // 'recursiveNEurdfIA:85' if jointType(i-1) == 1
        d = DynParInf_Jtype[static_cast<int>((static_cast<double>(b_i) + 1.0) -
                                             1.0) -
                            1];
        if (d == 1.0) {
          //  revolute joint
          // 'recursiveNEurdfIA:86' omega(:,i) = R_T(:,:,i-1)*(omega(:,i-1)) +
          // dq(i-1)*z0;
          b_Rz_tmp = dq[b_i - 1];
          for (i = 0; i < 3; i++) {
            Fi_inf[i] =
                ((R_T[i + 9 * (b_i - 1)] * omega[3 * (b_i - 1)] +
                  R_T[(i + 9 * (b_i - 1)) + 3] * omega[3 * (b_i - 1) + 1]) +
                 R_T[(i + 9 * (b_i - 1)) + 6] * omega[3 * (b_i - 1) + 2]) +
                b_Rz_tmp * static_cast<double>(b[i]);
          }
          omega[3 * b_i] = Fi_inf[0];
          omega[3 * b_i + 1] = Fi_inf[1];
          omega[3 * b_i + 2] = Fi_inf[2];
          //  dq(i-1) because it starts from i = 1
          // 'recursiveNEurdfIA:87' omegaDot(:,i) =
          // R_T(:,:,i-1)*(omegaDot(:,i-1)) +
          // Xcross(R_T(:,:,i-1)*omega(:,i-1))*dq(i-1)*z0 + ddq(i-1)*z0;
          for (i = 0; i < 3; i++) {
            Fi_inf[i] =
                (R_T[i + 9 * (b_i - 1)] * omega[3 * (b_i - 1)] +
                 R_T[(i + 9 * (b_i - 1)) + 3] * omega[3 * (b_i - 1) + 1]) +
                R_T[(i + 9 * (b_i - 1)) + 6] * omega[3 * (b_i - 1) + 2];
          }
          // UNTITLED7 Summary of this function goes here
          //    Detailed explanation goes here
          // 'Xcross:5' XMat = [0 -x(3) x(2);
          // 'Xcross:6'         x(3) 0 -x(1);
          // 'Xcross:7'         -x(2) x(1) 0];
          dv[0] = 0.0;
          dv[3] = -Fi_inf[2] * b_Rz_tmp;
          dv[6] = Fi_inf[1] * b_Rz_tmp;
          dv[1] = Fi_inf[2] * b_Rz_tmp;
          dv[4] = 0.0;
          dv[7] = -Fi_inf[0] * b_Rz_tmp;
          dv[2] = -Fi_inf[1] * b_Rz_tmp;
          dv[5] = Fi_inf[0] * b_Rz_tmp;
          dv[8] = 0.0;
          for (i = 0; i < 3; i++) {
            Fi_sup[i] = dv[i + 6];
            Fi_inf[i] =
                (R_T[i + 9 * (b_i - 1)] * omegaDot[3 * (b_i - 1)] +
                 R_T[(i + 9 * (b_i - 1)) + 3] * omegaDot[3 * (b_i - 1) + 1]) +
                R_T[(i + 9 * (b_i - 1)) + 6] * omegaDot[3 * (b_i - 1) + 2];
          }
          omegaDot[3 * b_i] = Fi_inf[0] + Fi_sup[0];
          omegaDot[3 * b_i + 1] = Fi_inf[1] + Fi_sup[1];
          omegaDot[3 * b_i + 2] = (Fi_inf[2] + Fi_sup[2]) + ddq[b_i - 1];
          // 'recursiveNEurdfIA:88' a_j(:,i) = R_T(:,:,i-1)*(a_j(:,i-1) +
          // Xcross(omegaDot(:,i-1))*r_link(:,i-1) +
          // Xcross(omega(:,i-1))*Xcross(omega(:,i-1))*r_link(:,i-1)); UNTITLED7
          // Summary of this function goes here
          //    Detailed explanation goes here
          // 'Xcross:5' XMat = [0 -x(3) x(2);
          // 'Xcross:6'         x(3) 0 -x(1);
          // 'Xcross:7'         -x(2) x(1) 0];
          // UNTITLED7 Summary of this function goes here
          //    Detailed explanation goes here
          // 'Xcross:5' XMat = [0 -x(3) x(2);
          // 'Xcross:6'         x(3) 0 -x(1);
          // 'Xcross:7'         -x(2) x(1) 0];
          // UNTITLED7 Summary of this function goes here
          //    Detailed explanation goes here
          // 'Xcross:5' XMat = [0 -x(3) x(2);
          // 'Xcross:6'         x(3) 0 -x(1);
          // 'Xcross:7'         -x(2) x(1) 0];
          dv[0] = 0.0;
          d = omegaDot[3 * (b_i - 1) + 2];
          dv[3] = -d;
          d1 = omegaDot[3 * (b_i - 1) + 1];
          dv[6] = d1;
          dv[1] = d;
          dv[4] = 0.0;
          d = omegaDot[3 * (b_i - 1)];
          dv[7] = -d;
          dv[2] = -d1;
          dv[5] = d;
          dv[8] = 0.0;
          c_Ry_tmp[0] = 0.0;
          b_Rz_tmp = omega[3 * (b_i - 1) + 2];
          c_Ry_tmp[3] = -b_Rz_tmp;
          Ry_tmp = omega[3 * (b_i - 1) + 1];
          c_Ry_tmp[6] = Ry_tmp;
          c_Ry_tmp[1] = b_Rz_tmp;
          c_Ry_tmp[4] = 0.0;
          b_Ry_tmp = omega[3 * (b_i - 1)];
          c_Ry_tmp[7] = -b_Ry_tmp;
          c_Ry_tmp[2] = -Ry_tmp;
          c_Ry_tmp[5] = b_Ry_tmp;
          c_Ry_tmp[8] = 0.0;
          c_Rz_tmp[0] = 0.0;
          c_Rz_tmp[3] = -b_Rz_tmp;
          c_Rz_tmp[6] = Ry_tmp;
          c_Rz_tmp[1] = b_Rz_tmp;
          c_Rz_tmp[4] = 0.0;
          c_Rz_tmp[7] = -b_Ry_tmp;
          c_Rz_tmp[2] = -Ry_tmp;
          c_Rz_tmp[5] = b_Ry_tmp;
          c_Rz_tmp[8] = 0.0;
          for (i = 0; i < 3; i++) {
            d = 0.0;
            d1 = 0.0;
            b_Rz_tmp = c_Ry_tmp[i];
            Ry_tmp = c_Ry_tmp[i + 3];
            b_Ry_tmp = c_Ry_tmp[i + 6];
            for (b_k = 0; b_k < 3; b_k++) {
              d_Rz_tmp = r_link[b_k + 3 * (b_i - 1)];
              d += dv[i + 3 * b_k] * d_Rz_tmp;
              d1 += ((b_Rz_tmp * c_Rz_tmp[3 * b_k] +
                      Ry_tmp * c_Rz_tmp[3 * b_k + 1]) +
                     b_Ry_tmp * c_Rz_tmp[3 * b_k + 2]) *
                    d_Rz_tmp;
            }
            Fi_inf[i] = (a_j[i + 3 * (b_i - 1)] + d) + d1;
          }
          d = Fi_inf[0];
          d1 = Fi_inf[1];
          b_Rz_tmp = Fi_inf[2];
          for (i = 0; i < 3; i++) {
            a_j[i + 3 * b_i] = 0.0;
            a_j[i + 3 * b_i] = a_j[i + 3 * b_i] + R_T[i + 9 * (b_i - 1)] * d;
            a_j[i + 3 * b_i] =
                a_j[i + 3 * b_i] + R_T[(i + 9 * (b_i - 1)) + 3] * d1;
            a_j[i + 3 * b_i] =
                a_j[i + 3 * b_i] + R_T[(i + 9 * (b_i - 1)) + 6] * b_Rz_tmp;
          }
          guard1 = true;
        } else if (d == 0.0) {
          // 'recursiveNEurdfIA:89' elseif jointType(i-1) == 0
          //  prismatic joint
          //  not implemented
          // 'recursiveNEurdfIA:91' u_out_inf = zeros(1,k-1);
          u_out_inf.set_size(1, static_cast<int>(k));
          loop_ub = static_cast<int>(k);
          for (i = 0; i < loop_ub; i++) {
            u_out_inf[i] = 0.0;
          }
          // ones(1,k-1)*NaN;
          // 'recursiveNEurdfIA:92' u_out_sup = zeros(1,k-1);
          u_out_sup.set_size(1, static_cast<int>(k));
          loop_ub = static_cast<int>(k);
          for (i = 0; i < loop_ub; i++) {
            u_out_sup[i] = 0.0;
          }
          // ones(1,k-1)*NaN;
          exitg2 = 1;
        } else {
          guard1 = true;
        }
      }
      if (guard1) {
        // 'recursiveNEurdfIA:96' [dwdw_inf, dwdw_sup] =
        // prodMatMatInt(Xcross(omegaDot(:,i)), Xcross(omegaDot(:,i)),
        // r_inf(:,i), r_sup(:,i)); UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        dv[0] = 0.0;
        dv[3] = -omegaDot[3 * b_i + 2];
        dv[6] = omegaDot[3 * b_i + 1];
        dv[1] = omegaDot[3 * b_i + 2];
        dv[4] = 0.0;
        dv[7] = -omegaDot[3 * b_i];
        dv[2] = -omegaDot[3 * b_i + 1];
        dv[5] = omegaDot[3 * b_i];
        dv[8] = 0.0;
        c_Ry_tmp[0] = 0.0;
        c_Ry_tmp[3] = -omegaDot[3 * b_i + 2];
        c_Ry_tmp[6] = omegaDot[3 * b_i + 1];
        c_Ry_tmp[1] = omegaDot[3 * b_i + 2];
        c_Ry_tmp[4] = 0.0;
        c_Ry_tmp[7] = -omegaDot[3 * b_i];
        c_Ry_tmp[2] = -omegaDot[3 * b_i + 1];
        c_Ry_tmp[5] = omegaDot[3 * b_i];
        c_Ry_tmp[8] = 0.0;
        prodMatMatInt(dv, c_Ry_tmp, *(double(*)[3]) & DynParInf_r[3 * b_i],
                      *(double(*)[3]) & DynParSup->r[3 * b_i], Fi_inf, Fi_sup);
        // 'recursiveNEurdfIA:97' [ww_inf, ww_sup] =
        // prodMatMatInt(Xcross(omega(:,i))*Xcross(omega(:,i)),
        // Xcross(omega(:,i))*Xcross(omega(:,i)), r_inf(:,i), r_sup(:,i));
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        dv[0] = 0.0;
        dv[3] = -omega[3 * b_i + 2];
        dv[6] = omega[3 * b_i + 1];
        dv[1] = omega[3 * b_i + 2];
        dv[4] = 0.0;
        dv[7] = -omega[3 * b_i];
        dv[2] = -omega[3 * b_i + 1];
        dv[5] = omega[3 * b_i];
        dv[8] = 0.0;
        c_Ry_tmp[0] = 0.0;
        c_Ry_tmp[3] = -omega[3 * b_i + 2];
        c_Ry_tmp[6] = omega[3 * b_i + 1];
        c_Ry_tmp[1] = omega[3 * b_i + 2];
        c_Ry_tmp[4] = 0.0;
        c_Ry_tmp[7] = -omega[3 * b_i];
        c_Ry_tmp[2] = -omega[3 * b_i + 1];
        c_Ry_tmp[5] = omega[3 * b_i];
        c_Ry_tmp[8] = 0.0;
        c_Rz_tmp[0] = 0.0;
        c_Rz_tmp[3] = -omega[3 * b_i + 2];
        c_Rz_tmp[6] = omega[3 * b_i + 1];
        c_Rz_tmp[1] = omega[3 * b_i + 2];
        c_Rz_tmp[4] = 0.0;
        c_Rz_tmp[7] = -omega[3 * b_i];
        c_Rz_tmp[2] = -omega[3 * b_i + 1];
        c_Rz_tmp[5] = omega[3 * b_i];
        c_Rz_tmp[8] = 0.0;
        e_Rz_tmp[0] = 0.0;
        e_Rz_tmp[3] = -omega[3 * b_i + 2];
        e_Rz_tmp[6] = omega[3 * b_i + 1];
        e_Rz_tmp[1] = omega[3 * b_i + 2];
        e_Rz_tmp[4] = 0.0;
        e_Rz_tmp[7] = -omega[3 * b_i];
        e_Rz_tmp[2] = -omega[3 * b_i + 1];
        e_Rz_tmp[5] = omega[3 * b_i];
        e_Rz_tmp[8] = 0.0;
        for (i = 0; i < 3; i++) {
          d = dv[i];
          d1 = dv[i + 3];
          b_Rz_tmp = dv[i + 6];
          for (b_k = 0; b_k < 3; b_k++) {
            xw_I_sup[i + 3 * b_k] =
                (d * c_Ry_tmp[3 * b_k] + d1 * c_Ry_tmp[3 * b_k + 1]) +
                b_Rz_tmp * c_Ry_tmp[3 * b_k + 2];
          }
          d = c_Rz_tmp[i];
          d1 = c_Rz_tmp[i + 3];
          b_Rz_tmp = c_Rz_tmp[i + 6];
          for (b_k = 0; b_k < 3; b_k++) {
            dv[i + 3 * b_k] =
                (d * e_Rz_tmp[3 * b_k] + d1 * e_Rz_tmp[3 * b_k + 1]) +
                b_Rz_tmp * e_Rz_tmp[3 * b_k + 2];
          }
        }
        prodMatMatInt(xw_I_sup, dv, *(double(*)[3]) & DynParInf_r[3 * b_i],
                      *(double(*)[3]) & DynParSup->r[3 * b_i], Rf_inf, Rf_sup);
        // 'recursiveNEurdfIA:98' a_c_inf(:,i) = a_j(:,i) + dwdw_inf + ww_inf;
        // 'recursiveNEurdfIA:99' a_c_sup(:,i) = a_j(:,i) + dwdw_sup + ww_sup;
        a_c_inf[3 * b_i] = (a_j[3 * b_i] + Fi_inf[0]) + Rf_inf[0];
        a_c_sup[3 * b_i] = (a_j[3 * b_i] + Fi_sup[0]) + Rf_sup[0];
        a_c_inf[3 * b_i + 1] = (a_j[3 * b_i + 1] + Fi_inf[1]) + Rf_inf[1];
        a_c_sup[3 * b_i + 1] = (a_j[3 * b_i + 1] + Fi_sup[1]) + Rf_sup[1];
        a_c_inf[3 * b_i + 2] = (a_j[3 * b_i + 2] + Fi_inf[2]) + Rf_inf[2];
        a_c_sup[3 * b_i + 2] = (a_j[3 * b_i + 2] + Fi_sup[2]) + Rf_sup[2];
        // a_c(:,i) = a_j(:,i) + Xcross(omegaDot(:,i))*r(:,i) +
        // Xcross(omega(:,i))*Xcross(omega(:,i))*r(:,i);
        b_i++;
      }
    } else {
      //  backward recursion
      // 'recursiveNEurdfIA:105' for i = k:-1:2
      b_i = 0;
      exitg2 = 2;
    }
  } while (exitg2 == 0);
  if (exitg2 != 1) {
    int exitg1;
    do {
      exitg1 = 0;
      if (b_i <=
          -static_cast<int>((-1.0 - static_cast<double>(k + 1U)) + 2.0) - 1) {
        double varargin_1_tmp[4];
        double Rmu_inf[3];
        double Rmu_sup[3];
        double x_Rf_rr_inf[3];
        double x_Rf_rr_sup[3];
        double b_varargin_1_tmp[2];
        double c_varargin_1_tmp[2];
        double d_varargin_1_tmp[2];
        double e_varargin_1_tmp[2];
        double f_varargin_1_tmp[2];
        double g_varargin_1_tmp[2];
        double c_i;
        c_i = static_cast<double>(k + 1U) + -static_cast<double>(b_i);
        //  k-1:-1:1
        // 'recursiveNEurdfIA:107' if i == k
        if (c_i == k + 1U) {
          // k-1 % force and torques at the end effector
          // 'recursiveNEurdfIA:108' f_inf(:,i+1) = zeros(3,1);
          // 'recursiveNEurdfIA:109' f_sup(:,i+1) = zeros(3,1);
          // 'recursiveNEurdfIA:110' mu_inf(:,i+1) = zeros(3,1);
          // 'recursiveNEurdfIA:111' mu_sup(:,i+1) = zeros(3,1);
          // 'recursiveNEurdfIA:112' R(:,:,i+1) = eye(3);
          loop_ub = static_cast<int>(c_i);
          for (i = 0; i < 3; i++) {
            f_inf[i + 3 * loop_ub] = 0.0;
            f_sup[i + 3 * loop_ub] = 0.0;
            mu_inf[i + 3 * loop_ub] = 0.0;
            mu_sup[i + 3 * loop_ub] = 0.0;
            R[3 * i + 9 * loop_ub] = 0.0;
            R[(3 * i + 9 * loop_ub) + 1] = 0.0;
            R[(3 * i + 9 * loop_ub) + 2] = 0.0;
          }
          loop_ub = static_cast<int>(c_i);
          R[9 * loop_ub] = 1.0;
          R[9 * loop_ub + 4] = 1.0;
          R[9 * loop_ub + 8] = 1.0;
        }
        // 'recursiveNEurdfIA:115' [Fi_inf, Fi_sup] =
        // prodMatScalInt(a_c_inf(:,i), a_c_sup(:,i), m_inf(i), m_sup(i));
        Ry_tmp = DynParInf_m[static_cast<int>(c_i) - 1];
        b_Ry_tmp = DynParSup->m[static_cast<int>(c_i) - 1];
        // UNTITLED3 Summary of this function goes here
        //    Detailed explanation goes here
        // 'prodMatScalInt:5' CInf = zeros(size(AInf));
        // 'prodMatScalInt:6' CSup = zeros(size(AInf));
        // 'prodMatScalInt:8' for i = 1:size(AInf,1)
        //  Fi = m(i)*a_c(:,i);
        // 'recursiveNEurdfIA:117' [Rf_inf, Rf_sup] = prodMatMatInt(R(:,:,i),
        // R(:,:,i), f_inf(:,i+1), f_sup(:,i+1));
        prodMatMatInt(
            *(double(*)[9]) & R[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[9]) & R[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[3]) &
                f_inf[3 * static_cast<int>(static_cast<unsigned int>(c_i))],
            *(double(*)[3]) &
                f_sup[3 * static_cast<int>(static_cast<unsigned int>(c_i))],
            Rf_inf, Rf_sup);
        // 'recursiveNEurdfIA:118' [f_inf(:,i), f_sup(:,i)] =
        // sumMatMatInt(Rf_inf, Rf_sup, Fi_inf, Fi_sup); UNTITLED3 Summary of
        // this function goes here
        //    Detailed explanation goes here
        // 'sumMatMatInt:5' CInf = zeros(size(AInf));
        // 'sumMatMatInt:6' CSup = zeros(size(AInf));
        // 'sumMatMatInt:8' for i = 1:size(AInf,1)
        //  f(:,i) = R(:,:,i)*f(:,i+1) + Fi;
        // 'recursiveNEurdfIA:121' [Rmu_inf, Rmu_sup] =
        // prodMatMatInt(R(:,:,i),R(:,:,i),mu_inf(:,i+1),mu_sup(:,i+1));
        prodMatMatInt(
            *(double(*)[9]) & R[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[9]) & R[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[3]) &
                mu_inf[3 * static_cast<int>(static_cast<unsigned int>(c_i))],
            *(double(*)[3]) &
                mu_sup[3 * static_cast<int>(static_cast<unsigned int>(c_i))],
            Rmu_inf, Rmu_sup);
        // 'recursiveNEurdfIA:122' [Rf_inf2, Rf_sup2] =
        // prodMatMatInt(R(:,:,i),R(:,:,i),f_inf(:,i+1),f_sup(:,i+1));
        // 'recursiveNEurdfIA:123' [r_rlink_inf, r_rlink_sup] =
        // subMatMatInt(r_inf(:,i),r_sup(:,i),r_link(:,i),r_link(:,i));
        // UNTITLED3 Summary of this function goes here
        //    Detailed explanation goes here
        // 'subMatMatInt:5' CInf = zeros(size(AInf));
        // 'subMatMatInt:6' CSup = zeros(size(AInf));
        // 'subMatMatInt:8' for i = 1:size(AInf,1)
        // 'prodMatScalInt:9' for j = 1:size(AInf,2)
        // 'prodMatScalInt:10' [CInf(i,j), CSup(i,j)] = prodInt(AInf(i,j),
        // ASup(i,j), bInf, bSup); UNTITLED3 Summary of this function goes here
        //    Detailed explanation goes here
        //  if bInf > bSup || aInf > aSup
        //      display('inf value greater than sup')
        //  end
        // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf, aSup*bSup]);
        Rz_tmp = a_c_inf[3 * (static_cast<int>(c_i) - 1)];
        varargin_1_tmp[0] = Rz_tmp * Ry_tmp;
        varargin_1_tmp[1] = Rz_tmp * b_Ry_tmp;
        Rz_tmp = a_c_sup[3 * (static_cast<int>(c_i) - 1)];
        varargin_1_tmp[2] = Rz_tmp * Ry_tmp;
        varargin_1_tmp[3] = Rz_tmp * b_Ry_tmp;
        // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
        // aSup*bSup]);
        //         if bInf > bSup || AInf(i,j) > ASup(i,j)
        //             display('inf value greater than sup')
        //         end
        // 'sumMatMatInt:9' for j = 1:size(AInf,2)
        // 'sumMatMatInt:10' CInf(i,j) = AInf(i,j) + BInf(i,j);
        // 'sumMatMatInt:11' CSup(i,j) = ASup(i,j) + BSup(i,j);
        //         if BInf(i,j) > BSup(i,j) || AInf(i,j) > ASup(i,j)
        //             display('inf value greater than sup')
        //         end
        f_inf[3 * (static_cast<int>(c_i) - 1)] =
            Rf_inf[0] + control_coder::coder::internal::minimum(varargin_1_tmp);
        f_sup[3 * (static_cast<int>(c_i) - 1)] =
            Rf_sup[0] + control_coder::coder::internal::maximum(varargin_1_tmp);
        // 'subMatMatInt:9' for j = 1:size(AInf,2)
        // 'subMatMatInt:10' CInf(i,j) = AInf(i,j) - BSup(i,j);
        b_Rz_tmp = r_link[3 * (static_cast<int>(c_i) - 1)];
        Rf_inf[0] = DynParInf_r[3 * (static_cast<int>(c_i) - 1)] - b_Rz_tmp;
        // 'subMatMatInt:11' CSup(i,j) = ASup(i,j) - BInf(i,j);
        Rf_sup[0] = DynParSup->r[3 * (static_cast<int>(c_i) - 1)] - b_Rz_tmp;
        // 'prodMatScalInt:9' for j = 1:size(AInf,2)
        // 'prodMatScalInt:10' [CInf(i,j), CSup(i,j)] = prodInt(AInf(i,j),
        // ASup(i,j), bInf, bSup); UNTITLED3 Summary of this function goes here
        //    Detailed explanation goes here
        //  if bInf > bSup || aInf > aSup
        //      display('inf value greater than sup')
        //  end
        // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf, aSup*bSup]);
        Rz_tmp = a_c_inf[3 * (static_cast<int>(c_i) - 1) + 1];
        varargin_1_tmp[0] = Rz_tmp * Ry_tmp;
        varargin_1_tmp[1] = Rz_tmp * b_Ry_tmp;
        Rz_tmp = a_c_sup[3 * (static_cast<int>(c_i) - 1) + 1];
        varargin_1_tmp[2] = Rz_tmp * Ry_tmp;
        varargin_1_tmp[3] = Rz_tmp * b_Ry_tmp;
        // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
        // aSup*bSup]);
        //         if bInf > bSup || AInf(i,j) > ASup(i,j)
        //             display('inf value greater than sup')
        //         end
        // 'sumMatMatInt:9' for j = 1:size(AInf,2)
        // 'sumMatMatInt:10' CInf(i,j) = AInf(i,j) + BInf(i,j);
        // 'sumMatMatInt:11' CSup(i,j) = ASup(i,j) + BSup(i,j);
        //         if BInf(i,j) > BSup(i,j) || AInf(i,j) > ASup(i,j)
        //             display('inf value greater than sup')
        //         end
        f_inf[3 * (static_cast<int>(c_i) - 1) + 1] =
            Rf_inf[1] + control_coder::coder::internal::minimum(varargin_1_tmp);
        f_sup[3 * (static_cast<int>(c_i) - 1) + 1] =
            Rf_sup[1] + control_coder::coder::internal::maximum(varargin_1_tmp);
        // 'subMatMatInt:9' for j = 1:size(AInf,2)
        // 'subMatMatInt:10' CInf(i,j) = AInf(i,j) - BSup(i,j);
        b_Rz_tmp = r_link[3 * (static_cast<int>(c_i) - 1) + 1];
        Rf_inf[1] = DynParInf_r[3 * (static_cast<int>(c_i) - 1) + 1] - b_Rz_tmp;
        // 'subMatMatInt:11' CSup(i,j) = ASup(i,j) - BInf(i,j);
        Rf_sup[1] =
            DynParSup->r[3 * (static_cast<int>(c_i) - 1) + 1] - b_Rz_tmp;
        // 'prodMatScalInt:9' for j = 1:size(AInf,2)
        // 'prodMatScalInt:10' [CInf(i,j), CSup(i,j)] = prodInt(AInf(i,j),
        // ASup(i,j), bInf, bSup); UNTITLED3 Summary of this function goes here
        //    Detailed explanation goes here
        //  if bInf > bSup || aInf > aSup
        //      display('inf value greater than sup')
        //  end
        // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf, aSup*bSup]);
        Rz_tmp = a_c_inf[3 * (static_cast<int>(c_i) - 1) + 2];
        varargin_1_tmp[0] = Rz_tmp * Ry_tmp;
        varargin_1_tmp[1] = Rz_tmp * b_Ry_tmp;
        Rz_tmp = a_c_sup[3 * (static_cast<int>(c_i) - 1) + 2];
        varargin_1_tmp[2] = Rz_tmp * Ry_tmp;
        varargin_1_tmp[3] = Rz_tmp * b_Ry_tmp;
        // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
        // aSup*bSup]);
        //         if bInf > bSup || AInf(i,j) > ASup(i,j)
        //             display('inf value greater than sup')
        //         end
        // 'sumMatMatInt:9' for j = 1:size(AInf,2)
        // 'sumMatMatInt:10' CInf(i,j) = AInf(i,j) + BInf(i,j);
        // 'sumMatMatInt:11' CSup(i,j) = ASup(i,j) + BSup(i,j);
        //         if BInf(i,j) > BSup(i,j) || AInf(i,j) > ASup(i,j)
        //             display('inf value greater than sup')
        //         end
        f_inf[3 * (static_cast<int>(c_i) - 1) + 2] =
            Rf_inf[2] + control_coder::coder::internal::minimum(varargin_1_tmp);
        f_sup[3 * (static_cast<int>(c_i) - 1) + 2] =
            Rf_sup[2] + control_coder::coder::internal::maximum(varargin_1_tmp);
        // 'subMatMatInt:9' for j = 1:size(AInf,2)
        // 'subMatMatInt:10' CInf(i,j) = AInf(i,j) - BSup(i,j);
        b_Rz_tmp = r_link[3 * (static_cast<int>(c_i) - 1) + 2];
        Rf_inf[2] = DynParInf_r[3 * (static_cast<int>(c_i) - 1) + 2] - b_Rz_tmp;
        // 'subMatMatInt:11' CSup(i,j) = ASup(i,j) - BInf(i,j);
        Rf_sup[2] =
            DynParSup->r[3 * (static_cast<int>(c_i) - 1) + 2] - b_Rz_tmp;
        prodMatMatInt(
            *(double(*)[9]) & R[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[9]) & R[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[3]) &
                f_inf[3 * static_cast<int>(static_cast<unsigned int>(c_i))],
            *(double(*)[3]) &
                f_sup[3 * static_cast<int>(static_cast<unsigned int>(c_i))],
            Fi_inf, Fi_sup);
        // 'recursiveNEurdfIA:124' [cross_Rf_inf, cross_Rf_sup] =
        // XcrossArithm(Rf_inf2,Rf_sup2); UNTITLED7 Summary of this function
        // goes here
        //    Detailed explanation goes here
        // 'XcrossArithm:5' XMat_inf = [0 min([-x_inf(3),-x_sup(3)])
        // min([x_inf(2),x_sup(2)]); 'XcrossArithm:6' min([x_inf(3),x_sup(3)]) 0
        // min([-x_inf(1),-x_sup(1)]); 'XcrossArithm:7'
        // min([-x_inf(2),-x_sup(2)]) min([x_inf(1),x_sup(1)]) 0];
        b_varargin_1_tmp[0] = -Fi_inf[2];
        b_varargin_1_tmp[1] = -Fi_sup[2];
        c_varargin_1_tmp[0] = Fi_inf[1];
        c_varargin_1_tmp[1] = Fi_sup[1];
        d_varargin_1_tmp[0] = Fi_inf[2];
        d_varargin_1_tmp[1] = Fi_sup[2];
        e_varargin_1_tmp[0] = -Fi_inf[0];
        e_varargin_1_tmp[1] = -Fi_sup[0];
        f_varargin_1_tmp[0] = -Fi_inf[1];
        f_varargin_1_tmp[1] = -Fi_sup[1];
        g_varargin_1_tmp[0] = Fi_inf[0];
        g_varargin_1_tmp[1] = Fi_sup[0];
        // 'XcrossArithm:9' XMat_sup = [0 max([-x_inf(3),-x_sup(3)])
        // max([x_inf(2),x_sup(2)]); 'XcrossArithm:10' max([x_inf(3),x_sup(3)])
        // 0 max([-x_inf(1),-x_sup(1)]); 'XcrossArithm:11'
        // max([-x_inf(2),-x_sup(2)]) max([x_inf(1),x_sup(1)]) 0];
        // 'recursiveNEurdfIA:125' [x_Rf_rr_inf, x_Rf_rr_sup] =
        // prodMatMatInt(cross_Rf_inf,cross_Rf_sup,r_rlink_inf,r_rlink_sup);
        dv[0] = 0.0;
        dv[3] = control_coder::coder::internal::b_minimum(b_varargin_1_tmp);
        dv[6] = control_coder::coder::internal::b_minimum(c_varargin_1_tmp);
        dv[1] = control_coder::coder::internal::b_minimum(d_varargin_1_tmp);
        dv[4] = 0.0;
        dv[7] = control_coder::coder::internal::b_minimum(e_varargin_1_tmp);
        dv[2] = control_coder::coder::internal::b_minimum(f_varargin_1_tmp);
        dv[5] = control_coder::coder::internal::b_minimum(g_varargin_1_tmp);
        dv[8] = 0.0;
        c_Ry_tmp[0] = 0.0;
        c_Ry_tmp[3] =
            control_coder::coder::internal::b_maximum(b_varargin_1_tmp);
        c_Ry_tmp[6] =
            control_coder::coder::internal::b_maximum(c_varargin_1_tmp);
        c_Ry_tmp[1] =
            control_coder::coder::internal::b_maximum(d_varargin_1_tmp);
        c_Ry_tmp[4] = 0.0;
        c_Ry_tmp[7] =
            control_coder::coder::internal::b_maximum(e_varargin_1_tmp);
        c_Ry_tmp[2] =
            control_coder::coder::internal::b_maximum(f_varargin_1_tmp);
        c_Ry_tmp[5] =
            control_coder::coder::internal::b_maximum(g_varargin_1_tmp);
        c_Ry_tmp[8] = 0.0;
        prodMatMatInt(dv, c_Ry_tmp, Rf_inf, Rf_sup, x_Rf_rr_inf, x_Rf_rr_sup);
        // 'recursiveNEurdfIA:126' [cross_f_inf, cross_f_sup] =
        // XcrossArithm(f_inf(:,i),f_sup(:,i)); UNTITLED7 Summary of this
        // function goes here
        //    Detailed explanation goes here
        // 'XcrossArithm:5' XMat_inf = [0 min([-x_inf(3),-x_sup(3)])
        // min([x_inf(2),x_sup(2)]); 'XcrossArithm:6' min([x_inf(3),x_sup(3)]) 0
        // min([-x_inf(1),-x_sup(1)]); 'XcrossArithm:7'
        // min([-x_inf(2),-x_sup(2)]) min([x_inf(1),x_sup(1)]) 0];
        Rz_tmp = f_inf[3 * (static_cast<int>(c_i) - 1) + 2];
        b_varargin_1_tmp[0] = -Rz_tmp;
        b_Rz_tmp = f_sup[3 * (static_cast<int>(c_i) - 1) + 2];
        b_varargin_1_tmp[1] = -b_Rz_tmp;
        Ry_tmp = f_inf[3 * (static_cast<int>(c_i) - 1) + 1];
        c_varargin_1_tmp[0] = Ry_tmp;
        b_Ry_tmp = f_sup[3 * (static_cast<int>(c_i) - 1) + 1];
        c_varargin_1_tmp[1] = b_Ry_tmp;
        d_varargin_1_tmp[0] = Rz_tmp;
        d_varargin_1_tmp[1] = b_Rz_tmp;
        Rz_tmp = f_inf[3 * (static_cast<int>(c_i) - 1)];
        e_varargin_1_tmp[0] = -Rz_tmp;
        b_Rz_tmp = f_sup[3 * (static_cast<int>(c_i) - 1)];
        e_varargin_1_tmp[1] = -b_Rz_tmp;
        f_varargin_1_tmp[0] = -Ry_tmp;
        f_varargin_1_tmp[1] = -b_Ry_tmp;
        g_varargin_1_tmp[0] = Rz_tmp;
        g_varargin_1_tmp[1] = b_Rz_tmp;
        // 'XcrossArithm:9' XMat_sup = [0 max([-x_inf(3),-x_sup(3)])
        // max([x_inf(2),x_sup(2)]); 'XcrossArithm:10' max([x_inf(3),x_sup(3)])
        // 0 max([-x_inf(1),-x_sup(1)]); 'XcrossArithm:11'
        // max([-x_inf(2),-x_sup(2)]) max([x_inf(1),x_sup(1)]) 0];
        // 'recursiveNEurdfIA:127' [x_f_r_inf, x_f_r_sup] =
        // prodMatMatInt(cross_f_inf,cross_f_sup,r_inf(:,i),r_sup(:,i));
        dv[0] = 0.0;
        dv[3] = control_coder::coder::internal::b_minimum(b_varargin_1_tmp);
        dv[6] = control_coder::coder::internal::b_minimum(c_varargin_1_tmp);
        dv[1] = control_coder::coder::internal::b_minimum(d_varargin_1_tmp);
        dv[4] = 0.0;
        dv[7] = control_coder::coder::internal::b_minimum(e_varargin_1_tmp);
        dv[2] = control_coder::coder::internal::b_minimum(f_varargin_1_tmp);
        dv[5] = control_coder::coder::internal::b_minimum(g_varargin_1_tmp);
        dv[8] = 0.0;
        c_Ry_tmp[0] = 0.0;
        c_Ry_tmp[3] =
            control_coder::coder::internal::b_maximum(b_varargin_1_tmp);
        c_Ry_tmp[6] =
            control_coder::coder::internal::b_maximum(c_varargin_1_tmp);
        c_Ry_tmp[1] =
            control_coder::coder::internal::b_maximum(d_varargin_1_tmp);
        c_Ry_tmp[4] = 0.0;
        c_Ry_tmp[7] =
            control_coder::coder::internal::b_maximum(e_varargin_1_tmp);
        c_Ry_tmp[2] =
            control_coder::coder::internal::b_maximum(f_varargin_1_tmp);
        c_Ry_tmp[5] =
            control_coder::coder::internal::b_maximum(g_varargin_1_tmp);
        c_Ry_tmp[8] = 0.0;
        prodMatMatInt(
            dv, c_Ry_tmp,
            *(double(*)[3]) & DynParInf_r[3 * (static_cast<int>(c_i) - 1)],
            *(double(*)[3]) & DynParSup->r[3 * (static_cast<int>(c_i) - 1)],
            Fi_inf, Fi_sup);
        // 'recursiveNEurdfIA:128' [I_dw_inf,I_dw_sup] =
        // prodMatMatInt(I_inf(:,:,i),I_sup(:,:,i),omegaDot(:,i),omegaDot(:,i));
        prodMatMatInt(
            *(double(*)[9]) & DynParInf_I[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[9]) & DynParSup->b_I[9 * (static_cast<int>(c_i) - 1)],
            *(double(*)[3]) & omegaDot[3 * (static_cast<int>(c_i) - 1)],
            *(double(*)[3]) & omegaDot[3 * (static_cast<int>(c_i) - 1)], Rf_inf,
            Rf_sup);
        // 'recursiveNEurdfIA:129' [xw_I_inf, xw_I_sup] = prodMatMatInt(
        // Xcross(omega(:,i)),Xcross(omega(:,i)),I_inf(:,:,i),I_sup(:,:,i));
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        c_Ry_tmp[0] = 0.0;
        b_Rz_tmp = omega[3 * (static_cast<int>(c_i) - 1) + 2];
        c_Ry_tmp[3] = -b_Rz_tmp;
        Ry_tmp = omega[3 * (static_cast<int>(c_i) - 1) + 1];
        c_Ry_tmp[6] = Ry_tmp;
        c_Ry_tmp[1] = b_Rz_tmp;
        c_Ry_tmp[4] = 0.0;
        b_Ry_tmp = omega[3 * (static_cast<int>(c_i) - 1)];
        c_Ry_tmp[7] = -b_Ry_tmp;
        c_Ry_tmp[2] = -Ry_tmp;
        c_Ry_tmp[5] = b_Ry_tmp;
        c_Ry_tmp[8] = 0.0;
        // UNTITLED7 Summary of this function goes here
        //    Detailed explanation goes here
        // 'Xcross:5' XMat = [0 -x(3) x(2);
        // 'Xcross:6'         x(3) 0 -x(1);
        // 'Xcross:7'         -x(2) x(1) 0];
        c_Rz_tmp[0] = 0.0;
        c_Rz_tmp[3] = -b_Rz_tmp;
        c_Rz_tmp[6] = Ry_tmp;
        c_Rz_tmp[1] = b_Rz_tmp;
        c_Rz_tmp[4] = 0.0;
        c_Rz_tmp[7] = -b_Ry_tmp;
        c_Rz_tmp[2] = -Ry_tmp;
        c_Rz_tmp[5] = b_Ry_tmp;
        c_Rz_tmp[8] = 0.0;
        // UNTITLED3 Summary of this function goes here
        //    Detailed explanation goes here
        // 'prodMatMatInt:5' CInf = zeros(size(AInf,1),size(BInf,2));
        // 'prodMatMatInt:6' CSup = zeros(size(AInf,1),size(BInf,2));
        // 'prodMatMatInt:8' temp1Inf = 0;
        //  inizialization for simulink
        // 'prodMatMatInt:9' temp1Sup = 0;
        //  inizialization for simulink
        // 'prodMatMatInt:11' for i = 1:size(AInf,1)
        for (d_i = 0; d_i < 3; d_i++) {
          // 'prodMatMatInt:12' for k = 1:size(BInf,2)
          d = c_Ry_tmp[d_i];
          Rz_tmp = c_Ry_tmp[d_i + 3];
          b_Rz_tmp = c_Ry_tmp[d_i + 6];
          d1 = c_Rz_tmp[d_i];
          Ry_tmp = c_Rz_tmp[d_i + 3];
          b_Ry_tmp = c_Rz_tmp[d_i + 6];
          for (b_k = 0; b_k < 3; b_k++) {
            // 'prodMatMatInt:13' for j = 1:size(AInf,2)
            // 'prodMatMatInt:14' if j == 1
            // 'prodMatMatInt:15' [temp1Inf, temp1Sup] =
            // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3
            // Summary of this function goes here
            //    Detailed explanation goes here
            //  if bInf > bSup || aInf > aSup
            //      display('inf value greater than sup')
            //  end
            // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
            // aSup*bSup]);
            d_Rz_tmp = DynParInf_I[3 * b_k + 9 * (static_cast<int>(c_i) - 1)];
            varargin_1_tmp[0] = d * d_Rz_tmp;
            Rx_tmp = DynParSup->b_I[3 * b_k + 9 * (static_cast<int>(c_i) - 1)];
            varargin_1_tmp[1] = d * Rx_tmp;
            varargin_1_tmp[2] = d1 * d_Rz_tmp;
            varargin_1_tmp[3] = d1 * Rx_tmp;
            temp1Inf = control_coder::coder::internal::minimum(varargin_1_tmp);
            // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
            // aSup*bSup]);
            temp1Sup = control_coder::coder::internal::maximum(varargin_1_tmp);
            //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
            //                     display('inf value greater than sup')
            //                 end
            // 'prodMatMatInt:14' if j == 1
            // 'prodMatMatInt:19' else
            // 'prodMatMatInt:20' [temp2Inf, temp2Sup] =
            // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3
            // Summary of this function goes here
            //    Detailed explanation goes here
            //  if bInf > bSup || aInf > aSup
            //      display('inf value greater than sup')
            //  end
            // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
            // aSup*bSup]);
            d_Rz_tmp =
                DynParInf_I[(3 * b_k + 9 * (static_cast<int>(c_i) - 1)) + 1];
            varargin_1_tmp[0] = Rz_tmp * d_Rz_tmp;
            Rx_tmp =
                DynParSup->b_I[(3 * b_k + 9 * (static_cast<int>(c_i) - 1)) + 1];
            varargin_1_tmp[1] = Rz_tmp * Rx_tmp;
            varargin_1_tmp[2] = Ry_tmp * d_Rz_tmp;
            varargin_1_tmp[3] = Ry_tmp * Rx_tmp;
            // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
            // aSup*bSup]); 'prodMatMatInt:21' [temp1Inf, temp1Sup] =
            // sumInt(temp1Inf,temp1Sup,temp2Inf,temp2Sup); UNTITLED3 Summary of
            // this function goes here
            //    Detailed explanation goes here
            //  if bInf > bSup || aInf > aSup
            //      display('inf value greater than sup')
            //  end
            // 'sumInt:9' cInf = aInf + bInf;
            temp1Inf += control_coder::coder::internal::minimum(varargin_1_tmp);
            // 'sumInt:10' cSup = aSup + bSup;
            temp1Sup += control_coder::coder::internal::maximum(varargin_1_tmp);
            //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
            //                     display('inf value greater than sup')
            //                 end
            // 'prodMatMatInt:14' if j == 1
            // 'prodMatMatInt:19' else
            // 'prodMatMatInt:20' [temp2Inf, temp2Sup] =
            // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3
            // Summary of this function goes here
            //    Detailed explanation goes here
            //  if bInf > bSup || aInf > aSup
            //      display('inf value greater than sup')
            //  end
            // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
            // aSup*bSup]);
            d_Rz_tmp =
                DynParInf_I[(3 * b_k + 9 * (static_cast<int>(c_i) - 1)) + 2];
            varargin_1_tmp[0] = b_Rz_tmp * d_Rz_tmp;
            Rx_tmp =
                DynParSup->b_I[(3 * b_k + 9 * (static_cast<int>(c_i) - 1)) + 2];
            varargin_1_tmp[1] = b_Rz_tmp * Rx_tmp;
            varargin_1_tmp[2] = b_Ry_tmp * d_Rz_tmp;
            varargin_1_tmp[3] = b_Ry_tmp * Rx_tmp;
            // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
            // aSup*bSup]); 'prodMatMatInt:21' [temp1Inf, temp1Sup] =
            // sumInt(temp1Inf,temp1Sup,temp2Inf,temp2Sup); UNTITLED3 Summary of
            // this function goes here
            //    Detailed explanation goes here
            //  if bInf > bSup || aInf > aSup
            //      display('inf value greater than sup')
            //  end
            // 'sumInt:9' cInf = aInf + bInf;
            temp1Inf += control_coder::coder::internal::minimum(varargin_1_tmp);
            // 'sumInt:10' cSup = aSup + bSup;
            temp1Sup += control_coder::coder::internal::maximum(varargin_1_tmp);
            //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
            //                     display('inf value greater than sup')
            //                 end
            // 'prodMatMatInt:27' CInf(i,k) = temp1Inf;
            loop_ub = d_i + 3 * b_k;
            e_Rz_tmp[loop_ub] = temp1Inf;
            // 'prodMatMatInt:28' CSup(i,k) = temp1Sup;
            xw_I_sup[loop_ub] = temp1Sup;
          }
        }
        double wIw_inf[3];
        double wIw_sup[3];
        // 'recursiveNEurdfIA:130' [wIw_inf, wIw_sup] =
        // prodMatMatInt(xw_I_inf,xw_I_sup,omega(:,i),omega(:,i));
        prodMatMatInt(e_Rz_tmp, xw_I_sup,
                      *(double(*)[3]) & omega[3 * (static_cast<int>(c_i) - 1)],
                      *(double(*)[3]) & omega[3 * (static_cast<int>(c_i) - 1)],
                      wIw_inf, wIw_sup);
        // 'recursiveNEurdfIA:132' mu_inf(:,i) = Rmu_inf + x_Rf_rr_inf -
        // x_f_r_sup + I_dw_inf + wIw_inf; 'recursiveNEurdfIA:133' mu_sup(:,i) =
        // Rmu_sup + x_Rf_rr_sup - x_f_r_inf + I_dw_sup + wIw_sup;
        mu_inf[3 * (static_cast<int>(c_i) - 1)] =
            (((Rmu_inf[0] + x_Rf_rr_inf[0]) - Fi_sup[0]) + Rf_inf[0]) +
            wIw_inf[0];
        mu_sup[3 * (static_cast<int>(c_i) - 1)] =
            (((Rmu_sup[0] + x_Rf_rr_sup[0]) - Fi_inf[0]) + Rf_sup[0]) +
            wIw_sup[0];
        mu_inf[3 * (static_cast<int>(c_i) - 1) + 1] =
            (((Rmu_inf[1] + x_Rf_rr_inf[1]) - Fi_sup[1]) + Rf_inf[1]) +
            wIw_inf[1];
        mu_sup[3 * (static_cast<int>(c_i) - 1) + 1] =
            (((Rmu_sup[1] + x_Rf_rr_sup[1]) - Fi_inf[1]) + Rf_sup[1]) +
            wIw_sup[1];
        mu_inf[3 * (static_cast<int>(c_i) - 1) + 2] =
            (((Rmu_inf[2] + x_Rf_rr_inf[2]) - Fi_sup[2]) + Rf_inf[2]) +
            wIw_inf[2];
        mu_sup[3 * (static_cast<int>(c_i) - 1) + 2] =
            (((Rmu_sup[2] + x_Rf_rr_sup[2]) - Fi_inf[2]) + Rf_sup[2]) +
            wIw_sup[2];
        //  mu(:,i) = R(:,:,i)*mu(:,i+1) +
        //  Xcross(R(:,:,i)*f(:,i+1))*(r(:,i)-r_link(:,i)) -
        //  Xcross(f(:,i))*r(:,i) + I(:,:,i)*omegaDot(:,i) +
        //  Xcross(omega(:,i))*I(:,:,i)*omega(:,i);
        // 'recursiveNEurdfIA:136' if jointType(i-1) == 1
        d = DynParInf_Jtype[static_cast<int>(c_i) - 2];
        if (d == 1.0) {
          double dv1[4];
          //  revolute joint
          // 'recursiveNEurdfIA:137' [mu_z_inf, mu_z_sup] =
          // prodMatMatInt(mu_inf(:,i)',mu_sup(:,i)',z0,z0); UNTITLED3 Summary
          // of this function goes here
          //    Detailed explanation goes here
          // 'prodMatMatInt:5' CInf = zeros(size(AInf,1),size(BInf,2));
          // 'prodMatMatInt:6' CSup = zeros(size(AInf,1),size(BInf,2));
          // 'prodMatMatInt:8' temp1Inf = 0;
          //  inizialization for simulink
          // 'prodMatMatInt:9' temp1Sup = 0;
          //  inizialization for simulink
          // 'prodMatMatInt:11' for i = 1:size(AInf,1)
          // 'prodMatMatInt:12' for k = 1:size(BInf,2)
          // 'prodMatMatInt:13' for j = 1:size(AInf,2)
          // 'prodMatMatInt:14' if j == 1
          // 'prodMatMatInt:15' [temp1Inf, temp1Sup] =
          // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3 Summary
          // of this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]);
          varargin_1_tmp[0] = 0.0;
          varargin_1_tmp[1] = 0.0;
          varargin_1_tmp[2] = 0.0;
          varargin_1_tmp[3] = 0.0;
          temp1Inf = control_coder::coder::internal::minimum(varargin_1_tmp);
          // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]);
          temp1Sup = control_coder::coder::internal::maximum(varargin_1_tmp);
          //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
          //                     display('inf value greater than sup')
          //                 end
          // 'prodMatMatInt:14' if j == 1
          // 'prodMatMatInt:19' else
          // 'prodMatMatInt:20' [temp2Inf, temp2Sup] =
          // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3 Summary
          // of this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]);
          varargin_1_tmp[0] = 0.0;
          varargin_1_tmp[1] = 0.0;
          varargin_1_tmp[2] = 0.0;
          varargin_1_tmp[3] = 0.0;
          // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]); 'prodMatMatInt:21' [temp1Inf, temp1Sup] =
          // sumInt(temp1Inf,temp1Sup,temp2Inf,temp2Sup); UNTITLED3 Summary of
          // this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'sumInt:9' cInf = aInf + bInf;
          temp1Inf += control_coder::coder::internal::minimum(varargin_1_tmp);
          // 'sumInt:10' cSup = aSup + bSup;
          temp1Sup += control_coder::coder::internal::maximum(varargin_1_tmp);
          //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
          //                     display('inf value greater than sup')
          //                 end
          // 'prodMatMatInt:14' if j == 1
          // 'prodMatMatInt:19' else
          // 'prodMatMatInt:20' [temp2Inf, temp2Sup] =
          // prodInt(AInf(i,j),ASup(i,j),BInf(j,k),BSup(j,k)); UNTITLED3 Summary
          // of this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]);
          Rz_tmp = mu_inf[3 * (static_cast<int>(c_i) - 1) + 2];
          varargin_1_tmp[0] = Rz_tmp;
          varargin_1_tmp[1] = Rz_tmp;
          Rz_tmp = mu_sup[3 * (static_cast<int>(c_i) - 1) + 2];
          varargin_1_tmp[2] = Rz_tmp;
          varargin_1_tmp[3] = Rz_tmp;
          // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]); 'prodMatMatInt:21' [temp1Inf, temp1Sup] =
          // sumInt(temp1Inf,temp1Sup,temp2Inf,temp2Sup); UNTITLED3 Summary of
          // this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'sumInt:9' cInf = aInf + bInf;
          temp1Inf += control_coder::coder::internal::minimum(varargin_1_tmp);
          // 'sumInt:10' cSup = aSup + bSup;
          temp1Sup += control_coder::coder::internal::maximum(varargin_1_tmp);
          //                 if BInf(j,k) > BSup(j,k) || AInf(i,j) > ASup(i,j)
          //                     display('inf value greater than sup')
          //                 end
          // 'prodMatMatInt:27' CInf(i,k) = temp1Inf;
          // 'prodMatMatInt:28' CSup(i,k) = temp1Sup;
          // 'recursiveNEurdfIA:138' [Fv_dq_inf, Fv_dq_sup] =
          // prodInt(Fv_inf(i-1),Fv_sup(i-1),dq(i-1),dq(i-1)); UNTITLED3 Summary
          // of this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]);
          b_Rz_tmp = dq[static_cast<int>(c_i) - 2];
          Rz_tmp = DynParInf_Fv[static_cast<int>(c_i) - 2] * b_Rz_tmp;
          varargin_1_tmp[0] = Rz_tmp;
          varargin_1_tmp[1] = Rz_tmp;
          Rz_tmp = DynParSup->Fv[static_cast<int>(c_i) - 2] * b_Rz_tmp;
          varargin_1_tmp[2] = Rz_tmp;
          varargin_1_tmp[3] = Rz_tmp;
          // 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]); 'recursiveNEurdfIA:139' [Fc_dq_inf, Fc_dq_sup] =
          // prodInt(Fc_inf(i-1),Fc_sup(i-1),sin(atan(1000*dq(i-1))),sin(atan(1000*dq(i-1))));
          // UNTITLED3 Summary of this function goes here
          //    Detailed explanation goes here
          //  if bInf > bSup || aInf > aSup
          //      display('inf value greater than sup')
          //  end
          // 'prodInt:9' cInf = min([aInf*bInf, aInf*bSup, aSup*bInf,
          // aSup*bSup]); 'prodInt:10' cSup = max([aInf*bInf, aInf*bSup,
          // aSup*bInf, aSup*bSup]); 'recursiveNEurdfIA:141' u_inf(i) = mu_z_inf
          // + Fv_dq_inf + Fc_dq_inf + ddq(i-1)*kr(i-1)^2*Im(i-1);
          dv1[0] = 0.0;
          dv1[1] = 0.0;
          dv1[2] = 0.0;
          dv1[3] = 0.0;
          b_Rz_tmp = DynParInf_kr[static_cast<int>(c_i) - 2];
          b_Rz_tmp = ddq[static_cast<int>(c_i) - 2] * (b_Rz_tmp * b_Rz_tmp) *
                     DynParInf_Im[static_cast<int>(c_i) - 2];
          u_inf[static_cast<int>(c_i) - 1] =
              ((temp1Inf +
                control_coder::coder::internal::minimum(varargin_1_tmp)) +
               control_coder::coder::internal::minimum(dv1)) +
              b_Rz_tmp;
          // 'recursiveNEurdfIA:142' u_sup(i) = mu_z_sup + Fv_dq_sup + Fc_dq_sup
          // + ddq(i-1)*kr(i-1)^2*Im(i-1);
          dv1[0] = 0.0;
          dv1[1] = 0.0;
          dv1[2] = 0.0;
          dv1[3] = 0.0;
          u_sup[static_cast<int>(c_i) - 1] =
              ((temp1Sup +
                control_coder::coder::internal::maximum(varargin_1_tmp)) +
               control_coder::coder::internal::maximum(dv1)) +
              b_Rz_tmp;
          //  u(i) = mu(:,i)'*z0 + Fv(i-1)*dq(i-1) +
          //  Fc(i-1)*sin(atan(1000*dq(i-1))) + ddq(i-1)*kr(i-1)^2*Im(i-1);
          b_i++;
        } else if (d == 0.0) {
          // 'recursiveNEurdfIA:145' elseif jointType(i-1) == 0
          //  prismatic joint
          //  not implemented
          // 'recursiveNEurdfIA:147' u_out_inf = zeros(1,k-1);
          u_out_inf.set_size(1, static_cast<int>(k));
          loop_ub = static_cast<int>(k);
          for (i = 0; i < loop_ub; i++) {
            u_out_inf[i] = 0.0;
          }
          // ones(1,k-1)*NaN;
          // 'recursiveNEurdfIA:148' u_out_sup = zeros(1,k-1);
          u_out_sup.set_size(1, static_cast<int>(k));
          loop_ub = static_cast<int>(k);
          for (i = 0; i < loop_ub; i++) {
            u_out_sup[i] = 0.0;
          }
          //  ones(1,k-1)*NaN;
          exitg1 = 1;
        } else {
          b_i++;
        }
      } else {
        // 'recursiveNEurdfIA:153' u_out_inf = u_inf(2:end);
        if (u_inf.size(1) < 2) {
          i = 0;
          b_k = 0;
        } else {
          i = 1;
          b_k = u_inf.size(1);
        }
        loop_ub = b_k - i;
        u_out_inf.set_size(1, loop_ub);
        for (b_k = 0; b_k < loop_ub; b_k++) {
          u_out_inf[b_k] = u_inf[i + b_k];
        }
        // 'recursiveNEurdfIA:154' u_out_sup = u_sup(2:end);
        if (u_sup.size(1) < 2) {
          i = 0;
          b_k = 0;
        } else {
          i = 1;
          b_k = u_sup.size(1);
        }
        loop_ub = b_k - i;
        u_out_sup.set_size(1, loop_ub);
        for (b_k = 0; b_k < loop_ub; b_k++) {
          u_out_sup[b_k] = u_sup[i + b_k];
        }
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
}

//
// File trailer for recursiveNEurdfIA.cpp
//
// [EOF]
//
