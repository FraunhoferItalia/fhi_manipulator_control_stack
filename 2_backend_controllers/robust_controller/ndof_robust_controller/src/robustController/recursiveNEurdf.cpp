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
#include "recursiveNEurdf.h"
#include "robustControllerURDF_data.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// function [ u_out ] = recursiveNEurdf( q,dq,ddq,KinPar,DynPar,g )
//
// UNTITLED18 Summary of this function goes here
//    Detailed explanation goes here
//
// Arguments    : const coder::array<double, 2U> &q
//                const coder::array<double, 2U> &dq
//                const coder::array<double, 2U> &ddq
//                const coder::array<double, 2U> &KinPar_xyz
//                const coder::array<double, 2U> &KinPar_rpy
//                const coder::array<double, 2U> &DynPar_m
//                const coder::array<double, 2U> &DynPar_r
//                const coder::array<double, 2U> &DynPar_kr
//                const coder::array<double, 2U> &DynPar_Im
//                const coder::array<double, 2U> &DynPar_Fv
//                const coder::array<double, 2U> &DynPar_Jtype
//                const coder::array<double, 3U> &DynPar_I
//                const double g[3]
//                coder::array<double, 2U> &u_out
// Return Type  : void
//
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
                     const double g[3], coder::array<double, 2U> &u_out)
{
  static const signed char b[3]{0, 0, 1};
  static const signed char b_iv[3]{1, 0, 0};
  coder::array<double, 3U> R;
  coder::array<double, 3U> R_T;
  coder::array<double, 2U> a_c;
  coder::array<double, 2U> a_j;
  coder::array<double, 2U> f;
  coder::array<double, 2U> mu;
  coder::array<double, 2U> omega;
  coder::array<double, 2U> omegaDot;
  coder::array<double, 2U> r_link;
  coder::array<double, 2U> u;
  double c_Ry_tmp[9];
  double dv[9];
  double e_Rz_tmp[9];
  double dv1[3];
  double y[3];
  double Ry_tmp;
  double b_Ry_tmp;
  double c_Rz_tmp;
  double d;
  double d1;
  double d_Rz_tmp;
  int b_i;
  int exitg2;
  int i;
  int i1;
  unsigned int k;
  int loop_ub;
  //  recursive Newton-Euler algorithm for URDF convention
  // 'recursiveNEurdf:7' jointType = DynPar.Jtype;
  // 'recursiveNEurdf:8' m = DynPar.m;
  // 'recursiveNEurdf:9' I = DynPar.I;
  // 'recursiveNEurdf:10' r = DynPar.r;
  // 'recursiveNEurdf:12' k = length(jointType) + 1;
  k = static_cast<unsigned int>(DynPar_Jtype.size(1));
  //  number of arms = number of joints + 1
  // 'recursiveNEurdf:14' xyz = KinPar.xyz;
  //  transposed
  // 'recursiveNEurdf:15' rpy = KinPar.rpy;
  //  transposed
  // 'recursiveNEurdf:17' kr = DynPar.kr;
  // 'recursiveNEurdf:18' Im = DynPar.Im;
  // 'recursiveNEurdf:19' Fv = DynPar.Fv;
  // 'recursiveNEurdf:20' Fc = DynPar.Fc;
  //  initialize parameters
  // 'recursiveNEurdf:24' omega       = zeros(3,k);
  omega.set_size(3, DynPar_Jtype.size(1) + 1);
  loop_ub = 3 * (DynPar_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    omega[i] = 0.0;
  }
  // 'recursiveNEurdf:25' omegaDot    = zeros(3,k);
  omegaDot.set_size(3, DynPar_Jtype.size(1) + 1);
  loop_ub = 3 * (DynPar_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    omegaDot[i] = 0.0;
  }
  // 'recursiveNEurdf:26' a_j         = zeros(3,k);
  a_j.set_size(3, DynPar_Jtype.size(1) + 1);
  loop_ub = 3 * (DynPar_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    a_j[i] = 0.0;
  }
  // 'recursiveNEurdf:27' a_c         = zeros(3,k);
  a_c.set_size(3, DynPar_Jtype.size(1) + 1);
  loop_ub = 3 * (DynPar_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    a_c[i] = 0.0;
  }
  // 'recursiveNEurdf:28' R           = zeros(3,3,k+1);
  R.set_size(3, 3,
             static_cast<int>(
                 (static_cast<double>(DynPar_Jtype.size(1)) + 1.0) + 1.0));
  loop_ub = 9 * static_cast<int>(
                    (static_cast<double>(DynPar_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    R[i] = 0.0;
  }
  //  from i to i-1
  // 'recursiveNEurdf:29' R_T         = zeros(3,3,k);
  R_T.set_size(3, 3, DynPar_Jtype.size(1) + 1);
  loop_ub = 9 * (DynPar_Jtype.size(1) + 1);
  for (i = 0; i < loop_ub; i++) {
    R_T[i] = 0.0;
  }
  //  from i-1 to i
  // 'recursiveNEurdf:30' z           = zeros(3,k);
  // 'recursiveNEurdf:31' u           = zeros(1,k);
  u.set_size(1, DynPar_Jtype.size(1) + 1);
  loop_ub = DynPar_Jtype.size(1) + 1;
  for (i = 0; i < loop_ub; i++) {
    u[i] = 0.0;
  }
  // 'recursiveNEurdf:32' z0          = [0;0;1];
  // 'recursiveNEurdf:33' f           = zeros(3,k+1);
  f.set_size(3, static_cast<int>(
                    (static_cast<double>(DynPar_Jtype.size(1)) + 1.0) + 1.0));
  loop_ub = 3 * static_cast<int>(
                    (static_cast<double>(DynPar_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    f[i] = 0.0;
  }
  //  forces (k+1 for end effector forces)
  // 'recursiveNEurdf:34' mu          = zeros(3,k+1);
  mu.set_size(3, static_cast<int>(
                     (static_cast<double>(DynPar_Jtype.size(1)) + 1.0) + 1.0));
  loop_ub = 3 * static_cast<int>(
                    (static_cast<double>(DynPar_Jtype.size(1)) + 1.0) + 1.0);
  for (i = 0; i < loop_ub; i++) {
    mu[i] = 0.0;
  }
  //  torques (k+1 for end effector torques)
  // 'recursiveNEurdf:36' r_link = zeros(3,k);
  r_link.set_size(3, DynPar_Jtype.size(1) + 1);
  // 'recursiveNEurdf:38' for i = 1:k
  i = DynPar_Jtype.size(1);
  for (b_i = 0; b_i <= i; b_i++) {
    // 'recursiveNEurdf:39' if i < k
    if (static_cast<unsigned int>(b_i) < k) {
      double g_Rz_tmp[16];
      double f_Rz_tmp[9];
      double Rx_tmp;
      double Rz_tmp;
      double b_Rx_tmp;
      double b_Rz_tmp;
      int i2;
      //  all but last link
      // 'recursiveNEurdf:41' RzCurrent = RzMat(q(i));
      // UNTITLED3 Summary of this function goes here
      //    Detailed explanation goes here
      // 'RzMat:5' Rz = [cos(x), -sin(x), 0 ,0;
      // 'RzMat:6'       sin(x), cos(x), 0, 0;
      // 'RzMat:7'       0, 0, 1, 0;
      // 'RzMat:8'       0, 0, 0, 1];
      Rz_tmp = std::sin(q[b_i]);
      b_Rz_tmp = std::cos(q[b_i]);
      // 'recursiveNEurdf:42' R(:,:,i) =
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
      b_Rx_tmp = std::cos(d);
      // 'rotMatURDF:10' Ry = [cos(p), 0, -sin(p); 0, 1, 0; sin(p), 0, cos(p)];
      d = KinPar_rpy[3 * b_i + 1];
      Ry_tmp = std::sin(d);
      b_Ry_tmp = std::cos(d);
      // 'rotMatURDF:11' Rz = [cos(y), -sin(y), 0; sin(y), cos(y), 0; 0, 0, 1];
      d = KinPar_rpy[3 * b_i + 2];
      c_Rz_tmp = std::sin(d);
      d_Rz_tmp = std::cos(d);
      // 'rotMatURDF:13' R = Rz * Ry * Rx;
      e_Rz_tmp[0] = d_Rz_tmp;
      e_Rz_tmp[3] = -c_Rz_tmp;
      e_Rz_tmp[6] = 0.0;
      e_Rz_tmp[1] = c_Rz_tmp;
      e_Rz_tmp[4] = d_Rz_tmp;
      e_Rz_tmp[7] = 0.0;
      c_Ry_tmp[0] = b_Ry_tmp;
      c_Ry_tmp[3] = 0.0;
      c_Ry_tmp[6] = -Ry_tmp;
      e_Rz_tmp[2] = 0.0;
      c_Ry_tmp[1] = 0.0;
      e_Rz_tmp[5] = 0.0;
      c_Ry_tmp[4] = 1.0;
      e_Rz_tmp[8] = 1.0;
      c_Ry_tmp[7] = 0.0;
      c_Ry_tmp[2] = Ry_tmp;
      c_Ry_tmp[5] = 0.0;
      c_Ry_tmp[8] = b_Ry_tmp;
      for (i1 = 0; i1 < 3; i1++) {
        d = e_Rz_tmp[i1];
        d1 = e_Rz_tmp[i1 + 3];
        loop_ub = static_cast<int>(e_Rz_tmp[i1 + 6]);
        for (i2 = 0; i2 < 3; i2++) {
          f_Rz_tmp[i1 + 3 * i2] =
              (d * c_Ry_tmp[3 * i2] + d1 * c_Ry_tmp[3 * i2 + 1]) +
              static_cast<double>(loop_ub) * c_Ry_tmp[3 * i2 + 2];
        }
        dv[3 * i1] = b_iv[i1];
      }
      dv[1] = 0.0;
      dv[4] = b_Rx_tmp;
      dv[7] = -Rx_tmp;
      dv[2] = 0.0;
      dv[5] = Rx_tmp;
      dv[8] = b_Rx_tmp;
      g_Rz_tmp[0] = b_Rz_tmp;
      g_Rz_tmp[4] = -Rz_tmp;
      g_Rz_tmp[8] = 0.0;
      g_Rz_tmp[12] = 0.0;
      g_Rz_tmp[1] = Rz_tmp;
      g_Rz_tmp[5] = b_Rz_tmp;
      g_Rz_tmp[9] = 0.0;
      g_Rz_tmp[13] = 0.0;
      g_Rz_tmp[2] = 0.0;
      g_Rz_tmp[3] = 0.0;
      g_Rz_tmp[6] = 0.0;
      g_Rz_tmp[7] = 0.0;
      g_Rz_tmp[10] = 1.0;
      g_Rz_tmp[11] = 0.0;
      g_Rz_tmp[14] = 0.0;
      g_Rz_tmp[15] = 1.0;
      // 'recursiveNEurdf:43' R_T(:,:,i) = R(:,:,i)';
      // 'recursiveNEurdf:45' z(:,i) = R_T(:,:,i)*z0;
      // 'recursiveNEurdf:47' r_link(:,i) = xyz(:,i);
      for (i1 = 0; i1 < 3; i1++) {
        d = f_Rz_tmp[i1];
        d1 = f_Rz_tmp[i1 + 3];
        Ry_tmp = f_Rz_tmp[i1 + 6];
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          e_Rz_tmp[i1 + 3 * loop_ub] =
              (d * dv[3 * loop_ub] + d1 * dv[3 * loop_ub + 1]) +
              Ry_tmp * dv[3 * loop_ub + 2];
        }
        d = e_Rz_tmp[i1];
        d1 = e_Rz_tmp[i1 + 3];
        Ry_tmp = e_Rz_tmp[i1 + 6];
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          i2 = loop_ub << 2;
          b_Ry_tmp = (d * g_Rz_tmp[i2] + d1 * g_Rz_tmp[i2 + 1]) +
                     Ry_tmp * g_Rz_tmp[i2 + 2];
          R[(i1 + 3 * loop_ub) + 9 * b_i] = b_Ry_tmp;
          R_T[(loop_ub + 3 * i1) + 9 * b_i] = b_Ry_tmp;
        }
        r_link[i1 + 3 * b_i] = KinPar_xyz[i1 + 3 * b_i];
      }
    } else {
      double f_Rz_tmp[9];
      double Rx_tmp;
      double Rz_tmp;
      double b_Rx_tmp;
      double b_Rz_tmp;
      // 'recursiveNEurdf:49' else
      //  last link
      // 'recursiveNEurdf:50' R(:,:,i) = rotMatURDF(rpy(:,i));
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
      b_Rx_tmp = std::cos(d);
      // 'rotMatURDF:10' Ry = [cos(p), 0, -sin(p); 0, 1, 0; sin(p), 0, cos(p)];
      d = KinPar_rpy[3 * b_i + 1];
      Ry_tmp = std::sin(d);
      b_Ry_tmp = std::cos(d);
      // 'rotMatURDF:11' Rz = [cos(y), -sin(y), 0; sin(y), cos(y), 0; 0, 0, 1];
      d = KinPar_rpy[3 * b_i + 2];
      Rz_tmp = std::sin(d);
      b_Rz_tmp = std::cos(d);
      // 'rotMatURDF:13' R = Rz * Ry * Rx;
      e_Rz_tmp[0] = b_Rz_tmp;
      e_Rz_tmp[3] = -Rz_tmp;
      e_Rz_tmp[6] = 0.0;
      e_Rz_tmp[1] = Rz_tmp;
      e_Rz_tmp[4] = b_Rz_tmp;
      e_Rz_tmp[7] = 0.0;
      c_Ry_tmp[0] = b_Ry_tmp;
      c_Ry_tmp[3] = 0.0;
      c_Ry_tmp[6] = -Ry_tmp;
      e_Rz_tmp[2] = 0.0;
      c_Ry_tmp[1] = 0.0;
      e_Rz_tmp[5] = 0.0;
      c_Ry_tmp[4] = 1.0;
      e_Rz_tmp[8] = 1.0;
      c_Ry_tmp[7] = 0.0;
      c_Ry_tmp[2] = Ry_tmp;
      c_Ry_tmp[5] = 0.0;
      c_Ry_tmp[8] = b_Ry_tmp;
      for (i1 = 0; i1 < 3; i1++) {
        d = e_Rz_tmp[i1];
        d1 = e_Rz_tmp[i1 + 3];
        loop_ub = static_cast<int>(e_Rz_tmp[i1 + 6]);
        for (int i2{0}; i2 < 3; i2++) {
          f_Rz_tmp[i1 + 3 * i2] =
              (d * c_Ry_tmp[3 * i2] + d1 * c_Ry_tmp[3 * i2 + 1]) +
              static_cast<double>(loop_ub) * c_Ry_tmp[3 * i2 + 2];
        }
        dv[3 * i1] = b_iv[i1];
      }
      dv[1] = 0.0;
      dv[4] = b_Rx_tmp;
      dv[7] = -Rx_tmp;
      dv[2] = 0.0;
      dv[5] = Rx_tmp;
      dv[8] = b_Rx_tmp;
      // 'recursiveNEurdf:51' R_T(:,:,i) = R(:,:,i)';
      // 'recursiveNEurdf:53' z(:,i) = R_T(:,:,i)*z0;
      // 'recursiveNEurdf:55' r_link(:,i) = xyz(:,i);
      for (i1 = 0; i1 < 3; i1++) {
        d = f_Rz_tmp[i1];
        d1 = f_Rz_tmp[i1 + 3];
        Ry_tmp = f_Rz_tmp[i1 + 6];
        for (loop_ub = 0; loop_ub < 3; loop_ub++) {
          b_Ry_tmp = (d * dv[3 * loop_ub] + d1 * dv[3 * loop_ub + 1]) +
                     Ry_tmp * dv[3 * loop_ub + 2];
          R[(i1 + 3 * loop_ub) + 9 * b_i] = b_Ry_tmp;
          R_T[(loop_ub + 3 * i1) + 9 * b_i] = b_Ry_tmp;
        }
        r_link[i1 + 3 * b_i] = KinPar_xyz[i1 + 3 * b_i];
      }
    }
    // 'recursiveNEurdf:58' z(:,i) = R_T(:,:,i)*z0;
  }
  // 'recursiveNEurdf:62' B = eye(3);
  // 'recursiveNEurdf:64' for i = 1:k
  b_i = -1;
  do {
    exitg2 = 0;
    if (b_i + 1 <= static_cast<int>(k)) {
      boolean_T guard1{false};
      //  initialization
      // 'recursiveNEurdf:66' if i == 1
      guard1 = false;
      if (static_cast<double>(b_i + 1) + 1.0 == 1.0) {
        // 'recursiveNEurdf:67' z(:,i) = z0;
        // 'recursiveNEurdf:68' a_j(:,i) = -B'*[g(1);g(2);g(3)];
        d = g[0];
        d1 = g[1];
        Ry_tmp = g[2];
        for (i = 0; i < 3; i++) {
          a_j[i] = (static_cast<double>(iv[i]) * d +
                    static_cast<double>(iv[i + 3]) * d1) +
                   static_cast<double>(iv[i + 6]) * Ry_tmp;
        }
        guard1 = true;
      } else {
        // 'recursiveNEurdf:70' else
        // 'recursiveNEurdf:71' if jointType(i-1) == 1
        d = DynPar_Jtype[static_cast<int>((static_cast<double>(b_i + 1) + 1.0) -
                                          1.0) -
                         1];
        if (d == 1.0) {
          //  revolute joint
          // 'recursiveNEurdf:72' omega(:,i) = R_T(:,:,i-1)*(omega(:,i-1)) +
          // dq(i-1)*z0;
          d_Rz_tmp = dq[b_i];
          for (i = 0; i < 3; i++) {
            y[i] = ((R_T[i + 9 * b_i] * omega[3 * b_i] +
                     R_T[(i + 9 * b_i) + 3] * omega[3 * b_i + 1]) +
                    R_T[(i + 9 * b_i) + 6] * omega[3 * b_i + 2]) +
                   d_Rz_tmp * static_cast<double>(b[i]);
          }
          omega[3 * (b_i + 1)] = y[0];
          omega[3 * (b_i + 1) + 1] = y[1];
          omega[3 * (b_i + 1) + 2] = y[2];
          //  dq(i-1) because it starts from i = 1
          // 'recursiveNEurdf:73' omegaDot(:,i) = R_T(:,:,i-1)*(omegaDot(:,i-1))
          // + Xcross(R_T(:,:,i-1)*omega(:,i-1))*dq(i-1)*z0 + ddq(i-1)*z0;
          for (i = 0; i < 3; i++) {
            y[i] = (R_T[i + 9 * b_i] * omega[3 * b_i] +
                    R_T[(i + 9 * b_i) + 3] * omega[3 * b_i + 1]) +
                   R_T[(i + 9 * b_i) + 6] * omega[3 * b_i + 2];
          }
          // UNTITLED7 Summary of this function goes here
          //    Detailed explanation goes here
          // 'Xcross:5' XMat = [0 -x(3) x(2);
          // 'Xcross:6'         x(3) 0 -x(1);
          // 'Xcross:7'         -x(2) x(1) 0];
          d_Rz_tmp = dq[b_i];
          dv[0] = 0.0;
          dv[3] = -y[2] * d_Rz_tmp;
          dv[6] = y[1] * d_Rz_tmp;
          dv[1] = y[2] * d_Rz_tmp;
          dv[4] = 0.0;
          dv[7] = -y[0] * d_Rz_tmp;
          dv[2] = -y[1] * d_Rz_tmp;
          dv[5] = y[0] * d_Rz_tmp;
          dv[8] = 0.0;
          for (i = 0; i < 3; i++) {
            dv1[i] = dv[i + 6];
            y[i] = (R_T[i + 9 * b_i] * omegaDot[3 * b_i] +
                    R_T[(i + 9 * b_i) + 3] * omegaDot[3 * b_i + 1]) +
                   R_T[(i + 9 * b_i) + 6] * omegaDot[3 * b_i + 2];
          }
          omegaDot[3 * (b_i + 1)] = y[0] + dv1[0];
          omegaDot[3 * (b_i + 1) + 1] = y[1] + dv1[1];
          omegaDot[3 * (b_i + 1) + 2] = (y[2] + dv1[2]) + ddq[b_i];
          // 'recursiveNEurdf:74' a_j(:,i) = R_T(:,:,i-1)*(a_j(:,i-1) +
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
          dv[3] = -omegaDot[3 * b_i + 2];
          dv[6] = omegaDot[3 * b_i + 1];
          dv[1] = omegaDot[3 * b_i + 2];
          dv[4] = 0.0;
          dv[7] = -omegaDot[3 * b_i];
          dv[2] = -omegaDot[3 * b_i + 1];
          dv[5] = omegaDot[3 * b_i];
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
            d = 0.0;
            d1 = 0.0;
            Ry_tmp = c_Ry_tmp[i];
            b_Ry_tmp = c_Ry_tmp[i + 3];
            d_Rz_tmp = c_Ry_tmp[i + 6];
            for (i1 = 0; i1 < 3; i1++) {
              d += dv[i + 3 * i1] * r_link[i1 + 3 * b_i];
              d1 += ((Ry_tmp * e_Rz_tmp[3 * i1] +
                      b_Ry_tmp * e_Rz_tmp[3 * i1 + 1]) +
                     d_Rz_tmp * e_Rz_tmp[3 * i1 + 2]) *
                    r_link[i1 + 3 * b_i];
            }
            y[i] = (a_j[i + 3 * b_i] + d) + d1;
          }
          d = y[0];
          d1 = y[1];
          Ry_tmp = y[2];
          for (i = 0; i < 3; i++) {
            a_j[i + 3 * (b_i + 1)] = 0.0;
            a_j[i + 3 * (b_i + 1)] =
                a_j[i + 3 * (b_i + 1)] + R_T[i + 9 * b_i] * d;
            a_j[i + 3 * (b_i + 1)] =
                a_j[i + 3 * (b_i + 1)] + R_T[(i + 9 * b_i) + 3] * d1;
            a_j[i + 3 * (b_i + 1)] =
                a_j[i + 3 * (b_i + 1)] + R_T[(i + 9 * b_i) + 6] * Ry_tmp;
          }
          guard1 = true;
        } else if (d == 0.0) {
          // 'recursiveNEurdf:75' elseif jointType(i-1) == 0
          //  prismatic joint
          //  not implemented
          // 'recursiveNEurdf:77' u_out = zeros(1,k-1);
          u_out.set_size(1, static_cast<int>(k));
          loop_ub = static_cast<int>(k);
          for (i = 0; i < loop_ub; i++) {
            u_out[i] = 0.0;
          }
          //  ones(1,k-1)*NaN;
          exitg2 = 1;
        } else {
          guard1 = true;
        }
      }
      if (guard1) {
        // 'recursiveNEurdf:81' a_c(:,i) = a_j(:,i) +
        // Xcross(omegaDot(:,i))*r(:,i) +
        // Xcross(omega(:,i))*Xcross(omega(:,i))*r(:,i); UNTITLED7 Summary of
        // this function goes here
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
        d = omegaDot[3 * (b_i + 1) + 2];
        dv[3] = -d;
        d1 = omegaDot[3 * (b_i + 1) + 1];
        dv[6] = d1;
        dv[1] = d;
        dv[4] = 0.0;
        d = omegaDot[3 * (b_i + 1)];
        dv[7] = -d;
        dv[2] = -d1;
        dv[5] = d;
        dv[8] = 0.0;
        c_Ry_tmp[0] = 0.0;
        d_Rz_tmp = omega[3 * (b_i + 1) + 2];
        c_Ry_tmp[3] = -d_Rz_tmp;
        c_Rz_tmp = omega[3 * (b_i + 1) + 1];
        c_Ry_tmp[6] = c_Rz_tmp;
        c_Ry_tmp[1] = d_Rz_tmp;
        c_Ry_tmp[4] = 0.0;
        Ry_tmp = omega[3 * (b_i + 1)];
        c_Ry_tmp[7] = -Ry_tmp;
        c_Ry_tmp[2] = -c_Rz_tmp;
        c_Ry_tmp[5] = Ry_tmp;
        c_Ry_tmp[8] = 0.0;
        e_Rz_tmp[0] = 0.0;
        e_Rz_tmp[3] = -d_Rz_tmp;
        e_Rz_tmp[6] = c_Rz_tmp;
        e_Rz_tmp[1] = d_Rz_tmp;
        e_Rz_tmp[4] = 0.0;
        e_Rz_tmp[7] = -Ry_tmp;
        e_Rz_tmp[2] = -c_Rz_tmp;
        e_Rz_tmp[5] = Ry_tmp;
        e_Rz_tmp[8] = 0.0;
        for (i = 0; i < 3; i++) {
          d = 0.0;
          d1 = 0.0;
          Ry_tmp = c_Ry_tmp[i];
          b_Ry_tmp = c_Ry_tmp[i + 3];
          d_Rz_tmp = c_Ry_tmp[i + 6];
          for (i1 = 0; i1 < 3; i1++) {
            c_Rz_tmp = DynPar_r[i1 + 3 * (b_i + 1)];
            d += dv[i + 3 * i1] * c_Rz_tmp;
            d1 +=
                ((Ry_tmp * e_Rz_tmp[3 * i1] + b_Ry_tmp * e_Rz_tmp[3 * i1 + 1]) +
                 d_Rz_tmp * e_Rz_tmp[3 * i1 + 2]) *
                c_Rz_tmp;
          }
          a_c[i + 3 * (b_i + 1)] = (a_j[i + 3 * (b_i + 1)] + d) + d1;
        }
        b_i++;
      }
    } else {
      //  backward recursion
      // 'recursiveNEurdf:86' for i = k:-1:2
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
        double b_R[3];
        double c_R[3];
        c_Rz_tmp = static_cast<double>(k + 1U) + -static_cast<double>(b_i);
        //  k-1:-1:1
        // 'recursiveNEurdf:88' if i == k
        if (c_Rz_tmp == k + 1U) {
          // k-1 % force and torques at the end effector
          // 'recursiveNEurdf:89' f(:,i+1) = zeros(3,1);
          // 'recursiveNEurdf:90' mu(:,i+1) = zeros(3,1);
          // 'recursiveNEurdf:91' R(:,:,i+1) = eye(3);
          loop_ub = static_cast<int>(c_Rz_tmp);
          for (i = 0; i < 3; i++) {
            f[i + 3 * loop_ub] = 0.0;
            mu[i + 3 * loop_ub] = 0.0;
            R[3 * i + 9 * loop_ub] = 0.0;
            R[(3 * i + 9 * loop_ub) + 1] = 0.0;
            R[(3 * i + 9 * loop_ub) + 2] = 0.0;
          }
          loop_ub = static_cast<int>(c_Rz_tmp);
          R[9 * loop_ub] = 1.0;
          R[9 * loop_ub + 4] = 1.0;
          R[9 * loop_ub + 8] = 1.0;
        }
        // 'recursiveNEurdf:94' Fi = m(i)*a_c(:,i);
        // 'recursiveNEurdf:95' f(:,i) = R(:,:,i)*f(:,i+1) + Fi;
        d_Rz_tmp = DynPar_m[static_cast<int>(c_Rz_tmp) - 1];
        for (i = 0; i < 3; i++) {
          b_R[i] = ((R[i + 9 * (static_cast<int>(c_Rz_tmp) - 1)] *
                         f[3 * static_cast<int>(c_Rz_tmp)] +
                     R[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 3] *
                         f[3 * static_cast<int>(c_Rz_tmp) + 1]) +
                    R[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 6] *
                        f[3 * static_cast<int>(c_Rz_tmp) + 2]) +
                   d_Rz_tmp * a_c[i + 3 * (static_cast<int>(c_Rz_tmp) - 1)];
        }
        f[3 * (static_cast<int>(c_Rz_tmp) - 1)] = b_R[0];
        f[3 * (static_cast<int>(c_Rz_tmp) - 1) + 1] = b_R[1];
        f[3 * (static_cast<int>(c_Rz_tmp) - 1) + 2] = b_R[2];
        // 'recursiveNEurdf:97' mu(:,i) = R(:,:,i)*mu(:,i+1) +
        // Xcross(R(:,:,i)*f(:,i+1))*(r(:,i)-r_link(:,i)) -
        // Xcross(f(:,i))*r(:,i) + I(:,:,i)*omegaDot(:,i) +
        // Xcross(omega(:,i))*I(:,:,i)*omega(:,i);
        for (i = 0; i < 3; i++) {
          y[i] = (R[i + 9 * (static_cast<int>(c_Rz_tmp) - 1)] *
                      f[3 * static_cast<int>(c_Rz_tmp)] +
                  R[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 3] *
                      f[3 * static_cast<int>(c_Rz_tmp) + 1]) +
                 R[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 6] *
                     f[3 * static_cast<int>(c_Rz_tmp) + 2];
        }
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
        dv[3] = -y[2];
        dv[6] = y[1];
        dv[1] = y[2];
        dv[4] = 0.0;
        dv[7] = -y[0];
        dv[2] = -y[1];
        dv[5] = y[0];
        dv[8] = 0.0;
        for (i = 0; i < 3; i++) {
          y[i] = DynPar_r[i + 3 * (static_cast<int>(c_Rz_tmp) - 1)] -
                 r_link[i + 3 * (static_cast<int>(c_Rz_tmp) - 1)];
          b_R[i] = (R[i + 9 * (static_cast<int>(c_Rz_tmp) - 1)] *
                        mu[3 * static_cast<int>(c_Rz_tmp)] +
                    R[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 3] *
                        mu[3 * static_cast<int>(c_Rz_tmp) + 1]) +
                   R[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 6] *
                       mu[3 * static_cast<int>(c_Rz_tmp) + 2];
        }
        d = y[0];
        d1 = y[1];
        Ry_tmp = y[2];
        for (i = 0; i < 3; i++) {
          dv1[i] = (dv[i] * d + dv[i + 3] * d1) + dv[i + 6] * Ry_tmp;
        }
        dv[0] = 0.0;
        d = f[3 * (static_cast<int>(c_Rz_tmp) - 1) + 2];
        dv[3] = -d;
        d1 = f[3 * (static_cast<int>(c_Rz_tmp) - 1) + 1];
        dv[6] = d1;
        dv[1] = d;
        dv[4] = 0.0;
        d = f[3 * (static_cast<int>(c_Rz_tmp) - 1)];
        dv[7] = -d;
        dv[2] = -d1;
        dv[5] = d;
        dv[8] = 0.0;
        for (i = 0; i < 3; i++) {
          y[i] = (DynPar_I[i + 9 * (static_cast<int>(c_Rz_tmp) - 1)] *
                      omegaDot[3 * (static_cast<int>(c_Rz_tmp) - 1)] +
                  DynPar_I[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 3] *
                      omegaDot[3 * (static_cast<int>(c_Rz_tmp) - 1) + 1]) +
                 DynPar_I[(i + 9 * (static_cast<int>(c_Rz_tmp) - 1)) + 6] *
                     omegaDot[3 * (static_cast<int>(c_Rz_tmp) - 1) + 2];
          c_R[i] =
              (b_R[i] + dv1[i]) -
              ((dv[i] * DynPar_r[3 * (static_cast<int>(c_Rz_tmp) - 1)] +
                dv[i + 3] *
                    DynPar_r[3 * (static_cast<int>(c_Rz_tmp) - 1) + 1]) +
               dv[i + 6] * DynPar_r[3 * (static_cast<int>(c_Rz_tmp) - 1) + 2]);
        }
        dv[0] = 0.0;
        d = omega[3 * (static_cast<int>(c_Rz_tmp) - 1) + 2];
        dv[3] = -d;
        d1 = omega[3 * (static_cast<int>(c_Rz_tmp) - 1) + 1];
        dv[6] = d1;
        dv[1] = d;
        dv[4] = 0.0;
        d = omega[3 * (static_cast<int>(c_Rz_tmp) - 1)];
        dv[7] = -d;
        dv[2] = -d1;
        dv[5] = d;
        dv[8] = 0.0;
        for (i = 0; i < 3; i++) {
          d = 0.0;
          d1 = dv[i];
          Ry_tmp = dv[i + 3];
          b_Ry_tmp = dv[i + 6];
          for (i1 = 0; i1 < 3; i1++) {
            d +=
                ((d1 * DynPar_I[3 * i1 + 9 * (static_cast<int>(c_Rz_tmp) - 1)] +
                  Ry_tmp *
                      DynPar_I[(3 * i1 + 9 * (static_cast<int>(c_Rz_tmp) - 1)) +
                               1]) +
                 b_Ry_tmp *
                     DynPar_I[(3 * i1 + 9 * (static_cast<int>(c_Rz_tmp) - 1)) +
                              2]) *
                omega[i1 + 3 * (static_cast<int>(c_Rz_tmp) - 1)];
          }
          mu[i + 3 * (static_cast<int>(c_Rz_tmp) - 1)] = (c_R[i] + y[i]) + d;
        }
        // 'recursiveNEurdf:98' if jointType(i-1) == 1
        d = DynPar_Jtype[static_cast<int>(c_Rz_tmp) - 2];
        if (d == 1.0) {
          //  revolute joint
          // 'recursiveNEurdf:99' u(i) = mu(:,i)'*z0 +
          // Fv(i-1)*dq(i-1)+Fc(i-1)*sin(atan(1000*dq(i-1))) +
          // ddq(i-1)*kr(i-1)^2*Im(i-1);
          d_Rz_tmp = DynPar_kr[static_cast<int>(c_Rz_tmp) - 2];
          u[static_cast<int>(c_Rz_tmp) - 1] =
              (mu[3 * (static_cast<int>(c_Rz_tmp) - 1) + 2] +
               DynPar_Fv[static_cast<int>(c_Rz_tmp) - 2] *
                   dq[static_cast<int>(c_Rz_tmp) - 2]) +
              ddq[static_cast<int>(c_Rz_tmp) - 2] * (d_Rz_tmp * d_Rz_tmp) *
                  DynPar_Im[static_cast<int>(c_Rz_tmp) - 2];
          //  + friq(dq[i]) + sigma[i]**2*Im[i,:,:]*ddq[i]
          b_i++;
        } else if (d == 0.0) {
          // 'recursiveNEurdf:100' elseif jointType(i-1) == 0
          //  prismatic joint
          //  not implemented
          // 'recursiveNEurdf:102' u_out = zeros(1,k-1);
          u_out.set_size(1, static_cast<int>(k));
          loop_ub = static_cast<int>(k);
          for (i = 0; i < loop_ub; i++) {
            u_out[i] = 0.0;
          }
          // ones(1,k-1)*NaN;
          exitg1 = 1;
        } else {
          b_i++;
        }
      } else {
        // 'recursiveNEurdf:107' u_out = u(2:end);
        if (u.size(1) < 2) {
          i = 0;
          i1 = 0;
        } else {
          i = 1;
          i1 = u.size(1);
        }
        loop_ub = i1 - i;
        u_out.set_size(1, loop_ub);
        for (i1 = 0; i1 < loop_ub; i1++) {
          u_out[i1] = u[i + i1];
        }
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
}

//
// File trailer for recursiveNEurdf.cpp
//
// [EOF]
//
