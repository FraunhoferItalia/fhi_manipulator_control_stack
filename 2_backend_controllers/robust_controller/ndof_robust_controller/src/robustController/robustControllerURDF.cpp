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
#include "robustControllerURDF.h"
#include "recursiveNEurdf.h"
#include "recursiveNEurdfIA.h"
#include "robustControllerURDF_types.h"
#include "coder_array.h"
#include <cmath>

// Function Declarations
static void binary_expand_op(coder::array<double, 1U> &in1,
                             const coder::array<double, 1U> &in2, double in3,
                             const coder::array<double, 1U> &in4, double in5,
                             const coder::array<double, 1U> &in6);

static void binary_expand_op(coder::array<double, 1U> &in1,
                             const coder::array<double, 1U> &in2,
                             const coder::array<double, 1U> &in3);

static void minus(coder::array<double, 1U> &in1,
                  const coder::array<double, 1U> &in2,
                  const coder::array<double, 1U> &in3);

// Function Definitions
//
// Arguments    : coder::array<double, 1U> &in1
//                const coder::array<double, 1U> &in2
//                double in3
//                const coder::array<double, 1U> &in4
//                double in5
//                const coder::array<double, 1U> &in6
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 1U> &in1,
                             const coder::array<double, 1U> &in2, double in3,
                             const coder::array<double, 1U> &in4, double in5,
                             const coder::array<double, 1U> &in6)
{
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in6.size(0) == 1) {
    if (in4.size(0) == 1) {
      i = in2.size(0);
    } else {
      i = in4.size(0);
    }
  } else {
    i = in6.size(0);
  }
  in1.set_size(i);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in4.size(0) != 1);
  stride_2_0 = (in6.size(0) != 1);
  if (in6.size(0) == 1) {
    if (in4.size(0) == 1) {
      loop_ub = in2.size(0);
    } else {
      loop_ub = in4.size(0);
    }
  } else {
    loop_ub = in6.size(0);
  }
  for (i = 0; i < loop_ub; i++) {
    in1[i] = (in2[i * stride_0_0] + in3 * in4[i * stride_1_0]) +
             in5 * in6[i * stride_2_0];
  }
}

//
// Arguments    : coder::array<double, 1U> &in1
//                const coder::array<double, 1U> &in2
//                const coder::array<double, 1U> &in3
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 1U> &in1,
                             const coder::array<double, 1U> &in2,
                             const coder::array<double, 1U> &in3)
{
  coder::array<double, 1U> b_in1;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in3.size(0) == 1) {
    if (in2.size(0) == 1) {
      i = in1.size(0);
    } else {
      i = in2.size(0);
    }
  } else {
    i = in3.size(0);
  }
  b_in1.set_size(i);
  stride_0_0 = (in1.size(0) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_2_0 = (in3.size(0) != 1);
  if (in3.size(0) == 1) {
    if (in2.size(0) == 1) {
      loop_ub = in1.size(0);
    } else {
      loop_ub = in2.size(0);
    }
  } else {
    loop_ub = in3.size(0);
  }
  for (i = 0; i < loop_ub; i++) {
    b_in1[i] =
        (in1[i * stride_0_0] + in2[i * stride_1_0]) + in3[i * stride_2_0];
  }
  in1.set_size(b_in1.size(0));
  loop_ub = b_in1.size(0);
  for (i = 0; i < loop_ub; i++) {
    in1[i] = b_in1[i];
  }
}

//
// Arguments    : coder::array<double, 1U> &in1
//                const coder::array<double, 1U> &in2
//                const coder::array<double, 1U> &in3
// Return Type  : void
//
static void minus(coder::array<double, 1U> &in1,
                  const coder::array<double, 1U> &in2,
                  const coder::array<double, 1U> &in3)
{
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in3.size(0);
  }
  in1.set_size(i);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  for (i = 0; i < loop_ub; i++) {
    in1[i] = in2[i * stride_0_0] - in3[i * stride_1_0];
  }
}

//
// function [tau, normRho, tauRobust, tauGravity] = robustControllerURDF(qRef,
// dqRef, ddqRef, q, dq, g, DynPar, KinPar, Kd, Kp, B, P, DynPar_inf,
// DynPar_sup, sampleT, gainTuneContr)
//
// Arguments    : const coder::array<double, 1U> &qRef
//                const coder::array<double, 1U> &dqRef
//                const coder::array<double, 1U> &ddqRef
//                const coder::array<double, 1U> &q
//                const coder::array<double, 1U> &dq
//                const double g[3]
//                const struct0_T *DynPar
//                const struct1_T *KinPar
//                double Kd
//                double Kp
//                const coder::array<double, 2U> &B
//                const coder::array<double, 2U> &P
//                const struct2_T *DynPar_inf
//                struct2_T *DynPar_sup
//                double sampleT
//                double gainTuneContr
//                coder::array<double, 1U> &tau
//                double *normRho
//                coder::array<double, 1U> &tauRobust
//                coder::array<double, 1U> &tauGravity
// Return Type  : void
//
void robustControllerURDF(
    const coder::array<double, 1U> &qRef, const coder::array<double, 1U> &dqRef,
    const coder::array<double, 1U> &ddqRef, const coder::array<double, 1U> &q,
    const coder::array<double, 1U> &dq, const double g[3],
    const struct0_T *DynPar, const struct1_T *KinPar, double Kd, double Kp,
    const coder::array<double, 2U> &B, const coder::array<double, 2U> &P,
    const struct2_T *DynPar_inf, struct2_T *DynPar_sup, double,
    double gainTuneContr, coder::array<double, 1U> &tau, double *normRho,
    coder::array<double, 1U> &tauRobust, coder::array<double, 1U> &tauGravity)
{
  coder::array<double, 2U> FcFW;
  coder::array<double, 2U> b_dq;
  coder::array<double, 2U> b_q;
  coder::array<double, 2U> b_y;
  coder::array<double, 2U> c_y;
  coder::array<double, 2U> result;
  coder::array<double, 2U> rho;
  coder::array<double, 2U> tau_sup;
  coder::array<double, 2U> varargin_1;
  coder::array<double, 1U> E;
  coder::array<double, 1U> dE;
  coder::array<double, 1U> y;
  double d;
  double scale;
  int i;
  int inner;
  int mc;
  int nx;
  // 'robustControllerURDF:4' FcFW = DynPar.Fc;
  FcFW.set_size(1, DynPar->Fc.size(1));
  nx = DynPar->Fc.size(1);
  for (i = 0; i < nx; i++) {
    FcFW[i] = DynPar->Fc[i];
  }
  unsigned int unnamed_idx_1;
  //  for feedforward
  //  to remove friction model in controller
  // 'robustControllerURDF:6' DynPar.Fc = zeros(size(DynPar.Fc));
  // 'robustControllerURDF:7' DynPar_inf.Fc = zeros(size(DynPar_inf.Fc));
  // 'robustControllerURDF:8' DynPar_sup.Fc = zeros(size(DynPar_sup.Fc));
  unnamed_idx_1 = static_cast<unsigned int>(DynPar_sup->Fc.size(1));
  DynPar_sup->Fc.set_size(1, DynPar_sup->Fc.size(1));
  nx = static_cast<int>(unnamed_idx_1);
  for (i = 0; i < nx; i++) {
    DynPar_sup->Fc[i] = 0.0;
  }
  //  ----
  //  FOR CONTROLLER PART 1
  // 'robustControllerURDF:14' dE = dqRef - dq;
  if (dqRef.size(0) == dq.size(0)) {
    dE.set_size(dqRef.size(0));
    nx = dqRef.size(0);
    for (i = 0; i < nx; i++) {
      dE[i] = dqRef[i] - dq[i];
    }
  } else {
    minus(dE, dqRef, dq);
  }
  // 'robustControllerURDF:15' E = qRef - q;
  if (qRef.size(0) == q.size(0)) {
    E.set_size(qRef.size(0));
    nx = qRef.size(0);
    for (i = 0; i < nx; i++) {
      E[i] = qRef[i] - q[i];
    }
  } else {
    minus(E, qRef, q);
  }
  // 'robustControllerURDF:17' y = ddqRef + Kd*dE + Kp*E;
  if (ddqRef.size(0) == 1) {
    i = dE.size(0);
  } else {
    i = ddqRef.size(0);
  }
  if ((ddqRef.size(0) == dE.size(0)) && (i == E.size(0))) {
    y.set_size(ddqRef.size(0));
    nx = ddqRef.size(0);
    for (i = 0; i < nx; i++) {
      y[i] = (ddqRef[i] + Kd * dE[i]) + Kp * E[i];
    }
  } else {
    binary_expand_op(y, ddqRef, Kd, dE, Kp, E);
  }
  // 'robustControllerURDF:19' tauKP =
  // recursiveNEurdf(q',dq',y',KinPar,DynPar,g)';
  b_q.set_size(1, q.size(0));
  nx = q.size(0);
  for (i = 0; i < nx; i++) {
    b_q[i] = q[i];
  }
  b_dq.set_size(1, dq.size(0));
  nx = dq.size(0);
  for (i = 0; i < nx; i++) {
    b_dq[i] = dq[i];
  }
  b_y.set_size(1, y.size(0));
  nx = y.size(0);
  for (i = 0; i < nx; i++) {
    b_y[i] = y[i];
  }
  recursiveNEurdf(b_q, b_dq, b_y, KinPar->xyz, KinPar->rpy, DynPar->m,
                  DynPar->r, DynPar->kr, DynPar->Im, DynPar->Fv, DynPar->Jtype,
                  DynPar->b_I, g, rho);
  tau.set_size(rho.size(1));
  nx = rho.size(1);
  for (i = 0; i < nx; i++) {
    tau[i] = rho[i];
  }
  //  transposed
  // 'robustControllerURDF:20' [tau_inf, tau_sup] =
  // recursiveNEurdfIA(q',dq',y',KinPar,DynPar_inf,DynPar_sup,g);
  b_q.set_size(1, q.size(0));
  nx = q.size(0);
  for (i = 0; i < nx; i++) {
    b_q[i] = q[i];
  }
  b_dq.set_size(1, dq.size(0));
  nx = dq.size(0);
  for (i = 0; i < nx; i++) {
    b_dq[i] = dq[i];
  }
  b_y.set_size(1, y.size(0));
  nx = y.size(0);
  for (i = 0; i < nx; i++) {
    b_y[i] = y[i];
  }
  recursiveNEurdfIA(b_q, b_dq, b_y, KinPar->xyz, KinPar->rpy, DynPar_inf->Jtype,
                    DynPar_inf->kr, DynPar_inf->Im, DynPar_inf->m,
                    DynPar_inf->b_I, DynPar_inf->Fv, DynPar_inf->r, DynPar_sup,
                    g, rho, tau_sup);
  // 'robustControllerURDF:21' [phi_inf, phi_sup] =
  // subMatMatInt(tau_inf,tau_sup,tauKP',tauKP'); UNTITLED3 Summary of this
  // function goes here
  //    Detailed explanation goes here
  // 'subMatMatInt:5' CInf = zeros(size(AInf));
  b_dq.set_size(1, rho.size(1));
  nx = rho.size(1);
  for (i = 0; i < nx; i++) {
    b_dq[i] = 0.0;
  }
  // 'subMatMatInt:6' CSup = zeros(size(AInf));
  b_y.set_size(1, rho.size(1));
  nx = rho.size(1);
  for (i = 0; i < nx; i++) {
    b_y[i] = 0.0;
  }
  // 'subMatMatInt:8' for i = 1:size(AInf,1)
  // 'subMatMatInt:9' for j = 1:size(AInf,2)
  i = rho.size(1);
  for (int j{0}; j < i; j++) {
    // 'subMatMatInt:10' CInf(i,j) = AInf(i,j) - BSup(i,j);
    b_dq[j] = rho[j] - tau[j];
    // 'subMatMatInt:11' CSup(i,j) = ASup(i,j) - BInf(i,j);
    b_y[j] = tau_sup[j] - tau[j];
  }
  // 'robustControllerURDF:22' rho = max(abs([phi_inf; phi_sup]));
  result.set_size(2, b_dq.size(1));
  nx = b_dq.size(1);
  for (i = 0; i < nx; i++) {
    result[2 * i] = b_dq[i];
  }
  nx = b_y.size(1);
  for (i = 0; i < nx; i++) {
    result[2 * i + 1] = b_y[i];
  }
  nx = result.size(1) << 1;
  varargin_1.set_size(2, result.size(1));
  for (int k{0}; k < nx; k++) {
    varargin_1[k] = std::abs(result[k]);
  }
  nx = varargin_1.size(1);
  rho.set_size(1, varargin_1.size(1));
  if (varargin_1.size(1) >= 1) {
    for (int j{0}; j < nx; j++) {
      d = varargin_1[2 * j];
      rho[j] = d;
      scale = varargin_1[2 * j + 1];
      if (d < scale) {
        rho[j] = scale;
      }
    }
  }
  // 'robustControllerURDF:23' normRho = norm(rho);
  if (rho.size(1) == 0) {
    *normRho = 0.0;
  } else {
    *normRho = 0.0;
    if (rho.size(1) == 1) {
      *normRho = rho[0];
    } else {
      scale = 3.3121686421112381E-170;
      nx = rho.size(1);
      for (int k{0}; k < nx; k++) {
        d = rho[k];
        if (d > scale) {
          double t;
          t = scale / d;
          *normRho = *normRho * t * t + 1.0;
          scale = d;
        } else {
          double t;
          t = d / scale;
          *normRho += t * t;
        }
      }
      *normRho = scale * std::sqrt(*normRho);
    }
  }
  // 'robustControllerURDF:25' z = B'*P*[E; dE];
  mc = B.size(1);
  inner = B.size(0);
  nx = P.size(1);
  c_y.set_size(B.size(1), P.size(1));
  for (int j{0}; j < nx; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * P.size(0);
    for (int b_i{0}; b_i < mc; b_i++) {
      c_y[coffset + b_i] = 0.0;
    }
    for (int k{0}; k < inner; k++) {
      scale = P[boffset + k];
      for (int b_i{0}; b_i < mc; b_i++) {
        i = coffset + b_i;
        c_y[i] = c_y[i] + B[b_i * B.size(0) + k] * scale;
      }
    }
  }
  y.set_size(E.size(0) + dE.size(0));
  nx = E.size(0);
  for (i = 0; i < nx; i++) {
    y[i] = E[i];
  }
  nx = dE.size(0);
  for (i = 0; i < nx; i++) {
    y[i + E.size(0)] = dE[i];
  }
  mc = c_y.size(0) - 1;
  inner = c_y.size(1);
  tauRobust.set_size(c_y.size(0));
  for (int b_i{0}; b_i <= mc; b_i++) {
    tauRobust[b_i] = 0.0;
  }
  for (int k{0}; k < inner; k++) {
    nx = k * c_y.size(0);
    for (int b_i{0}; b_i <= mc; b_i++) {
      tauRobust[b_i] = tauRobust[b_i] + c_y[nx + b_i] * y[k];
    }
  }
  //  CONTROLLER PART 2
  // 'robustControllerURDF:29' tauRobust = gainTuneContr * z;
  nx = tauRobust.size(0);
  for (i = 0; i < nx; i++) {
    tauRobust[i] = gainTuneContr * tauRobust[i];
  }
  // 'robustControllerURDF:31' tau = tauKP + tauRobust;
  //  for friction feedforward (TO BE SURE TO HAVE IT REMOVED IN recursiveNE
  //  function)
  // 'robustControllerURDF:35' tauFriction = zeros(length(DynPar.Jtype),1);
  y.set_size(DynPar->Jtype.size(1));
  nx = DynPar->Jtype.size(1);
  for (i = 0; i < nx; i++) {
    y[i] = 0.0;
  }
  // 'robustControllerURDF:36' for i = 1:length(DynPar.Jtype)
  i = DynPar->Jtype.size(1);
  for (int b_i{0}; b_i < i; b_i++) {
    // 'robustControllerURDF:37' tauFriction(i) =
    // FcFW(i)*sin(atan(1000*dqRef(i)));
    y[b_i] = FcFW[b_i] * std::sin(std::atan(1000.0 * dqRef[b_i]));
  }
  // 'robustControllerURDF:40' tau = tau + tauFriction;
  if (tau.size(0) == 1) {
    i = tauRobust.size(0);
  } else {
    i = tau.size(0);
  }
  if ((tau.size(0) == tauRobust.size(0)) && (i == y.size(0))) {
    nx = tau.size(0);
    for (i = 0; i < nx; i++) {
      tau[i] = (tau[i] + tauRobust[i]) + y[i];
    }
  } else {
    binary_expand_op(tau, tauRobust, y);
  }
  // 'robustControllerURDF:42' tauGravity =
  // recursiveNEurdf(q',zeros(size(q))',zeros(size(q))',KinPar,DynPar,g)';
  b_q.set_size(1, q.size(0));
  nx = q.size(0);
  for (i = 0; i < nx; i++) {
    b_q[i] = q[i];
  }
  b_dq.set_size(1, q.size(0));
  nx = q.size(0);
  for (i = 0; i < nx; i++) {
    b_dq[i] = 0.0;
  }
  b_y.set_size(1, q.size(0));
  nx = q.size(0);
  for (i = 0; i < nx; i++) {
    b_y[i] = 0.0;
  }
  recursiveNEurdf(b_q, b_dq, b_y, KinPar->xyz, KinPar->rpy, DynPar->m,
                  DynPar->r, DynPar->kr, DynPar->Im, DynPar->Fv, DynPar->Jtype,
                  DynPar->b_I, g, rho);
  tauGravity.set_size(rho.size(1));
  nx = rho.size(1);
  for (i = 0; i < nx; i++) {
    tauGravity[i] = rho[i];
  }
}

//
// File trailer for robustControllerURDF.cpp
//
// [EOF]
//
