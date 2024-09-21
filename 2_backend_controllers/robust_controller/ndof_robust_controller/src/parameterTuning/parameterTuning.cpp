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
#include "parameterTuning.h"
#include "combineVectorElements.h"
#include "fft.h"
#include "linspace.h"
#include "parameterTuning_data.h"
#include "parameterTuning_initialize.h"
#include "coder_array.h"
#include <cmath>

// Variable Definitions
static boolean_T gainTune_not_empty;

static boolean_T GainTuneOld_not_empty;

static boolean_T dqMat_not_empty;

static boolean_T torqueBuffer_not_empty;

static double counter;

static boolean_T meanPSDgain_not_empty;

// Function Declarations
static void binary_expand_op(coder::array<double, 2U> &in1, int in2,
                             const coder::array<double, 2U> &in3,
                             const coder::array<double, 2U> &in4);

static double rt_hypotd(double u0, double u1);

// Function Definitions
//
// Arguments    : coder::array<double, 2U> &in1
//                int in2
//                const coder::array<double, 2U> &in3
//                const coder::array<double, 2U> &in4
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 2U> &in1, int in2,
                             const coder::array<double, 2U> &in3,
                             const coder::array<double, 2U> &in4)
{
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  i = in3.size(0);
  stride_0_0 = (i != 1);
  stride_1_0 = (in4.size(1) != 1);
  if (in4.size(1) == 1) {
    loop_ub = i;
  } else {
    loop_ub = in4.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    in1[i + in1.size(0) * in2] =
        in3[i * stride_0_0 + in3.size(0) * in2] * in4[i * stride_1_0];
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

//
// Arguments    : void
// Return Type  : void
//
void GainTuneOld_not_empty_init()
{
  GainTuneOld_not_empty = false;
}

//
// Arguments    : void
// Return Type  : void
//
void gainTune_not_empty_init()
{
  gainTune_not_empty = false;
}

//
// function [gainTuneContr, avgPSDdata]  = parameterTuning(dq, normRho,
// robTorque, sampleT, dqThreshold, winSize, winPerc, PSDMaxUnique,
// PSDMinUnique, kpUnique, gainMaxUnique, freqLim, lambda1)
//
// Arguments    : const coder::array<double, 1U> &dq
//                double normRho
//                const coder::array<double, 1U> &robTorque
//                double sampleT
//                double dqThreshold
//                double winSize
//                double winPerc
//                double PSDMaxUnique
//                double PSDMinUnique
//                double kpUnique
//                double gainMaxUnique
//                double freqLim
//                double lambda1
//                double *gainTuneContr
//                coder::array<double, 1U> &avgPSDdata
// Return Type  : void
//
void parameterTuning(const coder::array<double, 1U> &dq, double normRho,
                     const coder::array<double, 1U> &robTorque, double sampleT,
                     double dqThreshold, double winSize, double winPerc,
                     double PSDMaxUnique, double PSDMinUnique, double kpUnique,
                     double gainMaxUnique, double freqLim, double lambda1,
                     double *gainTuneContr,
                     coder::array<double, 1U> &avgPSDdata)
{
  static coder::array<double, 2U> dqMat;
  static coder::array<double, 2U> meanPSDgain;
  static coder::array<double, 2U> torqueBuffer;
  static double GainTuneOld;
  static double gainTune;
  coder::array<creal_T, 2U> b_x;
  coder::array<creal_T, 2U> datafft;
  coder::array<double, 2U> dataHam;
  coder::array<double, 2U> hamWin;
  coder::array<double, 2U> r;
  coder::array<double, 2U> x;
  coder::array<boolean_T, 2U> c_x;
  double Fs;
  double PSDmean;
  double a;
  double d;
  double varargin_1;
  int i2;
  int i3;
  int k;
  int loop_ub;
  int nx;
  short i1;
  if (!isInitialized_parameterTuning) {
    parameterTuning_initialize();
  }
  // 'parameterTuning:4' PSDmean = (PSDMaxUnique + PSDMinUnique) / 2;
  PSDmean = (PSDMaxUnique + PSDMinUnique) / 2.0;
  //  winPerc = 7/8; % winPercentage, where the counter restarts
  //  threshold for update gain
  //  dqThreshold = 0.1; % [rad/s]
  // 'parameterTuning:12' if isempty(gainTune)
  if (!gainTune_not_empty) {
    // 'parameterTuning:13' gainTune = kpUnique;
    gainTune = kpUnique;
    gainTune_not_empty = true;
  }
  //  for smoothing
  // 'parameterTuning:17' if isempty(GainTuneOld)
  if (!GainTuneOld_not_empty) {
    // 'parameterTuning:18' GainTuneOld = kpUnique;
    GainTuneOld = kpUnique;
    GainTuneOld_not_empty = true;
  }
  //  ----
  // 'parameterTuning:23' if isempty(dqMat)
  if (!dqMat_not_empty) {
    short i;
    // 'parameterTuning:24' dqMat = zeros(int16(winSize),length(dq));
    d = std::round(winSize);
    if (d < 32768.0) {
      if (d >= -32768.0) {
        i = static_cast<short>(d);
        i1 = i;
      } else {
        i = MIN_int16_T;
        i1 = MIN_int16_T;
      }
    } else {
      i = MAX_int16_T;
      i1 = MAX_int16_T;
    }
    dqMat.set_size(static_cast<int>(i), dq.size(0));
    loop_ub = i1 * dq.size(0);
    for (i2 = 0; i2 < loop_ub; i2++) {
      dqMat[i2] = 0.0;
    }
    dqMat_not_empty = ((dqMat.size(0) != 0) && (dqMat.size(1) != 0));
  }
  //  -------- bound parameters ------
  // 'parameterTuning:29' gainUpperBound = gainMaxUnique;
  //  maximum gain to guarantee stability wrt sampling frequency
  // 'parameterTuning:30' gainLowerBound = (1 + normRho);
  //  minimum gain to guarantee robustness
  //  gainLowerBound = 0; % to test simple PID controller
  //  ----------------
  // 'parameterTuning:35' dqMat(1:end-1) = dqMat(2:end);
  i2 = dqMat.size(0) * dqMat.size(1);
  if (i2 < 2) {
    i3 = 0;
    i2 = 0;
  } else {
    i3 = 1;
  }
  loop_ub = i2 - i3;
  x.set_size(1, loop_ub);
  for (i2 = 0; i2 < loop_ub; i2++) {
    x[i2] = dqMat[i3 + i2];
  }
  loop_ub = x.size(1);
  for (i2 = 0; i2 < loop_ub; i2++) {
    dqMat[i2] = x[i2];
  }
  // 'parameterTuning:36' dqMat(end,:) = dq';
  i2 = dqMat.size(0) - 1;
  loop_ub = dq.size(0);
  for (i3 = 0; i3 < loop_ub; i3++) {
    dqMat[i2 + dqMat.size(0) * i3] = dq[i3];
  }
  //  ------------ for tuning ---------------
  //  hamming window
  // 'parameterTuning:41' hamWin = 0.54 -
  // 0.46*cos(2*pi/(winSize-1)*linspace(0,winSize-1,int16(winSize)));
  a = 6.2831853071795862 / (winSize - 1.0);
  d = std::round(winSize);
  if (d < 32768.0) {
    if (d >= -32768.0) {
      i1 = static_cast<short>(d);
    } else {
      i1 = MIN_int16_T;
    }
  } else {
    i1 = MAX_int16_T;
  }
  coderTuner::linspace(winSize - 1.0, i1, hamWin);
  hamWin.set_size(1, hamWin.size(1));
  loop_ub = hamWin.size(1) - 1;
  for (i2 = 0; i2 <= loop_ub; i2++) {
    hamWin[i2] = a * hamWin[i2];
  }
  nx = hamWin.size(1);
  for (k = 0; k < nx; k++) {
    hamWin[k] = std::cos(hamWin[k]);
  }
  hamWin.set_size(1, hamWin.size(1));
  loop_ub = hamWin.size(1) - 1;
  for (i2 = 0; i2 <= loop_ub; i2++) {
    hamWin[i2] = 0.54 - 0.46 * hamWin[i2];
  }
  // 'parameterTuning:44' if isempty(torqueBuffer)
  if (!torqueBuffer_not_empty) {
    // 'parameterTuning:45' torqueBuffer = zeros(int16(winSize),length(dq));
    d = std::round(winSize);
    if (d < 32768.0) {
      if (d >= -32768.0) {
        i1 = static_cast<short>(d);
      } else {
        i1 = MIN_int16_T;
      }
    } else {
      i1 = MAX_int16_T;
    }
    torqueBuffer.set_size(static_cast<int>(i1), dq.size(0));
    d = std::round(winSize);
    if (d < 32768.0) {
      if (d >= -32768.0) {
        i1 = static_cast<short>(d);
      } else {
        i1 = MIN_int16_T;
      }
    } else {
      i1 = MAX_int16_T;
    }
    loop_ub = i1 * dq.size(0);
    for (i2 = 0; i2 < loop_ub; i2++) {
      torqueBuffer[i2] = 0.0;
    }
    torqueBuffer_not_empty =
        ((torqueBuffer.size(0) != 0) && (torqueBuffer.size(1) != 0));
  }
  // 'parameterTuning:49' if isempty(counter)
  // 'parameterTuning:53' torqueBuffer(1:end-1) = torqueBuffer(2:end);
  i2 = torqueBuffer.size(0) * torqueBuffer.size(1);
  if (i2 < 2) {
    i3 = 0;
    i2 = 0;
  } else {
    i3 = 1;
  }
  loop_ub = i2 - i3;
  x.set_size(1, loop_ub);
  for (i2 = 0; i2 < loop_ub; i2++) {
    x[i2] = torqueBuffer[i3 + i2];
  }
  loop_ub = x.size(1);
  for (i2 = 0; i2 < loop_ub; i2++) {
    torqueBuffer[i2] = x[i2];
  }
  // 'parameterTuning:55' torqueBuffer(end,:) = robTorque';
  i2 = torqueBuffer.size(0) - 1;
  loop_ub = robTorque.size(0);
  for (i3 = 0; i3 < loop_ub; i3++) {
    torqueBuffer[i2 + torqueBuffer.size(0) * i3] = robTorque[i3];
  }
  // 'parameterTuning:57' Fs = 1/sampleT;
  Fs = 1.0 / sampleT;
  //  sampling frequency
  // 'parameterTuning:60' if isempty(meanPSDgain)
  if (!meanPSDgain_not_empty) {
    // 'parameterTuning:61' meanPSDgain = zeros(1,length(dq));
    meanPSDgain.set_size(1, dq.size(0));
    loop_ub = dq.size(0);
    for (i2 = 0; i2 < loop_ub; i2++) {
      meanPSDgain[i2] = 0.0;
    }
    meanPSDgain_not_empty = (meanPSDgain.size(1) != 0);
  }
  // 'parameterTuning:64' if counter == winSize
  if (counter == winSize) {
    int i4;
    int i5;
    int i6;
    // 'parameterTuning:66' GainTuneOld = gainTune;
    GainTuneOld = gainTune;
    // 'parameterTuning:68' counter = winSize*winPerc;
    counter = winSize * winPerc;
    // 'parameterTuning:70' dataHam = zeros(int16(winSize),length(dq));
    d = std::round(winSize);
    if (d < 32768.0) {
      if (d >= -32768.0) {
        i1 = static_cast<short>(d);
      } else {
        i1 = MIN_int16_T;
      }
    } else {
      i1 = MAX_int16_T;
    }
    dataHam.set_size(static_cast<int>(i1), dq.size(0));
    d = std::round(winSize);
    if (d < 32768.0) {
      if (d >= -32768.0) {
        i1 = static_cast<short>(d);
      } else {
        i1 = MIN_int16_T;
      }
    } else {
      i1 = MAX_int16_T;
    }
    loop_ub = i1 * dq.size(0);
    for (i2 = 0; i2 < loop_ub; i2++) {
      dataHam[i2] = 0.0;
    }
    // 'parameterTuning:71' for i = 1:length(dq)
    i2 = dq.size(0);
    for (nx = 0; nx < i2; nx++) {
      // 'parameterTuning:72' dataHam(:,i) = torqueBuffer(:,i).*hamWin';
      loop_ub = torqueBuffer.size(0);
      if (torqueBuffer.size(0) == hamWin.size(1)) {
        for (i3 = 0; i3 < loop_ub; i3++) {
          dataHam[i3 + dataHam.size(0) * nx] =
              torqueBuffer[i3 + torqueBuffer.size(0) * nx] * hamWin[i3];
        }
      } else {
        binary_expand_op(dataHam, nx, torqueBuffer, hamWin);
      }
    }
    // 'parameterTuning:75' datafft = fft(dataHam);
    coderTuner::fft(dataHam, datafft);
    // 'parameterTuning:77' datafftHalf = datafft(1:winSize/2+1,:);
    d = winSize / 2.0 + 1.0;
    if (d < 1.0) {
      loop_ub = 0;
    } else {
      loop_ub = static_cast<int>(d);
    }
    // 'parameterTuning:79' psddata = (1/(Fs*winSize)) * abs(datafftHalf).^2;
    a = 1.0 / (Fs * winSize);
    nx = datafft.size(1);
    b_x.set_size(loop_ub, datafft.size(1));
    for (i2 = 0; i2 < nx; i2++) {
      for (i3 = 0; i3 < loop_ub; i3++) {
        b_x[i3 + b_x.size(0) * i2] = datafft[i3 + datafft.size(0) * i2];
      }
    }
    nx = loop_ub * datafft.size(1);
    dataHam.set_size(loop_ub, datafft.size(1));
    for (k = 0; k < nx; k++) {
      dataHam[k] = rt_hypotd(b_x[k].re, b_x[k].im);
    }
    loop_ub = dataHam.size(0) * dataHam.size(1);
    for (i2 = 0; i2 < loop_ub; i2++) {
      varargin_1 = dataHam[i2];
      dataHam[i2] = varargin_1 * varargin_1;
    }
    loop_ub = dataHam.size(0) * dataHam.size(1);
    for (i2 = 0; i2 < loop_ub; i2++) {
      dataHam[i2] = a * dataHam[i2];
    }
    // 'parameterTuning:81' psddata(2:end-1,:) = 2*psddata(2:end-1,:);
    if (dataHam.size(0) - 1 < 2) {
      i2 = 0;
      i3 = 0;
      i4 = 0;
    } else {
      i2 = 1;
      i3 = dataHam.size(0) - 1;
      i4 = 1;
    }
    nx = dataHam.size(1) - 1;
    loop_ub = i3 - i2;
    r.set_size(loop_ub, dataHam.size(1));
    for (i3 = 0; i3 <= nx; i3++) {
      for (i5 = 0; i5 < loop_ub; i5++) {
        r[i5 + r.size(0) * i3] =
            2.0 * dataHam[(i2 + i5) + dataHam.size(0) * i3];
      }
    }
    loop_ub = r.size(1);
    for (i2 = 0; i2 < loop_ub; i2++) {
      nx = r.size(0);
      for (i3 = 0; i3 < nx; i3++) {
        dataHam[(i4 + i3) + dataHam.size(0) * i2] = r[i3 + r.size(0) * i2];
      }
    }
    //  freq = 0:Fs/winSize:Fs/2;
    //  discard zero frequency PSD
    // 'parameterTuning:87' psddatatr = psddata(2:end,:);
    if (dataHam.size(0) < 2) {
      i2 = -1;
      i3 = -2;
    } else {
      i2 = 0;
      i3 = dataHam.size(0) - 2;
    }
    // freqtr = freq(2:end);
    //  compute ratio
    // 'parameterTuning:93' indexMax = ceil(freqLim/Fs*winSize)-1;
    varargin_1 = std::ceil(freqLim / Fs * winSize);
    //  -1 becasue we remove 0 frequency
    //      freq = 0:500/1024:500/2;
    //      figure
    //      plot(freq(indexMax+1:end),psdkp2tr(indexMax:end,:))
    //      hold on
    //      plot(freq(indexMax+1:end),psdphiPtr(indexMax:end,:))
    // 'parameterTuning:100' meanPSDgain =
    // sum(psddatatr(indexMax:end,:),1)./length(psddatatr(indexMax:end,:));
    i3 -= i2;
    if (varargin_1 - 1.0 > i3 + 1) {
      i4 = -1;
      i3 = -1;
      i5 = -1;
      i6 = -1;
    } else {
      i4 = static_cast<int>(varargin_1 - 1.0) - 2;
      i5 = static_cast<int>(varargin_1 - 1.0) - 2;
      i6 = i3;
    }
    nx = i6 - i5;
    if ((nx == 0) || (dataHam.size(1) == 0)) {
      nx = 0;
    } else if (nx <= dataHam.size(1)) {
      nx = dataHam.size(1);
    }
    loop_ub = dataHam.size(1) - 1;
    for (i5 = 0; i5 <= loop_ub; i5++) {
      k = i3 - i4;
      for (i6 = 0; i6 < k; i6++) {
        dataHam[i6 + k * i5] =
            dataHam[(((i2 + i4) + i6) + dataHam.size(0) * i5) + 2];
      }
    }
    dataHam.set_size(i3 - i4, loop_ub + 1);
    coderTuner::combineVectorElements(dataHam, meanPSDgain);
    meanPSDgain.set_size(1, meanPSDgain.size(1));
    loop_ub = meanPSDgain.size(1) - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      meanPSDgain[i2] = meanPSDgain[i2] / static_cast<double>(nx);
    }
    meanPSDgain_not_empty = (meanPSDgain.size(1) != 0);
    // 'parameterTuning:102' if max(mean(abs(dqMat))) > dqThreshold
    nx = dqMat.size(0) * dqMat.size(1);
    dataHam.set_size(dqMat.size(0), dqMat.size(1));
    for (k = 0; k < nx; k++) {
      dataHam[k] = std::abs(dqMat[k]);
    }
    coderTuner::combineVectorElements(dataHam, x);
    x.set_size(1, x.size(1));
    loop_ub = x.size(1) - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      x[i2] = x[i2] / static_cast<double>(dataHam.size(0));
    }
    nx = x.size(1);
    if (x.size(1) <= 2) {
      if (x.size(1) == 1) {
        varargin_1 = x[0];
      } else if (x[0] < x[x.size(1) - 1]) {
        varargin_1 = x[x.size(1) - 1];
      } else {
        varargin_1 = x[0];
      }
    } else {
      varargin_1 = x[0];
      for (k = 2; k <= nx; k++) {
        d = x[k - 1];
        if (varargin_1 < d) {
          varargin_1 = d;
        }
      }
    }
    if (varargin_1 > dqThreshold) {
      boolean_T ex;
      // 'parameterTuning:103' if max(meanPSDgain > PSDMaxUnique)
      c_x.set_size(1, meanPSDgain.size(1));
      loop_ub = meanPSDgain.size(1);
      for (i2 = 0; i2 < loop_ub; i2++) {
        c_x[i2] = (meanPSDgain[i2] > PSDMaxUnique);
      }
      nx = c_x.size(1);
      ex = c_x[0];
      for (k = 2; k <= nx; k++) {
        ex = ((static_cast<int>(ex) < static_cast<int>(c_x[k - 1])) || ex);
      }
      if (ex) {
        // 'parameterTuning:104' gainTune =
        // gainTune*(1+lambda1*(PSDmean-max(meanPSDgain))./PSDmean);
        nx = meanPSDgain.size(1);
        if (meanPSDgain.size(1) <= 2) {
          if (meanPSDgain.size(1) == 1) {
            varargin_1 = meanPSDgain[0];
          } else if (meanPSDgain[0] < meanPSDgain[meanPSDgain.size(1) - 1]) {
            varargin_1 = meanPSDgain[meanPSDgain.size(1) - 1];
          } else {
            varargin_1 = meanPSDgain[0];
          }
        } else {
          varargin_1 = meanPSDgain[0];
          for (k = 2; k <= nx; k++) {
            d = meanPSDgain[k - 1];
            if (varargin_1 < d) {
              varargin_1 = d;
            }
          }
        }
        gainTune *= lambda1 * (PSDmean - varargin_1) / PSDmean + 1.0;
        //  gainTune = max([gainTune*(1-lambda1),gainLowerBound]);
        // 'parameterTuning:106' gainTune = max([gainTune,gainLowerBound]);
        if (gainTune < normRho + 1.0) {
          gainTune = normRho + 1.0;
        }
      }
      // 'parameterTuning:109' if max(meanPSDgain < PSDMinUnique)
      c_x.set_size(1, meanPSDgain.size(1));
      loop_ub = meanPSDgain.size(1);
      for (i2 = 0; i2 < loop_ub; i2++) {
        c_x[i2] = (meanPSDgain[i2] < PSDMinUnique);
      }
      nx = c_x.size(1);
      ex = c_x[0];
      for (k = 2; k <= nx; k++) {
        ex = ((static_cast<int>(ex) < static_cast<int>(c_x[k - 1])) || ex);
      }
      if (ex) {
        // 'parameterTuning:110' gainTune =
        // gainTune*(1+lambda1*(PSDmean-max(meanPSDgain))./PSDmean);
        nx = meanPSDgain.size(1);
        if (meanPSDgain.size(1) <= 2) {
          if (meanPSDgain.size(1) == 1) {
            varargin_1 = meanPSDgain[0];
          } else if (meanPSDgain[0] < meanPSDgain[meanPSDgain.size(1) - 1]) {
            varargin_1 = meanPSDgain[meanPSDgain.size(1) - 1];
          } else {
            varargin_1 = meanPSDgain[0];
          }
        } else {
          varargin_1 = meanPSDgain[0];
          for (k = 2; k <= nx; k++) {
            d = meanPSDgain[k - 1];
            if (varargin_1 < d) {
              varargin_1 = d;
            }
          }
        }
        gainTune *= lambda1 * (PSDmean - varargin_1) / PSDmean + 1.0;
        // 'parameterTuning:111' gainTune = min([gainTune,gainUpperBound]);
        if (gainTune > gainMaxUnique) {
          gainTune = gainMaxUnique;
        }
      }
    } else {
      boolean_T ex;
      // 'parameterTuning:114' else
      // 'parameterTuning:115' if max(meanPSDgain > PSDMaxUnique)
      c_x.set_size(1, meanPSDgain.size(1));
      loop_ub = meanPSDgain.size(1);
      for (i2 = 0; i2 < loop_ub; i2++) {
        c_x[i2] = (meanPSDgain[i2] > PSDMaxUnique);
      }
      nx = c_x.size(1);
      ex = c_x[0];
      for (k = 2; k <= nx; k++) {
        ex = ((static_cast<int>(ex) < static_cast<int>(c_x[k - 1])) || ex);
      }
      if (ex) {
        // 'parameterTuning:116' gainTune =
        // gainTune*(1+lambda1*(PSDmean-max(meanPSDgain))./PSDmean);
        nx = meanPSDgain.size(1);
        if (meanPSDgain.size(1) <= 2) {
          if (meanPSDgain.size(1) == 1) {
            varargin_1 = meanPSDgain[0];
          } else if (meanPSDgain[0] < meanPSDgain[meanPSDgain.size(1) - 1]) {
            varargin_1 = meanPSDgain[meanPSDgain.size(1) - 1];
          } else {
            varargin_1 = meanPSDgain[0];
          }
        } else {
          varargin_1 = meanPSDgain[0];
          for (k = 2; k <= nx; k++) {
            d = meanPSDgain[k - 1];
            if (varargin_1 < d) {
              varargin_1 = d;
            }
          }
        }
        gainTune *= lambda1 * (PSDmean - varargin_1) / PSDmean + 1.0;
        //  gainTune = max([gainTune*(1-lambda1),gainLowerBound]);
        // 'parameterTuning:118' gainTune = max([gainTune,gainLowerBound]);
        if (gainTune < normRho + 1.0) {
          gainTune = normRho + 1.0;
        }
      }
    }
  } else {
    // 'parameterTuning:123' else
    // 'parameterTuning:124' counter = counter + 1;
    counter++;
  }
  // 'parameterTuning:128' avgPSDdata = meanPSDgain';
  avgPSDdata.set_size(meanPSDgain.size(1));
  loop_ub = meanPSDgain.size(1);
  for (i2 = 0; i2 < loop_ub; i2++) {
    avgPSDdata[i2] = meanPSDgain[i2];
  }
  //  for smoothing
  // 'parameterTuning:131' gainTuneContr = GainTuneOld + (gainTune -
  // GainTuneOld)*(counter-winSize*winPerc)/((1-winPerc)*winSize+1);
  //  ----
  // 'parameterTuning:134' gainTuneContr = max([gainTuneContr,gainLowerBound]);
  varargin_1 = GainTuneOld + (gainTune - GainTuneOld) *
                                 (counter - winSize * winPerc) /
                                 ((1.0 - winPerc) * winSize + 1.0);
  if (varargin_1 < normRho + 1.0) {
    *gainTuneContr = normRho + 1.0;
  } else {
    *gainTuneContr = varargin_1;
  }
}

//
// function [gainTuneContr, avgPSDdata]  = parameterTuning(dq, normRho,
// robTorque, sampleT, dqThreshold, winSize, winPerc, PSDMaxUnique,
// PSDMinUnique, kpUnique, gainMaxUnique, freqLim, lambda1)
//
// Arguments    : void
// Return Type  : void
//
void parameterTuning_init()
{
  meanPSDgain_not_empty = false;
  torqueBuffer_not_empty = false;
  dqMat_not_empty = false;
  // 'parameterTuning:50' counter = 1;
  counter = 1.0;
}

//
// File trailer for parameterTuning.cpp
//
// [EOF]
//
