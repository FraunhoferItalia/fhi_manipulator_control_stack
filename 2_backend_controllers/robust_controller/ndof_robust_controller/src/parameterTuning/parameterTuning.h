//
// File: parameterTuning.h
//
// MATLAB Coder version            : 5.4
// C/C++ source code generated on  : 13-Feb-2023 16:08:56
//

#ifndef PARAMETERTUNING_H
#define PARAMETERTUNING_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
void GainTuneOld_not_empty_init();

void gainTune_not_empty_init();

extern void parameterTuning(const coder::array<double, 1U> &dq, double normRho,
                            const coder::array<double, 1U> &robTorque,
                            double sampleT, double dqThreshold, double winSize,
                            double winPerc, double PSDMaxUnique,
                            double PSDMinUnique, double kpUnique,
                            double gainMaxUnique, double freqLim,
                            double lambda1, double *gainTuneContr,
                            coder::array<double, 1U> &avgPSDdata);

void parameterTuning_init();

#endif
//
// File trailer for parameterTuning.h
//
// [EOF]
//
