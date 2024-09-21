#pragma once
#include "robustControllerURDF.h"
#include "robustControllerURDF_terminate.h"
#include "robustControllerURDF_types.h"
#include "coder_array.h"
#include "parameterTuning.h"
#include "parameterTuning_terminate.h"

// function declarations

// for both packages

static double argInit_real_T();

static coder::array<double, 1U> argInit_Unboundedx1_real_T();

// for parameter tuning

static void main_parameterTuning();

// for robust controller

static void argInit_1x3_real_T(double result[3]);

static coder::array<double, 2U> argInit_1xUnbounded_real_T();

static coder::array<double, 3U> argInit_3x3xUnbounded_real_T();

static coder::array<double, 2U> argInit_3xUnbounded_real_T();

static coder::array<double, 1U> argInit_Unboundedx1_real_T();

static coder::array<double, 2U> argInit_UnboundedxUnbounded_real_T();

static void argInit_struct0_T(struct0_T *result);

static struct1_T argInit_struct1_T();

static void argInit_struct2_T(struct2_T *result);

static void main_robustControllerURDF();