//
// File: FFTImplementationCallback.h
//
// MATLAB Coder version            : 5.4
// C/C++ source code generated on  : 13-Feb-2023 16:08:56
//

#ifndef FFTIMPLEMENTATIONCALLBACK_H
#define FFTIMPLEMENTATIONCALLBACK_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coderTuner {
namespace coder {
namespace internal {
namespace fft {
class FFTImplementationCallback {
public:
  static void r2br_r2dit_trig(const ::coder::array<double, 2U> &x,
                              int n1_unsigned,
                              const ::coder::array<double, 2U> &costab,
                              const ::coder::array<double, 2U> &sintab,
                              ::coder::array<creal_T, 2U> &y);
  static void dobluesteinfft(const ::coder::array<double, 2U> &x, int n2blue,
                             int nfft, const ::coder::array<double, 2U> &costab,
                             const ::coder::array<double, 2U> &sintab,
                             const ::coder::array<double, 2U> &sintabinv,
                             ::coder::array<creal_T, 2U> &y);

protected:
  static void get_half_twiddle_tables(const ::coder::array<double, 2U> &costab,
                                      const ::coder::array<double, 2U> &sintab,
                                      ::coder::array<double, 2U> &hcostab,
                                      ::coder::array<double, 2U> &hsintab);
  static void
  get_half_twiddle_tables(const ::coder::array<double, 2U> &costab,
                          const ::coder::array<double, 2U> &sintab,
                          const ::coder::array<double, 2U> &costabinv,
                          const ::coder::array<double, 2U> &sintabinv,
                          ::coder::array<double, 2U> &hcostab,
                          ::coder::array<double, 2U> &hsintab,
                          ::coder::array<double, 2U> &hcostabinv,
                          ::coder::array<double, 2U> &hsintabinv);
  static void r2br_r2dit_trig(const ::coder::array<creal_T, 1U> &x,
                              int n1_unsigned,
                              const ::coder::array<double, 2U> &costab,
                              const ::coder::array<double, 2U> &sintab,
                              ::coder::array<creal_T, 1U> &y);
  static void b_r2br_r2dit_trig(const ::coder::array<creal_T, 1U> &x,
                                int n1_unsigned,
                                const ::coder::array<double, 2U> &costab,
                                const ::coder::array<double, 2U> &sintab,
                                ::coder::array<creal_T, 1U> &y);
  static void doHalfLengthRadix2(const ::coder::array<double, 2U> &x,
                                 int xoffInit, ::coder::array<creal_T, 1U> &y,
                                 int unsigned_nRows,
                                 const ::coder::array<double, 2U> &costab,
                                 const ::coder::array<double, 2U> &sintab);
  static void
  doHalfLengthBluestein(const ::coder::array<double, 2U> &x, int xoffInit,
                        ::coder::array<creal_T, 1U> &y, int nrowsx, int nRows,
                        int nfft, const ::coder::array<creal_T, 1U> &wwc,
                        const ::coder::array<double, 2U> &costab,
                        const ::coder::array<double, 2U> &sintab,
                        const ::coder::array<double, 2U> &costabinv,
                        const ::coder::array<double, 2U> &sintabinv);
  static void getback_radix2_fft(::coder::array<creal_T, 1U> &y,
                                 const ::coder::array<creal_T, 1U> &reconVar1,
                                 const ::coder::array<creal_T, 1U> &reconVar2,
                                 const int wrapIndex_data[], int hnRows);
};

} // namespace fft
} // namespace internal
} // namespace coder
} // namespace coderTuner

#endif
//
// File trailer for FFTImplementationCallback.h
//
// [EOF]
//
