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
#include "FFTImplementationCallback.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const ::coder::array<creal_T, 1U> &x
//                int n1_unsigned
//                const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                ::coder::array<creal_T, 1U> &y
// Return Type  : void
//
namespace coderTuner {
namespace coder {
namespace internal {
namespace fft {
void FFTImplementationCallback::b_r2br_r2dit_trig(
    const ::coder::array<creal_T, 1U> &x, int n1_unsigned,
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab, ::coder::array<creal_T, 1U> &y)
{
  double temp_im;
  double temp_re;
  double temp_re_tmp;
  double twid_re;
  int i;
  int iDelta2;
  int iheight;
  int istart;
  int iy;
  int ju;
  int k;
  int nRowsD2;
  y.set_size(n1_unsigned);
  if (n1_unsigned > x.size(0)) {
    y.set_size(n1_unsigned);
    for (iDelta2 = 0; iDelta2 < n1_unsigned; iDelta2++) {
      y[iDelta2].re = 0.0;
      y[iDelta2].im = 0.0;
    }
  }
  iheight = x.size(0);
  if (iheight > n1_unsigned) {
    iheight = n1_unsigned;
  }
  istart = n1_unsigned - 2;
  nRowsD2 = n1_unsigned / 2;
  k = nRowsD2 / 2;
  iy = 0;
  ju = 0;
  for (i = 0; i <= iheight - 2; i++) {
    boolean_T tst;
    y[iy] = x[i];
    iDelta2 = n1_unsigned;
    tst = true;
    while (tst) {
      iDelta2 >>= 1;
      ju ^= iDelta2;
      tst = ((ju & iDelta2) == 0);
    }
    iy = ju;
  }
  y[iy] = x[iheight - 1];
  if (n1_unsigned > 1) {
    for (i = 0; i <= istart; i += 2) {
      temp_re_tmp = y[i + 1].re;
      temp_im = y[i + 1].im;
      temp_re = y[i].re;
      twid_re = y[i].im;
      y[i + 1].re = temp_re - temp_re_tmp;
      y[i + 1].im = twid_re - temp_im;
      y[i].re = temp_re + temp_re_tmp;
      y[i].im = twid_re + temp_im;
    }
  }
  iy = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_temp_re_tmp;
    for (i = 0; i < iheight; i += iDelta2) {
      b_temp_re_tmp = i + iy;
      temp_re = y[b_temp_re_tmp].re;
      temp_im = y[b_temp_re_tmp].im;
      y[b_temp_re_tmp].re = y[i].re - temp_re;
      y[b_temp_re_tmp].im = y[i].im - temp_im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }
    istart = 1;
    for (ju = k; ju < nRowsD2; ju += k) {
      double twid_im;
      int ihi;
      twid_re = costab[ju];
      twid_im = sintab[ju];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        b_temp_re_tmp = i + iy;
        temp_re_tmp = y[b_temp_re_tmp].im;
        temp_im = y[b_temp_re_tmp].re;
        temp_re = twid_re * temp_im - twid_im * temp_re_tmp;
        temp_im = twid_re * temp_re_tmp + twid_im * temp_im;
        y[b_temp_re_tmp].re = y[i].re - temp_re;
        y[b_temp_re_tmp].im = y[i].im - temp_im;
        y[i].re = y[i].re + temp_re;
        y[i].im = y[i].im + temp_im;
        i += iDelta2;
      }
      istart++;
    }
    k /= 2;
    iy = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iy;
  }
  if (y.size(0) > 1) {
    temp_im = 1.0 / static_cast<double>(y.size(0));
    iy = y.size(0);
    for (iDelta2 = 0; iDelta2 < iy; iDelta2++) {
      y[iDelta2].re = temp_im * y[iDelta2].re;
      y[iDelta2].im = temp_im * y[iDelta2].im;
    }
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                int xoffInit
//                ::coder::array<creal_T, 1U> &y
//                int nrowsx
//                int nRows
//                int nfft
//                const ::coder::array<creal_T, 1U> &wwc
//                const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                const ::coder::array<double, 2U> &costabinv
//                const ::coder::array<double, 2U> &sintabinv
// Return Type  : void
//
void FFTImplementationCallback::doHalfLengthBluestein(
    const ::coder::array<double, 2U> &x, int xoffInit,
    ::coder::array<creal_T, 1U> &y, int nrowsx, int nRows, int nfft,
    const ::coder::array<creal_T, 1U> &wwc,
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab,
    const ::coder::array<double, 2U> &costabinv,
    const ::coder::array<double, 2U> &sintabinv)
{
  ::coder::array<creal_T, 1U> b_y;
  ::coder::array<creal_T, 1U> fv;
  ::coder::array<creal_T, 1U> r;
  ::coder::array<creal_T, 1U> reconVar1;
  ::coder::array<creal_T, 1U> reconVar2;
  ::coder::array<creal_T, 1U> ytmp;
  ::coder::array<double, 2U> a__1;
  ::coder::array<double, 2U> b_costab;
  ::coder::array<double, 2U> b_sintab;
  ::coder::array<double, 2U> b_sintabinv;
  ::coder::array<double, 2U> costab1q;
  ::coder::array<double, 2U> costable;
  ::coder::array<double, 2U> hcostabinv;
  ::coder::array<double, 2U> hsintab;
  ::coder::array<double, 2U> hsintabinv;
  ::coder::array<double, 2U> sintable;
  cuint8_T y_data[32767];
  double b_im;
  double temp_im;
  double temp_re;
  double temp_re_tmp;
  double twid_im;
  double twid_re;
  int wrapIndex_data[16383];
  int hnRows;
  int i;
  int iheight;
  int ihi;
  int istart;
  int j;
  int ju;
  int k;
  int nRowsD2;
  int nd2;
  int nfft_tmp;
  boolean_T tst;
  hnRows = nRows / 2;
  if (hnRows > nrowsx) {
    if (hnRows - 1 >= 0) {
      std::memset(&y_data[0], 0, hnRows * sizeof(cuint8_T));
    }
  }
  ytmp.set_size(hnRows);
  for (ihi = 0; ihi < hnRows; ihi++) {
    ytmp[ihi].re = 0.0;
    ytmp[ihi].im = y_data[ihi].im;
  }
  if ((x.size(0) & 1) == 0) {
    tst = true;
    j = x.size(0);
  } else if (x.size(0) >= nRows) {
    tst = true;
    j = nRows;
  } else {
    tst = false;
    j = x.size(0) - 1;
  }
  if (j > nRows) {
    j = nRows;
  }
  nd2 = nRows << 1;
  temp_im = 6.2831853071795862 / static_cast<double>(nd2);
  istart = nd2 / 2 / 2;
  costab1q.set_size(1, istart + 1);
  costab1q[0] = 1.0;
  nd2 = istart / 2 - 1;
  for (k = 0; k <= nd2; k++) {
    costab1q[k + 1] = std::cos(temp_im * (static_cast<double>(k) + 1.0));
  }
  ihi = nd2 + 2;
  nd2 = istart - 1;
  for (k = ihi; k <= nd2; k++) {
    costab1q[k] = std::sin(temp_im * static_cast<double>(istart - k));
  }
  costab1q[istart] = 0.0;
  istart = costab1q.size(1) - 1;
  nd2 = (costab1q.size(1) - 1) << 1;
  b_costab.set_size(1, nd2 + 1);
  b_sintab.set_size(1, nd2 + 1);
  b_costab[0] = 1.0;
  b_sintab[0] = 0.0;
  b_sintabinv.set_size(1, nd2 + 1);
  for (k = 0; k < istart; k++) {
    b_sintabinv[k + 1] = costab1q[(istart - k) - 1];
  }
  ihi = costab1q.size(1);
  for (k = ihi; k <= nd2; k++) {
    b_sintabinv[k] = costab1q[k - istart];
  }
  for (k = 0; k < istart; k++) {
    b_costab[k + 1] = costab1q[k + 1];
    b_sintab[k + 1] = -costab1q[(istart - k) - 1];
  }
  ihi = costab1q.size(1);
  for (k = ihi; k <= nd2; k++) {
    b_costab[k] = -costab1q[nd2 - k];
    b_sintab[k] = -costab1q[k - istart];
  }
  costable.set_size(1, b_costab.size(1));
  nd2 = b_costab.size(1);
  for (ihi = 0; ihi < nd2; ihi++) {
    costable[ihi] = b_costab[ihi];
  }
  sintable.set_size(1, b_sintab.size(1));
  nd2 = b_sintab.size(1);
  for (ihi = 0; ihi < nd2; ihi++) {
    sintable[ihi] = b_sintab[ihi];
  }
  a__1.set_size(1, b_sintabinv.size(1));
  nd2 = b_sintabinv.size(1);
  for (ihi = 0; ihi < nd2; ihi++) {
    a__1[ihi] = b_sintabinv[ihi];
  }
  FFTImplementationCallback::get_half_twiddle_tables(costab, sintab, costabinv,
                                                     sintabinv, a__1, hsintab,
                                                     hcostabinv, hsintabinv);
  reconVar1.set_size(hnRows);
  reconVar2.set_size(hnRows);
  for (i = 0; i < hnRows; i++) {
    ihi = i << 1;
    temp_im = sintable[ihi];
    temp_re = costable[ihi];
    reconVar1[i].re = temp_im + 1.0;
    reconVar1[i].im = -temp_re;
    reconVar2[i].re = 1.0 - temp_im;
    reconVar2[i].im = temp_re;
    if (i + 1 != 1) {
      wrapIndex_data[i] = (hnRows - i) + 1;
    } else {
      wrapIndex_data[0] = 1;
    }
  }
  temp_im = static_cast<double>(j) / 2.0;
  ihi = static_cast<int>(temp_im);
  for (ju = 0; ju < ihi; ju++) {
    nd2 = (hnRows + ju) - 1;
    temp_re = wwc[nd2].re;
    twid_re = wwc[nd2].im;
    nd2 = xoffInit + (ju << 1);
    twid_im = x[nd2];
    b_im = x[nd2 + 1];
    ytmp[ju].re = temp_re * twid_im + twid_re * b_im;
    ytmp[ju].im = temp_re * b_im - twid_re * twid_im;
  }
  if (!tst) {
    nd2 = (hnRows + static_cast<int>(temp_im)) - 1;
    twid_im = x[xoffInit + (static_cast<int>(temp_im) << 1)];
    ytmp[static_cast<int>(temp_im)].re = wwc[nd2].re * twid_im;
    ytmp[static_cast<int>(temp_im)].im = 0.0 - wwc[nd2].im * twid_im;
    if (static_cast<int>(temp_im) + 2 <= hnRows) {
      ihi = static_cast<int>(static_cast<double>(j) / 2.0) + 2;
      for (i = ihi; i <= hnRows; i++) {
        ytmp[i - 1].re = 0.0;
        ytmp[i - 1].im = 0.0;
      }
    }
  } else if (static_cast<int>(temp_im) + 1 <= hnRows) {
    ihi = static_cast<int>(static_cast<double>(j) / 2.0) + 1;
    for (i = ihi; i <= hnRows; i++) {
      ytmp[i - 1].re = 0.0;
      ytmp[i - 1].im = 0.0;
    }
  }
  nfft_tmp = static_cast<int>(static_cast<double>(nfft) / 2.0);
  b_y.set_size(nfft_tmp);
  if (nfft_tmp > ytmp.size(0)) {
    b_y.set_size(nfft_tmp);
    for (ihi = 0; ihi < nfft_tmp; ihi++) {
      b_y[ihi].re = 0.0;
      b_y[ihi].im = 0.0;
    }
  }
  fv.set_size(b_y.size(0));
  nd2 = b_y.size(0);
  for (ihi = 0; ihi < nd2; ihi++) {
    fv[ihi] = b_y[ihi];
  }
  j = ytmp.size(0);
  if (j > nfft_tmp) {
    j = nfft_tmp;
  }
  iheight = nfft_tmp - 2;
  nRowsD2 = nfft_tmp / 2;
  k = nRowsD2 / 2;
  nd2 = 0;
  ju = 0;
  for (i = 0; i <= j - 2; i++) {
    fv[nd2] = ytmp[i];
    istart = nfft_tmp;
    tst = true;
    while (tst) {
      istart >>= 1;
      ju ^= istart;
      tst = ((ju & istart) == 0);
    }
    nd2 = ju;
  }
  fv[nd2] = ytmp[j - 1];
  b_y.set_size(fv.size(0));
  nd2 = fv.size(0);
  for (ihi = 0; ihi < nd2; ihi++) {
    b_y[ihi] = fv[ihi];
  }
  if (nfft_tmp > 1) {
    for (i = 0; i <= iheight; i += 2) {
      temp_re_tmp = b_y[i + 1].re;
      temp_re = b_y[i + 1].im;
      b_im = b_y[i].re;
      temp_im = b_y[i].im;
      b_y[i + 1].re = b_im - temp_re_tmp;
      b_y[i + 1].im = temp_im - temp_re;
      b_y[i].re = b_im + temp_re_tmp;
      b_y[i].im = temp_im + temp_re;
    }
  }
  nd2 = 2;
  ju = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_temp_re_tmp;
    for (i = 0; i < iheight; i += ju) {
      b_temp_re_tmp = i + nd2;
      temp_re = b_y[b_temp_re_tmp].re;
      temp_im = b_y[b_temp_re_tmp].im;
      b_y[b_temp_re_tmp].re = b_y[i].re - temp_re;
      b_y[b_temp_re_tmp].im = b_y[i].im - temp_im;
      b_y[i].re = b_y[i].re + temp_re;
      b_y[i].im = b_y[i].im + temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      twid_re = a__1[j];
      twid_im = hsintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        b_temp_re_tmp = i + nd2;
        temp_re_tmp = b_y[b_temp_re_tmp].im;
        temp_im = b_y[b_temp_re_tmp].re;
        temp_re = twid_re * temp_im - twid_im * temp_re_tmp;
        temp_im = twid_re * temp_re_tmp + twid_im * temp_im;
        b_y[b_temp_re_tmp].re = b_y[i].re - temp_re;
        b_y[b_temp_re_tmp].im = b_y[i].im - temp_im;
        b_y[i].re = b_y[i].re + temp_re;
        b_y[i].im = b_y[i].im + temp_im;
        i += ju;
      }
      istart++;
    }
    k /= 2;
    nd2 = ju;
    ju += ju;
    iheight -= nd2;
  }
  fv.set_size(b_y.size(0));
  nd2 = b_y.size(0);
  for (ihi = 0; ihi < nd2; ihi++) {
    fv[ihi] = b_y[ihi];
  }
  FFTImplementationCallback::r2br_r2dit_trig(wwc, nfft_tmp, a__1, hsintab, r);
  nd2 = fv.size(0);
  for (ihi = 0; ihi < nd2; ihi++) {
    b_im = fv[ihi].re;
    temp_im = r[ihi].im;
    temp_re = fv[ihi].im;
    twid_re = r[ihi].re;
    fv[ihi].re = b_im * twid_re - temp_re * temp_im;
    fv[ihi].im = b_im * temp_im + temp_re * twid_re;
  }
  FFTImplementationCallback::b_r2br_r2dit_trig(fv, nfft_tmp, hcostabinv,
                                               hsintabinv, r);
  fv.set_size(r.size(0));
  nd2 = r.size(0);
  for (ihi = 0; ihi < nd2; ihi++) {
    fv[ihi] = r[ihi];
  }
  ihi = wwc.size(0);
  for (k = hnRows; k <= ihi; k++) {
    temp_im = wwc[k - 1].re;
    temp_re = fv[k - 1].im;
    twid_re = wwc[k - 1].im;
    twid_im = fv[k - 1].re;
    nd2 = k - hnRows;
    ytmp[nd2].re = temp_im * twid_im + twid_re * temp_re;
    ytmp[nd2].im = temp_im * temp_re - twid_re * twid_im;
  }
  for (i = 0; i < hnRows; i++) {
    double b_ytmp_re_tmp;
    double ytmp_re_tmp;
    ihi = wrapIndex_data[i];
    temp_im = ytmp[i].re;
    b_im = reconVar1[i].im;
    temp_re = ytmp[i].im;
    twid_re = reconVar1[i].re;
    twid_im = ytmp[ihi - 1].re;
    temp_re_tmp = -ytmp[ihi - 1].im;
    ytmp_re_tmp = reconVar2[i].im;
    b_ytmp_re_tmp = reconVar2[i].re;
    y[i].re = 0.5 * ((temp_im * twid_re - temp_re * b_im) +
                     (twid_im * b_ytmp_re_tmp - temp_re_tmp * ytmp_re_tmp));
    y[i].im = 0.5 * ((temp_im * b_im + temp_re * twid_re) +
                     (twid_im * ytmp_re_tmp + temp_re_tmp * b_ytmp_re_tmp));
    ihi = hnRows + i;
    y[ihi].re = 0.5 * ((temp_im * b_ytmp_re_tmp - temp_re * ytmp_re_tmp) +
                       (twid_im * twid_re - temp_re_tmp * b_im));
    y[ihi].im = 0.5 * ((temp_im * ytmp_re_tmp + temp_re * b_ytmp_re_tmp) +
                       (twid_im * b_im + temp_re_tmp * twid_re));
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                int xoffInit
//                ::coder::array<creal_T, 1U> &y
//                int unsigned_nRows
//                const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
// Return Type  : void
//
void FFTImplementationCallback::doHalfLengthRadix2(
    const ::coder::array<double, 2U> &x, int xoffInit,
    ::coder::array<creal_T, 1U> &y, int unsigned_nRows,
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab)
{
  ::coder::array<creal_T, 1U> b_y;
  ::coder::array<creal_T, 1U> reconVar1;
  ::coder::array<creal_T, 1U> reconVar2;
  ::coder::array<double, 2U> hcostab;
  ::coder::array<double, 2U> hsintab;
  ::coder::array<int, 1U> bitrevIndex;
  double temp_im;
  double temp_re;
  double temp_re_tmp;
  double twid_re;
  double z;
  int wrapIndex_data[16383];
  int i;
  int iDelta;
  int istart;
  int iy;
  int j;
  int ju;
  int k;
  int nRows;
  int nRowsD2;
  int nRowsM2;
  boolean_T tst;
  nRows = unsigned_nRows / 2;
  istart = y.size(0);
  if (istart > nRows) {
    istart = nRows;
  }
  nRowsM2 = nRows - 2;
  nRowsD2 = nRows / 2;
  k = nRowsD2 / 2;
  FFTImplementationCallback::get_half_twiddle_tables(costab, sintab, hcostab,
                                                     hsintab);
  reconVar1.set_size(nRows);
  reconVar2.set_size(nRows);
  for (i = 0; i < nRows; i++) {
    temp_re = sintab[i];
    temp_im = costab[i];
    reconVar1[i].re = temp_re + 1.0;
    reconVar1[i].im = -temp_im;
    reconVar2[i].re = 1.0 - temp_re;
    reconVar2[i].im = temp_im;
    if (i + 1 != 1) {
      wrapIndex_data[i] = (nRows - i) + 1;
    } else {
      wrapIndex_data[0] = 1;
    }
  }
  z = static_cast<double>(unsigned_nRows) / 2.0;
  ju = 0;
  iy = 1;
  iDelta = static_cast<int>(z);
  bitrevIndex.set_size(iDelta);
  for (j = 0; j < iDelta; j++) {
    bitrevIndex[j] = 0;
  }
  for (j = 0; j <= istart - 2; j++) {
    bitrevIndex[j] = iy;
    iDelta = static_cast<int>(z);
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }
    iy = ju + 1;
  }
  bitrevIndex[istart - 1] = iy;
  if ((x.size(0) & 1) == 0) {
    tst = true;
    istart = x.size(0);
  } else if (x.size(0) >= unsigned_nRows) {
    tst = true;
    istart = unsigned_nRows;
  } else {
    tst = false;
    istart = x.size(0) - 1;
  }
  if (istart > unsigned_nRows) {
    istart = unsigned_nRows;
  }
  temp_re = static_cast<double>(istart) / 2.0;
  j = static_cast<int>(temp_re);
  for (i = 0; i < j; i++) {
    iDelta = xoffInit + (i << 1);
    y[bitrevIndex[i] - 1].re = x[iDelta];
    y[bitrevIndex[i] - 1].im = x[iDelta + 1];
  }
  if (!tst) {
    j = bitrevIndex[static_cast<int>(temp_re)] - 1;
    y[j].re = x[xoffInit + (static_cast<int>(temp_re) << 1)];
    y[j].im = 0.0;
  }
  b_y.set_size(y.size(0));
  iDelta = y.size(0);
  for (j = 0; j < iDelta; j++) {
    b_y[j] = y[j];
  }
  if (nRows > 1) {
    for (i = 0; i <= nRowsM2; i += 2) {
      temp_re_tmp = b_y[i + 1].re;
      temp_re = b_y[i + 1].im;
      temp_im = b_y[i].re;
      twid_re = b_y[i].im;
      b_y[i + 1].re = temp_im - temp_re_tmp;
      b_y[i + 1].im = twid_re - temp_re;
      b_y[i].re = temp_im + temp_re_tmp;
      b_y[i].im = twid_re + temp_re;
    }
  }
  iDelta = 2;
  iy = 4;
  ju = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (i = 0; i < ju; i += iy) {
      nRowsM2 = i + iDelta;
      temp_re = b_y[nRowsM2].re;
      temp_im = b_y[nRowsM2].im;
      b_y[nRowsM2].re = b_y[i].re - temp_re;
      b_y[nRowsM2].im = b_y[i].im - temp_im;
      b_y[i].re = b_y[i].re + temp_re;
      b_y[i].im = b_y[i].im + temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      double twid_im;
      twid_re = hcostab[j];
      twid_im = hsintab[j];
      i = istart;
      nRows = istart + ju;
      while (i < nRows) {
        nRowsM2 = i + iDelta;
        temp_re_tmp = b_y[nRowsM2].im;
        temp_im = b_y[nRowsM2].re;
        temp_re = twid_re * temp_im - twid_im * temp_re_tmp;
        temp_im = twid_re * temp_re_tmp + twid_im * temp_im;
        b_y[nRowsM2].re = b_y[i].re - temp_re;
        b_y[nRowsM2].im = b_y[i].im - temp_im;
        b_y[i].re = b_y[i].re + temp_re;
        b_y[i].im = b_y[i].im + temp_im;
        i += iy;
      }
      istart++;
    }
    k /= 2;
    iDelta = iy;
    iy += iy;
    ju -= iDelta;
  }
  y.set_size(b_y.size(0));
  iDelta = b_y.size(0);
  for (j = 0; j < iDelta; j++) {
    y[j] = b_y[j];
  }
  FFTImplementationCallback::getback_radix2_fft(
      y, reconVar1, reconVar2, wrapIndex_data, static_cast<int>(z));
}

//
// Arguments    : const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                ::coder::array<double, 2U> &hcostab
//                ::coder::array<double, 2U> &hsintab
// Return Type  : void
//
void FFTImplementationCallback::get_half_twiddle_tables(
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab,
    ::coder::array<double, 2U> &hcostab, ::coder::array<double, 2U> &hsintab)
{
  int hszCostab;
  hszCostab = costab.size(1) / 2;
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  for (int i{0}; i < hszCostab; i++) {
    int hcostab_tmp;
    hcostab_tmp = ((i + 1) << 1) - 2;
    hcostab[i] = costab[hcostab_tmp];
    hsintab[i] = sintab[hcostab_tmp];
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                const ::coder::array<double, 2U> &costabinv
//                const ::coder::array<double, 2U> &sintabinv
//                ::coder::array<double, 2U> &hcostab
//                ::coder::array<double, 2U> &hsintab
//                ::coder::array<double, 2U> &hcostabinv
//                ::coder::array<double, 2U> &hsintabinv
// Return Type  : void
//
void FFTImplementationCallback::get_half_twiddle_tables(
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab,
    const ::coder::array<double, 2U> &costabinv,
    const ::coder::array<double, 2U> &sintabinv,
    ::coder::array<double, 2U> &hcostab, ::coder::array<double, 2U> &hsintab,
    ::coder::array<double, 2U> &hcostabinv,
    ::coder::array<double, 2U> &hsintabinv)
{
  int hszCostab;
  hszCostab = costab.size(1) / 2;
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  hcostabinv.set_size(1, hszCostab);
  hsintabinv.set_size(1, hszCostab);
  for (int i{0}; i < hszCostab; i++) {
    int hcostab_tmp;
    hcostab_tmp = ((i + 1) << 1) - 2;
    hcostab[i] = costab[hcostab_tmp];
    hsintab[i] = sintab[hcostab_tmp];
    hcostabinv[i] = costabinv[hcostab_tmp];
    hsintabinv[i] = sintabinv[hcostab_tmp];
  }
}

//
// Arguments    : ::coder::array<creal_T, 1U> &y
//                const ::coder::array<creal_T, 1U> &reconVar1
//                const ::coder::array<creal_T, 1U> &reconVar2
//                const int wrapIndex_data[]
//                int hnRows
// Return Type  : void
//
void FFTImplementationCallback::getback_radix2_fft(
    ::coder::array<creal_T, 1U> &y,
    const ::coder::array<creal_T, 1U> &reconVar1,
    const ::coder::array<creal_T, 1U> &reconVar2, const int wrapIndex_data[],
    int hnRows)
{
  double b_temp1_re_tmp;
  double b_temp2_re_tmp;
  double b_y_re_tmp;
  double c_temp1_re_tmp;
  double c_temp2_re_tmp;
  double c_y_re_tmp;
  double d_y_re_tmp;
  double temp1_im_tmp;
  double temp1_re_tmp;
  double y_re_tmp;
  int b_i;
  int iterVar;
  iterVar = hnRows / 2;
  temp1_re_tmp = y[0].re;
  temp1_im_tmp = y[0].im;
  y[0].re =
      0.5 *
      ((temp1_re_tmp * reconVar1[0].re - temp1_im_tmp * reconVar1[0].im) +
       (temp1_re_tmp * reconVar2[0].re - -temp1_im_tmp * reconVar2[0].im));
  y[0].im =
      0.5 *
      ((temp1_re_tmp * reconVar1[0].im + temp1_im_tmp * reconVar1[0].re) +
       (temp1_re_tmp * reconVar2[0].im + -temp1_im_tmp * reconVar2[0].re));
  y[hnRows].re =
      0.5 *
      ((temp1_re_tmp * reconVar2[0].re - temp1_im_tmp * reconVar2[0].im) +
       (temp1_re_tmp * reconVar1[0].re - -temp1_im_tmp * reconVar1[0].im));
  y[hnRows].im =
      0.5 *
      ((temp1_re_tmp * reconVar2[0].im + temp1_im_tmp * reconVar2[0].re) +
       (temp1_re_tmp * reconVar1[0].im + -temp1_im_tmp * reconVar1[0].re));
  for (int i{2}; i <= iterVar; i++) {
    double temp2_im_tmp;
    double temp2_re_tmp;
    int i1;
    temp1_re_tmp = y[i - 1].re;
    temp1_im_tmp = y[i - 1].im;
    b_i = wrapIndex_data[i - 1];
    temp2_re_tmp = y[b_i - 1].re;
    temp2_im_tmp = y[b_i - 1].im;
    y_re_tmp = reconVar1[i - 1].im;
    b_y_re_tmp = reconVar1[i - 1].re;
    c_y_re_tmp = reconVar2[i - 1].im;
    d_y_re_tmp = reconVar2[i - 1].re;
    y[i - 1].re =
        0.5 * ((temp1_re_tmp * b_y_re_tmp - temp1_im_tmp * y_re_tmp) +
               (temp2_re_tmp * d_y_re_tmp - -temp2_im_tmp * c_y_re_tmp));
    y[i - 1].im =
        0.5 * ((temp1_re_tmp * y_re_tmp + temp1_im_tmp * b_y_re_tmp) +
               (temp2_re_tmp * c_y_re_tmp + -temp2_im_tmp * d_y_re_tmp));
    i1 = (hnRows + i) - 1;
    y[i1].re = 0.5 * ((temp1_re_tmp * d_y_re_tmp - temp1_im_tmp * c_y_re_tmp) +
                      (temp2_re_tmp * b_y_re_tmp - -temp2_im_tmp * y_re_tmp));
    y[i1].im = 0.5 * ((temp1_re_tmp * c_y_re_tmp + temp1_im_tmp * d_y_re_tmp) +
                      (temp2_re_tmp * y_re_tmp + -temp2_im_tmp * b_y_re_tmp));
    c_temp2_re_tmp = reconVar1[b_i - 1].im;
    b_temp2_re_tmp = reconVar1[b_i - 1].re;
    b_temp1_re_tmp = reconVar2[b_i - 1].im;
    c_temp1_re_tmp = reconVar2[b_i - 1].re;
    y[b_i - 1].re =
        0.5 *
        ((temp2_re_tmp * b_temp2_re_tmp - temp2_im_tmp * c_temp2_re_tmp) +
         (temp1_re_tmp * c_temp1_re_tmp - -temp1_im_tmp * b_temp1_re_tmp));
    y[b_i - 1].im =
        0.5 *
        ((temp2_re_tmp * c_temp2_re_tmp + temp2_im_tmp * b_temp2_re_tmp) +
         (temp1_re_tmp * b_temp1_re_tmp + -temp1_im_tmp * c_temp1_re_tmp));
    b_i = (b_i + hnRows) - 1;
    y[b_i].re =
        0.5 *
        ((temp2_re_tmp * c_temp1_re_tmp - temp2_im_tmp * b_temp1_re_tmp) +
         (temp1_re_tmp * b_temp2_re_tmp - -temp1_im_tmp * c_temp2_re_tmp));
    y[b_i].im =
        0.5 *
        ((temp2_re_tmp * b_temp1_re_tmp + temp2_im_tmp * c_temp1_re_tmp) +
         (temp1_re_tmp * c_temp2_re_tmp + -temp1_im_tmp * b_temp2_re_tmp));
  }
  if (iterVar != 0) {
    temp1_re_tmp = y[iterVar].re;
    temp1_im_tmp = y[iterVar].im;
    y_re_tmp = reconVar1[iterVar].im;
    b_y_re_tmp = reconVar1[iterVar].re;
    c_y_re_tmp = reconVar2[iterVar].im;
    d_y_re_tmp = reconVar2[iterVar].re;
    b_temp2_re_tmp = temp1_re_tmp * d_y_re_tmp;
    b_temp1_re_tmp = temp1_re_tmp * b_y_re_tmp;
    y[iterVar].re = 0.5 * ((b_temp1_re_tmp - temp1_im_tmp * y_re_tmp) +
                           (b_temp2_re_tmp - -temp1_im_tmp * c_y_re_tmp));
    c_temp1_re_tmp = temp1_re_tmp * c_y_re_tmp;
    c_temp2_re_tmp = temp1_re_tmp * y_re_tmp;
    y[iterVar].im = 0.5 * ((c_temp2_re_tmp + temp1_im_tmp * b_y_re_tmp) +
                           (c_temp1_re_tmp + -temp1_im_tmp * d_y_re_tmp));
    b_i = hnRows + iterVar;
    y[b_i].re = 0.5 * ((b_temp2_re_tmp - temp1_im_tmp * c_y_re_tmp) +
                       (b_temp1_re_tmp - -temp1_im_tmp * y_re_tmp));
    y[b_i].im = 0.5 * ((c_temp1_re_tmp + temp1_im_tmp * d_y_re_tmp) +
                       (c_temp2_re_tmp + -temp1_im_tmp * b_y_re_tmp));
  }
}

//
// Arguments    : const ::coder::array<creal_T, 1U> &x
//                int n1_unsigned
//                const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                ::coder::array<creal_T, 1U> &y
// Return Type  : void
//
void FFTImplementationCallback::r2br_r2dit_trig(
    const ::coder::array<creal_T, 1U> &x, int n1_unsigned,
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab, ::coder::array<creal_T, 1U> &y)
{
  double temp_im;
  double temp_re;
  double temp_re_tmp;
  double twid_re;
  int i;
  int iDelta2;
  int iheight;
  int iy;
  int ju;
  int k;
  int nRowsD2;
  y.set_size(n1_unsigned);
  if (n1_unsigned > x.size(0)) {
    y.set_size(n1_unsigned);
    for (iy = 0; iy < n1_unsigned; iy++) {
      y[iy].re = 0.0;
      y[iy].im = 0.0;
    }
  }
  iDelta2 = x.size(0);
  if (iDelta2 > n1_unsigned) {
    iDelta2 = n1_unsigned;
  }
  iheight = n1_unsigned - 2;
  nRowsD2 = n1_unsigned / 2;
  k = nRowsD2 / 2;
  iy = 0;
  ju = 0;
  for (i = 0; i <= iDelta2 - 2; i++) {
    boolean_T tst;
    y[iy] = x[i];
    iy = n1_unsigned;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }
    iy = ju;
  }
  y[iy] = x[iDelta2 - 1];
  if (n1_unsigned > 1) {
    for (i = 0; i <= iheight; i += 2) {
      temp_re_tmp = y[i + 1].re;
      temp_im = y[i + 1].im;
      temp_re = y[i].re;
      twid_re = y[i].im;
      y[i + 1].re = temp_re - temp_re_tmp;
      y[i + 1].im = twid_re - temp_im;
      y[i].re = temp_re + temp_re_tmp;
      y[i].im = twid_re + temp_im;
    }
  }
  iy = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_temp_re_tmp;
    for (i = 0; i < iheight; i += iDelta2) {
      b_temp_re_tmp = i + iy;
      temp_re = y[b_temp_re_tmp].re;
      temp_im = y[b_temp_re_tmp].im;
      y[b_temp_re_tmp].re = y[i].re - temp_re;
      y[b_temp_re_tmp].im = y[i].im - temp_im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }
    ju = 1;
    for (int j{k}; j < nRowsD2; j += k) {
      double twid_im;
      int ihi;
      twid_re = costab[j];
      twid_im = sintab[j];
      i = ju;
      ihi = ju + iheight;
      while (i < ihi) {
        b_temp_re_tmp = i + iy;
        temp_re_tmp = y[b_temp_re_tmp].im;
        temp_im = y[b_temp_re_tmp].re;
        temp_re = twid_re * temp_im - twid_im * temp_re_tmp;
        temp_im = twid_re * temp_re_tmp + twid_im * temp_im;
        y[b_temp_re_tmp].re = y[i].re - temp_re;
        y[b_temp_re_tmp].im = y[i].im - temp_im;
        y[i].re = y[i].re + temp_re;
        y[i].im = y[i].im + temp_im;
        i += iDelta2;
      }
      ju++;
    }
    k /= 2;
    iy = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iy;
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                int n2blue
//                int nfft
//                const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                const ::coder::array<double, 2U> &sintabinv
//                ::coder::array<creal_T, 2U> &y
// Return Type  : void
//
void FFTImplementationCallback::dobluesteinfft(
    const ::coder::array<double, 2U> &x, int n2blue, int nfft,
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab,
    const ::coder::array<double, 2U> &sintabinv, ::coder::array<creal_T, 2U> &y)
{
  ::coder::array<creal_T, 1U> b_fv;
  ::coder::array<creal_T, 1U> fv;
  ::coder::array<creal_T, 1U> r;
  ::coder::array<creal_T, 1U> wwc;
  double b_re_tmp;
  double c_re_tmp;
  double d_re_tmp;
  double re_tmp;
  int b_k;
  int b_y;
  int i;
  int i1;
  int minNrowsNx;
  int nInt2m1;
  int u0;
  int xoff;
  if ((nfft != 1) && ((nfft & 1) == 0)) {
    int nInt2;
    int nRows;
    int rt;
    nRows = nfft / 2;
    nInt2m1 = (nRows + nRows) - 1;
    wwc.set_size(nInt2m1);
    rt = 0;
    wwc[nRows - 1].re = 1.0;
    wwc[nRows - 1].im = 0.0;
    nInt2 = nRows << 1;
    for (int k{0}; k <= nRows - 2; k++) {
      double nt_im;
      double nt_re;
      b_y = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= b_y) {
        rt += b_y - nInt2;
      } else {
        rt += b_y;
      }
      nt_im = -3.1415926535897931 * static_cast<double>(rt) /
              static_cast<double>(nRows);
      if (nt_im == 0.0) {
        nt_re = 1.0;
        nt_im = 0.0;
      } else {
        nt_re = std::cos(nt_im);
        nt_im = std::sin(nt_im);
      }
      i = (nRows - k) - 2;
      wwc[i].re = nt_re;
      wwc[i].im = -nt_im;
    }
    i = nInt2m1 - 1;
    for (int k{i}; k >= nRows; k--) {
      wwc[k] = wwc[(nInt2m1 - k) - 1];
    }
  } else {
    int nInt2;
    int rt;
    nInt2m1 = (nfft + nfft) - 1;
    wwc.set_size(nInt2m1);
    rt = 0;
    wwc[nfft - 1].re = 1.0;
    wwc[nfft - 1].im = 0.0;
    nInt2 = nfft << 1;
    for (int k{0}; k <= nfft - 2; k++) {
      double nt_im;
      double nt_re;
      b_y = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= b_y) {
        rt += b_y - nInt2;
      } else {
        rt += b_y;
      }
      nt_im = -3.1415926535897931 * static_cast<double>(rt) /
              static_cast<double>(nfft);
      if (nt_im == 0.0) {
        nt_re = 1.0;
        nt_im = 0.0;
      } else {
        nt_re = std::cos(nt_im);
        nt_im = std::sin(nt_im);
      }
      i = (nfft - k) - 2;
      wwc[i].re = nt_re;
      wwc[i].im = -nt_im;
    }
    i = nInt2m1 - 1;
    for (int k{i}; k >= nfft; k--) {
      wwc[k] = wwc[(nInt2m1 - k) - 1];
    }
  }
  nInt2m1 = x.size(0);
  y.set_size(nfft, x.size(1));
  if (nfft > x.size(0)) {
    y.set_size(nfft, x.size(1));
    b_y = nfft * x.size(1);
    for (i = 0; i < b_y; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  b_y = x.size(1) - 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    fv, b_fv, r, xoff, i1, minNrowsNx, u0, b_k, re_tmp, b_re_tmp, c_re_tmp,    \
    d_re_tmp)

  for (int chan = 0; chan <= b_y; chan++) {
    xoff = chan * nInt2m1;
    r.set_size(nfft);
    if (nfft > x.size(0)) {
      r.set_size(nfft);
      for (i1 = 0; i1 < nfft; i1++) {
        r[i1].re = 0.0;
        r[i1].im = 0.0;
      }
    }
    if ((n2blue != 1) && ((nfft & 1) == 0)) {
      FFTImplementationCallback::doHalfLengthBluestein(
          x, xoff, r, x.size(0), nfft, n2blue, wwc, costab, sintab, costab,
          sintabinv);
    } else {
      minNrowsNx = x.size(0);
      if (nfft <= minNrowsNx) {
        minNrowsNx = nfft;
      }
      for (b_k = 0; b_k < minNrowsNx; b_k++) {
        u0 = (nfft + b_k) - 1;
        i1 = xoff + b_k;
        r[b_k].re = wwc[u0].re * x[i1];
        r[b_k].im = wwc[u0].im * -x[i1];
      }
      i1 = minNrowsNx + 1;
      for (b_k = i1; b_k <= nfft; b_k++) {
        r[b_k - 1].re = 0.0;
        r[b_k - 1].im = 0.0;
      }
      FFTImplementationCallback::r2br_r2dit_trig(r, n2blue, costab, sintab,
                                                 b_fv);
      FFTImplementationCallback::r2br_r2dit_trig(wwc, n2blue, costab, sintab,
                                                 fv);
      fv.set_size(b_fv.size(0));
      u0 = b_fv.size(0);
      for (i1 = 0; i1 < u0; i1++) {
        re_tmp = b_fv[i1].re;
        b_re_tmp = fv[i1].im;
        c_re_tmp = b_fv[i1].im;
        d_re_tmp = fv[i1].re;
        fv[i1].re = re_tmp * d_re_tmp - c_re_tmp * b_re_tmp;
        fv[i1].im = re_tmp * b_re_tmp + c_re_tmp * d_re_tmp;
      }
      FFTImplementationCallback::b_r2br_r2dit_trig(fv, n2blue, costab,
                                                   sintabinv, b_fv);
      i1 = wwc.size(0);
      for (b_k = nfft; b_k <= i1; b_k++) {
        re_tmp = wwc[b_k - 1].re;
        b_re_tmp = b_fv[b_k - 1].im;
        c_re_tmp = wwc[b_k - 1].im;
        d_re_tmp = b_fv[b_k - 1].re;
        u0 = b_k - nfft;
        r[u0].re = re_tmp * d_re_tmp + c_re_tmp * b_re_tmp;
        r[u0].im = re_tmp * b_re_tmp - c_re_tmp * d_re_tmp;
      }
    }
    u0 = r.size(0);
    for (i1 = 0; i1 < u0; i1++) {
      y[i1 + y.size(0) * chan] = r[i1];
    }
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                int n1_unsigned
//                const ::coder::array<double, 2U> &costab
//                const ::coder::array<double, 2U> &sintab
//                ::coder::array<creal_T, 2U> &y
// Return Type  : void
//
void FFTImplementationCallback::r2br_r2dit_trig(
    const ::coder::array<double, 2U> &x, int n1_unsigned,
    const ::coder::array<double, 2U> &costab,
    const ::coder::array<double, 2U> &sintab, ::coder::array<creal_T, 2U> &y)
{
  ::coder::array<creal_T, 1U> b_y;
  ::coder::array<creal_T, 1U> r;
  int i1;
  int loop_ub;
  int nrows;
  int u0;
  int xoff;
  nrows = x.size(0);
  y.set_size(n1_unsigned, x.size(1));
  if (n1_unsigned > x.size(0)) {
    y.set_size(n1_unsigned, x.size(1));
    loop_ub = n1_unsigned * x.size(1);
    for (int i{0}; i < loop_ub; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }
  loop_ub = x.size(1) - 1;
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
    b_y, r, xoff, i1, u0)

  for (int chan = 0; chan <= loop_ub; chan++) {
    xoff = chan * nrows;
    r.set_size(n1_unsigned);
    if (n1_unsigned > x.size(0)) {
      r.set_size(n1_unsigned);
      for (i1 = 0; i1 < n1_unsigned; i1++) {
        r[i1].re = 0.0;
        r[i1].im = 0.0;
      }
    }
    if (n1_unsigned != 1) {
      FFTImplementationCallback::doHalfLengthRadix2(x, xoff, r, n1_unsigned,
                                                    costab, sintab);
      u0 = r.size(0);
      for (i1 = 0; i1 < u0; i1++) {
        y[i1 + y.size(0) * chan] = r[i1];
      }
    } else {
      u0 = x.size(0);
      if (u0 > 1) {
        u0 = 1;
      }
      r[0].re = x[(xoff + u0) - 1];
      r[0].im = 0.0;
      b_y.set_size(r.size(0));
      u0 = r.size(0);
      for (i1 = 0; i1 < u0; i1++) {
        b_y[i1] = r[i1];
      }
      u0 = b_y.size(0);
      for (i1 = 0; i1 < u0; i1++) {
        y[i1 + y.size(0) * chan] = b_y[i1];
      }
    }
  }
}

} // namespace fft
} // namespace internal
} // namespace coder
} // namespace coderTuner

//
// File trailer for FFTImplementationCallback.cpp
//
// [EOF]
//
