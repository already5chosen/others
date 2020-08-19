#include <float.h>
#include <math.h>
#include "chol_internal_definitions1.h"

enum {SIMD_N=4};
typedef struct { double v[2]; } s_m128d;

static
void InnerLoop(complex_sm256d *rowR, complex_sm256d* rowC, unsigned r, unsigned nRows, unsigned c)
{
  unsigned nRowQuartets = nRows / 4;
  unsigned cHalf = c / 2;
  unsigned rHalf = (r / 2) + 1;
  while (nRowQuartets)
  {
    // calculate eight results
    // multiply-add two dual-rows by dual-column
    complex_sm256d *rowR0 = rowR;
    complex_sm256d *rowR2 = rowR + rHalf;

    s_m256d inp_re0x   = rowR0[cHalf].re;                     // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    s_m256d sum_re0011 = {{inp_re0x.v[0],0,0,inp_re0x.v[3]}}; // re(r+0,c+0) 0           0           re(r+1,c+1)
    s_m256d sum_re0110 = {{0,inp_re0x.v[1],inp_re0x.v[2],0}}; // 0           re(r+0,c+1) re(r+1,c+0) 0
    s_m256d inp_re2x   = rowR2[cHalf].re;                     // re(r+2,c+0) re(r+2,c+1) re(r+3,c+0) re(r+3,c+1)
    s_m256d sum_re2031 = {{inp_re2x.v[0],0,0,inp_re2x.v[3]}}; // re(r+2,c+0) 0           0           re(r+3,c+1)
    s_m256d sum_re2130 = {{0,inp_re2x.v[1],inp_re2x.v[2],0}}; // 0           re(r+2,c+1) re(r+3,c+0) 0
    s_m256d inp_im0x   = rowR0[cHalf].im;                     // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    s_m256d sum_im0011 = {{inp_im0x.v[0],0,0,inp_im0x.v[3]}}; // im(r+0,c+0) 0           0           im(r+1,c+1)
    s_m256d sum_im0110 = {{0,inp_im0x.v[1],inp_im0x.v[2],0}}; // 0           im(r+0,c+1) im(r+1,c+0) 0
    s_m256d inp_im2x   = rowR2[cHalf].im;                     // im(r+2,c+0) im(r+2,c+1) im(r+3,c+0) im(r+3,c+1)
    s_m256d sum_im2031 = {{inp_im2x.v[0],0,0,inp_im2x.v[3]}}; // im(r+2,c+0) 0           0           im(r+3,c+1)
    s_m256d sum_im2130 = {{0,inp_im2x.v[1],inp_im2x.v[2],0}}; // 0           im(r+2,c+1) im(r+3,c+0) 0

    unsigned kk = cHalf;
    s_m256d reC01 = rowC->re;  // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
    s_m256d imC01 = rowC->im;  // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
    s_m256d rval0 = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
    s_m256d rval2 = rowR2->re; // re(r+2,k*2+0) re(r+2,k*2+1) re(r+3,k*2+0) re(r+3,k*2+1)
    do { // multiply-add two dual-rows by dual-column
      //s_m256d reC01 = rowC->re;                                // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      //s_m256d imC01 = rowC->im;                                // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      s_m256d reC10 = {{reC01.v[2],reC01.v[3],reC01.v[0],reC01.v[1]}}; // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
      s_m256d imC10 = {{imC01.v[2],imC01.v[3],imC01.v[0],imC01.v[1]}}; // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
      //s_m256d rval0, rval2;

      //rval0 = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= rval0.v[k]*reC01.v[k];
      // 0: re(r+0,c+0) -= re(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= re(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= re(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: re(r+1,c+1) -= re(r+1,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0011.v[k] += rval0.v[k]*imC01.v[k];
      // 0: im(r+0,c+0) += re(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           += re(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           += re(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: im(r+1,c+1) += re(r+1,k*2+1)*im(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= rval0.v[k]*reC10.v[k];
      // 0: 0           -= re(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: re(r+0,c+1) -= re(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: re(r+1,c+0) -= re(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= re(r+1,k*2+1)*re(c+0,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] += rval0.v[k]*imC10.v[k];
      // 0: 0           += re(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: im(r+0,c+1) += re(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: im(r+1,c+0) += re(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           += re(r+1,k*2+1)*im(c+0,k*2+1)

      //rval2 = rowR2->re; // re(r+2,k*2+0) re(r+2,k*2+1) re(r+3,k*2+0) re(r+3,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re2031.v[k] -= rval2.v[k]*reC01.v[k];
      // 0: re(r+2,c+0) -= re(r+2,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= re(r+2,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= re(r+3,k*2+0)*re(c+1,k*2+0)
      // 3: re(r+3,c+1) -= re(r+3,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im2031.v[k] += rval2.v[k]*imC01.v[k];
      // 0: im(r+2,c+0) += re(r+2,k*2+0)*im(c+0,k*2+0)
      // 1: 0           += re(r+2,k*2+1)*im(c+0,k*2+1)
      // 2: 0           += re(r+3,k*2+0)*im(c+1,k*2+0)
      // 3: im(r+3,c+1) += re(r+3,k*2+1)*im(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re2130.v[k] -= rval2.v[k]*reC10.v[k];
      // 0: 0           -= re(r+2,k*2+0)*re(c+0,k*2+0)
      // 1: re(r+3,c+1) -= re(r+2,k*2+1)*re(c+0,k*2+1)
      // 2: re(r+2,c+0) -= re(r+3,k*2+0)*re(c+1,k*2+0)
      // 3: 0           -= re(r+3,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im2130.v[k] += rval2.v[k]*imC10.v[k];
      // 0: 0           += re(r+2,k*2+0)*im(c+0,k*2+0)
      // 1: im(r+2,c+1) += re(r+2,k*2+1)*im(c+0,k*2+1)
      // 2: im(r+3,c+0) += re(r+3,k*2+0)*im(c+1,k*2+0)
      // 3: 0           += re(r+3,k*2+1)*im(c+1,k*2+1)

      rval0 = rowR0->im; // im(r+0,k*2+0) im(r+0,k*2+1) im(r+1,k*2+0) im(r+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= rval0.v[k]*imC01.v[k];
      // 0: re(r+0,c+0) -= im(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: re(r+1,c+1) -= im(r+1,k*2+1)*im(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0011.v[k] -= rval0.v[k]*reC01.v[k];
      // 0: im(r+0,c+0) -= im(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: im(r+1,c+1) -= im(r+1,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= rval0.v[k]*imC10.v[k];
      // 0: 0           -= im(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: re(r+0,c+1) -= im(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: re(r+1,c+0) -= im(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*im(c+0,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] -= rval0.v[k]*reC10.v[k];
      // 0: 0           -= im(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: im(r+0,c+1) -= im(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: im(r+1,c+0) -= im(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*re(c+0,k*2+1)


      rval2 = rowR2->im; // im(r+2,k*2+0) im(r+2,k*2+1) im(r+3,k*2+0) im(r+3,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re2031.v[k] -= rval2.v[k]*imC01.v[k];
      // 0: re(r+2,c+0) -= im(r+2,k*2+0)*im(c+0,k*2+0)
      // 1: 0           -= im(r+2,k*2+1)*im(c+0,k*2+1)
      // 2: 0           -= im(r+3,k*2+0)*im(c+1,k*2+0)
      // 3: re(r+3,c+1) -= im(r+3,k*2+1)*im(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im2031.v[k] -= rval2.v[k]*reC01.v[k];
      // 0: im(r+2,c+0) -= im(r+2,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= im(r+2,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= im(r+3,k*2+0)*re(c+1,k*2+0)
      // 3: im(r+3,c+1) -= im(r+3,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re2130.v[k] -= rval2.v[k]*imC10.v[k];
      // 0: 0           -= im(r+2,k*2+0)*im(c+1,k*2+0)
      // 1: re(r+3,c+1) -= im(r+2,k*2+1)*im(c+1,k*2+1)
      // 2: re(r+2,c+0) -= im(r+3,k*2+0)*im(c+0,k*2+0)
      // 3: 0           -= im(r+3,k*2+1)*im(c+0,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im2130.v[k] -= rval2.v[k]*reC10.v[k];
      // 0: 0           -= im(r+2,k*2+0)*re(c+1,k*2+0)
      // 1: im(r+2,c+1) -= im(r+2,k*2+1)*re(c+1,k*2+1)
      // 2: im(r+3,c+0) -= im(r+3,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= im(r+3,k*2+1)*re(c+0,k*2+1)

      rowC += 1;
      rowR0 += 1;
      rowR2 += 1;
      reC01 = rowC->re;  // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      imC01 = rowC->im;  // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      rval0 = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
      rval2 = rowR2->re; // re(r+2,k*2+0) re(r+2,k*2+1) re(r+3,k*2+0) re(r+3,k*2+1)
    } while (--kk);

    {
    // sum up odd and even sub-sums
    s_m256d sum_reim0011 = {{
      sum_re0011.v[0]+sum_re0011.v[1], sum_im0011.v[0]+sum_im0011.v[1],
      sum_re0011.v[2]+sum_re0011.v[3], sum_im0011.v[2]+sum_im0011.v[3]}};// re(r+0,c+0) im(r+0,c+0) re(r+1,c+1) im(r+1,c+1)
    s_m256d sum_reim2031 = {{
      sum_re2031.v[0]+sum_re2031.v[1], sum_im2031.v[0]+sum_im2031.v[1],
      sum_re2031.v[2]+sum_re2031.v[3], sum_im2031.v[2]+sum_im2031.v[3]}};// re(r+2,c+0) im(r+2,c+0) re(r+3,c+1) im(r+3,c+1)
    s_m256d sum_reim0110 = {{
      sum_re0110.v[0]+sum_re0110.v[1], sum_im0110.v[0]+sum_im0110.v[1],
      sum_re0110.v[2]+sum_re0110.v[3], sum_im0110.v[2]+sum_im0110.v[3]}};// re(r+0,c+1) im(r+0,c+1) re(r+1,c+0) im(r+1,c+0)
    s_m256d sum_reim2130 = {{
      sum_re2130.v[0]+sum_re2130.v[1], sum_im2130.v[0]+sum_im2130.v[1],
      sum_re2130.v[2]+sum_re2130.v[3], sum_im2130.v[2]+sum_im2130.v[3]}};// re(r+2,c+1) im(r+2,c+1) re(r+3,c+0) im(r+3,c+0)
    // separate column c+0 from column c+1
    s_m256d sum_reim00 = {{sum_reim0011.v[0],sum_reim0011.v[1],sum_reim0110.v[2],sum_reim0110.v[3]}}; // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    s_m256d sum_reim01 = {{sum_reim0110.v[0],sum_reim0110.v[1],sum_reim0011.v[2],sum_reim0011.v[3]}}; // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)
    s_m256d sum_reim20 = {{sum_reim2031.v[0],sum_reim2031.v[1],sum_reim2130.v[2],sum_reim2130.v[3]}}; // re(r+2,c+0) im(r+2,c+0) re(r+3,c+0) im(r+3,c+0)
    s_m256d sum_reim21 = {{sum_reim2130.v[0],sum_reim2130.v[1],sum_reim2031.v[2],sum_reim2031.v[3]}}; // re(r+2,c+1) im(r+2,c+1) re(r+3,c+1) im(r+3,c+1)

    // calculate re/im(r+0:r+3, c+0)
    double aaInvSqrt0 = rowC->im.v[0];  // 1/re(c+0,c+0)
    s_m256d reim00, reim20;
    for (int k=0; k < SIMD_N; ++k) {
      reim00.v[k] = sum_reim00.v[k]*aaInvSqrt0; // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
      reim20.v[k] = sum_reim20.v[k]*aaInvSqrt0; // re(r+2,c+0) im(r+2,c+0) re(r+3,c+0) im(r+3,c+0)
    }
    // add elements (r+0:r+3,c+0) to dot product of elements (r+0:r+3,c+1)
    double reC10Last = rowC->re.v[2]; // re(c+1,c+0)
    for (int k=0; k < SIMD_N; ++k) sum_reim01.v[k] -= reim00.v[k]*reC10Last;
    // 0: re(r+0,c+1) -= re(r+0,c+0)*re(c+1,c+0)
    // 1: im(r+0,c+1) -= im(r+0,c+0)*re(c+1,c+0)
    // 2: re(r+1,c+1) -= re(r+1,c+0)*re(c+1,c+0)
    // 3: im(r+1,c+1) -= im(r+1,c+0)*re(c+1,c+0)
    for (int k=0; k < SIMD_N; ++k) sum_reim21.v[k] -= reim20.v[k]*reC10Last;
    // 0: re(r+2,c+1) -= re(r+2,c+0)*re(c+1,c+0)
    // 1: im(r+2,c+1) -= im(r+2,c+0)*re(c+1,c+0)
    // 2: re(r+3,c+1) -= re(r+3,c+0)*re(c+1,c+0)
    // 3: im(r+3,c+1) -= im(r+3,c+0)*re(c+1,c+0)
    {
    double imC10Last = rowC->im.v[2]; // im(c+1,c+0)
    sum_reim01.v[0] -= reim00.v[1]*imC10Last; // 0: re(r+0,c+1) -= im(r+0,c+0)*im(c+1,c+0)
    sum_reim01.v[1] += reim00.v[0]*imC10Last; // 1: im(r+0,c+1) += re(r+0,c+0)*im(c+1,c+0)
    sum_reim01.v[2] -= reim00.v[3]*imC10Last; // 2: re(r+1,c+1) -= im(r+1,c+0)*im(c+1,c+0)
    sum_reim01.v[3] += reim00.v[2]*imC10Last; // 3: im(r+1,c+1) += re(r+1,c+0)*im(c+1,c+0)

    sum_reim21.v[0] -= reim20.v[1]*imC10Last; // 0: re(r+2,c+1) -= im(r+2,c+0)*im(c+1,c+0)
    sum_reim21.v[1] += reim20.v[0]*imC10Last; // 1: im(r+2,c+1) += re(r+2,c+0)*im(c+1,c+0)
    sum_reim21.v[2] -= reim20.v[3]*imC10Last; // 2: re(r+3,c+1) -= im(r+3,c+0)*im(c+1,c+0)
    sum_reim21.v[3] += reim20.v[2]*imC10Last; // 3: im(r+3,c+1) += re(r+3,c+0)*im(c+1,c+0)
    {
    double aaInvSqrt1 = rowC->im.v[3];  // 1/re(c+1,c+1)
    s_m256d reim01, reim21;
    for (int k=0; k < SIMD_N; ++k) {
      reim01.v[k] = sum_reim01.v[k]*aaInvSqrt1; // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)
      reim21.v[k] = sum_reim21.v[k]*aaInvSqrt1; // re(r+2,c+1) im(r+2,c+1) re(r+3,c+1) im(r+3,c+1)
    }
    // store results
    s_m256d res01re = {{reim00.v[0], reim01.v[0],reim00.v[2], reim01.v[2]}}; // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    s_m256d res01im = {{reim00.v[1], reim01.v[1],reim00.v[3], reim01.v[3]}}; // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    rowR0->re = res01re; // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    rowR0->im = res01im; // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    s_m256d res23re = {{reim20.v[0], reim21.v[0],reim20.v[2], reim21.v[2]}}; // re(r+2,c+0) re(r+2,c+1) re(r+3,c+0) re(r+3,c+1)
    s_m256d res23im = {{reim20.v[1], reim21.v[1],reim20.v[3], reim21.v[3]}}; // im(r+2,c+0) im(r+2,c+1) im(r+3,c+0) im(r+3,c+1)
    rowR2->re = res23re; // re(r+2,c+0) re(r+2,c+1) re(r+3,c+0) re(r+3,c+1)
    rowR2->im = res23im; // im(r+2,c+0) im(r+2,c+1) im(r+3,c+0) im(r+3,c+1)
    }
    }
    }
    rowC  -= cHalf;
    rowR  += rHalf*2+1;
    rHalf += 2;
    --nRowQuartets;
  }

  if (0 != (nRows & 2))
  {
    // calculate four results
    // multiply-add one dual-rows by one dual-column
    complex_sm256d *rowR0 = rowR;

    s_m256d inp_re0x   = rowR0[cHalf].re;                     // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    s_m256d sum_re0011 = {{inp_re0x.v[0],0,0,inp_re0x.v[3]}}; // re(r+0,c+0) 0           0           re(r+1,c+1)
    s_m256d sum_re0110 = {{0,inp_re0x.v[1],inp_re0x.v[2],0}}; // 0           re(r+0,c+1) re(r+1,c+0) 0
    s_m256d inp_im0x   = rowR0[cHalf].im;                     // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    s_m256d sum_im0011 = {{inp_im0x.v[0],0,0,inp_im0x.v[3]}}; // im(r+0,c+0) 0           0           im(r+1,c+1)
    s_m256d sum_im0110 = {{0,inp_im0x.v[1],inp_im0x.v[2],0}}; // 0           im(r+0,c+1) im(r+1,c+0) 0

    unsigned kk = cHalf;
    do { // multiply-add two dual-rows by dual-column
      s_m256d reC01 = rowC->re;                                        // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      s_m256d imC01 = rowC->im;                                        // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      s_m256d reC10 = {{reC01.v[2],reC01.v[3],reC01.v[0],reC01.v[1]}}; // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
      s_m256d imC10 = {{imC01.v[2],imC01.v[3],imC01.v[0],imC01.v[1]}}; // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
      s_m256d rval;

      rval = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= rval.v[k]*reC01.v[k];
      // 0: re(r+0,c+0) -= re(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= re(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= re(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: re(r+1,c+1) -= re(r+1,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0011.v[k] += rval.v[k]*imC01.v[k];
      // 0: im(r+0,c+0) += re(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           += re(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           += re(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: im(r+1,c+1) += re(r+1,k*2+1)*im(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= rval.v[k]*reC10.v[k];
      // 0: 0           -= re(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: re(r+0,c+1) -= re(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: re(r+1,c+0) -= re(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= re(r+1,k*2+1)*re(c+0,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] += rval.v[k]*imC10.v[k];
      // 0: 0           += re(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: im(r+0,c+1) += re(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: im(r+1,c+0) += re(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           += re(r+1,k*2+1)*im(c+0,k*2+1)

      rval = rowR0->im; // im(r+0,k*2+0) im(r+0,k*2+1) im(r+1,k*2+0) im(r+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= rval.v[k]*imC01.v[k];
      // 0: re(r+0,c+0) -= im(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: re(r+1,c+1) -= im(r+1,k*2+1)*im(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0011.v[k] -= rval.v[k]*reC01.v[k];
      // 0: im(r+0,c+0) -= im(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: im(r+1,c+1) -= im(r+1,k*2+1)*re(c+1,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= rval.v[k]*imC10.v[k];
      // 0: 0           -= im(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: re(r+0,c+1) -= im(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: re(r+1,c+0) -= im(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*im(c+0,k*2+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] -= rval.v[k]*reC10.v[k];
      // 0: 0           -= im(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: im(r+0,c+1) -= im(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: im(r+1,c+0) -= im(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*re(c+0,k*2+1)

      rowC += 1;
      rowR0 += 1;
    } while (--kk);

    {
    // sum up odd and even sub-sums
    s_m256d sum_reim0011 = {{
      sum_re0011.v[0]+sum_re0011.v[1], sum_im0011.v[0]+sum_im0011.v[1],
      sum_re0011.v[2]+sum_re0011.v[3], sum_im0011.v[2]+sum_im0011.v[3]}};// re(r+0,c+0) im(r+0,c+0) re(r+1,c+1) im(r+1,c+1)
    s_m256d sum_reim0110 = {{
      sum_re0110.v[0]+sum_re0110.v[1], sum_im0110.v[0]+sum_im0110.v[1],
      sum_re0110.v[2]+sum_re0110.v[3], sum_im0110.v[2]+sum_im0110.v[3]}};// re(r+0,c+1) im(r+0,c+1) re(r+1,c+0) im(r+1,c+0)
    // separate column c+0 from column c+1
    s_m256d sum_reim00 = {{sum_reim0011.v[0],sum_reim0011.v[1],sum_reim0110.v[2],sum_reim0110.v[3]}}; // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    s_m256d sum_reim01 = {{sum_reim0110.v[0],sum_reim0110.v[1],sum_reim0011.v[2],sum_reim0011.v[3]}}; // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)

    // calculate re/im(r+0:r+1, c+0)
    double aaInvSqrt0 = rowC->im.v[0];  // 1/re(c+0,c+0)
    s_m256d reim00;
    for (int k=0; k < SIMD_N; ++k)
      reim00.v[k] = sum_reim00.v[k]*aaInvSqrt0; // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    {
    // add elements (r+0:r+3,c+0) to dot product of elements (r+0:r+3,c+1)
    double reC10Last = rowC->re.v[2]; // re(c+1,c+0)
    for (int k=0; k < SIMD_N; ++k) sum_reim01.v[k] -= reim00.v[k]*reC10Last;
    // 0: re(r+0,c+1) -= re(r+0,c+0)*re(c+1,c+0)
    // 1: im(r+0,c+1) -= im(r+0,c+0)*re(c+1,c+0)
    // 2: re(r+1,c+1) -= re(r+1,c+0)*re(c+1,c+0)
    // 3: im(r+1,c+1) -= im(r+1,c+0)*re(c+1,c+0)
    {
    double imC10Last = rowC->im.v[2]; // im(c+1,c+0)
    sum_reim01.v[0] -= reim00.v[1]*imC10Last; // 0: re(r+0,c+1) -= im(r+0,c+0)*im(c+1,c+0)
    sum_reim01.v[1] += reim00.v[0]*imC10Last; // 1: im(r+0,c+1) += re(r+0,c+0)*im(c+1,c+0)
    sum_reim01.v[2] -= reim00.v[3]*imC10Last; // 2: re(r+1,c+1) -= im(r+1,c+0)*im(c+1,c+0)
    sum_reim01.v[3] += reim00.v[2]*imC10Last; // 3: im(r+1,c+1) += re(r+1,c+0)*im(c+1,c+0)
    {
    double aaInvSqrt1 = rowC->im.v[3];  // 1/re(c+1,c+1)
    s_m256d reim01;
    for (int k=0; k < SIMD_N; ++k)
      reim01.v[k] = sum_reim01.v[k]*aaInvSqrt1; // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)
    // store results
    s_m256d res01re = {{reim00.v[0], reim01.v[0],reim00.v[2], reim01.v[2]}}; // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    s_m256d res01im = {{reim00.v[1], reim01.v[1],reim00.v[3], reim01.v[3]}}; // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    rowR0->re = res01re; // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    rowR0->im = res01im; // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    }
    }
    }
    }
    rowC  -= cHalf;
    rowR  += rHalf;
    rHalf += 1;
  }

  if (0 != (nRows & 1))
  { // process 2 columns of last row. Last row is non-interleaved
    // calculate two results
    // multiply-add one single row by one dual-column
    s_m256d *rowR0 = &rowR->re;

    s_m256d inp_reim   = rowR0[cHalf];                         // re(r,c+0) re(r,c+1) im(r,c+0) im(r,c+1)
    s_m256d sum_re0im1a = {{inp_reim.v[0],0,0,inp_reim.v[3]}}; // re(r,c+0) 0         0         im(r,c+1)
    s_m256d sum_im0re1ax= {{0,inp_reim.v[1],inp_reim.v[2],0}}; // 0         re(r,c+1) im(r,c+0) 0
    s_m256d sum_re0im1b = {{0}};
    s_m256d sum_im0re1b = {{0}};

    unsigned kk = cHalf;
    s_m256d sum_im0re1a = {{sum_im0re1ax.v[2],sum_im0re1ax.v[3],sum_im0re1ax.v[0],sum_im0re1ax.v[1]}}; // im(r,c+0) 0         0         re(r,c+1)
    do { // multiply-add two dual-rows by dual-column
      s_m256d reC01 = rowC->re;                                        // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      s_m256d imC01 = rowC->im;                                        // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      s_m256d reimR = *rowR0;                                          // re(r+0,k*2+0) re(r+0,k*2+1) im(r+0,k*2+0) im(r+0,k*2+1)
      s_m256d imreR = {{reimR.v[2],reimR.v[3],reimR.v[0],reimR.v[1]}}; // im(r+0,k*2+0) im(r+0,k*2+1) re(r+0,k*2+0) re(r+0,k*2+1)

      for (int k=0; k < SIMD_N; ++k) sum_re0im1a.v[k] -= reimR.v[k]*reC01.v[k];
      // 0: re(r,c+0) -= re(r,k*2+0)*re(c+0,k*2+0) : +re(r,c+0)
      // 1: 0         -= re(r,k*2+1)*re(c+0,k*2+1) : +re(r,c+0)
      // 2: 0         -= im(r,k*2+0)*re(c+1,k*2+0) : +im(r,c+1)
      // 3: im(r,c+1) -= im(r,k*2+1)*re(c+1,k*2+1) : +im(r,c+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0re1a.v[k] -= imreR.v[k]*reC01.v[k];
      // 0: 0         -= im(r,k*2+0)*re(c+0,k*2+0) : +im(r,c+0)
      // 1: 0         -= im(r,k*2+1)*re(c+0,k*2+1) : +im(r,c+0)
      // 2: 0         -= re(r,k*2+0)*re(c+1,k*2+0) : +re(r,c+1)
      // 3: 0         -= re(r,k*2+1)*re(c+1,k*2+1) : +re(r,c+1)
      for (int k=0; k < SIMD_N; ++k) sum_im0re1b.v[k] -= reimR.v[k]*imC01.v[k];
      // 0: 0         -= re(r,k*2+0)*im(c+0,k*2+0) : -im(r,c+0)
      // 1: 0         -= re(r,k*2+1)*im(c+0,k*2+1) : -im(r,c+0)
      // 2: 0         -= im(r,k*2+0)*im(c+1,k*2+0) : +re(r,c+1)
      // 3: 0         -= im(r,k*2+1)*im(c+1,k*2+1) : +re(r,c+1)
      for (int k=0; k < SIMD_N; ++k) sum_re0im1b.v[k] -= imreR.v[k]*imC01.v[k];
      // 0: 0         += im(r,k*2+0)*im(c+0,k*2+0) : +re(r,c+0)
      // 1: 0         += im(r,k*2+1)*im(c+0,k*2+1) : +re(r,c+0)
      // 2: 0         += re(r,k*2+0)*im(c+1,k*2+0) : -im(r,c+1)
      // 3: 0         += re(r,k*2+1)*im(c+1,k*2+1) : -im(r,c+1)

      rowC  += 1;
      rowR0 += 1;
    } while (--kk);

    {
    // sum up odd and even sub-sums
    s_m256d sum_imre_ax = {{
      sum_im0re1a.v[0]+sum_im0re1a.v[1], sum_re0im1a.v[0]+sum_re0im1a.v[1],
      sum_im0re1a.v[2]+sum_im0re1a.v[3], sum_re0im1a.v[2]+sum_re0im1a.v[3]}}; // +im(r,c+0) +re(r,c+0) +re(r,c+1) +im(r,c+1)
    s_m256d sum_imre_bx = {{
      sum_im0re1b.v[0]+sum_im0re1b.v[1], sum_re0im1b.v[0]+sum_re0im1b.v[1],
      sum_im0re1b.v[2]+sum_im0re1b.v[3], sum_re0im1b.v[2]+sum_re0im1b.v[3]}}; // -im(r,c+0) +re(r,c+0) +re(r,c+1) -im(r,c+1)

    // sum up partial sums
    s_m256d sum_imre_a = {{
      sum_imre_ax.v[0],sum_imre_ax.v[1],
      sum_imre_ax.v[3],sum_imre_ax.v[2]}}; // +im(r,c+0) +re(r,c+0) +im(r,c+1) +re(r,c+1)
    s_m256d sum_imre_b = {{
      sum_imre_bx.v[0],sum_imre_bx.v[1],
      sum_imre_bx.v[3],sum_imre_bx.v[2]}}; // -im(r,c+0) +re(r,c+0) -im(r,c+1) +re(r,c+1)
    {
      s_m256d sum_imre = {{
       sum_imre_a.v[0] - sum_imre_b.v[0],  // +im(r,c+0)
       sum_imre_a.v[1] + sum_imre_b.v[1],  // +re(r,c+0)
       sum_imre_a.v[2] - sum_imre_b.v[2],  // +im(r,c+1)
       sum_imre_a.v[3] + sum_imre_b.v[3]}};// +re(r,c+1)

    // calculate re/im(r, c+0)
    double aaInvSqrt0 = rowC->im.v[0];  // 1/re(c+0,c+0)
    s_m128d imre0 = {{sum_imre.v[0]*aaInvSqrt0,sum_imre.v[1]*aaInvSqrt0}}; // im(r,c+0) re(r,c+0)

    // add elements (r,c+0) to dot product of elements (r,c+1)
    double reC10Last = rowC->re.v[2]; // re(c+1,c+0)
    s_m128d sum_imre1 = {{sum_imre.v[2],sum_imre.v[3]}};         // im(r,c+1) re(r,c+1)
    sum_imre1.v[0] -= imre0.v[0]*reC10Last;
    sum_imre1.v[1] -= imre0.v[1]*reC10Last;
    {
    double imC10Last = rowC->im.v[2]; // im(c+1,c+0)
    s_m128d sum_reim1 = {{sum_imre1.v[1],sum_imre1.v[0]}}; // re(r,c+1) im(r,c+1)
    sum_reim1.v[0] -= imre0.v[0]*imC10Last;
    sum_reim1.v[1] += imre0.v[1]*imC10Last;
    {
    double aaInvSqrt1 = rowC->im.v[3];// 1/re(c+1,c+1)
    s_m128d reim1 = {{sum_reim1.v[0]*aaInvSqrt1,sum_reim1.v[1]*aaInvSqrt1}}; // re(r,c+1) im(r,c+1)
    // store results
    rowR0->v[0] = imre0.v[1]; // re(r,c+0)
    rowR0->v[1] = reim1.v[0]; // re(r,c+1)
    rowR0->v[2] = imre0.v[0]; // im(r,c+0)
    rowR0->v[3] = reim1.v[1]; // im(r,c+1)
    }
    }
    }
    }
  }
}

//
// row1st - always even
// return - 0 - matrix non positive defined, 1 - matrix positive defined
int chol_CrBa4_i2_Band_Factorize(complex_sm256d* triang, unsigned row1st, unsigned rowLast)
{
  const double DIAG_MIN   =  (double)(FLT_MIN);
  // FLT_MIN is not really special in context of double-precision calculations
  // but I wanted very small number that is still much bigger than DBL_MIN
  // and FLT_MIN looks good enough
  const double DIAG_SUBST =  1E40;

  int succ = 1;
  unsigned rHalf1st = row1st/2;
  complex_sm256d *rowDataset = &triang[rHalf1st*(rHalf1st+1)/2];

  // process Crout band
  // columns 0:1
  {
    double aaInvSqrt0 = ((double*)&triang[0].im)[0];
    double re10       = ((double*)&triang[0].re)[2];
    double im10       = ((double*)&triang[0].im)[2];
    double aaInvSqrt1 = ((double*)&triang[0].im)[3];

    complex_sm256d *rowR = rowDataset;
    unsigned rLast = rowLast & (-2);
    unsigned r;
    for (r = row1st; r < rLast; r += 2)
    {
      // multiply section of column 0 by 1/diagonal
      double reR00 = ((double*)&rowR->re)[0] * aaInvSqrt0;
      double imR00 = ((double*)&rowR->im)[0] * aaInvSqrt0;
      double reR10 = ((double*)&rowR->re)[2] * aaInvSqrt0;
      double imR10 = ((double*)&rowR->im)[2] * aaInvSqrt0;
      double reR01 = ((double*)&rowR->re)[1];
      double imR01 = ((double*)&rowR->im)[1];
      double reR11 = ((double*)&rowR->re)[3];
      double imR11 = ((double*)&rowR->im)[3];
      reR01 -= reR00 * re10;
      imR01 -= imR00 * re10;
      reR11 -= reR10 * re10;
      imR11 -= imR10 * re10;
      reR01 -= imR00 * im10;
      imR01 += reR00 * im10;
      reR11 -= imR10 * im10;
      imR11 += reR10 * im10;
      ((double*)&rowR->re)[0] = reR00;
      ((double*)&rowR->re)[1] = reR01 * aaInvSqrt1;
      ((double*)&rowR->re)[2] = reR10;
      ((double*)&rowR->re)[3] = reR11 * aaInvSqrt1;
      ((double*)&rowR->im)[0] = imR00;
      ((double*)&rowR->im)[1] = imR01 * aaInvSqrt1;
      ((double*)&rowR->im)[2] = imR10;
      ((double*)&rowR->im)[3] = imR11 * aaInvSqrt1;
      rowR += r/2 + 1;
    }
    if (rowLast & 1)
    { // non-interleaved last line
      double reR00 = rowR->re.v[0] * aaInvSqrt0;
      double imR00 = rowR->re.v[2] * aaInvSqrt0;
      double reR01 = rowR->re.v[1];
      double imR01 = rowR->re.v[3];
      reR01 -= reR00 * re10;
      imR01 -= imR00 * re10;
      reR01 -= imR00 * im10;
      imR01 += reR00 * im10;
      rowR->re.v[0] = reR00;
      rowR->re.v[1] = reR01 * aaInvSqrt1;
      rowR->re.v[2] = imR00;
      rowR->re.v[3] = imR01 * aaInvSqrt1;
    }
  }

  {
  complex_sm256d *rowCx2 = &triang[1];
  unsigned c;
  // rectangular part
  for (c = 2; c < row1st; c += 2)
  {
    // calculate section of two columns of result (row1st:rowLast-1,c:c+1)
    InnerLoop(
      rowDataset,
      rowCx2,
      row1st, rowLast - row1st, c);
    rowCx2 += c/2 + 1;
  }

  // triangular part
  for (c = row1st; c+1 < rowLast; c += 2)
  {  // calculate section of two columns (c:rowLast-1,c:c+1)
    unsigned nRows = rowLast - c;
    unsigned cHalf = c/2;
    if ((nRows & 2)==0)
    { // process two dual-rows
      complex_sm256d *rowR0 = rowCx2;
      complex_sm256d *rowR2 = rowR0 + cHalf + 1;
      s_m256d inp_re0x   = rowR0[cHalf].re;                     // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      s_m256d sum_re0011 = {{inp_re0x.v[0],0,0,inp_re0x.v[3]}}; // re(c+0,c+0) 0           0           re(c+1,c+1)
      s_m256d sum_re0110 = {{0,            0,inp_re0x.v[2],0}}; // 0           0           re(c+1,c+0) 0
      s_m256d inp_re2x   = rowR2[cHalf].re;                     // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      s_m256d sum_re2031 = {{inp_re2x.v[0],0,0,inp_re2x.v[3]}}; // re(c+2,c+0) 0           0           re(c+3,c+1)
      s_m256d sum_re2130 = {{0,inp_re2x.v[1],inp_re2x.v[2],0}}; // 0           re(c+2,c+1) re(c+3,c+0) 0
      s_m256d inp_im0x   = rowR0[cHalf].im;                     // x           x           im(c+1,c+0) x
      s_m256d sum_im0110 = {{0,            0,inp_im0x.v[2],0}}; // 0           0           im(c+1,c+0) 0
      s_m256d inp_im2x   = rowR2[cHalf].im;                     // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      s_m256d sum_im2031 = {{inp_im2x.v[0],0,0,inp_im2x.v[3]}}; // im(c+2,c+0) 0           0           im(c+3,c+1)
      s_m256d sum_im2130 = {{0,inp_im2x.v[1],inp_im2x.v[2],0}}; // 0           im(c+2,c+1) im(c+3,c+0) 0

      unsigned kk = cHalf;
      do { // multiply-add two dual-rows by the first dual-row of the two
        s_m256d reC01 = rowR0->re;                                       // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
        s_m256d imC01 = rowR0->im;                                       // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
        s_m256d reC10 = {{reC01.v[2],reC01.v[3],reC01.v[0],reC01.v[1]}}; // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
        s_m256d imC10 = {{imC01.v[2],imC01.v[3],imC01.v[0],imC01.v[1]}}; // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
        s_m256d rval;

        for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= reC01.v[k]*reC01.v[k];
        // 0: re(c+0,c+0) -= re(c+0,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= re(c+1,k*2+0)*re(c+1,k*2+0)
        // 3: re(c+1,c+1) -= re(c+1,k*2+1)*re(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= reC01.v[k]*reC10.v[k];
        // 0: 0           -= re(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: re(c+1,c+0) -= re(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= re(c+1,k*2+1)*re(c+0,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] += reC01.v[k]*imC10.v[k];
        // 0: 0           += re(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           += re(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: im(c+1,c+0) += re(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           += re(c+1,k*2+1)*im(c+0,k*2+1)

        rval = rowR2->re; // re(c+2,k*2+0) re(c+2,k*2+1) re(c+3,k*2+0) re(c+3,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re2031.v[k] -= rval.v[k]*reC01.v[k];
        // 0: re(c+2,c+0) -= re(c+2,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= re(c+2,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= re(c+3,k*2+0)*re(c+1,k*2+0)
        // 3: re(c+3,c+1) -= re(c+3,k*2+1)*re(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im2031.v[k] += rval.v[k]*imC01.v[k];
        // 0: im(c+2,c+0) += re(c+2,k*2+0)*im(c+0,k*2+0)
        // 1: 0           += re(c+2,k*2+1)*im(c+0,k*2+1)
        // 2: 0           += re(c+3,k*2+0)*im(c+1,k*2+0)
        // 3: im(c+3,c+1) += re(c+3,k*2+1)*im(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re2130.v[k] -= rval.v[k]*reC10.v[k];
        // 0: 0           -= re(c+2,k*2+0)*re(c+1,k*2+0)
        // 1: re(c+3,c+1) -= re(c+2,k*2+1)*re(c+1,k*2+1)
        // 2: re(c+2,c+0) -= re(c+3,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= re(c+3,k*2+1)*re(c+0,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im2130.v[k] += rval.v[k]*imC10.v[k];
        // 0: 0           += re(c+2,k*2+0)*im(c+1,k*2+0)
        // 1: im(c+2,c+1) += re(c+2,k*2+1)*im(c+1,k*2+1)
        // 2: im(c+3,c+0) += re(c+3,k*2+0)*im(c+0,k*2+0)
        // 3: 0           += re(c+3,k*2+1)*im(c+0,k*2+1)

        for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= imC01.v[k]*imC01.v[k];
        // 0: re(c+0,c+0) -= im(c+0,k*2+0)*im(c+0,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+0,k*2+1)
        // 2: 0           -= im(c+1,k*2+0)*im(c+1,k*2+0)
        // 3: re(c+1,c+1) -= im(c+1,k*2+1)*im(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= imC01.v[k]*imC10.v[k];
        // 0: 0           -= im(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: re(c+1,c+0) -= im(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*im(c+0,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] -= imC01.v[k]*reC10.v[k];
        // 0: 0           -= im(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: im(c+1,c+0) -= im(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*re(c+0,k*2+1)

        rval = rowR2->im; // im(c+2,k*2+0) im(c+2,k*2+1) im(c+3,k*2+0) im(c+3,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re2031.v[k] -= rval.v[k]*imC01.v[k];
        // 0: re(c+2,c+0) -= im(c+2,k*2+0)*im(c+0,k*2+0)
        // 1: 0           -= im(c+2,k*2+1)*im(c+0,k*2+1)
        // 2: 0           -= im(c+3,k*2+0)*im(c+1,k*2+0)
        // 3: re(c+3,c+1) -= im(c+3,k*2+1)*im(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im2031.v[k] -= rval.v[k]*reC01.v[k];
        // 0: im(c+2,c+0) -= im(c+2,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= im(c+2,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= im(c+3,k*2+0)*re(c+1,k*2+0)
        // 3: im(c+3,c+1) -= im(c+3,k*2+1)*re(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re2130.v[k] -= rval.v[k]*imC10.v[k];
        // 0: 0           -= im(c+2,k*2+0)*im(c+1,k*2+0)
        // 1: re(c+3,c+1) -= im(c+2,k*2+1)*im(c+1,k*2+1)
        // 2: re(c+2,c+0) -= im(c+3,k*2+0)*im(c+0,k*2+0)
        // 3: 0           -= im(c+3,k*2+1)*im(c+0,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im2130.v[k] -= rval.v[k]*reC10.v[k];
        // 0: 0           -= im(c+2,k*2+0)*re(c+1,k*2+0)
        // 1: im(c+2,c+1) -= im(c+2,k*2+1)*re(c+1,k*2+1)
        // 2: im(c+3,c+0) -= im(c+3,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= im(c+3,k*2+1)*re(c+0,k*2+1)

        rowR0 += 1;
        rowR2 += 1;
      } while (--kk);
      {
      // sum up odd and even sub-sums
      s_m256d sum_re0 = {{
        sum_re0011.v[0]+sum_re0011.v[1], sum_re0110.v[0]+sum_re0110.v[1],
        sum_re0011.v[2]+sum_re0011.v[3], sum_re0110.v[2]+sum_re0110.v[3]}};// re(c+0,c+0) x           re(c+1,c+1) re(c+1,c+0)
      s_m256d sum_im0 = {{
        sum_im0110.v[0]+sum_im0110.v[1], sum_im0110.v[0]+sum_im0110.v[1],
        sum_im0110.v[2]+sum_im0110.v[3], sum_im0110.v[2]+sum_im0110.v[3]}};// x           x           im(c+1,c+0) im(c+1,c+0)
      s_m256d sum_re2 = {{
        sum_re2031.v[0]+sum_re2031.v[1], sum_re2130.v[0]+sum_re2130.v[1],
        sum_re2031.v[2]+sum_re2031.v[3], sum_re2130.v[2]+sum_re2130.v[3]}};// re(c+2,c+0) re(c+2,c+1) re(c+3,c+1) re(c+3,c+0)
      s_m256d sum_im2 = {{
        sum_im2031.v[0]+sum_im2031.v[1], sum_im2130.v[0]+sum_im2130.v[1],
        sum_im2031.v[2]+sum_im2031.v[3], sum_im2130.v[2]+sum_im2130.v[3]}};// im(c+2,c+0) im(c+2,c+1) im(c+3,c+1) im(c+3,c+0)
      // store temporary results
      s_m256d sum_re0x = {{sum_re0.v[0],sum_re0.v[0],sum_re0.v[3],sum_re0.v[2]}}; // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      rowR0->re = sum_re0x;                                                       // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      rowR0->im = sum_im0;                                                        // x           x           im(c+1,c+0) im(c+1,c+0)
      s_m256d sum_re2x = {{sum_re2.v[0],sum_re2.v[1],sum_re2.v[3],sum_re2.v[2]}}; // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      rowR2->re = sum_re2x;                                                       // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      s_m256d sum_im2x = {{sum_im2.v[0],sum_im2.v[1],sum_im2.v[3],sum_im2.v[2]}}; // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      rowR2->im = sum_im2x;                                                       // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      }
    }
    else
    { // process one dual-row
      complex_sm256d *rowR0 = rowCx2;
      s_m256d inp_re0x   = rowR0[cHalf].re;                     // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      s_m256d sum_re0011 = {{inp_re0x.v[0],0,0,inp_re0x.v[3]}}; // re(c+0,c+0) 0           0           re(c+1,c+1)
      s_m256d sum_re0110 = {{0,            0,inp_re0x.v[2],0}}; // 0           0           re(c+1,c+0) 0
      s_m256d inp_im0x   = rowR0[cHalf].im;                     // x           x           im(c+1,c+0) x
      s_m256d sum_im0110 = {{0,            0,inp_im0x.v[2],0}}; // 0           0           im(c+1,c+0) 0

      unsigned kk = cHalf;
      do { // multiply-add two dual-rows by the first dual-row of the two
        s_m256d reC01 = rowR0->re;                                       // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
        s_m256d imC01 = rowR0->im;                                       // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
        s_m256d reC10 = {{reC01.v[2],reC01.v[3],reC01.v[0],reC01.v[1]}}; // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
        s_m256d imC10 = {{imC01.v[2],imC01.v[3],imC01.v[0],imC01.v[1]}}; // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)

        for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= reC01.v[k]*reC01.v[k];
        // 0: re(c+0,c+0) -= re(c+0,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= re(c+1,k*2+0)*re(c+1,k*2+0)
        // 3: re(c+1,c+1) -= re(c+1,k*2+1)*re(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= reC01.v[k]*reC10.v[k];
        // 0: 0           -= re(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: re(c+1,c+0) -= re(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= re(c+1,k*2+1)*re(c+0,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] += reC01.v[k]*imC10.v[k];
        // 0: 0           += re(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           += re(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: im(c+1,c+0) += re(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           += re(c+1,k*2+1)*im(c+0,k*2+1)

        for (int k=0; k < SIMD_N; ++k) sum_re0011.v[k] -= imC01.v[k]*imC01.v[k];
        // 0: re(c+0,c+0) -= im(c+0,k*2+0)*im(c+0,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+0,k*2+1)
        // 2: 0           -= im(c+1,k*2+0)*im(c+1,k*2+0)
        // 3: re(c+1,c+1) -= im(c+1,k*2+1)*im(c+1,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_re0110.v[k] -= imC01.v[k]*imC10.v[k];
        // 0: 0           -= im(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: re(c+1,c+0) -= im(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*im(c+0,k*2+1)
        for (int k=0; k < SIMD_N; ++k) sum_im0110.v[k] -= imC01.v[k]*reC10.v[k];
        // 0: 0           -= im(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: im(c+1,c+0) -= im(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*re(c+0,k*2+1)

        rowR0 += 1;
      } while (--kk);
      {
      // sum up odd and even sub-sums
      s_m256d sum_re0 = {{
        sum_re0011.v[0]+sum_re0011.v[1], sum_re0110.v[0]+sum_re0110.v[1],
        sum_re0011.v[2]+sum_re0011.v[3], sum_re0110.v[2]+sum_re0110.v[3]}};// re(c+0,c+0) x           re(c+1,c+1) re(c+1,c+0)
      s_m256d sum_im0 = {{
        sum_im0110.v[0]+sum_im0110.v[1], sum_im0110.v[0]+sum_im0110.v[1],
        sum_im0110.v[2]+sum_im0110.v[3], sum_im0110.v[2]+sum_im0110.v[3]}};// x           x           im(c+1,c+0) im(c+1,c+0)
      // store temporary results
      s_m256d sum_re0x = {{sum_re0.v[0],sum_re0.v[0],sum_re0.v[3],sum_re0.v[2]}}; // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      rowR0->re = sum_re0x;                                                       // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      rowR0->im = sum_im0;                                                        // x           x           im(c+1,c+0) im(c+1,c+0)
      }
    }
    {
    // calculate diagonal element (c,c)
    double aa = rowCx2[cHalf].re.v[0]; // current diagonal element (c,c)

    // check that we are positive defined
    //printf("%d %e\n", c, aa);
    if (aa < DIAG_MIN)
    {
      aa = DIAG_SUBST;
      succ = 0;
    }
    {
    double aaSqrt0    = sqrt(aa);
    double aaInvSqrt0 = 1.0 / aaSqrt0;
    rowCx2[cHalf].re.v[0] = aaSqrt0;
    rowCx2[cHalf].im.v[0] = aaInvSqrt0; // store inverse in image field of diagonal
    {
    double reR10 = rowCx2[cHalf].re.v[2] * aaInvSqrt0; // re(c+1,c)
    double imR10 = rowCx2[cHalf].im.v[2] * aaInvSqrt0; // im(c+1,c)
    rowCx2[cHalf].re.v[2] = reR10;
    rowCx2[cHalf].im.v[2] = imR10;

    // calculate diagonal element (c+1,c+1)
    aa = rowCx2[cHalf].re.v[3]; // re(c+1,c+1)
    aa -= reR10*reR10;
    aa -= imR10*imR10;

    // check that we are positive defined
    //printf("%d %e\n", c+1, aa);
    if (aa < DIAG_MIN)
    {
      aa = DIAG_SUBST;
      succ = 0;
    }
    {
    double aaSqrt1    = sqrt(aa);
    double aaInvSqrt1 = 1.0 / aaSqrt1;
    rowCx2[cHalf].re.v[3] = aaSqrt1;
    rowCx2[cHalf].im.v[3] = aaInvSqrt1; // store inverse in image field of diagonal
    }
    }
    }

    {
    unsigned r = c + 2;
    complex_sm256d *rowR0 = rowCx2 + cHalf + 1;
    if ((nRows & 2)==0)
    {
      // complete calculation of re/im(c+2:c+3,c+0:c+1)
      double aaInvSqrt0 = rowCx2[cHalf].im.v[0];       // 1/re(c+0,c+0)
      s_m256d reR2 = rowR0[cHalf].re;                  // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      s_m256d imR2 = rowR0[cHalf].im;                  // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      s_m256d reimR20 = {{
        reR2.v[0]*aaInvSqrt0,  imR2.v[0]*aaInvSqrt0,
        reR2.v[2]*aaInvSqrt0,  imR2.v[2]*aaInvSqrt0}}; // re(c+2,c+0) im(c+2,c+0) re(c+3,c+0) im(c+3,c+0)
      s_m256d reimR21 = {{
        reR2.v[1],  imR2.v[1],
        reR2.v[3],  imR2.v[3]}};                       // re(c+2,c+1) im(c+2,c+1) re(c+3,c+1) im(c+3,c+1)
      double reC = rowCx2[cHalf].re.v[2];              // re(c+1,c+0)
      for (int k=0; k < SIMD_N; ++k) reimR21.v[k] -= reimR20.v[k]*reC;
      {
      double imC = rowCx2[cHalf].im.v[2];              // im(c+1,c+0)
      reimR21.v[0] -= reimR20.v[1] * imC;
      reimR21.v[1] += reimR20.v[0] * imC;
      reimR21.v[2] -= reimR20.v[3] * imC;
      reimR21.v[3] += reimR20.v[2] * imC;
      }
      {
      double aaInvSqrt1 = rowCx2[cHalf].im.v[3];  // 1/re(c+1,c+1)
      for (int k=0; k < SIMD_N; ++k) reimR21.v[k] *= aaInvSqrt1;
      }
      s_m256d re2021 = {{reimR20.v[0],reimR21.v[0],reimR20.v[2],reimR21.v[2]}};
      s_m256d im2021 = {{reimR20.v[1],reimR21.v[1],reimR20.v[3],reimR21.v[3]}};
      rowR0[cHalf].re = re2021;
      rowR0[cHalf].im = im2021;

      rowR0 += cHalf + 2;
      r = c + 4;
      nRows -= 2;
    }
    nRows -= 2;

    if (nRows != 0)
      InnerLoop(rowR0, rowCx2, r, nRows, c);
    rowCx2 += cHalf + 1;
    }
    }
  }
  }

  return succ;
}
