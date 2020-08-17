#include <float.h>
#include <math.h>
#include <intrin.h>
#include <stdint.h>
#include "chol_internal_definitions.h"

#ifdef __AVX2__
#define MADD(   acc, x, y) acc = _mm256_fmadd_pd (x, y, acc)
#define MSUB(   acc, x, y) acc = _mm256_fnmadd_pd(x, y, acc)
#else
#define MADD(   acc, x, y) acc = _mm256_add_pd   (acc, _mm256_mul_pd(x, y))
#define MSUB(   acc, x, y) acc = _mm256_sub_pd   (acc, _mm256_mul_pd(x, y))
#endif
#define MADDSUB(acc, x, y) acc = _mm256_addsub_pd(acc, _mm256_mul_pd(x, y))


#ifdef __AVX2__
#define MADD128(   acc, x, y) acc = _mm_fmadd_pd (x, y, acc)
#define MSUB128(   acc, x, y) acc = _mm_fnmadd_pd(x, y, acc)
#else
#define MADD128(   acc, x, y) acc = _mm_add_pd   (acc, _mm_mul_pd(x, y))
#define MSUB128(   acc, x, y) acc = _mm_sub_pd   (acc, _mm_mul_pd(x, y))
#endif
#define MADDSUB128(acc, x, y) acc = _mm_addsub_pd(acc, _mm_mul_pd(x, y))

static
void InnerLoop(complex_m256d *rowR, complex_m256d* rowC, unsigned r, unsigned nRows, unsigned c)
{
  unsigned nRowQuartets = nRows / 4;
  unsigned cHalf = c / 2;
  unsigned rHalf = (r / 2) + 1;
  while (nRowQuartets)
  {
    // calculate eight results
    // multiply-add two dual-rows by dual-column
    complex_m256d *rowR0 = rowR;
    complex_m256d *rowR2 = rowR + rHalf;

    __m256d inp_re0x   = rowR0[cHalf].re;                                     // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    __m256d sum_re0011 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0x6); // re(r+0,c+0) 0           0           re(r+1,c+1)
    __m256d sum_re0110 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0x9); // 0           re(r+0,c+1) re(r+1,c+0) 0
    __m256d inp_re2x   = rowR2[cHalf].re;                                     // re(r+2,c+0) re(r+2,c+1) re(r+3,c+0) re(r+3,c+1)
    __m256d sum_re2031 = _mm256_blend_pd(inp_re2x, _mm256_setzero_pd(), 0x6); // re(r+2,c+0) 0           0           re(r+3,c+1)
    __m256d sum_re2130 = _mm256_blend_pd(inp_re2x, _mm256_setzero_pd(), 0x9); // 0           re(r+2,c+1) re(r+3,c+0) 0
    __m256d inp_im0x   = rowR0[cHalf].im;                                     // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    __m256d sum_im0011 = _mm256_blend_pd(inp_im0x, _mm256_setzero_pd(), 0x6); // im(r+0,c+0) 0           0           im(r+1,c+1)
    __m256d sum_im0110 = _mm256_blend_pd(inp_im0x, _mm256_setzero_pd(), 0x9); // 0           im(r+0,c+1) im(r+1,c+0) 0
    __m256d inp_im2x   = rowR2[cHalf].im;                                     // im(r+2,c+0) im(r+2,c+1) im(r+3,c+0) im(r+3,c+1)
    __m256d sum_im2031 = _mm256_blend_pd(inp_im2x, _mm256_setzero_pd(), 0x6); // im(r+2,c+0) 0           0           im(r+3,c+1)
    __m256d sum_im2130 = _mm256_blend_pd(inp_im2x, _mm256_setzero_pd(), 0x9); // 0           im(r+2,c+1) im(r+3,c+0) 0

    unsigned kk = cHalf;
    __m256d reC01 = rowC->re;  // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
    __m256d imC01 = rowC->im;  // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
    __m256d rval0 = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
    __m256d rval2 = rowR2->re; // re(r+2,k*2+0) re(r+2,k*2+1) re(r+3,k*2+0) re(r+3,k*2+1)
    do { // multiply-add two dual-rows by dual-column
      //__m256d reC01 = rowC->re;                                // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      //__m256d imC01 = rowC->im;                                // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      __m256d reC10 = _mm256_permute2f128_pd(reC01, reC01, 1); // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
      __m256d imC10 = _mm256_permute2f128_pd(imC01, imC01, 1); // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
      //__m256d rval0, rval2;

      _mm_prefetch((const char*)(rowC+3), _MM_HINT_T0);

      //rval0 = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
      MSUB(sum_re0011, rval0, reC01);
      // 0: re(r+0,c+0) -= re(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= re(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= re(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: re(r+1,c+1) -= re(r+1,k*2+1)*re(c+1,k*2+1)
      MADD(sum_im0011, rval0, imC01);
      // 0: im(r+0,c+0) += re(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           += re(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           += re(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: im(r+1,c+1) += re(r+1,k*2+1)*im(c+1,k*2+1)
      MSUB(sum_re0110, rval0, reC10);
      // 0: 0           -= re(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: re(r+0,c+1) -= re(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: re(r+1,c+0) -= re(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= re(r+1,k*2+1)*re(c+0,k*2+1)
      MADD(sum_im0110, rval0, imC10);
      // 0: 0           += re(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: im(r+0,c+1) += re(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: im(r+1,c+0) += re(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           += re(r+1,k*2+1)*im(c+0,k*2+1)

      //rval2 = rowR2->re; // re(r+2,k*2+0) re(r+2,k*2+1) re(r+3,k*2+0) re(r+3,k*2+1)
      MSUB(sum_re2031, rval2, reC01);
      // 0: re(r+2,c+0) -= re(r+2,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= re(r+2,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= re(r+3,k*2+0)*re(c+1,k*2+0)
      // 3: re(r+3,c+1) -= re(r+3,k*2+1)*re(c+1,k*2+1)
      MADD(sum_im2031, rval2, imC01);
      // 0: im(r+2,c+0) += re(r+2,k*2+0)*im(c+0,k*2+0)
      // 1: 0           += re(r+2,k*2+1)*im(c+0,k*2+1)
      // 2: 0           += re(r+3,k*2+0)*im(c+1,k*2+0)
      // 3: im(r+3,c+1) += re(r+3,k*2+1)*im(c+1,k*2+1)
      MSUB(sum_re2130, rval2, reC10);
      // 0: 0           -= re(r+2,k*2+0)*re(c+0,k*2+0)
      // 1: re(r+3,c+1) -= re(r+2,k*2+1)*re(c+0,k*2+1)
      // 2: re(r+2,c+0) -= re(r+3,k*2+0)*re(c+1,k*2+0)
      // 3: 0           -= re(r+3,k*2+1)*re(c+1,k*2+1)
      MADD(sum_im2130, rval2, imC10);
      // 0: 0           += re(r+2,k*2+0)*im(c+0,k*2+0)
      // 1: im(r+2,c+1) += re(r+2,k*2+1)*im(c+0,k*2+1)
      // 2: im(r+3,c+0) += re(r+3,k*2+0)*im(c+1,k*2+0)
      // 3: 0           += re(r+3,k*2+1)*im(c+1,k*2+1)

      rval0 = rowR0->im; // im(r+0,k*2+0) im(r+0,k*2+1) im(r+1,k*2+0) im(r+1,k*2+1)
      MSUB(sum_re0011, rval0, imC01);
      // 0: re(r+0,c+0) -= im(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: re(r+1,c+1) -= im(r+1,k*2+1)*im(c+1,k*2+1)
      MSUB(sum_im0011, rval0, reC01);
      // 0: im(r+0,c+0) -= im(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: im(r+1,c+1) -= im(r+1,k*2+1)*re(c+1,k*2+1)
      MSUB(sum_re0110, rval0, imC10);
      // 0: 0           -= im(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: re(r+0,c+1) -= im(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: re(r+1,c+0) -= im(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*im(c+0,k*2+1)
      MSUB(sum_im0110, rval0, reC10);
      // 0: 0           -= im(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: im(r+0,c+1) -= im(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: im(r+1,c+0) -= im(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*re(c+0,k*2+1)


      rval2 = rowR2->im; // im(r+2,k*2+0) im(r+2,k*2+1) im(r+3,k*2+0) im(r+3,k*2+1)
      MSUB(sum_re2031, rval2, imC01);
      // 0: re(r+2,c+0) -= im(r+2,k*2+0)*im(c+0,k*2+0)
      // 1: 0           -= im(r+2,k*2+1)*im(c+0,k*2+1)
      // 2: 0           -= im(r+3,k*2+0)*im(c+1,k*2+0)
      // 3: re(r+3,c+1) -= im(r+3,k*2+1)*im(c+1,k*2+1)
      MSUB(sum_im2031, rval2, reC01);
      // 0: im(r+2,c+0) -= im(r+2,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= im(r+2,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= im(r+3,k*2+0)*re(c+1,k*2+0)
      // 3: im(r+3,c+1) -= im(r+3,k*2+1)*re(c+1,k*2+1)
      MSUB(sum_re2130, rval2, imC10);
      // 0: 0           -= im(r+2,k*2+0)*im(c+1,k*2+0)
      // 1: re(r+3,c+1) -= im(r+2,k*2+1)*im(c+1,k*2+1)
      // 2: re(r+2,c+0) -= im(r+3,k*2+0)*im(c+0,k*2+0)
      // 3: 0           -= im(r+3,k*2+1)*im(c+0,k*2+1)
      MSUB(sum_im2130, rval2, reC10);
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
    __m256d sum_reim0011 = _mm256_hadd_pd(sum_re0011, sum_im0011); // re(r+0,c+0) im(r+0,c+0) re(r+1,c+1) im(r+1,c+1)
    __m256d sum_reim2031 = _mm256_hadd_pd(sum_re2031, sum_im2031); // re(r+2,c+0) im(r+2,c+0) re(r+3,c+1) im(r+3,c+1)
    __m256d sum_reim0110 = _mm256_hadd_pd(sum_re0110, sum_im0110); // re(r+0,c+1) im(r+0,c+1) re(r+1,c+0) im(r+1,c+0)
    __m256d sum_reim2130 = _mm256_hadd_pd(sum_re2130, sum_im2130); // re(r+2,c+1) im(r+2,c+1) re(r+3,c+0) im(r+3,c+0)
    // separate column c+0 from column c+1
    __m256d sum_reim00 = _mm256_blend_pd(sum_reim0011, sum_reim0110, 0xC); // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    __m256d sum_reim01 = _mm256_blend_pd(sum_reim0011, sum_reim0110, 0x3); // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)
    __m256d sum_reim20 = _mm256_blend_pd(sum_reim2031, sum_reim2130, 0xC); // re(r+2,c+0) im(r+2,c+0) re(r+3,c+0) im(r+3,c+0)
    __m256d sum_reim21 = _mm256_blend_pd(sum_reim2031, sum_reim2130, 0x3); // re(r+2,c+1) im(r+2,c+1) re(r+3,c+1) im(r+3,c+1)

    // calculate re/im(r+0:r+3, c+0)
    __m256d aaInvSqrt0 = _mm256_broadcast_sd((const double*)&rowC->im+0);  // 1/re(c+0,c+0)
    __m256d reim00 = _mm256_mul_pd(sum_reim00, aaInvSqrt0); // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    __m256d reim20 = _mm256_mul_pd(sum_reim20, aaInvSqrt0); // re(r+2,c+0) im(r+2,c+0) re(r+3,c+0) im(r+3,c+0)
    // add elements (r+0:r+3,c+0) to dot product of elements (r+0:r+3,c+1)
    __m256d reC10Last = _mm256_broadcast_sd((const double*)&rowC->re+2); // re(c+1,c+0)
    MSUB(sum_reim01, reim00, reC10Last);
    // 0: re(r+0,c+1) -= re(r+0,c+0)*re(c+1,c+0)
    // 1: im(r+0,c+1) -= im(r+0,c+0)*re(c+1,c+0)
    // 2: re(r+1,c+1) -= re(r+1,c+0)*re(c+1,c+0)
    // 3: im(r+1,c+1) -= im(r+1,c+0)*re(c+1,c+0)
    MSUB(sum_reim21, reim20, reC10Last);
    // 0: re(r+2,c+1) -= re(r+2,c+0)*re(c+1,c+0)
    // 1: im(r+2,c+1) -= im(r+2,c+0)*re(c+1,c+0)
    // 2: re(r+3,c+1) -= re(r+3,c+0)*re(c+1,c+0)
    // 3: im(r+3,c+1) -= im(r+3,c+0)*re(c+1,c+0)
    {
    __m256d imC10Last = _mm256_broadcast_sd((const double*)&rowC->im+2); // im(c+1,c+0)
    MADDSUB(sum_reim01, _mm256_shuffle_pd(reim00, reim00, 5), imC10Last);
    // 0: re(r+0,c+1) -= im(r+0,c+0)*im(c+1,c+0)
    // 1: im(r+0,c+1) += re(r+0,c+0)*im(c+1,c+0)
    // 2: re(r+1,c+1) -= im(r+1,c+0)*im(c+1,c+0)
    // 3: im(r+1,c+1) += re(r+1,c+0)*im(c+1,c+0)
    MADDSUB(sum_reim21, _mm256_shuffle_pd(reim20, reim20, 5), imC10Last);
    // 0: re(r+2,c+1) -= im(r+2,c+0)*im(c+1,c+0)
    // 1: im(r+2,c+1) += re(r+2,c+0)*im(c+1,c+0)
    // 2: re(r+3,c+1) -= im(r+3,c+0)*im(c+1,c+0)
    // 3: im(r+3,c+1) += re(r+3,c+0)*im(c+1,c+0)
    {
    __m256d aaInvSqrt1 = _mm256_broadcast_sd((const double*)&rowC->im+3);  // 1/re(c+1,c+1)
    __m256d reim01 = _mm256_mul_pd(sum_reim01, aaInvSqrt1); // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)
    __m256d reim21 = _mm256_mul_pd(sum_reim21, aaInvSqrt1); // re(r+2,c+1) im(r+2,c+1) re(r+3,c+1) im(r+3,c+1)
    // store results
    rowR0->re = _mm256_unpacklo_pd(reim00, reim01); // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    rowR0->im = _mm256_unpackhi_pd(reim00, reim01); // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    rowR2->re = _mm256_unpacklo_pd(reim20, reim21); // re(r+2,c+0) re(r+2,c+1) re(r+3,c+0) re(r+3,c+1)
    rowR2->im = _mm256_unpackhi_pd(reim20, reim21); // im(r+2,c+0) im(r+2,c+1) im(r+3,c+0) im(r+3,c+1)
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
    complex_m256d *rowR0 = rowR;

    __m256d inp_re0x   = rowR0[cHalf].re;                                     // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    __m256d sum_re0011 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0x6); // re(r+0,c+0) 0           0           re(r+1,c+1)
    __m256d sum_re0110 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0x9); // 0           re(r+0,c+1) re(r+1,c+0) 0
    __m256d inp_im0x   = rowR0[cHalf].im;                                     // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
    __m256d sum_im0011 = _mm256_blend_pd(inp_im0x, _mm256_setzero_pd(), 0x6); // im(r+0,c+0) 0           0           im(r+1,c+1)
    __m256d sum_im0110 = _mm256_blend_pd(inp_im0x, _mm256_setzero_pd(), 0x9); // 0           im(r+0,c+1) im(r+1,c+0) 0

    unsigned k = cHalf;
    do { // multiply-add two dual-rows by dual-column
      __m256d reC01 = rowC->re;                                // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      __m256d imC01 = rowC->im;                                // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      __m256d reC10 = _mm256_permute2f128_pd(reC01, reC01, 1); // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
      __m256d imC10 = _mm256_permute2f128_pd(imC01, imC01, 1); // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
      __m256d rval;

      rval = rowR0->re; // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
      MSUB(sum_re0011, rval, reC01);
      // 0: re(r+0,c+0) -= re(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= re(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= re(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: re(r+1,c+1) -= re(r+1,k*2+1)*re(c+1,k*2+1)
      MADD(sum_im0011, rval, imC01);
      // 0: im(r+0,c+0) += re(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           += re(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           += re(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: im(r+1,c+1) += re(r+1,k*2+1)*im(c+1,k*2+1)
      MSUB(sum_re0110, rval, reC10);
      // 0: 0           -= re(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: re(r+0,c+1) -= re(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: re(r+1,c+0) -= re(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= re(r+1,k*2+1)*re(c+0,k*2+1)
      MADD(sum_im0110, rval, imC10);
      // 0: 0           += re(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: im(r+0,c+1) += re(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: im(r+1,c+0) += re(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           += re(r+1,k*2+1)*im(c+0,k*2+1)

      rval = rowR0->im; // im(r+0,k*2+0) im(r+0,k*2+1) im(r+1,k*2+0) im(r+1,k*2+1)
      MSUB(sum_re0011, rval, imC01);
      // 0: re(r+0,c+0) -= im(r+0,k*2+0)*im(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*im(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*im(c+1,k*2+0)
      // 3: re(r+1,c+1) -= im(r+1,k*2+1)*im(c+1,k*2+1)
      MSUB(sum_im0011, rval, reC01);
      // 0: im(r+0,c+0) -= im(r+0,k*2+0)*re(c+0,k*2+0)
      // 1: 0           -= im(r+0,k*2+1)*re(c+0,k*2+1)
      // 2: 0           -= im(r+1,k*2+0)*re(c+1,k*2+0)
      // 3: im(r+1,c+1) -= im(r+1,k*2+1)*re(c+1,k*2+1)
      MSUB(sum_re0110, rval, imC10);
      // 0: 0           -= im(r+0,k*2+0)*im(c+1,k*2+0)
      // 1: re(r+0,c+1) -= im(r+0,k*2+1)*im(c+1,k*2+1)
      // 2: re(r+1,c+0) -= im(r+1,k*2+0)*im(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*im(c+0,k*2+1)
      MSUB(sum_im0110, rval, reC10);
      // 0: 0           -= im(r+0,k*2+0)*re(c+1,k*2+0)
      // 1: im(r+0,c+1) -= im(r+0,k*2+1)*re(c+1,k*2+1)
      // 2: im(r+1,c+0) -= im(r+1,k*2+0)*re(c+0,k*2+0)
      // 3: 0           -= im(r+1,k*2+1)*re(c+0,k*2+1)

      rowC += 1;
      rowR0 += 1;
    } while (--k);

    {
    // sum up odd and even sub-sums
    __m256d sum_reim0011 = _mm256_hadd_pd(sum_re0011, sum_im0011); // re(r+0,c+0) im(r+0,c+0) re(r+1,c+1) im(r+1,c+1)
    __m256d sum_reim0110 = _mm256_hadd_pd(sum_re0110, sum_im0110); // re(r+0,c+1) im(r+0,c+1) re(r+1,c+0) im(r+1,c+0)
    // separate column c+0 from column c+1
    __m256d sum_reim00 = _mm256_blend_pd(sum_reim0011, sum_reim0110, 0xC); // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    __m256d sum_reim01 = _mm256_blend_pd(sum_reim0011, sum_reim0110, 0x3); // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)

    // calculate re/im(r+0:r+1, c+0)
    __m256d aaInvSqrt0 = _mm256_broadcast_sd((const double*)&rowC->im+0);  // 1/re(c+0,c+0)
    __m256d reim00 = _mm256_mul_pd(sum_reim00, aaInvSqrt0); // re(r+0,c+0) im(r+0,c+0) re(r+1,c+0) im(r+1,c+0)
    {
    // add elements (r+0:r+3,c+0) to dot product of elements (r+0:r+3,c+1)
    __m256d reC10Last = _mm256_broadcast_sd((const double*)&rowC->re+2); // re(c+1,c+0)
    MSUB(sum_reim01, reim00, reC10Last);
    // 0: re(r+0,c+1) -= re(r+0,c+0)*re(c+1,c+0)
    // 1: im(r+0,c+1) -= im(r+0,c+0)*re(c+1,c+0)
    // 2: re(r+1,c+1) -= re(r+1,c+0)*re(c+1,c+0)
    // 3: im(r+1,c+1) -= im(r+1,c+0)*re(c+1,c+0)
    {
    __m256d imC10Last = _mm256_broadcast_sd((const double*)&rowC->im+2); // im(c+1,c+0)
    MADDSUB(sum_reim01, _mm256_shuffle_pd(reim00, reim00, 5), imC10Last);
    // 0: re(r+0,c+1) -= im(r+0,c+0)*im(c+1,c+0)
    // 1: im(r+0,c+1) += re(r+0,c+0)*im(c+1,c+0)
    // 2: re(r+1,c+1) -= im(r+1,c+0)*im(c+1,c+0)
    // 3: im(r+1,c+1) += re(r+1,c+0)*im(c+1,c+0)
    {
    __m256d aaInvSqrt1 = _mm256_broadcast_sd((const double*)&rowC->im+3);  // 1/re(c+1,c+1)
    __m256d reim01 = _mm256_mul_pd(sum_reim01, aaInvSqrt1); // re(r+0,c+1) im(r+0,c+1) re(r+1,c+1) im(r+1,c+1)
    // store results
    rowR0->re = _mm256_unpacklo_pd(reim00, reim01); // re(r+0,c+0) re(r+0,c+1) re(r+1,c+0) re(r+1,c+1)
    rowR0->im = _mm256_unpackhi_pd(reim00, reim01); // im(r+0,c+0) im(r+0,c+1) im(r+1,c+0) im(r+1,c+1)
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
    __m256d *rowR0 = &rowR->re;

    __m256d inp_reim   = rowR0[cHalf];                                         // re(r,c+0) re(r,c+1) im(r,c+0) im(r,c+1)
    __m256d sum_re0im1a = _mm256_blend_pd(inp_reim, _mm256_setzero_pd(), 0x6); // re(r,c+0) 0         0         im(r,c+1)
    __m256d sum_im0re1a = _mm256_blend_pd(inp_reim, _mm256_setzero_pd(), 0x9); // 0         re(r,c+1) im(r,c+0) 0
    __m256d sum_re0im1b = _mm256_setzero_pd();
    __m256d sum_im0re1b = _mm256_setzero_pd();

    unsigned k = cHalf;
    sum_im0re1a = _mm256_permute2f128_pd(sum_im0re1a, sum_im0re1a, 1);         // im(r,c+0) 0         0         re(r,c+1)
    do { // multiply-add two dual-rows by dual-column
      __m256d reC01 = rowC->re;                                // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      __m256d imC01 = rowC->im;                                // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
      __m256d reimR = *rowR0;                                  // re(r+0,k*2+0) re(r+0,k*2+1) im(r+0,k*2+0) im(r+0,k*2+1)
      __m256d imreR = _mm256_permute2f128_pd(reimR, reimR, 1); // im(r+0,k*2+0) im(r+0,k*2+1) re(r+0,k*2+0) re(r+0,k*2+1)

      MSUB(sum_re0im1a, reimR, reC01);
      // 0: re(r,c+0) -= re(r,k*2+0)*re(c+0,k*2+0) : +re(r,c+0)
      // 1: 0         -= re(r,k*2+1)*re(c+0,k*2+1) : +re(r,c+0)
      // 2: 0         -= im(r,k*2+0)*re(c+1,k*2+0) : +im(r,c+1)
      // 3: im(r,c+1) -= im(r,k*2+1)*re(c+1,k*2+1) : +im(r,c+1)
      MSUB(sum_im0re1a, imreR, reC01);
      // 0: 0         -= im(r,k*2+0)*re(c+0,k*2+0) : +im(r,c+0)
      // 1: 0         -= im(r,k*2+1)*re(c+0,k*2+1) : +im(r,c+0)
      // 2: 0         -= re(r,k*2+0)*re(c+1,k*2+0) : +re(r,c+1)
      // 3: 0         -= re(r,k*2+1)*re(c+1,k*2+1) : +re(r,c+1)
      MSUB(sum_im0re1b, reimR, imC01);
      // 0: 0         -= re(r,k*2+0)*im(c+0,k*2+0) : -im(r,c+0)
      // 1: 0         -= re(r,k*2+1)*im(c+0,k*2+1) : -im(r,c+0)
      // 2: 0         -= im(r,k*2+0)*im(c+1,k*2+0) : +re(r,c+1)
      // 3: 0         -= im(r,k*2+1)*im(c+1,k*2+1) : +re(r,c+1)
      MSUB(sum_re0im1b, imreR, imC01);
      // 0: 0         += im(r,k*2+0)*im(c+0,k*2+0) : +re(r,c+0)
      // 1: 0         += im(r,k*2+1)*im(c+0,k*2+1) : +re(r,c+0)
      // 2: 0         += re(r,k*2+0)*im(c+1,k*2+0) : -im(r,c+1)
      // 3: 0         += re(r,k*2+1)*im(c+1,k*2+1) : -im(r,c+1)

      rowC  += 1;
      rowR0 += 1;
    } while (--k);

    {
    // sum up odd and even sub-sums
    __m256d sum_imre_a = _mm256_hadd_pd(sum_im0re1a, sum_re0im1a); // +im(r,c+0) +re(r,c+0) +re(r,c+1) +im(r,c+1)
    __m256d sum_imre_b = _mm256_hadd_pd(sum_im0re1b, sum_re0im1b); // -im(r,c+0) +re(r,c+0) +re(r,c+1) -im(r,c+1)

    // sum up partial sums
    sum_imre_a = _mm256_shuffle_pd(sum_imre_a, sum_imre_a, 6);     // +im(r,c+0) +re(r,c+0) +im(r,c+1) +re(r,c+1)
    sum_imre_b = _mm256_shuffle_pd(sum_imre_b, sum_imre_b, 6);     // -im(r,c+0) +re(r,c+0) -im(r,c+1) +re(r,c+1)
    {
    __m256d sum_imre = _mm256_addsub_pd(sum_imre_a, sum_imre_b);   // +im(r,c+0) +re(r,c+0) +im(r,c+1) +re(r,c+1)

    // calculate re/im(r, c+0)
    __m128d aaInvSqrt0 = _mm_loaddup_pd((const double*)&rowC->im+0);  // 1/re(c+0,c+0)
    __m128d imre0 = _mm_mul_pd(_mm256_castpd256_pd128(sum_imre), aaInvSqrt0); // im(r,c+0) re(r,c+0)

    // add elements (r,c+0) to dot product of elements (r,c+1)
    __m128d reC10Last = _mm_loaddup_pd((const double*)&rowC->re+2); // re(c+1,c+0)
    __m128d sum_imre1 = _mm256_extractf128_pd(sum_imre, 1);         // im(r,c+1) re(r,c+1)
    sum_imre1 = _mm_sub_pd   (sum_imre1, _mm_mul_pd(imre0, reC10Last));
    {
    __m128d imC10Last = _mm_loaddup_pd((const double*)&rowC->im+2); // im(c+1,c+0)
    __m128d sum_reim1 = _mm_shuffle_pd(sum_imre1, sum_imre1, 1);    // re(r,c+1) im(r,c+1)
    sum_reim1 = _mm_addsub_pd(sum_reim1, _mm_mul_pd(imre0, imC10Last));
    {
    __m128d aaInvSqrt1 = _mm_loaddup_pd((const double*)&rowC->im+3);// 1/re(c+1,c+1)
    __m128d reim1 = _mm_mul_pd(sum_reim1, aaInvSqrt1);              // re(r,c+1) im(r,c+1)
    // store results
    ((__m128d*)rowR0)[0] = _mm_shuffle_pd(imre0, reim1, 1); // re(r,c+0) re(r,c+1)
    ((__m128d*)rowR0)[1] = _mm_shuffle_pd(imre0, reim1, 2); // im(r,c+0) im(r,c+1)
    }
    }
    }
    }
  }
}

//
// row1st - always even
// return - 0 - matrix non positive defined, 1 - matrix positive defined
int chol_CrBa4_i2_Band_Factorize(complex_m256d* triang, unsigned row1st, unsigned rowLast)
{
  const double DIAG_MIN   =  (double)(FLT_MIN);
  // FLT_MIN is not really special in context of double-precision calculations
  // but I wanted very small number that is still much bigger than DBL_MIN
  // and FLT_MIN looks good enough
  const double DIAG_SUBST =  1E40;

  int succ = 1;
  unsigned rHalf1st = row1st/2;
  complex_m256d *rowDataset = &triang[rHalf1st*(rHalf1st+1)/2];

  // process Crout band
  // columns 0:1
  {
    double aaInvSqrt0 = ((double*)&triang[0].im)[0];
    double re10       = ((double*)&triang[0].re)[2];
    double im10       = ((double*)&triang[0].im)[2];
    double aaInvSqrt1 = ((double*)&triang[0].im)[3];

    complex_m256d *rowR = rowDataset;
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
      __m128d* lastRowR = (__m128d*)(rowR);
      double reR00 = ((double*)lastRowR)[0] * aaInvSqrt0;
      double imR00 = ((double*)lastRowR)[2] * aaInvSqrt0;
      double reR01 = ((double*)lastRowR)[1];
      double imR01 = ((double*)lastRowR)[3];
      reR01 -= reR00 * re10;
      imR01 -= imR00 * re10;
      reR01 -= imR00 * im10;
      imR01 += reR00 * im10;
      ((double*)lastRowR)[0] = reR00;
      ((double*)lastRowR)[1] = reR01 * aaInvSqrt1;
      ((double*)lastRowR)[2] = imR00;
      ((double*)lastRowR)[3] = imR01 * aaInvSqrt1;
    }
  }

  {
  complex_m256d *rowCx2 = &triang[1];
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
      complex_m256d *rowR0 = rowCx2;
      complex_m256d *rowR2 = rowR0 + cHalf + 1;
      __m256d inp_re0x   = rowR0[cHalf].re;                                     // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      __m256d sum_re0011 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0x6); // re(c+0,c+0) 0           0           re(c+1,c+1)
      __m256d sum_re0110 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0xB); // 0           0           re(c+1,c+0) 0
      __m256d inp_re2x   = rowR2[cHalf].re;                                     // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      __m256d sum_re2031 = _mm256_blend_pd(inp_re2x, _mm256_setzero_pd(), 0x6); // re(c+2,c+0) 0           0           re(c+3,c+1)
      __m256d sum_re2130 = _mm256_blend_pd(inp_re2x, _mm256_setzero_pd(), 0x9); // 0           re(c+2,c+1) re(c+3,c+0) 0
      __m256d inp_im0x   = rowR0[cHalf].im;                                     // x           x           im(c+1,c+0) x
      __m256d sum_im0110 = _mm256_blend_pd(inp_im0x, _mm256_setzero_pd(), 0xB); // 0           0           im(c+1,c+0) 0
      __m256d inp_im2x   = rowR2[cHalf].im;                                     // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      __m256d sum_im2031 = _mm256_blend_pd(inp_im2x, _mm256_setzero_pd(), 0x6); // im(c+2,c+0) 0           0           im(c+3,c+1)
      __m256d sum_im2130 = _mm256_blend_pd(inp_im2x, _mm256_setzero_pd(), 0x9); // 0           im(c+2,c+1) im(c+3,c+0) 0

      unsigned k = cHalf;
      do { // multiply-add two dual-rows by the first dual-row of the two
        __m256d reC01 = rowR0->re;                               // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
        __m256d imC01 = rowR0->im;                               // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
        __m256d reC10 = _mm256_permute2f128_pd(reC01, reC01, 1); // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
        __m256d imC10 = _mm256_permute2f128_pd(imC01, imC01, 1); // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
        __m256d rval;

        MSUB(sum_re0011, reC01, reC01);
        // 0: re(c+0,c+0) -= re(c+0,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= re(c+1,k*2+0)*re(c+1,k*2+0)
        // 3: re(c+1,c+1) -= re(c+1,k*2+1)*re(c+1,k*2+1)
        MSUB(sum_re0110, reC01, reC10);
        // 0: 0           -= re(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: re(c+1,c+0) -= re(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= re(c+1,k*2+1)*re(c+0,k*2+1)
        MADD(sum_im0110, reC01, imC10);
        // 0: 0           += re(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           += re(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: im(c+1,c+0) += re(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           += re(c+1,k*2+1)*im(c+0,k*2+1)

        rval = rowR2->re; // re(c+2,k*2+0) re(c+2,k*2+1) re(c+3,k*2+0) re(c+3,k*2+1)
        MSUB(sum_re2031, rval, reC01);
        // 0: re(c+2,c+0) -= re(c+2,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= re(c+2,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= re(c+3,k*2+0)*re(c+1,k*2+0)
        // 3: re(c+3,c+1) -= re(c+3,k*2+1)*re(c+1,k*2+1)
        MADD(sum_im2031, rval, imC01);
        // 0: im(c+2,c+0) += re(c+2,k*2+0)*im(c+0,k*2+0)
        // 1: 0           += re(c+2,k*2+1)*im(c+0,k*2+1)
        // 2: 0           += re(c+3,k*2+0)*im(c+1,k*2+0)
        // 3: im(c+3,c+1) += re(c+3,k*2+1)*im(c+1,k*2+1)
        MSUB(sum_re2130, rval, reC10);
        // 0: 0           -= re(c+2,k*2+0)*re(c+1,k*2+0)
        // 1: re(c+3,c+1) -= re(c+2,k*2+1)*re(c+1,k*2+1)
        // 2: re(c+2,c+0) -= re(c+3,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= re(c+3,k*2+1)*re(c+0,k*2+1)
        MADD(sum_im2130, rval, imC10);
        // 0: 0           += re(c+2,k*2+0)*im(c+1,k*2+0)
        // 1: im(c+2,c+1) += re(c+2,k*2+1)*im(c+1,k*2+1)
        // 2: im(c+3,c+0) += re(c+3,k*2+0)*im(c+0,k*2+0)
        // 3: 0           += re(c+3,k*2+1)*im(c+0,k*2+1)

        MSUB(sum_re0011, imC01, imC01);
        // 0: re(c+0,c+0) -= im(c+0,k*2+0)*im(c+0,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+0,k*2+1)
        // 2: 0           -= im(c+1,k*2+0)*im(c+1,k*2+0)
        // 3: re(c+1,c+1) -= im(c+1,k*2+1)*im(c+1,k*2+1)
        MSUB(sum_re0110, imC01, imC10);
        // 0: 0           -= im(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: re(c+1,c+0) -= im(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*im(c+0,k*2+1)
        MSUB(sum_im0110, imC01, reC10);
        // 0: 0           -= im(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: im(c+1,c+0) -= im(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*re(c+0,k*2+1)

        rval = rowR2->im; // im(c+2,k*2+0) im(c+2,k*2+1) im(c+3,k*2+0) im(c+3,k*2+1)
        MSUB(sum_re2031, rval, imC01);
        // 0: re(c+2,c+0) -= im(c+2,k*2+0)*im(c+0,k*2+0)
        // 1: 0           -= im(c+2,k*2+1)*im(c+0,k*2+1)
        // 2: 0           -= im(c+3,k*2+0)*im(c+1,k*2+0)
        // 3: re(c+3,c+1) -= im(c+3,k*2+1)*im(c+1,k*2+1)
        MSUB(sum_im2031, rval, reC01);
        // 0: im(c+2,c+0) -= im(c+2,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= im(c+2,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= im(c+3,k*2+0)*re(c+1,k*2+0)
        // 3: im(c+3,c+1) -= im(c+3,k*2+1)*re(c+1,k*2+1)
        MSUB(sum_re2130, rval, imC10);
        // 0: 0           -= im(c+2,k*2+0)*im(c+1,k*2+0)
        // 1: re(c+3,c+1) -= im(c+2,k*2+1)*im(c+1,k*2+1)
        // 2: re(c+2,c+0) -= im(c+3,k*2+0)*im(c+0,k*2+0)
        // 3: 0           -= im(c+3,k*2+1)*im(c+0,k*2+1)
        MSUB(sum_im2130, rval, reC10);
        // 0: 0           -= im(c+2,k*2+0)*re(c+1,k*2+0)
        // 1: im(c+2,c+1) -= im(c+2,k*2+1)*re(c+1,k*2+1)
        // 2: im(c+3,c+0) -= im(c+3,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= im(c+3,k*2+1)*re(c+0,k*2+1)

        rowR0 += 1;
        rowR2 += 1;
      } while (--k);
      {
      // sum up odd and even sub-sums
      __m256d sum_re0 = _mm256_hadd_pd(sum_re0011, sum_re0110); // re(c+0,c+0) x           re(c+1,c+1) re(c+1,c+0)
      __m256d sum_im0 = _mm256_hadd_pd(sum_im0110, sum_im0110); // x           x           im(c+1,c+0) im(c+1,c+0)
      __m256d sum_re2 = _mm256_hadd_pd(sum_re2031, sum_re2130); // re(c+2,c+0) re(c+2,c+1) re(c+3,c+1) re(c+3,c+0)
      __m256d sum_im2 = _mm256_hadd_pd(sum_im2031, sum_im2130); // im(c+2,c+0) im(c+2,c+1) im(c+3,c+1) im(c+3,c+0)
      // store temporary results
      rowR0->re = _mm256_shuffle_pd(sum_re0, sum_re0, 4); // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      rowR0->im = sum_im0;                                // x           x           im(c+1,c+0) im(c+1,c+0)
      rowR2->re = _mm256_shuffle_pd(sum_re2, sum_re2, 6); // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      rowR2->im = _mm256_shuffle_pd(sum_im2, sum_im2, 6); // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      }
    }
    else
    { // process one dual-row
      complex_m256d *rowR0 = rowCx2;
      __m256d inp_re0x   = rowR0[cHalf].re;                                     // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      __m256d sum_re0011 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0x6); // re(c+0,c+0) 0           0           re(c+1,c+1)
      __m256d sum_re0110 = _mm256_blend_pd(inp_re0x, _mm256_setzero_pd(), 0xB); // 0           0           re(c+1,c+0) 0
      __m256d inp_im0x   = rowR0[cHalf].im;                                     // x           x           im(c+1,c+0) x
      __m256d sum_im0110 = _mm256_blend_pd(inp_im0x, _mm256_setzero_pd(), 0xB); // 0           0           im(c+1,c+0) 0

      unsigned k = cHalf;
      do { // multiply-add two dual-rows by the first dual-row of the two
        __m256d reC01 = rowR0->re;                               // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
        __m256d imC01 = rowR0->im;                               // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)
        __m256d reC10 = _mm256_permute2f128_pd(reC01, reC01, 1); // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
        __m256d imC10 = _mm256_permute2f128_pd(imC01, imC01, 1); // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)

        MSUB(sum_re0011, reC01, reC01);
        // 0: re(c+0,c+0) -= re(c+0,k*2+0)*re(c+0,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+0,k*2+1)
        // 2: 0           -= re(c+1,k*2+0)*re(c+1,k*2+0)
        // 3: re(c+1,c+1) -= re(c+1,k*2+1)*re(c+1,k*2+1)
        MSUB(sum_re0110, reC01, reC10);
        // 0: 0           -= re(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= re(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: re(c+1,c+0) -= re(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= re(c+1,k*2+1)*re(c+0,k*2+1)
        MADD(sum_im0110, reC01, imC10);
        // 0: 0           += re(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           += re(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: im(c+1,c+0) += re(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           += re(c+1,k*2+1)*im(c+0,k*2+1)

        MSUB(sum_re0011, imC01, imC01);
        // 0: re(c+0,c+0) -= im(c+0,k*2+0)*im(c+0,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+0,k*2+1)
        // 2: 0           -= im(c+1,k*2+0)*im(c+1,k*2+0)
        // 3: re(c+1,c+1) -= im(c+1,k*2+1)*im(c+1,k*2+1)
        MSUB(sum_re0110, imC01, imC10);
        // 0: 0           -= im(c+0,k*2+0)*im(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*im(c+1,k*2+1)
        // 2: re(c+1,c+0) -= im(c+1,k*2+0)*im(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*im(c+0,k*2+1)
        MSUB(sum_im0110, imC01, reC10);
        // 0: 0           -= im(c+0,k*2+0)*re(c+1,k*2+0)
        // 1: 0           -= im(c+0,k*2+1)*re(c+1,k*2+1)
        // 2: im(c+1,c+0) -= im(c+1,k*2+0)*re(c+0,k*2+0)
        // 3: 0           -= im(c+1,k*2+1)*re(c+0,k*2+1)

        rowR0 += 1;
      } while (--k);
      {
      // sum up odd and even sub-sums
      __m256d sum_re0 = _mm256_hadd_pd(sum_re0011, sum_re0110); // re(c+0,c+0) x           re(c+1,c+1) re(c+1,c+0)
      __m256d sum_im0 = _mm256_hadd_pd(sum_im0110, sum_im0110); // x           x           im(c+1,c+0) im(c+1,c+0)
      // store temporary results
      rowR0->re = _mm256_shuffle_pd(sum_re0, sum_re0, 4); // re(c+0,c+0) x           re(c+1,c+0) re(c+1,c+1)
      rowR0->im = sum_im0;                                // x           x           im(c+1,c+0) im(c+1,c+0)
      }
    }
    {
    // calculate diagonal element (c,c)
    double aa = ((double*)&rowCx2[cHalf].re)[0]; // current diagonal element (c,c)

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
    ((double*)&rowCx2[cHalf].re)[0] = aaSqrt0;
    ((double*)&rowCx2[cHalf].im)[0] = aaInvSqrt0; // store inverse in image field of diagonal
    {
    double reR10 = ((double*)&rowCx2[cHalf].re)[2] * aaInvSqrt0; // re(c+1,c)
    double imR10 = ((double*)&rowCx2[cHalf].im)[2] * aaInvSqrt0; // im(c+1,c)
    ((double*)&rowCx2[cHalf].re)[2] = reR10;
    ((double*)&rowCx2[cHalf].im)[2] = imR10;

    // calculate diagonal element (c+1,c+1)
    aa = ((double*)&rowCx2[cHalf].re)[3]; // re(c+1,c+1)
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
    ((double*)&rowCx2[cHalf].re)[3] = aaSqrt1;
    ((double*)&rowCx2[cHalf].im)[3] = aaInvSqrt1; // store inverse in image field of diagonal
    }
    }
    }

    {
    unsigned r = c + 2;
    complex_m256d *rowR0 = rowCx2 + cHalf + 1;
    if ((nRows & 2)==0)
    {
      // complete calculation of re/im(c+2:c+3,c+0:c+1)
      __m256d aaInvSqrt0_pd = _mm256_broadcast_sd(((double*)&rowCx2[cHalf].im) + 0);  // 1/re(c+0,c+0)
      __m256d reR2 = rowR0[cHalf].re;                                                 // re(c+2,c+0) re(c+2,c+1) re(c+3,c+0) re(c+3,c+1)
      __m256d imR2 = rowR0[cHalf].im;                                                 // im(c+2,c+0) im(c+2,c+1) im(c+3,c+0) im(c+3,c+1)
      __m256d reimR20 = _mm256_mul_pd(_mm256_unpacklo_pd(reR2, imR2), aaInvSqrt0_pd); // re(c+2,c+0) im(c+2,c+0) re(c+3,c+0) im(c+3,c+0)
      __m256d reimR21 = _mm256_unpackhi_pd(reR2, imR2);                               // re(c+2,c+1) im(c+2,c+1) re(c+3,c+1) im(c+3,c+1)
      __m256d reC = _mm256_broadcast_sd(((double*)&rowCx2[cHalf].re) + 2);            // re(c+1,c+0)
      MSUB(reimR21, reimR20, reC);
      {
      __m256d imC = _mm256_broadcast_sd(((double*)&rowCx2[cHalf].im) + 2);            // im(c+1,c+0)
      MADDSUB(reimR21, _mm256_shuffle_pd(reimR20, reimR20, 5), imC);
      }
      {
      __m256d aaInvSqrt1_pd = _mm256_broadcast_sd(((double*)&rowCx2[cHalf].im) + 3);  // 1/re(c+1,c+1)
      reimR21 = _mm256_mul_pd(reimR21, aaInvSqrt1_pd);
      }
      rowR0[cHalf].re = _mm256_unpacklo_pd(reimR20, reimR21);
      rowR0[cHalf].im = _mm256_unpackhi_pd(reimR20, reimR21);

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


void chol_CrBa4_i2_Band_ForwardSubstitute(
  __m256d*              result,
  const complex_m256d*  triang,
  unsigned row1st, unsigned rowLast)
{
  // solve L*y=vecB by forward substitution, result in result[]
  unsigned rHalf1st = row1st/2;
  const complex_m256d *row = &triang[rHalf1st*(rHalf1st+1)/2];
  unsigned rHalf;
  // process dual-rows
  unsigned rHalfLast = rowLast/2;
  for (rHalf = rHalf1st; rHalf != rHalfLast; rHalf += 1)
  {
    __m256d sumRexReIm = _mm256_setzero_pd();
    __m256d sumImxReIm = _mm256_setzero_pd();
    __m256d sumRexImRe = _mm256_setzero_pd();
    __m256d sumImxImRe = _mm256_setzero_pd();
    unsigned k;
    for (k = 0; k != rHalf; ++k)
    {
      __m256d reimRes = result[k];                                   // reRes(k*2+0)  reRes(k*2+1)  imRes(k*2+0)  imRes(k*2+1)
      __m256d imreRes = _mm256_permute2f128_pd(reimRes, reimRes, 1); // imRes(k*2+0)  imRes(k*2+1)  reRes(k*2+0)  reRes(k*2+1)
      __m256d reR = row->re;                                         // re(r+0,k*2+0) re(r+0,k*2+1) re(r+1,k*2+0) re(r+1,k*2+1)
      __m256d imR = row->im;                                         // im(r+0,k*2+0) im(r+0,k*2+1) im(r+1,k*2+0) im(r+1,k*2+1)

      MSUB(sumRexReIm, reR, reimRes);
      // 0: 0 -= re(r+0,k*2+0)*reRes(k*2+0) : reRes(r+0)
      // 1: 0 -= re(r+0,k*2+1)*reRes(k*2+1) : reRes(r+0)
      // 2: 0 -= re(r+1,k*2+0)*imRes(k*2+0) : imRes(r+1)
      // 3: 0 -= re(r+1,k*2+1)*imRes(k*2+1) : imRes(r+1)
      MSUB(sumImxReIm, imR, reimRes);
      // 0: 0 -= im(r+0,k*2+0)*reRes(k*2+0) : imRes(r+0)
      // 1: 0 -= im(r+0,k*2+1)*reRes(k*2+1) : imRes(r+0)
      // 2: 0 -= im(r+1,k*2+0)*imRes(k*2+0) :-reRes(r+1)
      // 3: 0 -= im(r+1,k*2+1)*imRes(k*2+1) :-reRes(r+1)
      MSUB(sumRexImRe, reR, imreRes);
      // 0: 0 -= re(r+0,k*2+0)*imRes(k*2+0) : imRes(r+0)
      // 1: 0 -= re(r+0,k*2+1)*imRes(k*2+1) : imRes(r+0)
      // 2: 0 -= re(r+1,k*2+0)*reRes(k*2+0) : reRes(r+1)
      // 3: 0 -= re(r+1,k*2+1)*reRes(k*2+1) : reRes(r+1)
      MSUB(sumImxImRe, imR, imreRes);
      // 0: 0 -= im(r+0,k*2+0)*imRes(k*2+0) :-reRes(r+0)
      // 1: 0 -= im(r+0,k*2+1)*imRes(k*2+1) :-reRes(r+0)
      // 2: 0 -= im(r+1,k*2+0)*reRes(k*2+0) : imRes(r+1)
      // 3: 0 -= im(r+1,k*2+1)*reRes(k*2+1) : imRes(r+1)
      row += 1;
    }
    {
    // sum up partial sums
    __m256d sum_p = _mm256_hadd_pd(sumRexReIm, sumRexImRe);      //  reRes(r+0) imRes(r+0)  imRes(r+1)  reRes(r+1)
    __m256d sum_n = _mm256_hadd_pd(sumImxImRe, sumImxReIm);      // -reRes(r+0) imRes(r+0)  imRes(r+1) -reRes(r+1)

    __m256d sum_re = _mm256_sub_pd(sum_p, sum_n);                //  reRes(r+0) x           x           reRes(r+1)
    __m256d sum_im = _mm256_add_pd(sum_p, sum_n);                //  x          imRes(r+0)  imRes(r+1)  x

    __m256d sum0 = _mm256_permute2f128_pd(sum_re, sum_im, 0x20); //  reRes(r+0) x           x           imRes(r+0)
    __m256d sum1 = _mm256_permute2f128_pd(sum_re, sum_im, 0x31); //  x          reRes(r+1)  imRes(r+1)  x
    __m256d reim = _mm256_shuffle_pd(sum0, sum1, 6);             //  reRes(r+0) reRes(r+1)  imRes(r+0)  imRes(r+1)

    reim = _mm256_add_pd(reim, result[rHalf]);                   //  reRes(r+0) reRes(r+1)  imRes(r+0)  imRes(r+1)
    {
    __m256d aaInvSqrt = _mm256_broadcast_sd((const double*)&row->im+0);  // 1/re(r+0,r+0)
    __m256d reimRes0  = _mm256_mul_pd(reim, aaInvSqrt);          //  reRes(r+0) x           imRes(r+0)  x
    __m256d imR10     = _mm256_broadcast_sd((const double*)&row->im+2);  // im(r+1,r+0)
    __m256d reR10     = _mm256_broadcast_sd((const double*)&row->re+2);  // re(r+1,r+0)

    sum_n = _mm256_mul_pd(imR10, reimRes0);                      //  imRes(r+1) x          -reRes(r+1)  x
    // 0: = im(r+1,r+0)*reRes(r+0) : imRes(r+1)
    // 1: = im(r+1,r+0)*x          : x
    // 2: = im(r+1,r+0)*imRes(r+0) :-reRes(r+1)
    // 3: = im(r+1,r+0)*x          : x
    sum_p = _mm256_mul_pd(reR10, reimRes0);                      //  reRes(r+1) x           imRes(r+1)  x
    // 0: = re(r+1,r+0)*reRes(r+0) : reRes(r+1)
    // 1: = re(r+1,r+0)*x          : x
    // 2: = re(r+1,r+0)*imRes(r+0) : imRes(r+1)
    // 3: = re(r+1,r+0)*x          : x
    sum_n  = _mm256_permute2f128_pd(sum_n, sum_n, 1);            // -reRes(r+1) x           imRes(r+1)  x
    sum_re = _mm256_sub_pd(sum_p, sum_n);                        //  reRes(r+1) x           x           x
    sum_im = _mm256_add_pd(sum_p, sum_n);                        //  x          x           imRes(r+1)  x
    {
    __m256d reimRes1 = _mm256_blend_pd(sum_re, sum_im, 0xC);     //  reRes(r+1) x           imRes(r+1)  x
    reimRes1 = _mm256_unpacklo_pd(reimRes1, reimRes1);           //  reRes(r+1) reRes(r+1)  imRes(r+1)  imRes(r+1)
    reimRes1 = _mm256_sub_pd(reim, reimRes1);                    //  x          reRes(r+1)  x           imRes(r+1)
    aaInvSqrt = _mm256_broadcast_sd((const double*)&row->im+3);  // 1/re(r+1,r+1)
    reimRes1  = _mm256_mul_pd(reimRes1, aaInvSqrt);              //  x          reRes(r+1)  x           imRes(r+1)
    result[rHalf] = _mm256_blend_pd(reimRes0, reimRes1, 0xA);    //  reRes(r+0) reRes(r+1)  imRes(r+0)  imRes(r+1)
    }
    }
    }
    row += 1;
  }
}

void chol_CrBa4_i2_BackSubstitute(
  __m256d*              result,
  const complex_m256d*  triang,
  unsigned n)
{
  // solve L'*x=res by column-wise back substitution:
  // algorithm structure that is similar to gaussian elimination
  // This algorithm do more memory accesses than the "natural" substitution,
  // but it avoids column-wise access to the lower triangular matrix.
  // Result in res[]
  unsigned cHalf1st = n/2;
  const complex_m256d *col = &triang[cHalf1st*(cHalf1st+1)/2];
  unsigned cHalf;
  for (cHalf = cHalf1st; cHalf != 0; )
  {
    col -= cHalf;
    --cHalf;
    {
    __m256d resC = result[cHalf];                                         // reRes(c+0) reRes(c+1) imRes(c+0) imRes(c+1)
    __m256d invD1 = _mm256_broadcast_sd((const double*)&col[cHalf].im+3); // 1/re(c+1,c+1)
    __m256d resC1 = _mm256_mul_pd(resC, invD1);                           // x          reRes(c+1) x          imRes(c+1)
    __m256d re10  = _mm256_broadcast_sd((const double*)&col[cHalf].re+2); // re(c+1,c+0)
    __m256d im10  = _mm256_broadcast_sd((const double*)&col[cHalf].im+2); // im(c+1,c+0)
    __m256d negMask = _mm256_castsi256_pd(
               _mm256_setr_epi64x(INT64_MIN, INT64_MIN, 0, 0));           // -           -           +           +
    resC1 = _mm256_unpackhi_pd(resC1, resC1);                             // reRes(c+1)  reRes(c+1)  imRes(c+1)  imRes(c+1)
    MSUB(resC, resC1, re10);
    // 0: reRes(c+0) -= re(c+1,c+0)*reRes(c+1) : reRes(c+0)
    // 1: x          -= re(c+1,c+0)*reRes(c+1) : x
    // 2: imRes(c+0) -= re(c+1,c+0)*imRes(c+1) : imRes(c+0)
    // 3: x          -= re(c+1,c+0)*imRes(c+1) : x
    {
    __m256d resC1i = _mm256_xor_pd(resC1, negMask);                       // -reRes(c+1) -reRes(c+1)  imRes(c+1)  imRes(c+1)
    resC1i =  _mm256_permute2f128_pd(resC1i, resC1i, 1);                  //  imRes(c+1)  imRes(c+1) -reRes(c+1) -reRes(c+1)
    MSUB(resC, resC1i, im10);
    // 0: reRes(c+0) -= im(c+1,c+0)*imRes(c+1) : reRes(c+0)
    // 1: x          -= im(c+1,c+0)*imRes(c+1) : x
    // 2: imRes(c+0) += im(c+1,c+0)*reRes(c+1) : imRes(c+0)
    // 3: x          += im(c+1,c+0)*reRes(c+1) : x
    }
    {
    __m256d invD0 = _mm256_broadcast_sd((const double*)&col[cHalf].im+0); // 1/re(c+0,c+0)
    resC = _mm256_mul_pd(resC, invD0);            // reRes(c+0) x          imRes(c+0) x
    resC = _mm256_blend_pd(resC, resC1, 0xA);     // reRes(c+0) reRes(c+1) imRes(c+0) imRes(c+1)
    result[cHalf] = resC;
    }
    {
    __m256d re0im1ResC = _mm256_shuffle_pd(resC, resC, 0xC); //  reRes(c+0)  reRes(c+0)  imRes(c+1)  imRes(c+1)
    __m256d re1im0ResC = _mm256_shuffle_pd(resC, resC, 0x3); //  reRes(c+1)  reRes(c+1)  imRes(c+0)  imRes(c+0)

    resC = _mm256_xor_pd(resC, negMask);                     // -reRes(c+1) -reRes(c+1)  imRes(c+0)  imRes(c+0)
    resC = _mm256_permute2f128_pd(resC, resC, 1);            //  imRes(c+0)  imRes(c+1) -reRes(c+0) -reRes(c+1)
    {
    __m256d im0re1ResC = _mm256_shuffle_pd(resC, resC, 0xC); //  imRes(c+0)  imRes(c+0) -reRes(c+1) -reRes(c+1)
    __m256d im1re0ResC = _mm256_shuffle_pd(resC, resC, 0x3); //  imRes(c+1)  imRes(c+1) -reRes(c+0) -reRes(c+0)

    unsigned k;
    for (k = 0; k != cHalf; ++k)
    {
      __m256d res = result[k]; // reRes(k*2+0)  reRes(k*2+1)  imRes(k*2+0)  imRes(k*2+1)
      __m256d re = col[k].re;  // re(c+0,k*2+0) re(c+0,k*2+1) re(c+1,k*2+0) re(c+1,k*2+1)
      __m256d im = col[k].im;  // im(c+0,k*2+0) im(c+0,k*2+1) im(c+1,k*2+0) im(c+1,k*2+1)

      MSUB(res, re0im1ResC, re);
      // 0: -= re(c+0,k*2+0)*reRes(c+0) : reRes(k*2+0)
      // 1: -= re(c+0,k*2+1)*reRes(c+0) : reRes(k*2+1)
      // 2: -= re(c+1,k*2+0)*imRes(c+1) : imRes(k*2+0)
      // 3: -= re(c+1,k*2+1)*imRes(c+1) : imRes(k*2+1)

      MSUB(res, im0re1ResC, im);
      // 0: -= im(c+0,k*2+0)*imRes(c+0) : reRes(k*2+0)
      // 1: -= im(c+0,k*2+1)*imRes(c+0) : reRes(k*2+1)
      // 2: += im(c+1,k*2+0)*reRes(c+1) : imRes(k*2+0)
      // 3: += im(c+1,k*2+1)*reRes(c+1) : imRes(k*2+1)

      re = _mm256_permute2f128_pd(re, re, 1);  // re(c+1,k*2+0) re(c+1,k*2+1) re(c+0,k*2+0) re(c+0,k*2+1)
      MSUB(res, re1im0ResC, re);
      // 0: -= re(c+1,k*2+0)*reRes(c+1) : reRes(k*2+0)
      // 1: -= re(c+1,k*2+1)*reRes(c+1) : reRes(k*2+1)
      // 2: += re(c+0,k*2+0)*imRes(c+0) : imRes(k*2+0)
      // 3: += re(c+0,k*2+1)*imRes(c+0) : imRes(k*2+1)

      im = _mm256_permute2f128_pd(im, im, 1);  // im(c+1,k*2+0) im(c+1,k*2+1) im(c+0,k*2+0) im(c+0,k*2+1)
      MSUB(res, im1re0ResC, im);
      // 0: -= im(c+1,k*2+0)*imRes(c+1) : reRes(k*2+0)
      // 1: -= im(c+1,k*2+1)*imRes(c+1) : reRes(k*2+1)
      // 2: += im(c+0,k*2+0)*reRes(c+0) : imRes(k*2+0)
      // 3: += im(c+0,k*2+1)*reRes(c+0) : imRes(k*2+1)

      result[k] = res;
    }
    }
    }
    }
  }
}
