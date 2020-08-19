#include <stdint.h>
#include <intrin.h>
#include "chol_internal_definitions2.h"

#ifdef use__AVX__

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

#endif /* use_AVX__ */
