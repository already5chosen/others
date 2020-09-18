#include <complex>
#include <cfloat>
#include <cstring>
#include <immintrin.h>
#include "chol.h"

#define CHOL_DBG 0

#if CHOL_DBG
bool chol_scalar(std::complex<double> *dst, const std::complex<double> *src, int n);
void PrintLowerTriangle(const std::complex<double> *src, unsigned n);
void PrintLowerTriangle(std::complex<double> *dst, const __m256d* triang, unsigned n);
void PrintInternalTriangle(const std::complex<double> *src, unsigned n);
#endif

#ifdef __AVX2__
#define MADD(   acc, x, y) acc = _mm256_fmadd_pd (x, y, acc)
#define MSUB(   acc, x, y) acc = _mm256_fnmadd_pd(x, y, acc)
#define MADDSUB(acc, x, y) acc = _mm256_addsub_pd(acc, _mm256_mul_pd(x, y))
#else
#define MADD(   acc, x, y) acc = _mm256_add_pd   (acc, _mm256_mul_pd(x, y))
#define MSUB(   acc, x, y) acc = _mm256_sub_pd   (acc, _mm256_mul_pd(x, y))
#define MADDSUB(acc, x, y) acc = _mm256_addsub_pd(acc, _mm256_mul_pd(x, y))
#endif


#ifdef __AVX2__
#define MADD128(   acc, x, y) acc = _mm_fmadd_pd (x, y, acc)
#define MSUB128(   acc, x, y) acc = _mm_fnmadd_pd(x, y, acc)
#define MADDSUB128(acc, x, y) acc = _mm_addsub_pd(acc, _mm_mul_pd(x, y))
#else
#define MADD128(   acc, x, y) acc = _mm_add_pd   (acc, _mm_mul_pd(x, y))
#define MSUB128(   acc, x, y) acc = _mm_sub_pd   (acc, _mm_mul_pd(x, y))
#define MADDSUB128(acc, x, y) acc = _mm_addsub_pd(acc, _mm_mul_pd(x, y))
#endif

namespace {

inline
void PackUpperTriangle_R(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // transpose - read rows, store column
  #if 0
  if ((n & 1)==0) { // N=2*k, pad odd rows of triangle
    for (unsigned row = 0; row < n; ++row) {
      std::complex<double>* pOut = &triang[row];
      for (unsigned col = 0; col < row+1; ++col) {
        *pOut = src[col];
        pOut += (n-col) & unsigned(-2); // 0:n, 1:n-2, 2:n-2, 3:n-3, ...
      }
      src += n;
    }
  } else {
    for (unsigned row = 0; row < n; ++row) {
      std::complex<double>* pOut = &triang[row+1];
      for (unsigned col = 0; col < row+1; ++col) {
        *pOut = src[col];
        pOut += (n-col) & unsigned(-2); // 0:n-1, 1:n-1, 2:n-3, 3:n-3, ...
      }
      src += n;
    }
  }
  #endif
  for (unsigned row = 0; row < n; ++row) {
    std::complex<double>* pOut = &triang[row+ (n&1)];
    for (unsigned col = 0; col < row+1; ++col) {
      *pOut = src[col];
      pOut += (n-col) & unsigned(-2); // 0:n-1, 1:n-1, 2:n-3, 3:n-3, ...
    }
    src += n;
  }

  auto pZero = triang;
  unsigned step = n*2;
  if ((n & 1)==0) {
    pZero += n;
    step  -= 2;
  }
  for (;step >= 4; step -= 4) {
    *pZero = 0;
    pZero += step;
  }
}

inline
void PackUpperTriangle_C(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  const unsigned hn = n / 2;
  if ((n & 1)==1) { // N=2*k+1
    triang[0] = 0;
    triang[1] = src[0];
    for (int col = 0; col < int(hn)*2; ++col)
      triang[col+2] = src[col+1];
    triang += n+1;
    src    += n+1;
  }
  for (unsigned hr = hn; ; --hr) {
    auto y0 = &triang[0];
    auto y1 = &triang[hr*2];
    auto x0 = &src[0];
    auto x1 = &src[n];
    y0[0] = x0[0];
    y0[1] = x0[1];
    y1[0] = 0;
    y1[1] = x1[1];
    if (hr == 1)
      break;
    for (int col = 0; col < int(hr-1)*2; ++col) {
      y0[col+2] = x0[col+2];
      y1[col+2] = x1[col+2];
    }
    triang += hr*4;
    src    += n*2+2;
  }
}

inline
void PackUpperTriangle_P(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  if ((n & 1)==0) { // N=2*k+1
    memcpy(triang, src, sizeof(*src)*n);
    triang += n;
    src    += n;
  }
  for (unsigned hr = (n-1)/2; hr > 0; --hr) {
    *triang++ = 0;
    memcpy(triang, src, sizeof(*src)*(hr*4+1));
    triang += hr*4+1;
    src    += hr*4+1;
  }
  *triang++ = 0;
  *triang++ = *src++;
}

inline
void PackUpperTriangle_Q(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // transpose - read rows, store column
  for (unsigned row = 0; row < n; ++row) {
    std::complex<double>* pOut = &triang[row+ (n&1)];
    for (unsigned col = 0; col < row+1; ++col) {
      *pOut = *src++;
      pOut += (n-col) & unsigned(-2); // 0:n-1, 1:n-1, 2:n-3, 3:n-3, ...
    }
  }

  auto pZero = triang;
  unsigned step = n*2;
  if ((n & 1)==0) {
    pZero += n;
    step  -= 2;
  }
  for (;step >= 4; step -= 4) {
    *pZero = 0;
    pZero += step;
  }
}

inline
void PackUpperTriangle(std::complex<double>* triang, const std::complex<double> *src, unsigned n, int srcLayout)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  if (n == 0)
    return;
  switch (srcLayout)
  {
  case 'C':
    PackUpperTriangle_C(triang, src, n);
    break;
  case 'P':
    PackUpperTriangle_P(triang, src, n);
    break;
  case 'Q':
    PackUpperTriangle_Q(triang, src, n);
    break;
  case 'R':
  default:
    PackUpperTriangle_R(triang, src, n);
    break;
  }
}

inline
void UnpackUpperTriangle(std::complex<double> *dst, const std::complex<double>* triang, unsigned n)
{ // unpack upper triangle with diagonal
  // The result in triang is a conjugate of proper upper triangle,
  // so copying it into lower triangle of dst matrix is what the doctor ordered
  // transpose - read rows, store column
  for (unsigned row = 0; row < n; ++row) {
    std::complex<double>* pOut = dst;
    triang += (n-row) & 1;      // skip padding
    *pOut = (*triang++).real(); // zeroize imagery part of diagonal
    pOut += n;
    for (unsigned col = 1; col < n-row; ++col) {
      *pOut = *triang++;
      pOut += n;
    }
    dst += n+1; // to next diagonal element
  }
}

const double DIAG_MIN   =  double(FLT_MIN);
// FLT_MIN is not really special in context of double-precision calculations
// but I wanted very small number that is still much bigger than DBL_MIN
// and FLT_MIN looks good enough
const double DIAG_SUBST =  1E40;

inline bool chol_FactorizeAndSolveFwd(std::complex<double>* triang, unsigned N, std::complex<double>* result, const bool solveFwd)
{
  bool succ = true;
  if ((N & 1) != 0) { // special handling for the first row of matrix with odd number of elements
    // process top row
    auto x0 = &triang[1];
    double aa = x0[0].real(); // diagonal element
    // check that we are positive definite
    // printf("%u %e\n", rlen, aa);
    if (aa < DIAG_MIN) {
      aa = DIAG_SUBST;
      succ = false;
    }
    double aaSqrt    = sqrt(aa);
    double aaInvSqrt = 1.0/aaSqrt;
    x0[0].real(aaSqrt);
    x0[0].imag(aaInvSqrt);

    unsigned hlen=N/2;
    if (hlen==0)
      return succ; // x0 was last (and only) row

    x0 += 1;
    if (!solveFwd) {
      for (int i = 0; i < int(hlen*2); ++i)
        x0[i] *= aaInvSqrt;
    } else { // combine multiplication by invSqrt with Forward propagation
      auto r = result[0] * aaInvSqrt;
      result[0] = r;
      result += 1;
      #if 1
      for (int i = 0; i < int(hlen*2); ++i) {
        auto x0val = x0[i] * aaInvSqrt;
        x0[i]      = x0val;
        result[i] -= x0val*r;
      }
      #else
      __m256d r_re = _mm256_set1_pd(r.real());
      __m256d r_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(r.imag()));
      __m256d vaInvSqrt = _mm256_set1_pd(aaInvSqrt);
      int cnt = int(hlen), idx = 0;
      #ifdef __clang__
      #pragma unroll 1
      #endif
      do {
        __m256d vx0 = _mm256_loadu_pd((const double*)&x0[idx]);
        __m256d vr  = _mm256_loadu_pd((const double*)&result[idx]);
        vx0 = _mm256_mul_pd(vx0, vaInvSqrt);
        MSUB(vr,  vx0, r_re);
        _mm256_storeu_pd((double*)&x0[idx], vx0);
        vx0 = _mm256_permute_pd(vx0, 5);
        MSUB(vr,  vx0, r_im);
        _mm256_storeu_pd((double*)&result[idx], vr);

        idx += 2;
      } while (--cnt);
      #endif
    }
    // subtract outer product of top row from lower part of the matrix
    // process two output rows together
    triang += N+1;
    auto y = triang;
    do {
      auto y0 = &y[0];
      auto y1 = &y[hlen*2];
      #if 0
      auto f00 = x0[0];
      auto f01 = x0[1];
      for (int c = 0; c < int(hlen)*2; ++c) {
          y0[c] -= x0[c]*conj(f00);
          y1[c] -= x0[c]*conj(f01);
      }
      #else
      __m256d f00_re = _mm256_set1_pd(x0[0].real());
      __m256d f01_re = _mm256_set1_pd(x0[1].real());
      __m256d f00_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(x0[0].imag()));
      __m256d f01_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(x0[1].imag()));
      int cnt = int(hlen), idx = 0;
      #ifdef __clang__
      #pragma unroll 1
      #endif
      do {
        __m256d vx0 = _mm256_loadu_pd((const double*)&x0[idx]);
        __m256d vy0 = _mm256_loadu_pd((const double*)&y0[idx]);
        __m256d vy1 = _mm256_loadu_pd((const double*)&y1[idx]);
        MSUB(vy0, vx0, f00_re);
        MSUB(vy1, vx0, f01_re);
        vx0 = _mm256_permute_pd(vx0, 5);
        MADD(vy0, vx0, f00_im);
        MADD(vy1, vx0, f01_im);
        _mm256_storeu_pd((double*)&y0[idx], vy0);
        _mm256_storeu_pd((double*)&y1[idx], vy1);
        idx += 2;
      } while (--cnt);
      #endif
      x0 += 2;
      y  += hlen*4;
      --hlen;
    } while (hlen > 0);

  }

  int ltMask = 0;
  for (unsigned rhlen = N / 2; ; rhlen-=1) {
    // process top row
    auto x0 = &triang[0];
    auto x1 = &triang[rhlen*2];

    double aa0 = x0[0].real(); // diagonal element
    auto f     = x0[1];
    // process next after top row
    // subtract outer product of top row from next after top row
    double aa1 = x1[1].real()*aa0 - norm(f); // diagonal element
    #if 0
    // check that we are positive definite
    // printf("%u %e %e\n", rlen, aa0. aa1);
    double aa[] = { aa0, aa1 }, aaSqrt[2], aaInvSqrt[2];
    for (int k=0; k < 2; ++k) {
      double v = aa[k];
      if (v < DIAG_MIN) {
        v = DIAG_SUBST;
        succ = false;
      }
      aaSqrt[k]    = v = sqrt(v);
      aaInvSqrt[k] = 1.0/v;
    }
    #else
    __m128d aa = _mm_setr_pd(aa0, aa1);
    __m128d lt = _mm_cmplt_pd(aa, _mm_set1_pd(DIAG_MIN));
    aa = _mm_blendv_pd(aa, _mm_set1_pd(DIAG_SUBST), lt);
    ltMask |= _mm_movemask_pd(lt);
    __m128d aaSqrt = _mm_sqrt_pd(aa);
    __m128d aaInvSqrt = _mm_div_pd(_mm_set1_pd(1.0),aaSqrt);
    #endif
    double aa0InvSqrt = aaInvSqrt[0];
    double aa1InvSqrt = aaInvSqrt[1]*aaSqrt[0];
    x0[0].real(aaSqrt[0]);
    x0[0].imag(aa0InvSqrt);
    x1[1].real(aaSqrt[1]*aa0InvSqrt);
    x1[1].imag(aa1InvSqrt);
    x0[1] = (f *= aa0InvSqrt);

    if (solveFwd) {
      auto r0 = result[0] * aa0InvSqrt;
      auto r1 = (result[1] - f*r0) * aa1InvSqrt;
      result[0] = r0;
      result[1] = r1;
    }

    if (rhlen<=1)
      break; // x1 was last row

    x0 += 2;
    x1 += 2;
    if (!solveFwd) {
      f *= aa0InvSqrt;
      for (int c = 0; c < int(rhlen-1)*2; ++c) {
        auto ax0 = x0[c];
        auto ax1 = x1[c];
        x0[c] = ax0 * aa0InvSqrt;
        x1[c] = (ax1 - ax0*conj(f))*aa1InvSqrt;
      }
    } else {
      #if 0
      f *= aa0InvSqrt;
      for (int c = 0; c < int(rhlen-1)*2; ++c) {
        auto ax0 = x0[c];
        auto ax1 = x1[c];
        x0[c] = ax0 * aa0InvSqrt;
        x1[c] = (ax1 - ax0*conj(f))*aa1InvSqrt;
      }
      // Forward propagation of pair of results
      auto r0 = result[0];
      auto r1 = result[1];
      result += 2;
      for (int c = 0; c < int(rhlen-1)*2; ++c) {
        result[c] -= x0[c]*r0 + x1[c]*r1;
      }
      #else
      // combine handling of pair of top rows with forward propagation of pair of results
      auto r0 = result[0];
      auto r1 = result[1];
      __m256d f_re  = _mm256_set1_pd(f.real());
      __m256d r0_re = _mm256_set1_pd(r0.real());
      __m256d r1_re = _mm256_set1_pd(r1.real());
      __m256d f_im  = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(f.imag()));
      __m256d r0_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(r0.imag()));
      __m256d r1_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(r1.imag()));
      __m256d va0InvSqrt = _mm256_set1_pd(aa0InvSqrt);
      __m256d va1InvSqrt = _mm256_set1_pd(aa1InvSqrt);
      result += 2;
      int cnt = int(rhlen-1), idx = 0;
      #ifdef __clang__
      #pragma unroll 1
      #endif
      do {
        __m256d vx0 = _mm256_loadu_pd((const double*)&x0[idx]);
        __m256d vx1 = _mm256_loadu_pd((const double*)&x1[idx]);
        __m256d vr  = _mm256_loadu_pd((const double*)&result[idx]);
        vx0 = _mm256_mul_pd(vx0, va0InvSqrt);
        MSUB(vx1, vx0, f_re);
        MSUB(vr,  vx0, r0_re);
        _mm256_storeu_pd((double*)&x0[idx], vx0);

        vx0 = _mm256_permute_pd(vx0, 5);
        MADD(vx1, vx0, f_im);
        MSUB(vr,  vx0, r0_im);
        vx1 = _mm256_mul_pd(vx1, va1InvSqrt);
        _mm256_storeu_pd((double*)&x1[idx], vx1);

        MSUB(vr,  vx1, r1_re);
        vx1 = _mm256_permute_pd(vx1, 5);
        MSUB(vr,  vx1, r1_im);
        _mm256_storeu_pd((double*)&result[idx], vr);

        idx += 2;
      } while (--cnt);
      #endif
    }

    // subtract outer product of two top rows from lower part of the matrix
    // process two output rows together
    triang += rhlen*4;
    auto y = triang;
    for (unsigned chlen = rhlen-1; chlen > 0; y += chlen*4, chlen -= 1) {
      #if 0
      auto f00 = x0[0];
      auto f10 = x1[0];
      auto f01 = x0[1];
      auto f11 = x1[1];
      auto y0 = &y[0];
      auto y1 = &y[chlen*2];
      for (int c = 0; c < int(chlen)*2; ++c) {
        y0[c] = y0[c] - x0[c]*conj(f00) - x1[c]*conj(f10);
        y1[c] = y1[c] - x0[c]*conj(f01) - x1[c]*conj(f11);
      }
      #else
      __m256d f00_re = _mm256_set1_pd(x0[0].real());
      __m256d f10_re = _mm256_set1_pd(x1[0].real());
      __m256d f01_re = _mm256_set1_pd(x0[1].real());
      __m256d f11_re = _mm256_set1_pd(x1[1].real());
      __m256d f00_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(x0[0].imag()));
      __m256d f10_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(x1[0].imag()));
      __m256d f01_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(x0[1].imag()));
      __m256d f11_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(x1[1].imag()));
      auto y0 = &y[0];
      auto y1 = &y[chlen*2];
      int cnt = int(chlen), idx = 0;
      #ifdef __clang__
      #pragma unroll 1
      #endif
      do {
        __m256d vx0 = _mm256_loadu_pd((const double*)&x0[idx]);
        __m256d vy0 = _mm256_loadu_pd((const double*)&y0[idx]);
        __m256d vx1 = _mm256_loadu_pd((const double*)&x1[idx]);
        __m256d vy1 = _mm256_loadu_pd((const double*)&y1[idx]);
        MSUB(vy0, vx0, f00_re);
        MSUB(vy1, vx0, f01_re);
        MSUB(vy0, vx1, f10_re);
        MSUB(vy1, vx1, f11_re);
        vx0 = _mm256_permute_pd(vx0, 5);
        vx1 = _mm256_permute_pd(vx1, 5);
        MADD(vy0, vx0, f00_im);
        MADD(vy1, vx0, f01_im);
        MADD(vy0, vx1, f10_im);
        MADD(vy1, vx1, f11_im);
        _mm256_storeu_pd((double*)&y0[idx], vy0);
        _mm256_storeu_pd((double*)&y1[idx], vy1);
        idx += 2;
      } while (--cnt);
      #endif
      x0 += 2;
      x1 += 2;
    }
  }

  if (ltMask)
    succ = false;
  return succ;
}

bool chol_Factorize(std::complex<double>* triang, unsigned N)
{
  return chol_FactorizeAndSolveFwd(triang, N, NULL, false);
}

bool chol_Factorize(std::complex<double>* triang, unsigned N, std::complex<double>* result)
{
  return chol_FactorizeAndSolveFwd(triang, N, result, true);
}

void chol_SolveFwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
{
  // x = R' \ conj(B);
  // R' is lower triangle - solve by forward propagation (caxpy-like)
  // x = B;
  // for r=1:N
    // x(r) = x(r)/R(r,r);
    // x(r+1:N) -= R(r, r+1:N).'*x(r);
  // end
  if ((N & 1) != 0) {    // special handling for the first row of matrix with odd number of elements
    auto y = &triang[1]; // point to diag element
    auto xr = x[0] * y[0].imag(); // imag() of diag element contains inverse of real() part
    x[0] = xr;
    if (N <= 1)
      return;
    x += 1;
    y += 1;
    unsigned hlen = N / 2;
    for (int c = 0; c < int(hlen)*2; ++c)
      x[c] -= y[c]*xr;
    triang += N+1;
  }

  // process two rows per iteration
  for (unsigned rhlen = N/2; ; --rhlen) {
    auto y0 = &triang[0];
    auto y1 = &triang[rhlen*2];

    auto xr0 = x[0] * y0[0].imag(); // imag() of diag element contains inverse of real() part
    auto xr1 = (x[1] - y0[1]*xr0) * y1[1].imag(); // odd row is padded

    x[0] = xr0;
    x[1] = xr1;
    if (rhlen == 1)
      break;

    y0 += 2;
    y1 += 2;
    x  += 2;
    #if 0
    for (int c = 0; c < int(rhlen-1)*2; ++c)
      x[c] = x[c] - y0[c]*xr0 - y1[c]*xr1;
    #else
    __m256d xr0_re = _mm256_set1_pd(xr0.real());
    __m256d xr1_re = _mm256_set1_pd(xr1.real());
    __m256d xr0_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(xr0.imag()));
    __m256d xr1_im = _mm256_addsub_pd(_mm256_setzero_pd(), _mm256_set1_pd(xr1.imag()));
    int cnt = int(rhlen-1), idx = 0;
    #ifdef __clang__
    #pragma unroll 1
    #endif
    do {
      __m256d vy0 = _mm256_loadu_pd((const double*)&y0[idx]);
      __m256d vy1 = _mm256_loadu_pd((const double*)&y1[idx]);
      __m256d vx  = _mm256_loadu_pd((const double*)&x[idx]);
      MSUB(vx, vy0, xr0_re);
      MSUB(vx, vy1, xr1_re);
      vy0 = _mm256_permute_pd(vy0, 5);
      vy1 = _mm256_permute_pd(vy1, 5);
      MSUB(vx, vy0, xr0_im);
      MSUB(vx, vy1, xr1_im);
      _mm256_storeu_pd((double*)&x[idx], vx);
      idx += 2;
    } while (--cnt);
    #endif

    triang += rhlen*4;
  }
}

void chol_SolveBwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
{
  // x = conj(R) \ x;
  // R is upper triangle - solve by backward substitution (dot-like)
  // for r=N:-1:1
   // x(r) = x(r) - sum(x(r+1:N) .* conj(R(r,r+1:N)(:)) );
   // x(r) = x(r)/R(r,r);
  // end
  triang += (N+1)*(N+1)/2; // point past end
  x += N;                  // point past end

  // process two rows per iteration
  const unsigned hlen = N/2;
  for (unsigned rhlen = 1; rhlen <= hlen; ++rhlen) {
    triang -= rhlen*4; // point to diag element
    auto y0 = &triang[0];
    auto y1 = &triang[rhlen*2];
    x -= 2;
    auto acc0 = x[0];
    auto acc1 = x[1];
    #if 0
    for (int c = 2; c < int(rhlen)*2; ++c) {
      acc0 -= x[c] * conj(y0[c]);
      acc1 -= x[c] * conj(y1[c]);
    }
    #else
    if (rhlen > 1) {
      __m256d vacc0_re = _mm256_setr_pd(acc0.real(),0,0,0);
      __m256d vacc1_re = _mm256_setr_pd(acc1.real(),0,0,0);
      __m256d vacc0_im = _mm256_setr_pd(acc0.imag(),0,0,0);
      __m256d vacc1_im = _mm256_setr_pd(acc1.imag(),0,0,0);
      int cnt = int(rhlen-1), idx = 2;
      #ifdef __clang__
      #pragma unroll 1
      #endif
      do {
        __m256d vy0 = _mm256_loadu_pd((const double*)&y0[idx]);
        __m256d vy1 = _mm256_loadu_pd((const double*)&y1[idx]);
        __m256d vx  = _mm256_loadu_pd((const double*)&x[idx]);
        MSUB(vacc0_re, vy0, vx);
        MSUB(vacc1_re, vy1, vx);
        vx = _mm256_permute_pd(vx, 5);
        MSUB(vacc0_im, vy0, vx);
        MSUB(vacc1_im, vy1, vx);
        idx += 2;
      } while (--cnt);
      __m256d vacc01_re = _mm256_hadd_pd(vacc0_re, vacc1_re);
      __m256d vacc01_im = _mm256_hsub_pd(vacc0_im, vacc1_im);
      __m256d vacc0 = _mm256_unpacklo_pd(vacc01_re, vacc01_im);
      __m256d vacc1 = _mm256_unpackhi_pd(vacc01_re, vacc01_im);
      __m256d vacc01a = _mm256_permute2f128_pd(vacc0, vacc1, 0x21);
      __m256d vacc01b = _mm256_blend_pd(vacc0, vacc1, 0xC);
      __m256d vacc01 = _mm256_add_pd(vacc01a, vacc01b);
      acc0.real(vacc01[0]);
      acc0.imag(vacc01[1]);
      acc1.real(vacc01[2]);
      acc1.imag(vacc01[3]);
    }
    #endif
    acc1 *= y1[1].imag(); // imag() of diag element contains inverse of it's real()
    acc0 -= acc1*conj(y0[1]);
    acc0 *= y0[0].imag();
    x[0] = acc0;
    x[1] = acc1;
  }

  if ((N & 1) != 0) { // special handling for the first row of matrix with odd number of elements
    triang -= hlen*2;
    #ifdef GCC10_WORKAROUND
    std::complex<double> acc = 0;
    for (int c = 0; c < int(hlen)*2; ++c)
      acc += x[c] * conj(triang[c]);
    acc = x[-1] - acc;
    #else
    auto acc = x[-1];
    for (int c = 0; c < int(hlen)*2; ++c)
      acc -= x[c] * conj(triang[c]);
    #endif
    x[-1] = acc * triang[-1].imag(); // imag() of diag element contains inverse of it's real()
  }
}

};

unsigned chol_getWorkBufferSize(int n)
{
  return (n+1)*(n+3)/2*sizeof(std::complex<double>);
}

// chol - Perform Cholesky decomposition of complex Hermitian matrix
// Decomposition uses basic Cholesky algorithm (outer product) without pivoting.
// Inner loop processes two input and two output columns
// Parameters:
// src - input matrix. See srcLayout for description of src layout.
//       In 'R' and 'C' layouts only lower triangle (including diagonal) of input matrix is used in calculation, values in upper triangle are ignored
// dst - output matrix stored in row-wise order
//       only lower triangle  (including diagonal) of output matrix guaranteed to be correct,
//       values in upper triangle are not guaranteed to be zeroed
// n   - number of rows/columns in matrix a. supported range [1..16383]
// workBuffer - temporary buffer.
//       The size of buffer is returned by chol_getWorkBufferSize
//       For best performance buffer has to be aligned of 64-byte boundary. However  for compatibility with future implementations it's recommended to align it on 128-byte boundary.
// srcLayout - source matrix layout. Possible values:
//       'R' - Input matrix stored in row-wise order. This is default.
//       'C' - Input matrix stored in column-wise order, as in Fortran/Matlab
//       'P' - Input matrix stored in packed column-wise order, i.e. only lower triangle present.
//             This layout is identical to result of Matlab's aa(tril(true(n))), where aa is full matrix.
//       'Q' - Input matrix stored in packed row-wise order, i.e. only lower triangle present.
// Return value:
// true = success, false - input matrix is not positive definite
bool chol(std::complex<double> *dst, const std::complex<double> *src, int n, void* workBuffer, int srcLayout)
{
  const int N_MIN = 1;
  const int N_MAX = 16383;
  if (n < N_MIN || n > N_MAX)
    return false;

#if CHOL_DBG
  PrintLowerTriangle(src, n);
  chol_scalar(dst, src, n);
  printf("ref:\n");
  PrintLowerTriangle(dst, n);
#endif

  std::complex<double>* packedResult = static_cast<std::complex<double>*>(workBuffer) + (n & 1);
  std::complex<double>* triang = packedResult + n;
  PackUpperTriangle(triang, src, n, srcLayout);
  bool succ = chol_Factorize(triang, n);
  UnpackUpperTriangle(dst, triang, n);

#if CHOL_DBG
  printf("internal:\n");
  PrintInternalTriangle(triang, n);
  printf("res:\n");
  PrintLowerTriangle(dst, n);
  return false;
#endif

  return succ;
}

// chol_solver - solve complex Hermitian matrix by means of Cholesky decomposition
// Decomposition uses basic Cholesky algorithm (outer product) without pivoting.
// Inner loop processes two input and two output columns
// Parameters:
// src    - input matrix. See srcLayout for description of src layout.
//          In 'R' and 'C' layouts only lower triangle (including diagonal) of input matrix is used in calculation, values in upper triangle are ignored
// vecB   - complex right-hand column vector
// result - complex solution column vector
//          It is legal for vecB and result to pont to the same location
// n      - number of rows/columns in matrix a, also number of elements (rows) in vectors vecB and result
// workBuffer - temporary buffer.
//          The size of buffer is returned by chol_getWorkBufferSize
//          For best performance buffer has to be aligned of 64-byte boundary. However  for compatibility with future implementations it's recommended to align it on 128-byte boundary.
// srcLayout - source matrix layout. Possible values:
//       'R' - Input matrix stored in row-wise order. This is default.
//       'C' - Input matrix stored in column-wise order, as in Fortran/Matlab
//       'P' - Input matrix stored in packed column-wise order, i.e. only lower triangle present.
//             This layout is identical to result of Matlab's aa(tril(true(n))), where aa is full matrix.
//       'Q' - Input matrix stored in packed row-wise order, i.e. only lower triangle present.
// Return value:
// true = success, false - input matrix is not positive definite
bool chol_solver(std::complex<double> *result, const std::complex<double> *src, const std::complex<double> *vecB, int n, void* workBuffer, int srcLayout)
{
  const int N_MIN = 1;
  const int N_MAX = 16383;
  if (n < N_MIN || n > N_MAX)
    return false;

#if CHOL_DBG
  PrintLowerTriangle(src, n);
  std::complex<double>* dst = new std::complex<double>[n*n];
  chol_scalar(dst, src, n);
  PrintLowerTriangle(dst, n);
  delete [] dst;
#endif

  std::complex<double>* packedResult = static_cast<std::complex<double>*>(workBuffer) + (n & 1);
  std::complex<double>* triang = packedResult + n;
  PackUpperTriangle(triang, src, n, srcLayout);
  memcpy(packedResult, vecB, sizeof(*vecB)*n); // PackVecB(packedResult, vecB, n);
  bool succ = chol_Factorize(triang, n, packedResult);
  if (succ) {
    chol_SolveBwd(packedResult, n, triang);
    memcpy(result, packedResult, sizeof(*result)*n);
  }
  return succ;
}

// chol_trif_solver - solve a system of linear equations A*result=vecB with a complex
// Hermitian positive definite matrix src utilizing ready triangular factor computed by
// chol or chol_solver.
// Parameters:
// vecB   - complex right-hand column vector
// result - complex solution column vector
//          It is legal for vecB and result to pont to the same location
// n      - number of rows/columns in matrix a, also number of elements (rows) in vectors vecB and result
// workBuffer - temporary buffer containing triangular factor for matrix A.
//          workBuffer should be prepared by earlier call to chol() or chol_solver() routines.
// Return value: none
// true = success, false - input matrix is not positive definite
void chol_trif_solver(std::complex<double> *result, const std::complex<double> *vecB, int n, void* workBuffer)
{
  const int N_MIN = 1;
  const int N_MAX = 16383;
  if (n < N_MIN || n > N_MAX)
    return;

  std::complex<double>* packedResult = static_cast<std::complex<double>*>(workBuffer) + (n & 1);
  memcpy(packedResult, vecB, sizeof(*vecB)*n); // PackVecB(packedResult, vecB, n);
  const std::complex<double>* triang = packedResult + n;
  chol_SolveFwd(packedResult, n, triang);
  chol_SolveBwd(packedResult, n, triang);
  memcpy(result, packedResult, sizeof(*result)*n);
}

#if CHOL_DBG
// chol - Perform Cholesky decomposition of complex Hermitian matrix
// Decomposition uses combination of Cholesky-Crout and Cholesky-Banachiewicz algorithms without pivoting.
// Inner loop calculates two subsequent elements of the same column
// Parameters:
// src - input matrix stored in row-wise order
//       only lower triangle (including diagonal) of input matrix is used in calculation, values in upper triangle are ignored
// dst - output matrix stored in row-wise order
//       only lower triangle  (including diagonal) of output matrix guaranteed to be correct,
//       values in upper triangle are not guaranteed to be zeroed
// n - number of rows/columns in matrix src
// Return value:
// true = success, false - input matrix is not positive definite
bool chol_scalar(std::complex<double> *dst, const std::complex<double> *src, int n)
{
  const int N_MIN = 1;
  const int N_MAX = 128;
  const int SIMD_FACTOR = 4;
  if (n < N_MIN || n > N_MAX)
    return false;

  // Copy lower triangle from source to destination.
  // Zero upper triangle - not necessary, but let's do it anyway
  const std::complex<double> *srcRow = src;
  std::complex<double> *dstRow = dst;
  for (int r = 0; r < n; ++r)
  {
    for (int col = 0; col <= r; ++col)
      dstRow[col] = srcRow[col];
    for (int col = r+1; col < n; ++col)
      dstRow[col] = 0;
    srcRow += n;
    dstRow += n;
  }

  const int MAX_DATASET_SIZE  = 24*1024;
  const int MAX_DATASET_NELEM = MAX_DATASET_SIZE/(2*sizeof(double));

  // Cholesky-Banachiewicz decomposition
  double invDiag[N_MAX];
  std::complex<double> *rowDataset = dst;

  { // dataset1st = 0
    // Calculate diagonal element
    double aa = rowDataset[0].real(); // current diagonal element

    // check that we are positive defined
    // FLT_MIN is not really special in context of double-precision calculations
    // but I wanted very small number that is still much bigger than DBL_MIN
    // and FLT_MIN looks good enough
    //printf("%d %e\n", r, aa);
    if (aa < double(FLT_MIN))
      return false;
    double aaSqrt = sqrt(aa);
    double aaInvSqrt = 1.0 / aaSqrt;
    rowDataset[0] = aaSqrt;
    invDiag[0] = aaInvSqrt;
    rowDataset += n;
  }

  for (int dataset1st = 1; dataset1st < n; )
  {
    // calculate height of the Crout band
    int datasetLast = dataset1st;
    int datasetNelem = 0;
    for (;datasetLast < n; datasetLast += 2)
    {
      datasetNelem += datasetLast*2+3;
      if (datasetNelem > MAX_DATASET_NELEM)
        break;
    }
    if (datasetLast == dataset1st)
      datasetLast = dataset1st + 2;
    if (datasetLast > n)
      datasetLast = n;

    std::complex<double> *rowC = dst;
    { // c=0
      // calculate section of column 0 of result
      std::complex<double> *rowR = rowDataset;
      double aaInvSqrt = invDiag[0];
      for (int r = dataset1st; r < datasetLast; ++r)
      {
        rowR[0] *= aaInvSqrt;
        rowR += n;
      }
      rowC += n;
    }

    for (int c = 1; c < datasetLast; ++c)
    {
      if (c >= dataset1st)
      {
        // Calculate diagonal element
        double diagSum = 0;
        for (int k = 0; k < c; ++k)
          diagSum += rowDataset[k].real()*rowDataset[k].real()+rowDataset[k].imag()*rowDataset[k].imag();
        double aa = rowDataset[c].real() - diagSum; // current diagonal element

        // check that we are positive defined
        // FLT_MIN is not really special in context of double-precision calculations
        // but I wanted very small number that is still much bigger than DBL_MIN
        // and FLT_MIN looks good enough
        //printf("%3d %e %d:%d\n", c, aa, dataset1st, datasetLast);
        if (aa < double(FLT_MIN))
          return false;
        double aaSqrt = sqrt(aa);
        double aaInvSqrt = 1.0 / aaSqrt;
        rowDataset[c] = aaSqrt;
        invDiag[c] = aaInvSqrt;
        rowDataset += n;
        ++dataset1st;
      }

      // calculate section of column c of result
      std::complex<double> *rowR0 = rowDataset;
      double aaInvSqrt = invDiag[c];
      bool startAtEven = (dataset1st & 1)==0;
      for (unsigned nRows = datasetLast - dataset1st; nRows != 0; )
      {
        if (startAtEven || nRows == 1)
        {
          std::complex<double> sum = 0;
          for (int k = 0; k < c; ++k)
            sum += rowR0[k]*conj(rowC[k]);
          rowR0[c] = (rowR0[c] - sum)*aaInvSqrt;
          rowR0 += n;
          nRows -= 1;
          startAtEven = false;
        }
        unsigned nRowPairs = nRows / 2;
        while (nRowPairs != 0)
        {
          std::complex<double> *rowR1 = rowR0 + n;
          std::complex<double> sum0 = 0;
          std::complex<double> sum1 = 0;
          for (int k = 0; k < c; ++k)
          {
            sum0 += rowR0[k]*conj(rowC[k]);
            sum1 += rowR1[k]*conj(rowC[k]);
          }
          rowR0[c] = (rowR0[c] - sum0)*aaInvSqrt;
          rowR1[c] = (rowR1[c] - sum1)*aaInvSqrt;
          rowR0 = rowR1 + n;
          --nRowPairs;
        }
        nRows &= 1;
      }
      rowC += n;
    }
  }

  return true;
}

void PrintLowerTriangle(const std::complex<double> *src, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
  {
    for (unsigned j=0; j <= i; ++j)
      printf("(%9f %9f) ", src[n*i+j].real(), src[n*i+j].imag());
    printf("\n");
  }
}

void PrintInternalTriangle(const std::complex<double> *src, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
  {
    for (unsigned j=0; j < n-i; ++j)
      printf("(%9f %9f) ", src[j].real(), src[j].imag());
    printf("\n");
    src += n-i;
  }
}

//void PrintLowerTriangle(std::complex<double> *dst, const complex_m256d* triang, unsigned n)
//{
//  UnpackLowerTriangle(dst, triang, n);
//  PrintLowerTriangle(dst, n);
//}

#endif
