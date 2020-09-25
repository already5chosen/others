#include <complex>
#include <cfloat>
#include <cstring>
//#include <intrin.h>
//#include <stdio.h>
#include "chol.h"

#define CHOL_DBG 0

#if CHOL_DBG
bool chol_scalar(std::complex<double> *dst, const std::complex<double> *src, int n);
void chol_SolveFwd_scalar(std::complex<double> *x, unsigned N, const std::complex<double>* mat, const std::complex<double> *vecB);
void chol_SolveBwd_scalar(std::complex<double> *x, unsigned N, const std::complex<double>* mat);
void PrintLowerTriangle(const std::complex<double> *src, unsigned n);
void PrintInternalTriangle(const double* src, unsigned n);
void PrintResult(const std::complex<double> *src, unsigned n);
#endif

namespace {

inline
void PackUpperTriangle_R(double* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // transpose - read rows, store column
  if (n & 3) {
    auto pZero = triang;
    for (unsigned row = 0; row < n; ++row) {
      memset(pZero, 0, sizeof(*pZero)*8);
      pZero += ((n-row+3) & unsigned(-4))*2;
    }
  }
  const unsigned xi0 = (0-n) & 3;
  for (unsigned row = 0; row < n; ++row) {
    const unsigned xi = xi0 + row;
    auto pOut = &triang[(xi/4)*8+(xi%4)];
    for (unsigned col = 0; col < row+1; ++col) {
      pOut[0] = src[col].real();
      pOut[4] = src[col].imag();
      pOut += ((n-col+2) & unsigned(-4))*2;
    }
    src += n;
  }
}

inline
void PackUpperTriangle_C(double* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  int src_i = 0;
  for (unsigned r = n; r > 0; --r) {
    int qr = (r + 3)/4;
    int pad_i = qr*4-r;
    if (src_i < pad_i) {
      for (int k = 0; k < 4; ++k) {
        triang[k+0] = k < pad_i ? 0 : src[src_i+k-pad_i].real();
        triang[k+4] = k < pad_i ? 0 : src[src_i+k-pad_i].imag();
      }
      triang += 8;
      qr     -= 1;
      pad_i  -= 4;
    }
    auto pSrc = &src[src_i-pad_i];
    src_i += n+1; // to next diagonal
    for (int i = 0; i < qr; ++i) {
      auto x0 = pSrc[0];
      auto x1 = pSrc[1];
      auto x2 = pSrc[2];
      auto x3 = pSrc[3];
      triang[0*4+0] = x0.real();
      triang[0*4+1] = x1.real();
      triang[0*4+2] = x2.real();
      triang[0*4+3] = x3.real();
      triang[1*4+0] = x0.imag();
      triang[1*4+1] = x1.imag();
      triang[1*4+2] = x2.imag();
      triang[1*4+3] = x3.imag();
      pSrc   += 4;
      triang += 8;
    }
  }
}

inline
void PackUpperTriangle_P(double* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  int src_i = 0;
  for (unsigned r = n; r > 0; --r) {
    int qr = (r + 3)/4;
    int pad_i = qr*4-r;
    if (src_i < pad_i) {
      for (int k = 0; k < 4; ++k) {
        triang[k+0] = k < pad_i ? 0 : src[src_i+k-pad_i].real();
        triang[k+4] = k < pad_i ? 0 : src[src_i+k-pad_i].imag();
      }
      triang += 8;
      qr     -= 1;
      pad_i  -= 4;
    }
    auto pSrc = &src[src_i-pad_i];
    src_i += r;
    for (int i = 0; i < qr; ++i) {
      auto x0 = pSrc[0];
      auto x1 = pSrc[1];
      auto x2 = pSrc[2];
      auto x3 = pSrc[3];
      triang[0*4+0] = x0.real();
      triang[0*4+1] = x1.real();
      triang[0*4+2] = x2.real();
      triang[0*4+3] = x3.real();
      triang[1*4+0] = x0.imag();
      triang[1*4+1] = x1.imag();
      triang[1*4+2] = x2.imag();
      triang[1*4+3] = x3.imag();
      pSrc   += 4;
      triang += 8;
    }
  }
}

inline
void PackUpperTriangle_Q(double* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // transpose - read rows, store column
  if (n & 3) {
    auto pZero = triang;
    for (unsigned row = 0; row < n; ++row) {
      memset(pZero, 0, sizeof(*pZero)*8);
      pZero += ((n-row+3) & unsigned(-4))*2;
    }
  }
  const unsigned xi0 = (0-n) & 3;
  for (unsigned row = 0; row < n; ++row) {
    const unsigned xi = xi0 + row;
    auto pOut = &triang[(xi/4)*8+(xi%4)];
    for (unsigned col = 0; col < row+1; ++col) {
      pOut[0] = src[col].real();
      pOut[4] = src[col].imag();
      pOut += ((n-col+2) & unsigned(-4))*2;
    }
    src += row+1;
  }
}

inline
void PackUpperTriangle(double* triang, const std::complex<double> *src, unsigned n, int srcLayout)
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
void UnpackUpperTriangle(std::complex<double> *dst, const double* triang, unsigned n)
{ // unpack upper triangle with diagonal
  // The result in triang is a conjugate of proper upper triangle,
  // so copying it into lower triangle of dst matrix is what the doctor ordered
  // transpose - read rows, store column
  for (unsigned rlen = n; rlen > 0; --rlen) {
    auto pOut = dst;
    if (rlen & 3) {
      for (int k = (0-rlen) & 3; k < 4; ++k) {
        *pOut = std::complex<double>(triang[k], triang[k+4]);
        pOut += n;
      }
      triang += 8;
    }
    for (unsigned rqlen = rlen / 4; rqlen > 0; --rqlen) {
      for (int k = 0; k < 4; ++k) {
        *pOut = std::complex<double>(triang[k], triang[k+4]);
        pOut += n;
      }
      triang += 8;
    }
    dst->imag(0); // zeroize imagery part of diagonal
    dst += n+1;   // to next diagonal element
  }
}

inline
void PackVecB(double* __restrict dst, const std::complex<double> *src, unsigned n)
{ // organize vecB in processing-friendly interleaved format
  int nr = int(n & 3);
  if (nr) {
    memset(dst, 0, sizeof(double)*8);
    int xi0 = 4-nr;
    for (int k = 0; k < nr; ++k) {
      dst[xi0+k+0] = src[k].real();
      dst[xi0+k+4] = src[k].imag();
    }
    src += nr;
    dst += 8;
  }
  for (int c = 0; c < int(n/4); ++c) {
    for (int k = 0; k < 4; ++k) {
      dst[k+0] = src[k].real();
      dst[k+4] = src[k].imag();
    }
    src += 4;
    dst += 8;
  }
}

inline
void UnpackResult(std::complex<double>* __restrict dst, const double* src, unsigned n)
{ // translate result from internal interleaved format to standard complex<double>
  int nr = int(n & 3);
  if (nr) {
    int xi0 = 4-nr;
    for (int k = 0; k < nr; ++k)
      dst[k] = std::complex<double>(src[xi0+k+0], src[xi0+k+4]);
    src += 8;
    dst += nr;
  }
  for (int c = 0; c < int(n/4); ++c) {
    for (int k = 0; k < 4; ++k)
      dst[k] = std::complex<double>(src[k+0], src[k+4]);
    src += 8;
    dst += 4;
  }
}


const double DIAG_MIN   =  double(FLT_MIN);
// FLT_MIN is not really special in context of double-precision calculations
// but I wanted very small number that is still much bigger than DBL_MIN
// and FLT_MIN looks good enough
const double DIAG_SUBST =  1E40;

inline bool chol_FactorizeAndSolveFwd(double* triang, unsigned N, double* result, const bool solveFwd)
{
  bool succ = true;
  if ((N & 1) != 0) { // special handling for the first row of matrix with odd number of elements
    // process top row
    int xi = (0-N) & 3;
    double aa = triang[4*0+xi]; // diagonal element
    // check that we are positive definite
    // printf("%u %e\n", rlen, aa);
    if (aa < DIAG_MIN) {
      aa = DIAG_SUBST;
      succ = false;
    }
    double aaSqrt    = sqrt(aa);
    double aaInvSqrt = 1.0/aaSqrt;
    triang[4*0+xi] = aaSqrt;
    triang[4*1+xi] = aaInvSqrt;
    double r0_re, r0_im;
    if (solveFwd) {
      result[4*0+xi] = r0_re = result[4*0+xi] * aaInvSqrt;
      result[4*1+xi] = r0_im = result[4*1+xi] * aaInvSqrt;
    }

    if (N==1)
      return succ;

    int qlen=(N+3)/4;
    if (!solveFwd) {
      for (int i = 0; i < qlen*8; ++i)
        triang[i] *= aaInvSqrt;
    } else { // combine multiplication by invSqrt with Forward propagation
      for (int i = 0; i < qlen; ++i) {
        for (int k = 0; k < 4; ++k) {
          auto x_re = triang[i*8+4*0+k] * aaInvSqrt;
          auto x_im = triang[i*8+4*1+k] * aaInvSqrt;
          auto r_re = result[i*8+4*0+k];
          auto r_im = result[i*8+4*1+k];
          triang[i*8+4*0+k] = x_re;
          triang[i*8+4*1+k] = x_im;
          result[i*8+4*0+k] = r_re - x_re*r0_re + x_im*r0_im;
          result[i*8+4*1+k] = r_im - x_re*r0_im - x_im*r0_re;
        }
      }
      result[4*0+xi] = r0_re;
      result[4*1+xi] = r0_im;
      result += ((xi + 2) & -4)*2;
    }
    triang[4*0+xi] = aaSqrt;
    triang[4*1+xi] = aaInvSqrt;

    // subtract outer product of top row from lower part of the matrix
    // process two output rows together
    auto x = triang;
    triang += qlen*8;
    auto y = triang;
    xi += 1;
    unsigned hlen = N/2;
    do {
      auto x0 = &x[(xi & -4)*2];
      auto f00_re = x0[4*0+(xi&3)+0];
      auto f00_im = x0[4*1+(xi&3)+0];
      auto f01_re = x0[4*0+(xi&3)+1];
      auto f01_im = x0[4*1+(xi&3)+1];
      int  qlen = (hlen+1)/2;
      auto y0 = &y[0];
      auto y1 = &y[qlen*8];
      for (int c = 0; c < qlen; ++c) {
        // y0[c] -= x0[c]*conj(f00);
        // y1[c] -= x0[c]*conj(f01);
        for (int k = 0; k < 4; ++k) {
          auto x_re = x0[c*8+4*0+k];
          auto x_im = x0[c*8+4*1+k];
          auto y0_re = y0[c*8+4*0+k];
          auto y0_im = y0[c*8+4*1+k];
          auto y1_re = y1[c*8+4*0+k];
          auto y1_im = y1[c*8+4*1+k];
          y0[c*8+4*0+k] = y0_re - x_re * f00_re - x_im * f00_im;
          y0[c*8+4*1+k] = y0_im + x_re * f00_im - x_im * f00_re;
          y1[c*8+4*0+k] = y1_re - x_re * f01_re - x_im * f01_im;
          y1[c*8+4*1+k] = y1_im + x_re * f01_im - x_im * f01_re;
        }
      }
      xi += 2;
      y  += qlen*4*4;
      --hlen;
    } while (hlen > 0);
  }

  for (unsigned rhlen = N/2; ; rhlen-=1) {
    // process top row
    int xi = rhlen & 1;
    int rlen = ((rhlen+1)/2)*4;
    auto x0 = &triang[0];
    auto x1 = &triang[rlen*2];

    double aa0 = x0[4*0+xi*2+0]; // diagonal element
    auto f_re  = x0[4*0+xi*2+1];
    auto f_im  = x0[4*1+xi*2+1];
    // process next after top row
    double aa1 = x1[4*0+xi*2+1]; // diagonal element
    // subtract outer product of top row from next after top row
    aa1 = aa1*aa0 - f_re*f_re - f_im*f_im;
    // check that we are positive definite
    // printf("%u %e %e\n", rhlen, aa0, aa1);
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
    double aa0InvSqrt = aaInvSqrt[0];
    double aa1InvSqrt = aaInvSqrt[1]*aaSqrt[0];
    x0[4*0+xi*2+0] = aaSqrt[0];
    x0[4*1+xi*2+0] = aa0InvSqrt;
    x1[4*0+xi*2+1] = aaSqrt[1]*aa0InvSqrt;
    x1[4*1+xi*2+1] = aa1InvSqrt;
    x0[4*0+xi*2+1] = (f_re *= aa0InvSqrt);
    x0[4*1+xi*2+1] = (f_im *= aa0InvSqrt);

    double r0_re, r0_im, r1_re, r1_im;
    if (solveFwd) {
      // auto r0 = result[0] * aa0InvSqrt;
      // auto r1 = (result[1] - f*r0) * aa1InvSqrt;
      result[4*0+xi*2+0] = r0_re = result[4*0+xi*2+0] * aa0InvSqrt;
      result[4*1+xi*2+0] = r0_im = result[4*1+xi*2+0] * aa0InvSqrt;
      result[4*0+xi*2+1] = r1_re = (result[4*0+xi*2+1] - f_re*r0_re + f_im*r0_im) * aa1InvSqrt;
      result[4*1+xi*2+1] = r1_im = (result[4*1+xi*2+1] - f_re*r0_im - f_im*r0_re) * aa1InvSqrt;
    }

    if (rhlen<=1)
      break; // x1 was last row

    x0 += xi*8;
    x1 += xi*8;
    int cqlen = rhlen / 2;
    for (int c = 0; c < cqlen; ++c) {
      // auto ax0 = x0[c] * aa0InvSqrt;
      // auto ax1 = x1[c];
      // x0[c] = ax0;
      // x1[c] = (ax1 - ax0*conj(f))*aa1InvSqrt;
      for (int k = 0; k < 4; ++k) {
        auto x0_re = x0[c*8+4*0+k] * aa0InvSqrt;
        auto x0_im = x0[c*8+4*1+k] * aa0InvSqrt;
        auto x1_re = x1[c*8+4*0+k];
        auto x1_im = x1[c*8+4*1+k];
        x0[c*8+4*0+k] = x0_re;
        x0[c*8+4*1+k] = x0_im;
        x1[c*8+4*0+k] = (x1_re - x0_re*f_re - x0_im*f_im) * aa1InvSqrt;
        x1[c*8+4*1+k] = (x1_im + x0_re*f_im - x0_im*f_re) * aa1InvSqrt;
      }
    }
    triang[4*0+xi*2+0] = aaSqrt[0];
    triang[4*1+xi*2+0] = aa0InvSqrt;
    triang[4*0+xi*2+1] = f_re;
    triang[4*1+xi*2+1] = f_im;
    triang += rlen*2;
    triang[4*0+xi*2+1] = aaSqrt[1]*aa0InvSqrt;
    triang[4*1+xi*2+1] = aa1InvSqrt;
    triang += rlen*2;
    if (solveFwd) {
      // Forward propagation of pair of results
      auto pr = &result[xi*8];
      for (int c = 0; c < cqlen; ++c) {
        // result[c] -= x0[c]*r0 + x1[c]*r1;
        for (int k = 0; k < 4; ++k) {
          auto x0_re = x0[c*8+4*0+k];
          auto x0_im = x0[c*8+4*1+k];
          auto x1_re = x1[c*8+4*0+k];
          auto x1_im = x1[c*8+4*1+k];
          auto r_re  = pr[c*8+4*0+k];
          auto r_im  = pr[c*8+4*1+k];
          pr[c*8+4*0+k] = r_re - x0_re*r0_re + x0_im*r0_im - x1_re*r1_re + x1_im*r1_im;
          pr[c*8+4*1+k] = r_im - x0_re*r0_im - x0_im*r0_re - x1_re*r1_im - x1_im*r1_re;
        }
      }
      result[4*0+xi*2+0] = r0_re;
      result[4*0+xi*2+1] = r1_re;
      result[4*1+xi*2+0] = r0_im;
      result[4*1+xi*2+1] = r1_im;
      result = pr;
    }

    // subtract outer product of two top rows from lower part of the matrix
    // process two output rows together
    auto y = triang;
    for (unsigned chlen = rhlen-1; chlen > 0; chlen -= 1) {
      int xi = chlen & 1;
      auto f00_re = x0[4*0+xi*2+0];
      auto f10_re = x1[4*0+xi*2+0];
      auto f01_re = x0[4*0+xi*2+1];
      auto f11_re = x1[4*0+xi*2+1];
      auto f00_im = x0[4*1+xi*2+0];
      auto f10_im = x1[4*1+xi*2+0];
      auto f01_im = x0[4*1+xi*2+1];
      auto f11_im = x1[4*1+xi*2+1];
      int cqlen = (chlen+1)/2;
      auto y0 = &y[0];
      auto y1 = &y[cqlen*8];
      for (int c = 0; c < cqlen; ++c) {
        // y0[c] = y0[c] - x0[c]*conj(f00) - x1[c]*conj(f10);
        // y1[c] = y1[c] - x0[c]*conj(f01) - x1[c]*conj(f11);
        for (int k = 0; k < 4; ++k) {
          auto x0_re = x0[c*8+4*0+k];
          auto x0_im = x0[c*8+4*1+k];
          auto y0_re = y0[c*8+4*0+k];
          auto y0_im = y0[c*8+4*1+k];
          auto y1_re = y1[c*8+4*0+k];
          auto y1_im = y1[c*8+4*1+k];
          y0_re = y0_re - x0_re * f00_re - x0_im * f00_im;
          y0_im = y0_im + x0_re * f00_im - x0_im * f00_re;
          y1_re = y1_re - x0_re * f01_re - x0_im * f01_im;
          y1_im = y1_im + x0_re * f01_im - x0_im * f01_re;
          auto x1_re = x1[c*8+4*0+k];
          auto x1_im = x1[c*8+4*1+k];
          y0_re = y0_re - x1_re * f10_re - x1_im * f10_im;
          y0_im = y0_im + x1_re * f10_im - x1_im * f10_re;
          y1_re = y1_re - x1_re * f11_re - x1_im * f11_im;
          y1_im = y1_im + x1_re * f11_im - x1_im * f11_re;
          y0[c*8+4*0+k] = y0_re;
          y0[c*8+4*1+k] = y0_im;
          y1[c*8+4*0+k] = y1_re;
          y1[c*8+4*1+k] = y1_im;
        }
      }
      y += cqlen*4*4;
      x0 += xi * 8;
      x1 += xi * 8;
    }
  }

  return succ;
}

bool chol_Factorize(double* triang, unsigned N)
{
  return chol_FactorizeAndSolveFwd(triang, N, NULL, false);
}

bool chol_Factorize(double* triang, unsigned N, double* result)
{
  return chol_FactorizeAndSolveFwd(triang, N, result, true);
}

void chol_SolveFwd(double *x, unsigned N, const double* triang)
{
  // x = R' \ conj(B);
  // R' is lower triangle - solve by forward propagation (caxpy-like)
  // x = B;
  // for r=1:N
    // x(r) = x(r)/R(r,r);
    // x(r+1:N) -= R(r, r+1:N).'*x(r);
  // end
  if ((N & 1) != 0) {    // special handling for the first row of matrix with odd number of elements
    int xi = ~N & 2;
    auto aaInvSqrt = triang[xi+1+4]; // imag() of diag element contains inverse of its real()
    auto x0 = &x[xi+1];
    auto xr_re = x0[0] * aaInvSqrt;
    auto xr_im = x0[4] * aaInvSqrt;
    int qlen = (N+1)/4;
    if (qlen != 0) {
      x += xi*4;
      auto y = &triang[xi*4];
      for (int c = 0; c < qlen; ++c) {
        // x[c] -= y[c]*xr;
        for (int k = 0; k < 4; ++k) {
          auto x_re = x[c*8+k+0];
          auto x_im = x[c*8+k+4];
          auto y_re = y[c*8+k+0];
          auto y_im = y[c*8+k+4];
          x[c*8+k+0] = x_re - y_re*xr_re + y_im*xr_im;
          x[c*8+k+4] = x_im - y_re*xr_im - y_im*xr_re;
        }
      }
    }
    x0[0] = xr_re;
    x0[4] = xr_im;
    if (qlen == 0)
      return;
    triang += ((N+3)/4)*8;
  }

  // process two rows per iteration
  for (unsigned rlen = N & (-2); ; rlen -= 2) {
    int xi   = rlen & 2;
    int ylen = (rlen + 2) & (-4);
    auto y0 = &triang[xi];
    auto y1 = &y0[ylen*2];
    auto aa0InvSqrt = y0[0+4]; // imag() of diag element contains inverse of its real()
    auto x0 = &x[xi];
    auto xr0_re = x0[0+0] * aa0InvSqrt;
    auto xr0_im = x0[0+4] * aa0InvSqrt;
    auto y01_re = y0[1+0];
    auto y01_im = y0[1+4];
    auto aa1InvSqrt = y1[1+4]; // imag() of diag element contains inverse of its real()
    auto xr1_re = (x0[1+0] - y01_re*xr0_re + y01_im*xr0_im) * aa1InvSqrt;
    auto xr1_im = (x0[1+4] - y01_re*xr0_im - y01_im*xr0_re) * aa1InvSqrt;

    int cqlen = rlen/4;
    if (cqlen != 0) {
      x += xi*4;
      y0 = &triang[xi*4];
      y1 = &y0[ylen*2];
      for (int c = 0; c < cqlen; ++c) {
        // x[c] = x[c] - y0[c]*xr0 - y1[c]*xr1;
        for (int k = 0; k < 4; ++k) {
          auto x_re = x[c*8+k+0];
          auto x_im = x[c*8+k+4];
          auto y0_re = y0[c*8+k+0];
          auto y0_im = y0[c*8+k+4];
          auto y1_re = y1[c*8+k+0];
          auto y1_im = y1[c*8+k+4];
          x[c*8+k+0] = x_re - y0_re*xr0_re + y0_im*xr0_im - y1_re*xr1_re + y1_im*xr1_im;
          x[c*8+k+4] = x_im - y0_re*xr0_im - y0_im*xr0_re - y1_re*xr1_im - y1_im*xr1_re;
        }
      }
    }
    x0[0+0] = xr0_re;
    x0[1+0] = xr1_re;
    x0[0+4] = xr0_im;
    x0[1+4] = xr1_im;
    if (cqlen == 0)
      break;
    triang += ylen*4;
  }
}

void chol_SolveBwd(double *x, unsigned N, const double* triang)
{
  // x = conj(R) \ x;
  // R is upper triangle - solve by backward substitution (dot-like)
  // for r=N:-1:1
   // x(r) = x(r) - sum(x(r+1:N) .* conj(R(r,r+1:N)(:)) );
   // x(r) = x(r)/R(r,r);
  // end
  triang += ((N+2)*(N+2)/8)*8; // point past end
  x += ((N-1)/4)*8;            // point to last bi-qaud

  // process two rows per iteration
  const unsigned hlen = N/2;
  for (unsigned rhlen = 1; rhlen <= hlen; ++rhlen) {
    int ylen = ((rhlen+1)/2)*4;
    triang -= ylen*4; // point to diag element
    int xi = rhlen & 1;
    auto x0 = &x[xi*2];
    auto acc0_re = x0[0+0]; auto acc1_re = x0[1+0];
    auto acc0_im = x0[0+4]; auto acc1_im = x0[1+4];
    x0[0+0]=x0[1+0] = 0;
    x0[0+4]=x0[1+4] = 0;
    x += xi*8;
    auto y0 = &triang[xi*8];
    auto y1 = &y0[ylen*2];
    int cqlen = rhlen/2;
    for (int c = 0; c < cqlen; ++c) {
      // acc0 -= x[c] * conj(y0[c]);
      // acc1 -= x[c] * conj(y1[c]);
      for (int k = 0; k < 4; ++k) {
        acc0_re = acc0_re - x[c*8+k+0]*y0[c*8+k+0] - x[c*8+k+4]*y0[c*8+k+4];
        acc1_re = acc1_re - x[c*8+k+0]*y1[c*8+k+0] - x[c*8+k+4]*y1[c*8+k+4];
        acc0_im = acc0_im + x[c*8+k+0]*y0[c*8+k+4] - x[c*8+k+4]*y0[c*8+k+0];
        acc1_im = acc1_im + x[c*8+k+0]*y1[c*8+k+4] - x[c*8+k+4]*y1[c*8+k+0];
      }
    }
    y0 = &triang[xi*2];
    y1 = &y0[ylen*2];
    auto aa1InvSqrt = y1[1+4]; // imag() of diag element contains inverse of its real()
    acc1_re *= aa1InvSqrt;
    acc1_im *= aa1InvSqrt;
    // acc0 -= acc1*conj(y0[1]);
    acc0_re = acc0_re - acc1_re*y0[1+0] - acc1_im*y0[1+4];
    acc0_im = acc0_im + acc1_re*y0[1+4] - acc1_im*y0[1+0];
    auto aa0InvSqrt = y0[0+4]; // imag() of diag element contains inverse of its real()
    acc0_re *= aa0InvSqrt;
    acc0_im *= aa0InvSqrt;

    x0[0+0] = acc0_re; x0[1+0] = acc1_re;
    x0[0+4] = acc0_im; x0[1+4] = acc1_im;
    x -= 8;
  }

  if ((N & 1) != 0) { // special handling for the first row of matrix with odd number of elements
    triang -= (N/4+1)*8;
    int xi = ~N & 2;
    auto x0 = &x[xi];
    auto acc_re = x0[1+0];
    auto acc_im = x0[1+4];
    x0[0+0]=x0[1+0]=0;
    x0[0+4]=x0[1+4]=0;
    x += xi*4;
    auto y = &triang[xi*4];
    int cqlen = (N+1)/4;
    for (int c = 0; c < cqlen; ++c) {
      // acc -= x[c] * conj(y[c]);
      for (int k = 0; k < 4; ++k) {
        acc_re = acc_re - x[c*8+k+0]*y[c*8+k+0] - x[c*8+k+4]*y[c*8+k+4];
        acc_im = acc_im + x[c*8+k+0]*y[c*8+k+4] - x[c*8+k+4]*y[c*8+k+0];
      }
    }
    auto aaInvSqrt = triang[xi+1+4]; // imag() of diag element contains inverse of its real()
    x0[1+0] = acc_re * aaInvSqrt;
    x0[1+4] = acc_im * aaInvSqrt;
  }
}

};

unsigned chol_getWorkBufferSize(int n)
{
  return ((n+2)*(n+2)/8 + (n+3)/4)*4*2*sizeof(double);
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

  double* packedResult = static_cast<double*>(workBuffer);
  double* triang = packedResult + (unsigned(n+3)/4)*4*2;
  PackUpperTriangle(triang, src, n, srcLayout);
#if CHOL_DBG
  printf("internal src:\n");
  PrintInternalTriangle(triang, n);
#endif
  bool succ = chol_Factorize(triang, n);
  UnpackUpperTriangle(dst, triang, n);

#if CHOL_DBG
  printf("internal:\n");
  PrintInternalTriangle(triang, n);
  printf("res:\n");
  PrintLowerTriangle(dst, n);
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
  printf("src tri:\n");
  PrintLowerTriangle(src, n);
  std::complex<double>* ref_tri = new std::complex<double>[n*n];
  chol_scalar(ref_tri, src, n);
  printf("ref tri:\n");
  PrintLowerTriangle(ref_tri, n);
  std::complex<double>* ref_x = new std::complex<double>[n];
  chol_SolveFwd_scalar(ref_x, n, ref_tri, vecB);
#endif

  double* packedResult = static_cast<double*>(workBuffer);
  double* triang = packedResult + (unsigned(n+3)/4)*4*2;
  PackUpperTriangle(triang, src, n, srcLayout);
  PackVecB(packedResult, vecB, n);
  bool succ = chol_Factorize(triang, n, packedResult);
  if (succ) {
#if CHOL_DBG
    printf("ref x:\n"); PrintResult(ref_x, n);
    UnpackResult(result, packedResult, n);
    printf("res x:\n"); PrintResult(result, n);
    chol_SolveBwd_scalar(ref_x, n, ref_tri);
#endif

    chol_SolveBwd(packedResult, n, triang);
    UnpackResult(result, packedResult, n);

#if CHOL_DBG
    printf("ref result:\n"); PrintResult(ref_x, n);
    printf("result:\n");     PrintResult(result, n);
    delete [] ref_tri;
    delete [] ref_x;
    //return false;
#endif
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

  double* packedResult = static_cast<double*>(workBuffer);
  PackVecB(packedResult, vecB, n);
  const double* triang = packedResult + (unsigned(n+3)/4)*4*2;
  chol_SolveFwd(packedResult, n, triang);
  chol_SolveBwd(packedResult, n, triang);
  UnpackResult(result, packedResult, n);
}

#if CHOL_DBG
// chol - Perform Cholesky decomposition of complex Hermitian matrix
// Decomposition uses combination of Cholesky-Crout and Cholesky-Banachiewicz algorithms without pivoting.
// Inner loop calculates two subsequent elements of the same column
// Parameters:
// src - input matrix stored in row-wise order
//       only lower triangle (including diagonal) of input matrix is used in calculation, values in upper triangle are ignored
// dst - output matrix stored in row-wise order
//       only lower triangle (including diagonal) of output matrix guaranteed to be correct,
//       values in upper triangle are not guaranteed to be zeroed
// n - number of rows/columns in matrix src
// Return value:
// true = success, false - input matrix is not positive definite
bool chol_scalar(std::complex<double> *dst, const std::complex<double> *src, int n)
{
  const int N_MIN = 1;
  const int N_MAX = 128;
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

void chol_SolveFwd_scalar(std::complex<double> *x, unsigned N, const std::complex<double>* mat, const std::complex<double> *vecB)
{
  for (unsigned i = 0; i < N; ++i)
    x[i] = vecB[i];
  for (unsigned r = 0; r < N; ++r) {
    auto x0 = x[r];
    for (unsigned c = 0; c < r; ++c)
      x0 -= x[c]*mat[c];
    x[r] = x0 / mat[r].real();
    mat += N;
  }
}

void chol_SolveBwd_scalar(std::complex<double> *x, unsigned N, const std::complex<double>* mat)
{
  mat += N*(N-1);
  for (unsigned rlen = N; rlen > 0; --rlen) {
    auto x0 = (x[rlen-1] /= mat[rlen-1].real());
    for (unsigned c = 0; c < rlen-1; ++c)
      x[c] -= x0*conj(mat[c]);
    mat -= N;
  }
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

void PrintResult(const std::complex<double> *src, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    printf("(%9f %9f) ", src[i].real(), src[i].imag());
  printf("\n");
}

void PrintInternalTriangle(const double* src, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
  {
    unsigned qlen = (n-i+3)/4;
    for (unsigned j=0; j < qlen; ++j)
      for (unsigned k=0; k < 4; ++k)
        printf("(%9f %9f) ", src[j*8+k], src[j*8+k+4]);
    printf("\n");
    src += qlen*8;
  }
}

#endif
