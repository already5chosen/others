#include <complex>
#include <cfloat>
#include <cstring>
#include <intrin.h>
//#include <stdio.h>
#include "chol.h"

#define CHOL_DBG 0

#if CHOL_DBG
bool chol_scalar(std::complex<double> *dst, const std::complex<double> *src, int n);
void PrintLowerTriangle(const std::complex<double> *src, unsigned n);
void PrintLowerTriangle(std::complex<double> *dst, const __m256d* triang, unsigned n);
void PrintInternalTriangle(const std::complex<double> *src, unsigned n);
#endif

namespace {

inline
void PackUpperTriangle_R(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // transpose - read rows, store column
  for (unsigned row = 0; row < n; ++row) {
    std::complex<double>* pOut = &triang[row];
    for (unsigned col = 0; col < row+1; ++col) {
      *pOut = src[col];
      pOut += n-col-1;
    }
    src += n;
  }
}

inline
void PackUpperTriangle_C(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  for (unsigned row = 0; row < n; ++row) {
    for (unsigned col = row; col < n; ++col)
      *triang++ = src[col];
    src += n;
  }
}

inline
void PackUpperTriangle_P(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  memcpy(triang, src, sizeof(*src)*n*(n+1)/2);
}

inline
void PackUpperTriangle_Q(std::complex<double>* triang, const std::complex<double> *src, unsigned n)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
  // transpose - read rows, store column
  for (unsigned row = 0; row < n; ++row) {
    std::complex<double>* pOut = &triang[row];
    for (unsigned col = 0; col < row+1; ++col) {
      *pOut = *src++;
      pOut += n-col-1;
    }
  }
}

inline
void PackUpperTriangle(std::complex<double>* triang, const std::complex<double> *src, unsigned n, int srcLayout)
{ // pack upper triangle with diagonal
  // since the input is taken from lower triangle, it is a conjugate of what is referenced in classic algorithm
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
  // for (unsigned i = 0; i < n*n; ++i) dst[i]=42;
  for (unsigned row = 0; row < n; ++row) {
    std::complex<double>* pOut = dst;
    *pOut = (*triang++).real(); // zeroize imagery part of diagonal
    pOut += n;
    for (unsigned col = 1; col < n-row; ++col) {
      *pOut = *triang++;
      pOut += n;
    }
    dst += n+1; // to next diaganal element
  }
}

// inline
// void PackVecB(__m256d* dst, const std::complex<double> *src, unsigned n)
// { // organize vecB in processing-friendly interleaved format
  // unsigned wMax = n/2;
  // for (unsigned w = 0; w != wMax; ++w)
  // {
    // __m128d w0 = *reinterpret_cast<const __m128d*>(&src[w*2+0]);
    // __m128d w1 = *reinterpret_cast<const __m128d*>(&src[w*2+1]);
    // reinterpret_cast<__m128d*>(&dst[w])[0] = _mm_unpacklo_pd(w0, w1); // re(w*2+0)  re(w*2+1)
    // reinterpret_cast<__m128d*>(&dst[w])[1] = _mm_unpackhi_pd(w0, w1); // im(w*2+0)  im(w*2+1)
  // }
  // if ((n & 1) != 0)
  // { // last word
    // __m128d w0 = *reinterpret_cast<const __m128d*>(&src[n-1]);
    // reinterpret_cast<__m128d*>(&dst[wMax])[0] = _mm_unpacklo_pd(w0, _mm_setzero_pd()); // re(w*2+0)  0
    // reinterpret_cast<__m128d*>(&dst[wMax])[1] = _mm_unpackhi_pd(w0, _mm_setzero_pd()); // im(w*2+0)  0
  // }
// }

// inline
// void UnpackResult(std::complex<double> *dst_, const __m256d* src, unsigned n)
// { // translate result from internal interleaved format to standard complex<double>
  // unsigned wMax = n/2;
  // __m128d* dst = reinterpret_cast<__m128d*>(dst_);
  // for (unsigned w = 0; w != wMax; ++w)
  // {
    // __m128d re = reinterpret_cast<const __m128d*>(&src[w])[0]; // re(w*2+0)  re(w*2+1)
    // __m128d im = reinterpret_cast<const __m128d*>(&src[w])[1]; // im(w*2+0)  im(w*2+1)
    // dst[w*2+0] = _mm_unpacklo_pd(re, im); // re(w*2+0)  im(w*2+0)
    // dst[w*2+1] = _mm_unpackhi_pd(re, im); // re(w*2+1)  im(w*2+1)
  // }
  // if ((n & 1) != 0)
  // { // last word
    // __m128d re = reinterpret_cast<const __m128d*>(&src[wMax])[0]; // re(n-1) x
    // __m128d im = reinterpret_cast<const __m128d*>(&src[wMax])[1]; // im(n-1) x
    // dst[n-1] = _mm_unpacklo_pd(re, im); // re(n-1)  im(n-1)
  // }
// }

const double DIAG_MIN   =  double(FLT_MIN);
// FLT_MIN is not really special in context of double-precision calculations
// but I wanted very small number that is still much bigger than DBL_MIN
// and FLT_MIN looks good enough
const double DIAG_SUBST =  1E40;

// inline
// int Process0to1(std::complex<double>* triang)
// {
  // // Cholesky-Banachiewicz/Crout decomposition
  // int succ = 1;
  // // c,r=0:1
  // // calculate diagonal element (0,0)
  // double aa = reinterpret_cast<double*>(&triang[0].re)[0]; // current diagonal element

  // // check that we are positive defined
  // //printf("%d %e\n", 0, aa);
  // if (aa < DIAG_MIN)
  // {
    // aa = DIAG_SUBST;
    // succ = 0;
  // }

  // double aaSqrt0    = sqrt(aa);
  // double aaInvSqrt0 = 1.0/aaSqrt0;

  // reinterpret_cast<double*>(&triang[0].re)[0] = aaSqrt0;
  // reinterpret_cast<double*>(&triang[0].im)[0] = aaInvSqrt0;

  // // calculate non-diagonal element (1,0)
  // double re10 = reinterpret_cast<double*>(&triang[0].re)[2] * aaInvSqrt0;
  // double im10 = reinterpret_cast<double*>(&triang[0].im)[2] * aaInvSqrt0;

  // reinterpret_cast<double*>(&triang[0].re)[2] = re10;
  // reinterpret_cast<double*>(&triang[0].im)[2] = im10;

  // // calculate diagonal element (1,1)
  // aa = reinterpret_cast<double*>(&triang[0].re)[3] - re10*re10 - im10*im10; // current diagonal element

  // // check that we are positive defined
  // //printf("%d %e\n", 1, aa);
  // if (aa < DIAG_MIN)
  // {
    // aa = DIAG_SUBST;
    // succ = 0;
  // }

  // double aaSqrt1 = sqrt(aa);
  // double aaInvSqrt1 = 1.0 / aaSqrt1;

  // reinterpret_cast<double*>(&triang[0].re)[3] = aaSqrt1;
  // reinterpret_cast<double*>(&triang[0].im)[3] = aaInvSqrt1;

  // return succ;
// }

// inline
// int ProcessLowerRight(complex_m256d* triang, unsigned n)
// {
  // // diagonal of last line, non-interleaved
  // __m256d re = _mm256_setzero_pd();
  // unsigned dotlen = unsigned(n-1); // always even
  // unsigned nPairs = dotlen/2;
  // __m256d* rowLast = &triang[nPairs*(nPairs+1)/2].re;
  // for (unsigned k = 0; k != nPairs; ++k)
  // {
    // __m256d reim = rowLast[k];
    // re = _mm256_add_pd(re, _mm256_mul_pd(reim, reim));
  // }
  // __m128d re_h = _mm256_extractf128_pd(re, 1);
  // __m128d re_l = _mm256_castpd256_pd128(re);
  // __m128d re128 = _mm_add_pd(re_h, re_l);
  // re128 = _mm_hadd_pd(re128, re128);

  // double aa = reinterpret_cast<double*>(&rowLast[nPairs])[0];
  // aa -= reinterpret_cast<double*>(&re128)[0]; // current diagonal element
  // // check that we are positive defined
  // //printf("%d %e\n", n-1, aa);
  // int succ = 1;
  // if (aa < DIAG_MIN)
  // {
    // aa = DIAG_SUBST;
    // succ = 0;
  // }
  // double aaSqrt = sqrt(aa);
  // reinterpret_cast<double*>(&rowLast[nPairs])[0] = aaSqrt;
  // reinterpret_cast<double*>(&rowLast[nPairs])[2] = 1/aaSqrt;
  // return succ;
// }

// inline
  // void ForwardSubstituteLowerRight(
  // __m256d*              result,
  // const complex_m256d*  triang,
  // unsigned n)
// {
  // // solve L*y=vecB by forward substitution, process last non-interleaved row
  // unsigned rHalf = n/2;
  // const __m256d* row = &triang[rHalf*(rHalf+1)/2].re;
  // __m256d sumRe = _mm256_setzero_pd();
  // __m256d sumIm = _mm256_setzero_pd();
  // for (unsigned k = 0; k != rHalf; ++k)
  // {
    // __m256d reimRes = *result;                                     // reRes(k*2+0)  reRes(k*2+1)  imRes(k*2+0)  imRes(k*2+1)
    // __m256d imreRes = _mm256_permute2f128_pd(reimRes, reimRes, 1); // imRes(k*2+0)  imRes(k*2+1)  reRes(k*2+0)  reRes(k*2+1)
    // __m256d reim = *row;                                           // re(r+0,k*2+0) re(r+0,k*2+1) im(r+0,k*2+0) im(r+0,k*2+1)

    // MSUB(sumRe, reim, reimRes);
    // // 0: 0 -= re(r+0,k*2+0)*reRes(k*2+0) : reRes(r+0)
    // // 1: 0 -= re(r+0,k*2+1)*reRes(k*2+1) : reRes(r+0)
    // // 2: 0 -= im(r+0,k*2+0)*imRes(k*2+0) :-reRes(r+0)
    // // 3: 0 -= im(r+0,k*2+1)*imRes(k*2+1) :-reRes(r+0)
    // MSUB(sumIm, reim, imreRes);
    // // 0: 0 -= re(r+0,k*2+0)*imRes(k*2+0) : imRes(r+0)
    // // 1: 0 -= re(r+0,k*2+1)*imRes(k*2+1) : imRes(r+0)
    // // 2: 0 -= im(r+0,k*2+0)*reRes(k*2+0) : imRes(r+0)
    // // 3: 0 -= im(r+0,k*2+1)*reRes(k*2+1) : imRes(r+0)
    // row    += 1;
    // result += 1;
  // }
  // __m256d sumReim = _mm256_hadd_pd(sumRe, sumIm);         //  reRes(r+0)  imRes(r+0) -reRes(r+0) imRes(r+0)
  // __m128d sum_p = _mm256_castpd256_pd128(sumReim);        //  reRes(r+0)  imRes(r+0)
  // __m128d sum_n = _mm256_extractf128_pd(sumReim, 1);      // -reRes(r+0)  imRes(r+0)
  // __m128d sum = _mm_addsub_pd(sum_p, sum_n);              //  reRes(r+0)  imRes(r+0)
  // __m128d reRes = reinterpret_cast<__m128d*>(result)[0];  //  reRes(r+0,r+0) x
  // __m128d imRes = reinterpret_cast<__m128d*>(result)[1];  //  imRes(r+0,r+0) x
  // __m128d reimRes = _mm_unpacklo_pd(reRes, imRes);        //  reRes(r+0)  imRes(r+0)
  // reimRes = _mm_add_pd(reimRes, sum);                     //  reRes(r+0)  imRes(r+0)
  // reimRes = _mm_mul_pd(reimRes,                           //  reRes(r+0)  imRes(r+0)
    // _mm_loaddup_pd(reinterpret_cast<const double*>(row)+2)); // 1/re(r+0,r+0)

  // __m128d imimRes = _mm_unpackhi_pd(reimRes, reimRes);    //  imRes(r+0)  imRes(r+0)
  // reinterpret_cast<double*>(result)[0] = *reinterpret_cast<double*>(&reimRes);
  // reinterpret_cast<double*>(result)[2] = *reinterpret_cast<double*>(&imimRes);
// }


// inline
  // void BackSubstituteLowerRight(
  // __m256d*              result,
  // const complex_m256d*  triang,
  // unsigned n)
// {
  // // solve L'*x=res by column-wise back substitution:
  // // algorithm structure that is similar to gaussian elimination
  // // This algorithm do more memory accesses than the "natural" substitution,
  // // but it avoids column-wise access to the lower triangular matrix.
  // // Process last non-interleaved row
  // unsigned cHalf = n/2;
  // const __m256d *col = &triang[cHalf*(cHalf+1)/2].re;
  // __m256d reImC = result[cHalf];      // reRes(c) x imRes(c) x
  // __m256d invD = _mm256_broadcast_sd(reinterpret_cast<const double*>(&col[cHalf])+2); // 1/re(c,c)
  // reImC = _mm256_mul_pd(reImC, invD); // reRes(c) x imRes(c) x
  // result[cHalf] = reImC;

  // reImC = _mm256_unpacklo_pd(reImC, reImC);      // reRes(c) reRes(c) imRes(c) imRes(c)
  // __m128d imC = _mm256_extractf128_pd(reImC, 1); // imRes(c) imRes(c)
  // __m128d reC = _mm256_castpd256_pd128(reImC);   // reRes(c) reRes(c)
  // for (unsigned k = 0; k != cHalf; ++k)
  // {
    // __m128d re = reinterpret_cast<const __m128d*>(&col[k])[0]; // re(k*2+0)    re(k*2+1)
    // __m128d im = reinterpret_cast<const __m128d*>(&col[k])[1]; // im(k*2+0)    im(k*2+1)
    // __m128d reRes = reinterpret_cast<__m128d*>(&result[k])[0]; // reRes(k*2+0) reRes(k*2+1)
    // __m128d imRes = reinterpret_cast<__m128d*>(&result[k])[1]; // reRes(k*2+0) reRes(k*2+1)

    // MSUB128(reRes, re, reC); // reRes(k*2+0:k*2+1) -= re(k*2+0:k*2+1)*reRes(c)
    // MSUB128(imRes, re, imC); // imRes(k*2+0:k*2+1) -= re(k*2+0:k*2+1)*imRes(c)
    // MSUB128(reRes, im, imC); // reRes(k*2+0:k*2+1) -= im(k*2+0:k*2+1)*imRes(c)
    // MADD128(imRes, im, reC); // imRes(k*2+0:k*2+1) += im(k*2+0:k*2+1)*reRes(c)

    // reinterpret_cast<__m128d*>(&result[k])[0] = reRes;
    // reinterpret_cast<__m128d*>(&result[k])[1] = imRes;
  // }
// #if 0
  // unsigned cHalf = n/2;
  // const __m256d *col = &triang[cHalf*(cHalf+1)/2].re;
  // double invD = col[cHalf].m256d_f64[2];
  // double re = result[cHalf].m256d_f64[0] * invD;
  // double im = result[cHalf].m256d_f64[2] * invD;
  // result[cHalf].m256d_f64[0] = re;
  // result[cHalf].m256d_f64[2] = im;
  // for (unsigned rHalf = 0; rHalf != cHalf; ++rHalf)
  // {
    // result[rHalf].m256d_f64[0] -= re*col[rHalf].m256d_f64[0] + im*col[rHalf].m256d_f64[2];
    // result[rHalf].m256d_f64[1] -= re*col[rHalf].m256d_f64[1] + im*col[rHalf].m256d_f64[3];
    // result[rHalf].m256d_f64[2] -= im*col[rHalf].m256d_f64[0] - re*col[rHalf].m256d_f64[2];
    // result[rHalf].m256d_f64[3] -= im*col[rHalf].m256d_f64[1] - re*col[rHalf].m256d_f64[3];
  // }
// #endif
// }

bool chol_Factorize(std::complex<double>* triang, unsigned N)
{
  bool succ = true;
  for (unsigned rlen = N; ; --rlen) {
    // process top row
    double aa = (*triang).real(); // diaganal element
    // check that we are positive defined
    // printf("%u %e\n", rlen, aa);
    if (aa < DIAG_MIN)
    {
      aa = DIAG_SUBST;
      succ = false;
    }
    double aaSqrt    = sqrt(aa);
    double aaInvSqrt = 1.0/aaSqrt;
    (*triang).real(aaSqrt);
    (*triang).imag(aaInvSqrt);

    if (rlen<=1)
      break;

    for (unsigned i = 1; i < rlen; ++i)
      triang[i] *= aaInvSqrt;

    // subtract outer product of top row from lower part of the matrix
    auto y = &triang[rlen];
    auto x = &triang[1];
    for (unsigned clen = rlen-1; clen > 0; --clen) {
      auto f = conj(*x);
      y[0].real(y[0].real() - ((*x)*f).real());
      for (unsigned c = 1; c < clen; ++c)
        y[c] -= x[c]*f;
      x += 1;
      y += clen;
    }
    triang += rlen;
  }
  return succ;
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
 for (unsigned rlen = N; ; --rlen) {
   auto xr = x[0] * triang[0].imag(); // imag() of diag element contains inverse of it's real()
   x[0] = xr;
   if (rlen <= 1)
     break;
   for (unsigned c = 1; c < rlen; ++c)
     x[c] -= triang[c]*xr;
   x += 1;
   triang += rlen;
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
 triang += N*(N+1)/2; // point past end
 x += N;              // point past end
 for (unsigned rlen = 1; rlen <= N; ++rlen) {
   triang -= rlen; // point to diag element
   x -= 1;
   auto acc = x[0];
   for (unsigned c = 1; c < rlen; ++c)
     acc -= x[c] * conj(triang[c]);
   x[0] = acc * triang[0].imag(); // imag() of diag element contains inverse of it's real()
 }
}

void chol_Solve(std::complex<double> *result, const std::complex<double> *vecB, unsigned N, const std::complex<double>* triang)
{
  memcpy(result, vecB, sizeof(*vecB)*N);
  chol_SolveFwd(result, N, triang);
  chol_SolveBwd(result, N, triang);
}

};

unsigned chol_getWorkBufferSize(int n)
{
  return n*(n+1)*(sizeof(std::complex<double>)/2);
  // unsigned nh = (n - 1)/2 + 1;
  // unsigned nq = (n - 1)/4 + 1;
  // return (nh*(nh+1)*2+nq*4+1)*sizeof(std::complex<double>);
}

// chol - Perform Cholesky decomposition of complex Hermitian matrix
// Decomposition uses basic Cholesky algorithm (outer product) without pivoting.
// Inner loop processes a single column
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

  std::complex<double>* packedResult = static_cast<std::complex<double>*>(workBuffer);
  std::complex<double>* triang = packedResult;
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
// Inner loop processes a single column
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

  std::complex<double>* packedResult = static_cast<std::complex<double>*>(workBuffer);
  std::complex<double>* triang = packedResult;
  // PackVecB(packedResult, vecB, n);
  PackUpperTriangle(triang, src, n, srcLayout);
  bool succ = chol_Factorize(triang, n);
  if (succ)
    chol_Solve(result, vecB, n, triang);
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

  const std::complex<double>* packedResult = static_cast<const std::complex<double>*>(workBuffer);
  const std::complex<double>* triang = packedResult;
  chol_Solve(result, vecB, n, triang);

  // PackVecB(packedResult, vecB, n);

  // const int MAX_DATASET_SIZE  = 24*1024;
  // const int MAX_DATASET_NELEM = MAX_DATASET_SIZE/(sizeof(complex_m256d));
  // for (unsigned dataset1st = 2; dataset1st < unsigned(n); )
  // {
    // // calculate height of the Crout band
    // int datasetNelem = 0;
    // unsigned datasetLast;
    // for (datasetLast = dataset1st; datasetLast < unsigned(n); datasetLast += 4)
    // {
      // datasetNelem += datasetLast+3; // two dual-rows
      // if (datasetNelem > MAX_DATASET_NELEM)
        // break;
    // }
    // if (datasetLast == dataset1st)
      // datasetLast = dataset1st + 4;
    // if (datasetLast + 4 > unsigned(n))
      // datasetLast = n;

    // chol_CrBa4_i2_Band_ForwardSubstitute(packedResult, triang,
                           // dataset1st == 2 ? 0 : dataset1st, datasetLast);
    // dataset1st = datasetLast;
  // }

  // if ((n & 1) == 1)
  // {
    // ForwardSubstituteLowerRight(packedResult, triang, n);
    // BackSubstituteLowerRight(packedResult, triang, n);
  // }

  // chol_CrBa4_i2_BackSubstitute(packedResult, triang, n);
  // UnpackResult(result, packedResult, n);
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
