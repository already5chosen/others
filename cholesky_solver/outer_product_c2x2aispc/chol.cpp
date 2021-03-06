#include <complex>
#include <cfloat>
#include <cstring>
#include <intrin.h>
//#include <stdio.h>
#include "chol.h"
#include "chol_FactorizeAndSolveFwd.h"
#include "chol_SolveBwd.h"
#include "chol_SolveFwd.h"

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
    dst += n+1; // to next diaganal element
  }
}

bool chol_Factorize(std::complex<double>* triang, unsigned N)
{
  // return chol_FactorizeAndSolveFwd(triang, N, NULL, false);
  return ispc::chol_Factorize(reinterpret_cast<ispc::DComplex*>(triang), N);
}

bool chol_Factorize(std::complex<double>* triang, unsigned N, std::complex<double>* result)
{
#if 0
  bool succ = ispc::chol_Factorize(reinterpret_cast<double*>(triang), N);
  ispc::chol_SolveFwd(reinterpret_cast<double*>(result), N, reinterpret_cast<const double*>(triang));
  return succ;
#else
  return ispc::chol_FactorizeAndSolveFwd(
    reinterpret_cast<ispc::DComplex*>(triang),
    N,
    reinterpret_cast<ispc::DComplex*>(result));
#endif
}

void chol_SolveFwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
{
  #if 1
  ispc::chol_SolveFwd(reinterpret_cast<ispc::DComplex*>(x), N, reinterpret_cast<const ispc::DComplex*>(triang));
  #else
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
    for (int c = 0; c < int(rhlen-1)*2; ++c)
      x[c] = x[c] - y0[c]*xr0 - y1[c]*xr1;

    triang += rhlen*4;
  }
  #endif
}

void chol_SolveBwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
{
  #if 1
  ispc::chol_SolveBwd(reinterpret_cast<ispc::DComplex*>(x), N, reinterpret_cast<const ispc::DComplex*>(triang));
  #else
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
    auto acc0 = x[-2];
    auto acc1 = x[-1];
    for (int c = 0; c < int(rhlen-1)*2; ++c) {
      acc0 -= x[c] * conj(y0[2+c]);
      acc1 -= x[c] * conj(y1[2+c]);
    }
    acc1 *= y1[1].imag(); // imag() of diag element contains inverse of it's real()
    acc0 -= acc1*conj(y0[1]);
    acc0 *= y0[0].imag();
    x[-2] = acc0;
    x[-1] = acc1;
    x -= 2;
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
  #endif
}
//void chol_SolveBwd(std::complex<double> *x, unsigned N, const std::complex<double>* triang)
//{
//  std::complex<double> *x1 = new std::complex<double>[N];
//  memcpy(x1, x, sizeof(*x1)*N);
//  ispc::chol_SolveBwd(reinterpret_cast<double*>(x), N, reinterpret_cast<const double*>(triang));
//  chol_SolveBwd0(x1, N, triang);
//  for (int i = 0; i < N; ++i)
//    printf("%d (%18f %18f) (%18f %18f)\n", i, x[i].real(), x[i].imag(), x1[i].real(), x1[i].imag());
//  delete [] x1;
//}

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
