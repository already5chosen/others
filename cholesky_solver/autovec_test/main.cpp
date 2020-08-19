#include <cstdio>
#include <cstring>
#include <complex>
#include <random>
#include <algorithm>
#include <vector>
#include <chrono>
#include <thread>

#include "chol.h"

const int MIN_N = 5;
const int MAX_N = 128;

static void GenerateHermitianMatrix(std::complex<double> *mat, int n, std::mt19937_64& gen);
static void GenerateRandomVector(std::complex<double> *vec, int n, std::mt19937_64& gen);
static void Transpose(std::complex<double>*dst, const std::complex<double>* src, int n);
static void ExtractLowerTriangleColumnwise(std::complex<double>*dst, const std::complex<double>* src, int n);
static void ExtractLowerTriangleRowwise(std::complex<double>*dst, const std::complex<double>* src, int n);
static char* alignUp(char* p, unsigned al);

static double validate_chol(
  const std::complex<double> *matA,
  const std::complex<double> *matL,
  int n);
static double validate_lineq_solution(
  const std::complex<double> *matA,
  const std::complex<double> *vecB,
  const std::complex<double> *vecX,
  int n);

int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "Usage:\n"
      "chol matrix-size [-r | -c | -p ] [-t]\n"
      " -r  - row-wise layout of input matrix, as standard C/C++. It is default\n"
      " -c  - column-wise layout of input matrix, like in Fortran/Matlab\n"
      " -p  - packed column-wise layout of input matrix, lower triangle only\n"
      " -q  - packed row-wise layout of input matrix, lower triangle only\n"
      " -t  - profile combination of chol_solver() followed by chol_trif_solver()\n"
      );
    return 1;
  }

  char* nstr = argv[1];
  int n = strtol(nstr, NULL, 0);
  if (n < MIN_N || n > MAX_N)
  {
    fprintf(stderr
      ,"matrix-size parameter '%s' out of range. Please specify number in range [%d..%d]\n"
      ,nstr, MIN_N, MAX_N);
    return 1;
  }

  int matrixFormat = 'R';
  bool profileTrifSolver = false;
  for (int argI= 2; argI < argz; ++argI)
  {
    if      (0==_stricmp(argv[argI], "-R")) matrixFormat = 'R';
    else if (0==_stricmp(argv[argI], "-C")) matrixFormat = 'C';
    else if (0==_stricmp(argv[argI], "-P")) matrixFormat = 'P';
    else if (0==_stricmp(argv[argI], "-Q")) matrixFormat = 'Q';
    else if (0==_stricmp(argv[argI], "-T")) profileTrifSolver = true;
    else
    {
      fprintf(stderr, "Unknown option %s\n", argv[argI]);
      return 1;
    }
  }

  const int MAX_TEST_BUF_SIZE = 2000000;
  const int matSz             = n * n;
  const int N_TEST_BUFFERS    = MAX_TEST_BUF_SIZE / matSz;


  std::vector<char> workBufferV(chol_getWorkBufferSize(n)+128);
  void* workBuffer = (void*)alignUp(workBufferV.data(), 128);

  std::mt19937_64 gen;
  gen.seed(1);

  std::vector<std::complex<double>> matA(matSz*N_TEST_BUFFERS);
  for (int i = 0; i < N_TEST_BUFFERS; ++i)
    GenerateHermitianMatrix(&matA.at(i*matSz), n, gen);

  std::vector<std::complex<double>> vecB(n*N_TEST_BUFFERS);
  for (int i = 0; i < N_TEST_BUFFERS; ++i)
    GenerateRandomVector(&vecB.at(i*n), n, gen);

  std::vector<std::complex<double>> vecX(n*N_TEST_BUFFERS);
  std::vector<std::complex<double>> matL(n*n);
  std::vector<std::complex<double>> matT(n*n); // for 'C' and 'P' input formats
  double decompMaxErr = 0;
  double solverMaxErr = 0;
  double solver1MaxErr= 0;
  int errIt = -1;
  for (int i = 0; i < N_TEST_BUFFERS; ++i)
  {
    const std::complex<double>* pInpMat = &matA.at(i*matSz);
    switch (matrixFormat)
    {
    case 'C':
      Transpose(&matT.at(0), &matA.at(i*matSz), n);
      pInpMat = &matT.at(0);
      break;

    case 'P':
      ExtractLowerTriangleColumnwise(&matT.at(0), &matA.at(i*matSz), n);
      pInpMat = &matT.at(0);
      break;

    case 'Q':
      ExtractLowerTriangleRowwise(&matT.at(0), &matA.at(i*matSz), n);
      pInpMat = &matT.at(0);
      break;

    case 'R':
    default:
      break;
    }

    if (!chol(&matL.at(0), pInpMat, n, workBuffer, matrixFormat))
    {
      printf("N=%4d. Decomposition #%d failed\n", n, i);
      return 1;
    }
    if (!chol_solver(&vecX.at(i*n), pInpMat, &vecB.at(i*n), n, workBuffer, matrixFormat))
    {
      printf("N=%4d. Decomposition failed in solver\n", n);
      return 1;
    }
    int i1 = (i + 3) % N_TEST_BUFFERS;
    chol_trif_solver(&vecX.at(i1*n), &vecB.at(i1*n), n, workBuffer);

    decompMaxErr = std::max(decompMaxErr, validate_chol(&matA.at(i*matSz), &matL.at(0), n));
    solverMaxErr = std::max(solverMaxErr, validate_lineq_solution(&matA.at(i*matSz), &vecB.at(i*n), &vecX.at(i*n), n));
    solver1MaxErr= std::max(solver1MaxErr,validate_lineq_solution(&matA.at(i*matSz), &vecB.at(i1*n),&vecX.at(i1*n),n));
    if (decompMaxErr > 1e-8 || solverMaxErr > 1e-8 || solver1MaxErr > 1e-8)
    {
      errIt = i;
      break;
    }
  }
  printf("Layout='%c'. N=%4d. max. err: decomposition %.3e, solver %.3e. ", matrixFormat, n, decompMaxErr, std::max(solverMaxErr, solver1MaxErr));
  if (errIt >= 0)
  {
    printf(
      "Failed at test buffer #%d.\n"
      "max. err: decomposition %.3e, solver %.3e. 2nd solver %.3e.\n",
      errIt,
      decompMaxErr, solverMaxErr, solver1MaxErr);
    return -1;
  }

  unsigned inpMatSz = matSz;
  switch (matrixFormat)
  {
  case 'C':
    for (int i = 0; i < N_TEST_BUFFERS; ++i)
    {
      Transpose(&matT.at(0), &matA.at(i*matSz), n);
      std::copy(matT.begin()+0, matT.begin()+matSz, matA.begin()+matSz*i);
    }
    break;

  case 'P':
    inpMatSz = n * (n+1) / 2;
    for (int i = 0; i < N_TEST_BUFFERS; ++i)
    {
      ExtractLowerTriangleColumnwise(&matT.at(0), &matA.at(i*matSz), n);
      std::copy(matT.begin()+0, matT.begin()+inpMatSz, matA.begin()+inpMatSz*i);
    }
    break;

  case 'Q':
    inpMatSz = n * (n+1) / 2;
    for (int i = 0; i < N_TEST_BUFFERS; ++i)
    {
      ExtractLowerTriangleRowwise(&matT.at(0), &matA.at(i*matSz), n);
      std::copy(matT.begin()+0, matT.begin()+inpMatSz, matA.begin()+inpMatSz*i);
    }
    break;

  case 'R':
  default:
    break;
  }

  const int N_REPS = 117;
  std::chrono::duration<int64_t, std::nano> dt[N_REPS];
  std::vector<std::complex<double>> vecX1(n);
  std::this_thread::sleep_for(std::chrono::milliseconds(20));
  for (int rep = 0; rep < N_REPS; ++rep)
  {
    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < N_TEST_BUFFERS; ++it)
    {
      if (!chol_solver(&vecX.at(it*n), &matA.at(it*inpMatSz), &vecB.at(it*n), n, workBuffer, matrixFormat))
      {
        printf("N=%4d. Decomposition failed in solver. rep %d, it %d.\n", n, rep, it);
        return 1;
      }
      if (profileTrifSolver)
        chol_trif_solver(&vecX1.at(0), &vecX.at(it*n), n, workBuffer);
    }
    auto t1 = std::chrono::steady_clock::now();
    dt[rep] = t1 - t0;
  }

  std::nth_element(&dt[0], &dt[N_REPS/2], &dt[N_REPS]); // median
  auto ddt = std::chrono::duration_cast<std::chrono::duration<double>>(dt[N_REPS/2]);
  double tmPerTest = ddt.count()/N_TEST_BUFFERS;
  double deconvFlops = n*(n+1)*(n+2)/6*4.0; // *4 because matrix is complex
  double substFlops  = n*(n+1)*4.0;         // *4 because matrix is complex
  double flops = deconvFlops+substFlops;
  if (profileTrifSolver)
    flops += substFlops;
  printf("N=%4d. %.1f usec. %.3f GMADD/s\n", n, tmPerTest*1E6, flops/(tmPerTest*1E9));

#if 0
  unsigned wsz = chol_getWorkBufferSize(n);
  for (unsigned i=0; i < wsz; ++i)
    if (static_cast<char*>(workBuffer)[wsz-1-i] != '.')
    {
      printf("workBuffer size=%u bytes. %u unused.\n", wsz, i);
      break;
    }
#endif

  return 0;
}

static char* alignUp(char* p, unsigned al)
{
  // align work buffer on 2**N boundary
  // the code below is non-portable, but o.k. for a test
  uintptr_t offs = reinterpret_cast<uintptr_t>(p);
  return p + ((0 - offs) & (al-1));
}

static void GenerateRandomVector(std::complex<double> *vec, int n, std::mt19937_64& gen)
{
  std::uniform_real_distribution<double> uniDistr(-1,1);
  for (int i = 0; i < n; ++i)
  {
    vec[i].real(uniDistr(gen));
    vec[i].imag(uniDistr(gen));
  }
}

static void GenerateHermitianMatrix(std::complex<double> *mat, int n, std::mt19937_64& gen)
{
  // generate pseudorandom matrix r
  int rlen = n + 1;
  int rsize = rlen * n;
  std::complex<double>* r = new std::complex<double>[rsize];
  std::normal_distribution<double> normalDistr;

  for (int i = 0; i < rsize; ++i)
  {
    r[i].real(normalDistr(gen));
    r[i].imag(normalDistr(gen));
  }
  //for (int i = 0; i < rsize; ++i)
  //  printf("%15.7f %15.7f\n", r[i].real(), r[i].imag());

  for (int row = 0; row < n; ++row)
  {
    std::complex<double>* a = &r[row*rlen];
    for (int col = row; col < n; ++col)
    {
      std::complex<double>* b = &r[col*rlen];
      std::complex<double> sum = 0;
      for (int i = 0; i < rlen; ++i)
        sum += a[i] * conj(b[i]);
      if (col == row)
        sum.imag(0);
      mat[row*n+col] = conj(sum);
      mat[col*n+row] = sum;
    }
  }
  delete [] r;
#if 0
  for (int row = 0; row < n; ++row)
  {
    for (int col = 0; col < n; ++col)
      printf("(%8.4f %8.4f)", mat[row*n+col].real(), mat[row*n+col].imag());
    printf("\n");
  }
#endif
}

static void Transpose(std::complex<double>* dst, const std::complex<double>* src, int n)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      *dst++ = src[j*n+i];
}

static void ExtractLowerTriangleColumnwise(std::complex<double>*dst, const std::complex<double>* src, int n)
{
  for (int i = 0; i < n; ++i)
    for (int j = i; j < n; ++j)
      *dst++ = src[j*n+i];
}

static void ExtractLowerTriangleRowwise(std::complex<double>*dst, const std::complex<double>* src, int n)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j)
      *dst++ = src[i*n+j];
}

// cMat - Lower triangle matrix
static double validate_chol(const std::complex<double> *aMat, const std::complex<double> *cMat, int n)
{
  double maxErr = 0;
  for (int row = 0; row < n; ++row)
  {
    const std::complex<double>* a = &cMat[row*n];
    for (int col = 0; col <= row; ++col)
    {
      const std::complex<double>* b = &cMat[col*n];
      // multiply rows
      std::complex<double> sum = 0;
      for (int i = 0; i <= col; ++i)
        sum += a[i] * conj(b[i]);
      sum -= aMat[row*n+col];
      double err = (sum * conj(sum)).real();
      if (err > maxErr)
        maxErr = err;
    }
  }
  return sqrt(maxErr);
}

static double validate_lineq_solution(
  const std::complex<double> *matA,
  const std::complex<double> *vecB,
  const std::complex<double> *vecX,
  int n)
{
  double maxErr = 0;
  const std::complex<double>* row = matA;
  for (int r = 0; r < n; ++r)
  {
    std::complex<double> sum = vecB[r];
    for (int c = 0; c < n; ++c)
      sum -= vecX[c]*row[c]; // multiply row by result
    double err = (sum * conj(sum)).real();
    if (err > maxErr)
      maxErr = err;
    row += n;
  }
  return sqrt(maxErr);
}
