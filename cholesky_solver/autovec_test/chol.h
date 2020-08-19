#include <complex>

unsigned chol_getWorkBufferSize(int n);

// chol - Perform Cholesky decomposition of complex Hermitian matrix
// Decomposition uses combination of Cholesky-Crout and Cholesky-Banachiewicz algorithms without pivoting.
// Inner loop calculates two subsequent elements of the same column
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
bool chol(std::complex<double> *dst, const std::complex<double> *src, int n, void* workBuffer, int srcLayout);

// chol_solver - solve complex Hermitian matrix by means of Cholesky decomposition
// Decomposition uses combination of Cholesky-Crout and Cholesky-Banachiewicz algorithms without pivoting.
// Inner loop calculates two subsequent elements of the same column
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
bool chol_solver(std::complex<double> *result, const std::complex<double> *src, const std::complex<double> *vecB, int n, void* workBuffer, int srcLayout);

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
void chol_trif_solver(std::complex<double> *result, const std::complex<double> *vecB, int n, void* workBuffer);
