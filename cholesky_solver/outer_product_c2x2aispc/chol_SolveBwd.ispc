#include "isps_complex.ispc"

export void chol_SolveBwd(
  uniform       DComplex x[],
  uniform       int      N,
  uniform const DComplex triang[])
{
  // x = conj(R) \ x;
  // R is upper triangle - solve by backward substitution (dot-like)
  // for r=N:-1:1
   // x(r) = x(r) - sum(x(r+1:N) .* conj(R(r,r+1:N)(:)) );
   // x(r) = x(r)/R(r,r);
  // end

  triang += ((N+1)*(N+1)/2); // point past end
  x += N;                    // point past end

  // process two rows per iteration
  uniform const int hlen = N/2;
  for (uniform unsigned int rhlen = 1; rhlen <= hlen; ++rhlen) {
    triang -= rhlen*4; // point to diag element
    uniform const DComplex* uniform y0 = &triang[0];
    uniform const DComplex* uniform y1 = &triang[rhlen*2];
    double zeros = 0;
    double acc0_re = insert(zeros, 0, x[-2].re);
    double acc0_im = insert(zeros, 0, x[-2].im);
    double acc1_re = insert(zeros, 0, x[-1].re);
    double acc1_im = insert(zeros, 0, x[-1].im);
    y0 += 2;
    y1 += 2;
    uniform double* uniform dx = &x[0].re;
    uniform const double* uniform dy0 = &y0[0].re;
    uniform const double* uniform dy1 = &y1[0].re;
    foreach (c = 0 ... (rhlen-1)*2*2) {
     // acc0 -= x[c] * conj(y0[c]);
     // acc1 -= x[c] * conj(y1[c]);
     double vx  = dx[c];
     double vy0 = dy0[c];
     double vy1 = dy1[c];
     double vxx = shuffle(vx, programIndex ^1);
     acc0_re -= vx *vy0;
     acc0_im -= vxx*vy0;
     acc1_re -= vx *vy1;
     acc1_im -= vxx*vy1;
    }
    y0 -= 2;
    y1 -= 2;

    // reduction
    uniform DComplex x1 = { reduce_add(acc1_re), 0 };
    uniform DComplex x0 = { reduce_add(acc0_re), 0 };
    for (uniform int k = 0; k < programCount; k += 2) {
      x0.im += extract(acc0_im, k+0);
      x0.im -= extract(acc0_im, k+1);
      x1.im += extract(acc1_im, k+0);
      x1.im -= extract(acc1_im, k+1);
    }

    x1 = x1 * y1[1].im;                  // .im of diag element contains inverse of it's real()
    x0 = (x0 - x1*conj(y0[1]))* y0[0].im;

    x -= 2;
    x[0] = x0;
    x[1] = x1;
  }

  if ((N & 1) != 0) { // special handling for the first row of matrix with odd number of elements
    triang -= hlen*2;
    double zeros = 0;
    double acc_re = insert(zeros, 0, x[-1].re);
    double acc_im = insert(zeros, 0, x[-1].im);
    uniform double* uniform dx = &x[0].re;
    uniform const double* uniform dy = &triang[0].re;
    foreach (c = 0 ... hlen*2*2) {
      // acc -= x[c] * conj(triang[c]);
      double vx  = dx[c];
      double vy  = dy[c];
      double vxx = shuffle(vx, programIndex ^1);
      acc_re -= vx*vy;
      acc_im -= vxx*vy;
    }
    uniform DComplex xx = {reduce_add(acc_re), 0};
    for (uniform int k = 0; k < programCount; k += 2) {
      xx.im += extract(acc_im, k+0);
      xx.im -= extract(acc_im, k+1);
    }
    x[-1] = xx * triang[-1].im; // .imag() of diag element contains inverse of it's real()
  }
}

