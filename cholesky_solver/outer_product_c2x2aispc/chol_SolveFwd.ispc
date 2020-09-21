#include "isps_complex.ispc"

inline double minplus() { // return varying double containing {-1,+1,-1+1, ... }
  double ret;
  const uniform double<2> mp = {-1, 1};
  for (uniform int k=0; k < programCount; ++k)
    ret = insert(ret, k, mp[k %2]);
  return ret;
}

export void chol_SolveFwd(
  uniform       DComplex     x[], // complex
  uniform       unsigned int N,
  uniform const DComplex     triang[]) // complex
{
  // x = R' \ conj(B);
  // R' is lower triangle - solve by forward propagation (caxpy-like)
  // x = B;
  // for r=1:N
    // x(r) = x(r)/R(r,r);
    // x(r+1:N) -= R(r, r+1:N).'*x(r);
  // end
  if ((N & 1) != 0) {    // special handling for the first row of matrix with odd number of elements
    const DComplex* uniform y = &triang[1]; // point to diag element
    uniform DComplex xr = x[0] * y[0].im;
    x[0] = xr;
    if (N <= 1)
      return;
    x += 1;
    y += 1;
    uniform unsigned int hlen = N / 2;
    const double xr_imv = xr.im*minplus();
    uniform double* uniform dx = &x[0].re;
    const uniform double* uniform dy = &y[0].re;
    foreach (c = 0 ... (int)hlen*2*2) {
      double vx = dx[c];
      double vy = dy[c];
      double vyx = shuffle(vy, programIndex ^1); // swap re/im
      vx = vx - vy*xr.re - vyx*xr_imv;
      dx[c] = vx;
    }
    triang += N+1;
  }

  // process two rows per iteration
  for (uniform int rhlen = N/2; rhlen > 0; --rhlen) {
    const DComplex* uniform y0 = &triang[0];
    const DComplex* uniform y1 = &triang[rhlen*2];

    uniform double invY00 = y0[0].im; // .im of diag element contains inverse of real() part
    uniform double invY11 = y1[1].im; // odd row is padded
    
    uniform DComplex xr0 = x[0] * invY00;
    // xr1 = (x[1] - y0[1]*xr0) * invY11;
    uniform DComplex xr1 = (x[1] - y0[1]*xr0)*invY11;

    x[0] = xr0;
    x[1] = xr1;
    if (rhlen > 1) {
      y0 += 2;
      y1 += 2;
      x  += 2;
      const double xr0_imv = xr0.im*minplus();
      const double xr1_imv = xr1.im*minplus();
      uniform double* uniform dx = &x[0].re;
      const uniform double* uniform dy0 = &y0[0].re;
      const uniform double* uniform dy1 = &y1[0].re;
      foreach (c = 0 ... (int)(rhlen-1)*2*2) {
        // x[c] = x[c] - y0[c]*xr0 - y1[c]*xr1;
        double vx  = dx[c];
        double vy0 = dy0[c];
        double vy1 = dy1[c];
        double vy0x = shuffle(vy0, programIndex ^1); // swap re/im
        double vy1x = shuffle(vy1, programIndex ^1); // swap re/im
        vx = vx - vy0*xr0.re - vy0x*xr0_imv;
        vx = vx - vy1*xr1.re - vy1x*xr1_imv;
        dx[c] = vx;
      }    
      triang += rhlen*4;
    }
  }
}

