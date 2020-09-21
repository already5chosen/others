#include "isps_complex.ispc"

inline double minplus() {
  double ret;
  const uniform double<2> mp = {-1, 1};
  for (uniform int k=0; k < programCount; ++k)
    ret = insert(ret, k, mp[k %2]);
  return ret;
}

inline void caxpy1x2tri(uniform DComplex triang[], uniform int hlen) // hlen > 0
{
  uniform       DComplex* uniform dst = &triang[hlen*2+2];
  uniform const DComplex* uniform src = &triang[2];
  do {
    // process 2 rows at time
    uniform const double f0re = src[0].re;
    const double f0im = src[0].im*minplus();
    const double f1re = src[1].re;
    const double f1im = src[1].im*minplus();
    const double* uniform px = &src[0].re;
    uniform double* uniform py0 = &dst[0].re;
    uniform double* uniform py1 = &dst[hlen*2].re;
    foreach (index = 0 ... hlen*4) {
      double x  = px[index];
      double y0 = py0[index];
      double y1 = py1[index];
      double xx = shuffle(x, programIndex ^1);
      y0 = y0 - x*f0re + xx*f0im;
      y1 = y1 - x*f1re + xx*f1im;
      py0[index] = y0;
      py1[index] = y1;
    }
    src += 2;
    dst += hlen*4; // 2 * 2 rows
    hlen -= 1;
  } while (hlen > 0);
}

inline void caxpy2x2tri(uniform DComplex triang[], uniform int xhlen) // xhlen > 1
{
  uniform int xlen = xhlen*2;
  uniform const DComplex* uniform px  = &triang[2];
  uniform       DComplex* uniform dst = &triang[xlen*2];
  uniform int hlen = xhlen - 1;
  do {
    // process 2 rows at time
    const double f00re = px[0].re;
    const double f00im = px[0].im*minplus();
    const double f01re = px[1].re;
    const double f01im = px[1].im*minplus();
    const double f10re = px[xlen+0].re;
    const double f10im = px[xlen+0].im*minplus();
    const double f11re = px[xlen+1].re;
    const double f11im = px[xlen+1].im*minplus();
    uniform const double* uniform dpx0 = &px[0].re;
    uniform const double* uniform dpx1 = &px[xlen].re;
    uniform double* uniform py0 = &dst[0].re;
    uniform double* uniform py1 = &dst[hlen*2].re;
    foreach (index = 0 ... hlen*4) {
      double x0 = dpx0[index];
      double y0 = py0[index];
      double y1 = py1[index];
      double xx0= shuffle(x0, programIndex ^1);
      y0 = y0 - x0*f00re + xx0*f00im;
      y1 = y1 - x0*f01re + xx0*f01im;
      double x1 = dpx1[index];
      double xx1= shuffle(x1, programIndex ^1);
      y0 = y0 - x1*f10re + xx1*f10im;
      y1 = y1 - x1*f11re + xx1*f11im;
      py0[index] = y0;
      py1[index] = y1;
    }
    px += 2;
    dst += hlen*4; // 2 * 2 rows
    hlen -= 1;
  } while (hlen > 0);
}

inline uniform bool chol_FactorizeAndSolveFwdEx(
  uniform DComplex     triang[],
  uniform unsigned int N,
  uniform DComplex     result[],
  uniform const bool   solve)
{
  uniform const double DIAG_MIN   = 0x1p-128;
  uniform const double DIAG_SUBST = 0x1p+130;
  uniform bool succ = true;
  if ((N & 1) != 0) { // special handling for the first row of matrix with odd number of elements
    // process top row
    uniform DComplex* uniform x0 = &triang[1];
    uniform double  aa = x0[0].re; // diagonal element
    // check that we are positive definite
    if (aa < DIAG_MIN) {
      aa = DIAG_SUBST;
      succ = false;
    }
    uniform double aaSqrt    = sqrt(aa);
    uniform double aaInvSqrt = 1.0/aaSqrt;
    x0[0].re = aaSqrt;
    x0[0].im = aaInvSqrt;
    if (solve)
      result[0] = result[0] * aaInvSqrt;

    uniform unsigned int hlen=N/2;
    if (hlen==0)
      return succ; // x0 was last (and only) row

    x0 += 1;
    if (solve) {
      // combine multiplication by invSqrt with Forward propagation
      uniform DComplex r = result[0];
      result += 1;
      const double r_imv = r.im*minplus();
      uniform double* uniform dresult = &result[0].re;
      uniform double* uniform dx0     = &x0[0].re;
      foreach (i = 0 ... hlen*2*2) {
        // result[i] -= x0val*r;
        double xv = dx0[i] * aaInvSqrt;
        double rv = dresult[i];
        double xvx = shuffle(xv, programIndex ^1); // swap re and im
        rv = rv - xv * r.re - xvx * r_imv;
        dx0[i]     = xv;
        dresult[i] = rv;
      }
    } else {
      // multiply top row by invSqrt
      uniform double* uniform dx0 = &x0[0].re;
      foreach (i = 0 ... hlen*2*2) {
        dx0[i] = dx0[i] * aaInvSqrt;
      }
    }

    // subtract outer product of top row from lower part of the matrix
    // process two output rows together
    caxpy1x2tri(triang, hlen);
    triang += N+1;
  }

  for (uniform unsigned int rhlen = N/2; ; --rhlen) {
    // process top row
    uniform DComplex* uniform x0 = &triang[0];
    uniform DComplex* uniform x1 = &triang[rhlen*2];

    uniform double aa0 = x0[0].re; // diagonal element
    uniform DComplex f = x0[1];
    // process next after top row
    // subtract outer product of top row from next after top row
    uniform double aa1 = x1[1].re*aa0 - norm(f); // diagonal element
    // check that we are positive definite
    uniform double<2> aa = { aa0, aa1 }, aaSqrt, aaInvSqrt;
    foreach (k=0 ... 2) {
      double v = aa[k];
      if (v < DIAG_MIN) {
        v = DIAG_SUBST;
        succ = false;
      }
      aaSqrt[k]    = v = sqrt(v);
      aaInvSqrt[k] = 1.0/v;
    }
    uniform double aa0InvSqrt = aaInvSqrt[0];
    uniform double aa1InvSqrt = aaInvSqrt[1]*aaSqrt[0];
    x0[0].re = aaSqrt[0];
    x0[0].im = aa0InvSqrt;
    x1[1].re = aaSqrt[1]*aa0InvSqrt;
    x1[1].im = aa1InvSqrt;
    f = f * aa0InvSqrt;
    x0[1] = f;

    uniform DComplex r0, r1;
    if (solve) {
      result[0] = r0 = result[0] * aa0InvSqrt;
      result[1] = r1 = (result[1] - f*r0) * aa1InvSqrt;
      result += 2;
    }

    if (rhlen<=1)
      break; // x1 was last row

    x0 += 2;
    x1 += 2;

    if (solve) {
      // Combine handling of pair of top rows with forward propagation of pair of results
      const double f_imv  = f.im  * minplus();
      const double r0_imv = r0.im * minplus();
      const double r1_imv = r1.im * minplus();
      uniform double* uniform dx0 = &x0[0].re;
      uniform double* uniform dx1 = &x1[0].re;
      uniform double* uniform dresult = &result[0].re;
      foreach (c = 0 ... (int)(rhlen-1)*2*2) {
        double vx0 = dx0[c];
        double vx1 = dx1[c];
        double vr  = dresult[c];
        vx0 *= aa0InvSqrt;
        double vx0x = shuffle(vx0, programIndex ^ 1); // swap re and im
        vx1 = vx1 - vx0*f.re  + vx0x * f_imv;
        vr  = vr  - vx0*r0.re - vx0x * r0_imv;
        vx1 *= aa1InvSqrt;
        double vx1x = shuffle(vx1, programIndex ^ 1); // swap re and im
        vr  = vr  - vx1*r1.re - vx1x * r1_imv;
        dx0[c]     = vx0;
        dx1[c]     = vx1;
        dresult[c] = vr;
      }
    } else {
      // Handle pair of top rows
      const double f_imv = f.im*minplus();
      uniform double* uniform dx0 = &x0[0].re;
      uniform double* uniform dx1 = &x1[0].re;
      foreach (c = 0 ... (int)(rhlen-1)*2*2) {
        double vx0 = dx0[c];
        double vx1 = dx1[c];
        vx0 *= aa0InvSqrt;
        double vx0x = shuffle(vx0, programIndex ^ 1); // swap re and im
        vx1 = vx1 - vx0*f.re + vx0x * f_imv;
        vx1 *= aa1InvSqrt;
        dx0[c]     = vx0;
        dx1[c]     = vx1;
      }
    }

    // subtract outer product of two top rows from lower part of the matrix
    // process two output rows together
    caxpy2x2tri(triang, rhlen);
    triang += rhlen*4;
  }

  return succ;
}

export uniform bool chol_FactorizeAndSolveFwd(
  uniform DComplex     triang[],
  uniform unsigned int N,
  uniform DComplex     result[])
{
  return chol_FactorizeAndSolveFwdEx(triang, N, result, true);
}

export uniform bool chol_Factorize(
  uniform DComplex     triang[],
  uniform unsigned int N)
{
  return chol_FactorizeAndSolveFwdEx(triang, N, NULL, false);
}
