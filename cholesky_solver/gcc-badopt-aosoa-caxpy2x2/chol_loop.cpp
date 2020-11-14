void chol_FactorizeMainLoop(double* y, double* x0, double* x1, unsigned rhlen)
{
  // subtract outer product of two top rows from lower part of the matrix
  // process two output rows together
  for (unsigned chlen = rhlen-1; chlen > 0; chlen -= 1) {
    int xi = chlen & 1;
    double f00_re = x0[4*0+xi*2+0];
    double f10_re = x1[4*0+xi*2+0];
    double f01_re = x0[4*0+xi*2+1];
    double f11_re = x1[4*0+xi*2+1];
    double f00_im = x0[4*1+xi*2+0];
    double f10_im = x1[4*1+xi*2+0];
    double f01_im = x0[4*1+xi*2+1];
    double f11_im = x1[4*1+xi*2+1];
    int cqlen = (chlen+1)/2;
    double* y0 = &y[0];
    double* y1 = &y[cqlen*8];
    for (int c = 0; c < cqlen; ++c) {
      // y0[c] = y0[c] - x0[c]*conj(f00) - x1[c]*conj(f10);
      // y1[c] = y1[c] - x0[c]*conj(f01) - x1[c]*conj(f11);
      for (int k = 0; k < 4; ++k) {
        double x0_re = x0[c*8+4*0+k];
        double x0_im = x0[c*8+4*1+k];
        double y0_re = y0[c*8+4*0+k];
        double y0_im = y0[c*8+4*1+k];
        double y1_re = y1[c*8+4*0+k];
        double y1_im = y1[c*8+4*1+k];
        y0_re = y0_re - x0_re * f00_re - x0_im * f00_im;
        y0_im = y0_im + x0_re * f00_im - x0_im * f00_re;
        y1_re = y1_re - x0_re * f01_re - x0_im * f01_im;
        y1_im = y1_im + x0_re * f01_im - x0_im * f01_re;
        double x1_re = x1[c*8+4*0+k];
        double x1_im = x1[c*8+4*1+k];
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
