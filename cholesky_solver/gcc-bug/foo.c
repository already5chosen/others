#include <complex.h>

void foo(double complex *x, int N, const double complex* triang)
{
  triang += N*(N+1)/2; // point past end
  x += N;              // point past end
  for (int rlen = 1; rlen <= N; ++rlen) {
   triang -= rlen; // point to diag element
   x -= 1;
   double complex acc = x[0];
   for (int c = 1; c < rlen; ++c)
     acc -= x[c] * triang[c];
   x[0] = acc * creal(triang[0]);
   // double complex acc = x[0];
   // for (int c = 1; c < rlen; ++c)
     // acc -= x[c] * triang[c];
   // x[0] = acc * creal(triang[0]);
  }
}

