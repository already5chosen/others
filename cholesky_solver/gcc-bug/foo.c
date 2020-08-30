#include <complex.h>

void foo(double complex *x, int N, const double complex* tri)
{
  for (int r = 0; r < N; ++r) {
   double complex acc = x[r];
   for (int c = 0; c < r; ++c)
     acc -= x[c] * tri[c];
   x[r] = acc;
   tri += r;
  }
}

