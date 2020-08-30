#include <complex.h>

double complex foo(double complex acc, const double complex *x, const double complex* y, int N)
{
  for (int c = 0; c < N; ++c)
    acc -= x[c] * y[c];
  return acc;
}
