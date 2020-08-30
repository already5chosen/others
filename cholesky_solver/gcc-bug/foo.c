#include <complex.h>

void foo(double complex *x, int N, const double complex* y)
{
  double complex acc = x[N];
  for (int c = 0; c < N; ++c)
    acc -= x[c] * y[c];
  x[N] = acc;
}
