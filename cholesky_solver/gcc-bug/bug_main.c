#include <complex.h>
#include <stdio.h>

static const double complex y[]  = { 1, 2, };
static const double complex x0[] = {1, 3, 0,};

void foo1(double complex *x, int N, const double complex* y);
void foo2(double complex *x, int N, const double complex* y);

int main(void)
{
  enum { N = sizeof(x0)/sizeof(x0[0]) };
  double complex x[2][N];
  for (int i = 0; i < N; ++i)
    x[0][i] = x[1][i] = x0[i];

  foo1(x[0], N-1, y); // reference
  foo2(x[1], N-1, y); // buggy

  for (int i = 0; i < N; ++i)
    printf(" %+.0f%+.0fi %+.0f%+.0fi\n"
      ,creal(x[0][i]), cimag(x[0][i])
      ,creal(x[1][i]), cimag(x[1][i])
    );
  return 0;
}

