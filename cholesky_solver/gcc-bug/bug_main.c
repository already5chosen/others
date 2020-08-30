#include <complex.h>
#include <stdio.h>
#include <string.h>

static const double complex triang[] = {
  0,
  1, 2,
};

static const double complex b[] = {
  1,
  3,
  0,
};

void foo1(double complex *x, int N, const double complex* triang);
void foo2(double complex *x, int N, const double complex* triang);

int main(void)
{
  enum { N = sizeof(b)/sizeof(b[0]) };
  double complex x[2][N];
  memcpy(x[0], b, sizeof(b));
  memcpy(x[1], b, sizeof(b));
  foo1(x[0], N, triang); // reference
  foo2(x[1], N, triang); // buggy
  for (int i = 0; i < N; ++i)
    printf(" %+.0f%+.0fi %+.0f%+.0fi\n"
      ,creal(x[0][i]), cimag(x[0][i])
      ,creal(x[1][i]), cimag(x[1][i])
    );
  return 0;
}

