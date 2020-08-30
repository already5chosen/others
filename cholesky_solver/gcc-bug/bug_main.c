#include <complex.h>
#include <stdio.h>


double complex foo1(double complex acc, const double complex *x, const double complex* y, int N);
double complex foo2(double complex acc, const double complex *x, const double complex* y, int N);

int main(void)
{
  static const double complex y[] = { 1, 2, };
  static const double complex x[] = { 1, 3, };

  double complex ref = foo1(0, x, y, 2); // reference
  double complex res = foo2(0, x, y, 2); // buggy

  printf(" %+.0f%+.0fi %+.0f%+.0fi\n"
    ,creal(ref), cimag(ref)
    ,creal(res), cimag(res)
  );
  return 0;
}

