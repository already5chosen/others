#include <complex.h>
#include <stdio.h>
#include <string.h>

static const double complex triang[] = {
  3.60272769429001827e-01,  1.52194066772070302e+00+I* -1.29392026331265364e+00,  2.80434829472096282e-02+I* -1.25774443360784272e+00,  2.74124543353924555e-01+I*  1.60702518617406076e+00, -1.22222056891831882e-01+I*  1.68066150976443107e+00,
  3.34999380495743626e-01,  1.93915228990529237e-01+I* -1.20073922014596546e+00, -1.11873472798761076e+00+I* -5.14628835678236030e-01, -1.72706196972782400e+00+I* -1.53092492686058179e+00,
  3.03597651587843076e-01, -7.38213788193278764e-01+I* -8.78211708649634937e-02, -6.97378247097579784e-01+I* -8.32802773091502968e-01,
  9.37972701732749270e-01, -9.35873236893126487e-02+I*  2.92247627881069283e+00,
  6.13985388478296534e-01,
};

static const double complex b[] = {
 -2.18064505683800092e-01+I* -2.87175979465980821e-02,  
  1.38552944621677138e-01+I* -1.65599337865485058e-01, 
 -2.55819174296737950e-01+I*  1.91520814800752670e-01,  
  5.75103187838959640e-01+I* -6.73040426368472278e-01,  
  1.50872629467316521e+00+I*  1.25206237895507333e+00,
};

void foo1(double complex *x, int N, const double complex* triang);
void foo2(double complex *x, int N, const double complex* triang);

int main(void)
{
  double complex x[2][5];
  memcpy(x[0], b, sizeof(b));
  memcpy(x[1], b, sizeof(b));
  foo1(x[0], 5, triang);
  foo2(x[1], 5, triang);
  for (int i = 0; i < 5; ++i)
    printf("%25.17e+I*%25.17e   %25.17e+I*%25.17e\n"
      ,creal(x[0][i]), cimag(x[0][i])
      ,creal(x[1][i]), cimag(x[1][i])
    );
  return 0;
}

