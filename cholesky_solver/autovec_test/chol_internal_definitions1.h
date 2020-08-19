typedef struct { double v[4]; } s_m256d;
typedef struct {
  s_m256d re, im;
} complex_sm256d;

#ifdef __cplusplus
extern "C"
#endif
int chol_CrBa4_i2_Band_Factorize(complex_sm256d* triang, unsigned row1st, unsigned rowLast);
