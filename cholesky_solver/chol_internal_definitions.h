#ifdef __AVX__
// macro __AVX__ is not available in VS2010, so we use our own use__AVX__
// However on newer versions of VS we wont to utilize automatic recognition of /arch:avx 
 #ifndef use__AVX__
  #define use__AVX__
 #endif
#endif

#ifdef CBF_SIMD256
 #ifndef use__AVX__
  #define use__AVX__
 #endif
#endif

#ifdef CBF_SIMD128
 #ifdef use__AVX__
  #undef use__AVX__
 #endif
#endif

#ifdef use__AVX__

typedef struct {
  __m256d re, im;
} complex_m256d;

#ifdef __cplusplus
extern "C"
#endif
int chol_CrBa4_i2_Band_Factorize(complex_m256d* triang, unsigned row1st, unsigned rowLast);

#ifdef __cplusplus
extern "C"
#endif
void chol_CrBa4_i2_Band_ForwardSubstitute(
  __m256d*              result, 
  const complex_m256d*  triang, 
  unsigned row1st, unsigned rowLast);

#ifdef __cplusplus
extern "C"
#endif
void chol_CrBa4_i2_BackSubstitute(
  __m256d*              result, 
  const complex_m256d*  triang, 
  unsigned n);

#else

typedef struct {
  __m128d re, im;
} complex_m128d;

#ifdef __cplusplus
extern "C"
#endif
int chol_CrBa4_i2_Band_Factorize(complex_m128d* triang, unsigned row1st, unsigned rowLast);

#ifdef __cplusplus
extern "C"
#endif
void chol_CrBa4_i2_Band_ForwardSubstitute(
  __m128d*             result, 
  const complex_m128d* triang, 
  unsigned row1st, unsigned rowLast);

#ifdef __cplusplus
extern "C"
#endif
void chol_CrBa4_i2_BackSubstitute(
  __m128d*             result, 
  const complex_m128d* triang, 
  unsigned n);

#endif