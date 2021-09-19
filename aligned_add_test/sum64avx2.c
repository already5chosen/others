#include <stdint.h>
#include <string.h>
#include <intrin.h>

uint64_t sum64(char* buf, int incr)
{
  __m256i acc0=_mm256_setzero_si256();
  __m256i acc1=_mm256_setzero_si256();
  __m256i acc2=_mm256_setzero_si256();
  __m256i acc3=_mm256_setzero_si256();
  for (int i = 0; i < 4000; ++i) {
    char* p = buf;
    for (int k = 0; k < 125; ++k) {
      __m256i x0 = _mm256_loadu_si256((__m256i*)&p[0*32]);
      __m256i x1 = _mm256_loadu_si256((__m256i*)&p[1*32]);
      __m256i x2 = _mm256_loadu_si256((__m256i*)&p[2*32]);
      __m256i x3 = _mm256_loadu_si256((__m256i*)&p[3*32]);
      acc0 = _mm256_add_epi64(acc0, x0);
      acc1 = _mm256_add_epi64(acc1, x1);
      acc2 = _mm256_add_epi64(acc2, x2);
      acc3 = _mm256_add_epi64(acc3, x3);
      p += 32*4;
    }
    buf += incr;
  }
  acc0 = _mm256_add_epi64(acc0, acc1);
  acc2 = _mm256_add_epi64(acc2, acc3);
  acc0 = _mm256_add_epi64(acc0, acc2);
  return
    _mm256_extract_epi64(acc0, 0) +
    _mm256_extract_epi64(acc0, 1) +
    _mm256_extract_epi64(acc0, 2) +
    _mm256_extract_epi64(acc0, 3) ;
}