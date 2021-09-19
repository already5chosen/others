#include <stdint.h>
#include <string.h>
#include <intrin.h>

uint64_t sum64(char* buf, int incr)
{
  __m128i acc0=_mm_setzero_si128();
  __m128i acc1=_mm_setzero_si128();
  __m128i acc2=_mm_setzero_si128();
  __m128i acc3=_mm_setzero_si128();
  for (int i = 0; i < 4000; ++i) {
    char* p = buf;
    for (int k = 0; k < 250; ++k) {
      __m128i x0 = _mm_loadu_si128((__m128i*)&p[0*16]);
      __m128i x1 = _mm_loadu_si128((__m128i*)&p[1*16]);
      __m128i x2 = _mm_loadu_si128((__m128i*)&p[2*16]);
      __m128i x3 = _mm_loadu_si128((__m128i*)&p[3*16]);
      acc0 = _mm_add_epi64(acc0, x0);
      acc1 = _mm_add_epi64(acc1, x1);
      acc2 = _mm_add_epi64(acc2, x2);
      acc3 = _mm_add_epi64(acc3, x3);
      p += 16*4;
    }
    buf += incr;
  }
  acc0 = _mm_add_epi64(acc0, acc1);
  acc2 = _mm_add_epi64(acc2, acc3);
  acc0 = _mm_add_epi64(acc0, acc2);
  return 
    _mm_extract_epi64(acc0, 0) +
    _mm_extract_epi64(acc0, 1) ;
}