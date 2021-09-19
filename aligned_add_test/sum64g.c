#include <stdint.h>
#include <string.h>

uint64_t sum64(char* buf, int incr)
{
  uint64_t acc0=0;
  uint64_t acc1=0;
  uint64_t acc2=0;
  uint64_t acc3=0;
  for (int i = 0; i < 4000; ++i) {
    char* p = buf;
    for (int k = 0; k < 500; ++k) {
      uint64_t x0, x1, x2, x3;
      memcpy(&x0, &p[0*8], sizeof(x0));
      memcpy(&x1, &p[1*8], sizeof(x1));
      memcpy(&x2, &p[2*8], sizeof(x2));
      memcpy(&x3, &p[3*8], sizeof(x3));
      acc0 += x0;
      acc1 += x1;
      acc2 += x2;
      acc3 += x3;
      p += 32;
    }
    buf += incr;
  }
  return acc0+acc1+acc2+acc3;
}