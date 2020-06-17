#include <stdint.h>

volatile unsigned memset32_c_dummy;
void memset32_c(uint32_t* dst, uint32_t val, unsigned nKbytes, unsigned n_it)
{
  do {
    uint32_t* p = dst;
    for (unsigned i = 0; i < nKbytes; ++i) {
      for (int wi = 0; wi < 1024/4; ++wi)
        p[wi] = val;
      p += 1024/4;
    }
    memset32_c_dummy = n_it;
  } while (--n_it);
}