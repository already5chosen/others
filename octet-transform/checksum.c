#include <stdlib.h>
#include <stdint.h>

uint64_t Checksum(const uint32_t *input, size_t quarterlen)
{
  uint64_t sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
  do {
    sum0 += input[0];
    sum1 += input[1];
    sum2 += input[2];
    sum3 += input[3];
    input += 4;
  } while (--quarterlen);
  return sum0+sum1*3+sum2*5+sum3*7;
}
