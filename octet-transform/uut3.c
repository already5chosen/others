#include <stdlib.h>
#include <string.h>
#include <stdint.h>

void uut_ByteConversion(size_t length, const uint8_t *lut, const uint8_t *input, uint8_t *output)
{
  if (length > 0) {
    do {
      uint64_t inpx;
      memcpy(&inpx, input, sizeof(inpx));
      size_t inp0 = inpx & 255;
      size_t inp1 = (inpx >> (8*1)) & 255;
      size_t inp2 = (inpx >> (8*2)) & 255;
      size_t inp3 = (inpx >> (8*3)) & 255;
      size_t inp4 = (inpx >> (8*4)) & 255;
      size_t inp5 = (inpx >> (8*5)) & 255;
      size_t inp6 = (inpx >> (8*6)) & 255;
      size_t inp7 = (inpx >> (8*7)) & 255;
      input += 8;
      asm volatile("" ::: "memory");
      uint32_t out0 = lut[inp0];
      uint32_t out1 = lut[inp1];
      uint32_t out2 = lut[inp2];
      uint32_t out3 = lut[inp3];
      uint32_t out4 = lut[inp4];
      uint32_t out5 = lut[inp5];
      uint32_t out6 = lut[inp6];
      uint32_t out7 = lut[inp7];

      out0 |= out1 << 8;
      out2 |= out3 << 8;
      out4 |= out5 << 8;
      out6 |= out7 << 8;

      out0 |= out2 << 16;
      out4 |= out6 << 16;

      uint64_t outx = (uint64_t)out0 | ((uint64_t)out4 << 32);

      memcpy(output, &outx, sizeof(outx));
      output += 8;
    } while (--length);
  }
}
