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
      uint8_t out0 = lut[inp0];
      uint8_t out1 = lut[inp1];
      uint8_t out2 = lut[inp2];
      uint8_t out3 = lut[inp3];
      uint8_t out4 = lut[inp4];
      uint8_t out5 = lut[inp5];
      uint8_t out6 = lut[inp6];
      uint8_t out7 = lut[inp7];

      output[0] = out0;
      output[1] = out1;
      output[2] = out2;
      output[3] = out3;
      output[4] = out4;
      output[5] = out5;
      output[6] = out6;
      output[7] = out7;
      output += 8;
    } while (--length);
  }
}
