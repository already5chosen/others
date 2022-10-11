#include <stdlib.h>
#include <stdint.h>

void ref_ByteConversion(size_t length, const uint8_t *lut, const uint8_t *input, uint8_t *output)
{
  for (size_t i=0; i<length; i++)
    output[i] = lut[input[i]];
}
