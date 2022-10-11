#include <cstdint>
#include <cstdio>
#include <random>
#include <memory>
#include <vector>
#include <algorithm>
#include <chrono>

extern "C" {
void ref_ByteConversion(size_t length, const uint8_t *lut, const uint8_t *input, uint8_t *output);
void uut_ByteConversion(size_t length, const uint8_t *lut, const uint8_t *input, uint8_t *output);
uint64_t Checksum(const uint32_t *input, size_t quarterlen);
}


int main(int , char** )
{
  const int BUF_LEN = 1024*12;   // for both input and output to fit comfortably into 32K L1D
  const int REP_PER_ITER = 5000; // for individual test to take order of 1 msec
  const int N_ITER  = 201;

  // prepare aligned buffers
  const int BUF_ALIGNMENT = 64;
  const int BUF_MEM_SZ = BUF_LEN + BUF_ALIGNMENT;

  std::size_t align_space;
  std::vector<char> lutmem(256+BUF_ALIGNMENT), inpmem(BUF_MEM_SZ), outmem1(BUF_MEM_SZ), outmem2(BUF_MEM_SZ);

  void* lutbuf = lutmem.data();
  align_space = 256+BUF_ALIGNMENT;
  std::align(BUF_ALIGNMENT, 1, lutbuf, align_space);

  void* inpbuf = inpmem.data();
  align_space = BUF_MEM_SZ;
  std::align(BUF_ALIGNMENT, 1, inpbuf, align_space);

  void* outbuf1 = outmem1.data();
  align_space = BUF_MEM_SZ;
  std::align(BUF_ALIGNMENT, 1, outbuf1, align_space);

  void* outbuf2 = outmem2.data();
  align_space = BUF_MEM_SZ;
  std::align(BUF_ALIGNMENT, 1, outbuf2, align_space);

  // fill input and lut with random data
  std::mt19937_64 rndGen(42);
  for (int i = 0; i < 256; ++i)
    ((uint8_t*)lutbuf)[i] = rndGen();
  for (int i = 0; i < BUF_LEN; ++i)
    ((uint8_t*)inpbuf)[i] = rndGen();

  std::chrono::steady_clock::now();
  std::vector<long long>  ref_dt(N_ITER), uut_dt(N_ITER);
  for (int it = 0; it < N_ITER; ++it) {
    auto t0 = std::chrono::steady_clock::now();
    uint64_t ref_cs = 0;
    for (int i = 0; i < REP_PER_ITER; ++i) {
      ref_ByteConversion(BUF_LEN, (uint8_t*)lutbuf, (uint8_t*)inpbuf, (uint8_t*)outbuf1);
      ref_cs += Checksum((uint32_t*)outbuf1, BUF_LEN/16);
    }
    auto t1 = std::chrono::steady_clock::now();
    uint64_t uut_cs = 0;
    for (int i = 0; i < REP_PER_ITER; ++i) {
      uut_ByteConversion(BUF_LEN/8, (uint8_t*)lutbuf, (uint8_t*)inpbuf, (uint8_t*)outbuf1);
      uut_cs += Checksum((uint32_t*)outbuf1, BUF_LEN/16);
    }
    auto t2 = std::chrono::steady_clock::now();

    std::chrono::duration<long long, std::nano> dur1 = t1 - t0;
    ref_dt[it] = dur1.count();
    std::chrono::duration<long long, std::nano> dur2 = t2 - t1;
    uut_dt[it] = dur2.count();

    if (ref_cs != uut_cs) {
      printf("Checksum mismatch.\n");
      return 1;
    }
  }

  std::nth_element(ref_dt.data(), ref_dt.data() + N_ITER/2, ref_dt.data() + N_ITER); // median
  long long ref_dtMed = ref_dt[N_ITER/2];

  std::nth_element(uut_dt.data(), uut_dt.data() + N_ITER/2, uut_dt.data() + N_ITER); // median
  long long uut_dtMed = uut_dt[N_ITER/2];

  printf(
    "ref %12.6f msec. %12.3f MB/s\n"
    "uut %12.6f msec. %12.3f MB/s\n"
    ,ref_dtMed * 1e-6, BUF_LEN*REP_PER_ITER*1e3/ref_dtMed
    ,uut_dtMed * 1e-6, BUF_LEN*REP_PER_ITER*1e3/uut_dtMed
  );

  return 0;
}
