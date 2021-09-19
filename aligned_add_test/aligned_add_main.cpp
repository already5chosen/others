#include <cstdio>
#include <cstdint>
#include <algorithm>
#include <chrono>

extern "C" uint64_t sum64(char* buf, int incr);

static char test_buf[32*501] __attribute__((aligned (64)));

int main()
{
  const int NITER = 19;
  uint64_t dummy = 0;

  int64_t dt[32][NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int offs = 0; offs < 32; ++offs) {
      std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
      dummy += sum64(&test_buf[offs], 0);
      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
      dt[offs][it] = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    }
  }
  if (dummy==42)
    printf("Blue moon.\n");
    
  for (int offs = 0; offs < 32; ++offs) {
    std::nth_element(&dt[offs][0], &dt[offs][NITER/2],&dt[offs][NITER]);
    printf("%2d %5lld\n", offs, (long long)dt[offs][NITER/2]);
  }
  return 0;
}