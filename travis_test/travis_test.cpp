#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <memory>
#include <chrono>
#ifdef WINDOWS_LARGE_PAGES
#include <windows.h>
extern "C" int RightRite(void);
#endif

extern "C" void memset32_c(uint32_t* dst, uint32_t val, unsigned nKbytes, unsigned n_it);

enum {
  nSizes = 32,
  nIter  = 7,
  alignment = 1<<12,
};
int main(int argz, char** argv)
{
  if (argz < 3) {
    fprintf(stderr, "Usage:\ntravis_test fill0 fill1\n");
    return 1;
  }
  uint32_t fill[2];
  for (int arg_i = 1; arg_i < 3; ++arg_i) {
    char* arg = argv[arg_i];
    char* endp;
    fill[arg_i-1] = strtoul(arg, &endp, 0);
    if (endp==arg) {
      fprintf(stderr, "Argument %s is not a number\n", arg);
      return 1;
    }
  }
  printf("fill0=%u, fill1=%u\n", fill[0], fill[1]);

  const int nSamples = nSizes*nIter*2;
  int test_plan[nSamples];
  for (int i = 0; i < nSamples; ++i)
    test_plan[i] = i;
  srand(42);
  std::random_shuffle(test_plan, test_plan+nSamples);

  int sizes[nSizes]; // in KB
  for (int i = 0; i < nSizes; ++i) {
    unsigned sz = 1u << (i/2);
    sizes[i] = (i % 2 == 0) ? 2*sz : 3*sz;
  }

  const unsigned MAX_SZ_KB = sizes[nSizes-1];
  const unsigned RUN_SZ_KB = MAX_SZ_KB*6;
  const size_t MAX_SZ = size_t(MAX_SZ_KB)*1024;

#ifdef WINDOWS_LARGE_PAGES
  RightRite();
  DWORD allocationType = MEM_RESERVE+MEM_COMMIT+MEM_LARGE_PAGES;
  uint32_t* abuf = reinterpret_cast<uint32_t*>(::VirtualAlloc(NULL, MAX_SZ, allocationType, PAGE_READWRITE));
  if (NULL==abuf) {
    fprintf(stderr, "VirtualAlloc(NULL,%llu,0x%x,0x%x) failed. Error %d\n"
      , MAX_SZ, allocationType, PAGE_READWRITE, ::GetLastError());
    return 1;
  }
#else
  std::vector<uint32_t> buf((MAX_SZ+alignment)/4);
  void* al_ptr    = buf.data();
  size_t al_space = MAX_SZ+alignment;
  uint32_t* abuf = reinterpret_cast<uint32_t*>(std::align(alignment, 1, al_ptr, al_space));
#endif

  int64_t res[2][nSizes][nIter];
  for (int samp_i = 0; samp_i < nSamples; ++samp_i) {
    int tst  = test_plan[samp_i];
    int f_i  = tst % 2; tst /= 2;
    int it_i = tst % nIter; tst /= nIter;
    int sz_i = tst;
    int nKbytes = sizes[sz_i];
    int n_it = RUN_SZ_KB / nKbytes;
    memset32_c(abuf, fill[f_i], nKbytes, 1); // pre-fill buffer

    using namespace std::chrono;
    auto t0 = steady_clock::now();
    memset32_c(abuf, fill[f_i], nKbytes, n_it);
    auto t1 = steady_clock::now();
    auto dt = duration_cast<microseconds>(t1 - t0);
    res[f_i][sz_i][it_i] = dt.count();
  }

#if 0
  for (int sz_i = 0; sz_i < nSizes; ++sz_i) {
    printf("%7d KB:", sizes[sz_i]);
    for (int f_i = 0; f_i < 2; ++f_i) {
      printf(" fill=%u", fill[f_i]);
      for (int it_i = 0; it_i < nIter; ++it_i) {
        printf(" %6lld", res[f_i][sz_i][it_i]);
      }
    }
    printf("\n");
  }
#endif

  for (int sz_i = 0; sz_i < nSizes; ++sz_i) {
    printf("%7d KB:", sizes[sz_i]);
    for (int f_i = 0; f_i < 2; ++f_i) {
      std::nth_element(&res[f_i][sz_i][0],&res[f_i][sz_i][nIter/2],&res[f_i][sz_i][nIter]);
      int64_t usec = res[f_i][sz_i][nIter/2]; // median
      // printf(" %6lld", usec);
      printf(" %8.1f", (1024*double(RUN_SZ_KB))/usec);
    }
    printf("\n");
  }

  return 0;
}
