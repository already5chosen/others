#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <thread>

#include "report_t.h"

enum { NTHREADS_MAX=64, N_INC = 100000000 };
extern "C" report_t fadd_loop(volatile long long* pCnt, long long maxN);
extern "C" int SetHighPriorityAndAffinity(unsigned i);
static void test(int nThreads);

int main(int argz, char** argv)
{
  if (argz < 2) {
    fprintf(stderr,
      "Usage:\n"
      "%s nThreads\n"
      ,argv[0]
      );
    return 1;
  }
  long nThreads = strtol(argv[1], 0, 0);
  if (nThreads < 1 || nThreads > (long)(NTHREADS_MAX)) {
    fprintf(stderr, "nThreads = %s out of range [1:%d]\n", argv[1], NTHREADS_MAX);
    return 1;
  }

  for (int it=0; it < 5; ++it)
    test(nThreads);
  return 0;
}

typedef struct {
  volatile long long c;
  volatile long      begC, endC;
  std::chrono::steady_clock::time_point      t0, t1;
} common_data_t;

typedef struct {
  report_t       r;
  common_data_t* common;
} thread_data_t;

static void TestThreadProc(thread_data_t* threadData)
{
  if (!SetHighPriorityAndAffinity(threadData->r.c[0]))
    return;

  common_data_t* common = threadData->common;
  long begRem = __atomic_add_fetch(&common->begC, -1L, __ATOMIC_RELAXED);
  if (begRem > 0) {
    while (common->begC > 0); // wait for start of remaining threads
  } else {
    common->t0 = std::chrono::steady_clock::now(); // time stamp
  }

  threadData->r = fadd_loop(&common->c, N_INC);

  if (__atomic_add_fetch(&common->endC, -1L, __ATOMIC_RELAXED)==0) {
    // last thread
    common->t1 = std::chrono::steady_clock::now(); // time stamp
  }
}

static void test(int nThreads)
{
  common_data_t commonData = {0};
  commonData.begC = nThreads;
  commonData.endC = nThreads;

  thread_data_t threadsData[NTHREADS_MAX];
  std::thread threads[NTHREADS_MAX];
  for (int ti = 0; ti < nThreads; ++ti) {
    threadsData[ti].common = &commonData;
    threadsData[ti].r.c[0] = ti;
    threads[ti] = std::thread(TestThreadProc, &threadsData[ti]);
  }
  for (int ti = 0; ti < nThreads; ++ti) {
    threads[ti].join();
  }

  // double sec = (double)(commonData.t1 - commonData.t0);
  using namespace std::chrono;
  double sec = duration_cast<duration<double>>(commonData.t1 - commonData.t0).count();
  printf("%.6f %.1f  ", sec, sec*1e9/N_INC);
  long long mn=threadsData[0].r.c[0];
  long long mx=mn;
  long long totNBounce = 0;
  for (int ti = 0; ti < nThreads; ++ti) {
    long long ni = threadsData[ti].r.c[0];
    long long nBounce = threadsData[ti].r.c[1];
    totNBounce += nBounce;
    double bounceR = 0;
    if (ni > 0) {
      bounceR = (double)nBounce/ni;
    }
    printf(" %lld %5.3f", ni, bounceR);
    if (mn > ni) mn = ni;
    if (mx < ni) mx = ni;
  }
  double bounceR = (double)totNBounce/(N_INC-1);
  printf("   %.4f %5.3f\n", (double)((mx-mn)*nThreads)/(N_INC-1), bounceR);
  fflush(stdout);
}