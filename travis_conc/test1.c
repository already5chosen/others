#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

enum { NTHREADS_MAX=128, N_INC = 100000000 };
long long fadd_loop(volatile long long* pCnt, long long maxN);
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
  LARGE_INTEGER      t0, t1;
  HANDLE             hEv;
} common_data_t;

typedef struct {
  long long      ni;
  common_data_t* common;
} thread_data_t;

static DWORD WINAPI TestThreadProc(void* pPrm)
{
  thread_data_t* threadData = pPrm;
  common_data_t* common = threadData->common;
  long begRem = InterlockedDecrement(&common->begC);
  if (begRem > 0) {
    while (common->begC > 0); // wait for start of remaining threads
  } else {
    QueryPerformanceCounter(&common->t0); // time stamp
  }

  threadData->ni = fadd_loop(&common->c, N_INC);

  if (InterlockedDecrement(&common->endC)==0) {
    // last thread
    QueryPerformanceCounter(&common->t1); // time stamp
    SetEvent(common->hEv); // tell to test that we're done
  }

  return 0;
}

static void test(int nThreads)
{
  common_data_t commonData = {0};
  commonData.begC = nThreads;
  commonData.endC = nThreads;
  commonData.hEv = CreateEvent(NULL, FALSE, FALSE, NULL);
  if (commonData.hEv==NULL) {
    fprintf(stderr, "CreateEvent() failed. Error %lu\n", GetLastError());
    exit(1);
  }
  thread_data_t threadsData[NTHREADS_MAX];
  for (int ti = 0; ti < nThreads; ++ti) {
    threadsData[ti].common = &commonData;
    if (!QueueUserWorkItem(TestThreadProc, &threadsData[ti], 0)) {
      fprintf(stderr, "QueueUserWorkItem() at ti=%d failed. Error %lu\n", ti, GetLastError());
      exit(1);
    }
  }
  WaitForSingleObject(commonData.hEv, INFINITE);
  CloseHandle(commonData.hEv);

  LARGE_INTEGER fr;
  QueryPerformanceFrequency(&fr);
  double sec = (double)(commonData.t1.QuadPart - commonData.t0.QuadPart)/fr.QuadPart;
  printf("%.6f %.1f  ", sec, sec*1e9/N_INC);
  long long mn=threadsData[0].ni;
  long long mx=mn;
  for (int ti = 0; ti < nThreads; ++ti) {
    long long ni = threadsData[ti].ni;
    printf(" %lld", ni);
    if (mn > ni) mn = ni;
    if (mx < ni) mx = ni;
  }
  printf(" %.4f\n", (double)((mx-mn)*nThreads)/N_INC);
  fflush(stdout);
}