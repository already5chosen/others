#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_SRWL

typedef SRWLOCK   mtx_t;
#define PRG_NAME  "srw_lock_test"
static __inline void InitializeMtx(mtx_t* pMtx) { InitializeSRWLock(pMtx); }
static __inline void AcquireMtx(mtx_t* pMtx)    { AcquireSRWLockExclusive(pMtx); }
static __inline void ReleaseMtx(mtx_t* pMtx)    { ReleaseSRWLockExclusive(pMtx); }

#else

typedef CRITICAL_SECTION mtx_t;
#define PRG_NAME  "crit_sec_test"
static __inline void InitializeMtx(mtx_t* pMtx) { InitializeCriticalSection(pMtx); }
// static __inline void InitializeMtx(mtx_t* pMtx) { InitializeCriticalSectionAndSpinCount(pMtx, 1000); }
static __inline void AcquireMtx(mtx_t* pMtx)    { EnterCriticalSection(pMtx); }
static __inline void ReleaseMtx(mtx_t* pMtx)    { LeaveCriticalSection(pMtx); }

#endif

static mtx_t* mtx;
static int RunTest(int n_obj);

int main(int argz, char** argv)
{
  if (argz < 2) {
    printf(
      "Usage:\n"
      "%s n-objects [n-iter]\n"
      , PRG_NAME
    );
    return 1;
  }

  long arg_val = strtol(argv[1], 0, 0);
  if (arg_val < 100 || arg_val > 100000000) {
    printf("Parameter n-objects = '%s' out of range. Please specify number in range [100:100,000,000]\n", argv[1]);
    return 1;
  }
  int n_obj = (int)arg_val;

  int n_iter = 1;
  if (argz > 2) {
    arg_val = strtol(argv[2], 0, 0);
    if (arg_val < 1 || arg_val > 1000) {
      printf("Parameter n-iter = '%s' out of range. Please specify number in range [1:1000]\n", argv[2]);
      return 1;
    }
    n_iter = (int)arg_val;
  }

  printf("Running test for %d objects. %d iterations\n", n_obj, n_iter);

  printf("Allocating %zd*%d = %zu bytes...", sizeof(mtx_t), n_obj, sizeof(mtx_t)*n_obj);
  DWORD t0 = GetTickCount();
  mtx_t* p = calloc(sizeof(*p), n_obj);
  if (p == NULL) {
    printf("fail.\n");
    return 1;
  }
  DWORD t1 = GetTickCount();
  printf("done. %u msec.\n", t1 - t0);

  printf("Initializing objects...");
  t0 = GetTickCount();
  for (int i = 0; i < n_obj; ++i)
    InitializeMtx(&p[i]);
  t1 = GetTickCount();
  printf("done. %u msec.\n", t1 - t0);

  DWORD handleCount;
  if (GetProcessHandleCount(GetCurrentProcess(), &handleCount))
     printf("%u handles\n", handleCount);

  printf("Priming objects...");
  t0 = GetTickCount();
  for (int i = 0; i < n_obj; ++i) {
    AcquireMtx(&p[i]);
    ReleaseMtx(&p[i]);
  }
  t1 = GetTickCount();
  printf("done. %u msec.\n", t1 - t0);
  if (GetProcessHandleCount(GetCurrentProcess(), &handleCount))
     printf("%u handles\n", handleCount);

  mtx = p;
  int ret = 0;
  for (int it = 0; it < n_iter && ret == 0; ++it) {
    ret = RunTest(n_obj);
    if (GetProcessHandleCount(GetCurrentProcess(), &handleCount))
       printf("%u handles\n", handleCount);
  }

#ifndef USE_SRWL
  printf("Destroying objects...");
  t0 = GetTickCount();
  for (int i = 0; i < n_obj; ++i)
    DeleteCriticalSection(&p[i]);
  t1 = GetTickCount();
  printf("done. %u msec.\n", t1 - t0);
  if (GetProcessHandleCount(GetCurrentProcess(), &handleCount))
     printf("%u handles\n", handleCount);
#endif

  printf("Freeing...");
  t0 = GetTickCount();
  free(p);
  t1 = GetTickCount();
  printf("done. %u msec.\n", t1 - t0);

  return ret;
}


enum {
  STATE_INIT,
  STATE_IDLE_THREAD_RUNNING,
  STATE_BEGIN,
  STATE_ACQUIRE_ENA,
  STATE_ACQUIRE,
  STATE_AT_IDLE_THREAD,
};
static volatile int state;
static volatile int st_n_steps;
static volatile int test_step;
static volatile unsigned long long thr1_tsc_sum1, thr1_tsc_sum2, thr2_tsc_sum;

static DWORD WINAPI Thread1(void* dummy)
{
  (void)dummy;
  state = STATE_BEGIN;

  const int n_steps = st_n_steps;
  mtx_t* p = mtx;
  unsigned long long sum1 = 0, sum2 = 0;
  for (;;) {
    while (state != STATE_ACQUIRE_ENA);
    // state == STATE_ACQUIRE_ENA
    unsigned long long tsc1 = __rdtsc();
    int step = test_step;
    state = STATE_ACQUIRE; // tell to thread2 to trigger release of mutex after we are blocked
    if (step >= n_steps)
      break;
    AcquireMtx(&p[step]); // we will be blocked and then preempted
    unsigned long long tsc2 = __rdtsc();
    state = STATE_BEGIN;  // tell to main thread that it can continue
    ReleaseMtx(&p[step]);
    sum1 += tsc1;
    sum2 += tsc2;
  }

  thr1_tsc_sum1 = sum1;
  thr1_tsc_sum2 = sum2;
  return 0;
}

static DWORD WINAPI Thread2(void* dummy)
{
  (void)dummy;
  state = STATE_IDLE_THREAD_RUNNING;

  const int n_steps = st_n_steps;
  unsigned long long sum = 0;
  for (;;) {
    while (state != STATE_ACQUIRE);
    // state == STATE_ACQUIRE
    // we come here only after Thread1 is blocked by mutex
    unsigned long long tsc = __rdtsc();
    int step = test_step;
    state = STATE_AT_IDLE_THREAD; // tell to main thread to release mutex
    if (step >= n_steps)
      break;
    sum += tsc;
  }

  thr2_tsc_sum = sum;
  return 0;
}

static HANDLE CreateTestThread(LPTHREAD_START_ROUTINE pStartAddress, int prio)
{
  HANDLE h = CreateThread(NULL, 0, pStartAddress, NULL, CREATE_SUSPENDED, NULL);
  if (h != NULL) {
    if (SetThreadAffinityMask(h, 1)) {
      if (SetThreadPriority(h, prio)) {
        return h; // success
      } else {
        printf("SetThreadPriority(%d) failed. Error %d.", prio, GetLastError());
      }
    } else {
      printf("SetThreadAffinityMask(1) failed. Error %d.", GetLastError());
    }
    CloseHandle(h);
  } else {
    printf("CreateThread(%p) failed. Error %d.", pStartAddress, GetLastError());
  }
  return NULL; // fail
}

static int RunTest(int n_obj)
{
  if (!SetThreadAffinityMask(GetCurrentThread(), (DWORD)~1)) {
    printf("main thread: SetThreadAffinityMask() failed. Error %d.", GetLastError());
    return 1;
  }

  st_n_steps = n_obj;
  test_step = 0;
  state = STATE_INIT;
  int err = 1;
  printf("Creating threads...");
  DWORD t0 = GetTickCount();
  HANDLE hThr1 = CreateTestThread(Thread1, THREAD_PRIORITY_TIME_CRITICAL);
  if (hThr1 != NULL) {
    HANDLE hThr2 = CreateTestThread(Thread2, THREAD_PRIORITY_IDLE);
    if (hThr2 != NULL) {
      DWORD t1 = GetTickCount();
      printf("done. %u msec. ", t1 - t0);
      // both threads forced to run on core 0

      printf("Threads initialization...");
      ResumeThread(hThr2);
      while (state != STATE_IDLE_THREAD_RUNNING) ;
      ResumeThread(hThr1);
      while (state != STATE_BEGIN) ;
      DWORD t2 = GetTickCount();
      printf("done. %u msec.\n", t2 - t1);

      // run test loop
      unsigned long long tsc_sum[8] = {0};
      printf("Running test...");
      t0 = GetTickCount();
      mtx_t* p = mtx;
      for (int step = 0; step < n_obj; ++step) {
        test_step = step;
        tsc_sum[0] += __rdtsc();
        AcquireMtx(&p[step]);
        tsc_sum[1] += __rdtsc();
        state = STATE_ACQUIRE_ENA;
        while (state != STATE_AT_IDLE_THREAD) ;
        // state == STATE_AT_IDLE_THREAD
        tsc_sum[4] += __rdtsc();
        ReleaseMtx(&p[step]);
        tsc_sum[5] += __rdtsc();
        while (state != STATE_BEGIN) ;
        tsc_sum[7] += __rdtsc();
      }
      test_step = n_obj;
      state = STATE_ACQUIRE_ENA;
      t1 = GetTickCount();
      printf("done. %u msec. ", t1 - t0);

      printf("Waiting for threads termination...");
      t0 = GetTickCount();
      WaitForSingleObject(hThr2, INFINITE);
      t1 = GetTickCount();
      WaitForSingleObject(hThr1, INFINITE);
      t2 = GetTickCount();
      printf("done. %u+%u msec.\n"
        , t1 - t0
        , t2 - t1
      );

      tsc_sum[2] = thr1_tsc_sum1;
      tsc_sum[3] = thr2_tsc_sum;
      tsc_sum[6] = thr1_tsc_sum2;

      static const char* cntr_names[7] = {
        "Fast Acquire (at main)",
        "main->thread1",
        "Slow Acquire Entrance + wakeup of idle thread",
        "idle thread->main",
        "Slow Release (at main)",
        "Slow Acquire Exit == wakeup of thread1",
        "thread1->main",
      };
      for (int i = 0; i < 7; ++i)  {
        unsigned long long d = tsc_sum[i+1]-tsc_sum[i];
        printf(" %9llu_%06llu %s\n", d /1000000, d % 1000000, cntr_names[i]);
      }
      // printf("\n");

      CloseHandle(hThr2);
      err = 0; // success
    }
    CloseHandle(hThr1);
  }
  return err;
}
