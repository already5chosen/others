#include <windows.h>
#include <stdio.h>

int SetHighPriorityAndAffinity(unsigned i)
{
  HANDLE hThr = GetCurrentThread();
  SetThreadPriority(hThr, THREAD_PRIORITY_HIGHEST);
  if (0==SetThreadAffinityMask(hThr, 1ull << (i % 64))) {
    fprintf(stderr, "SetThreadAffinityMask() at i=%d failed. Error %lu\n", i, GetLastError());
    return 0;
  }
  return 1;
}
