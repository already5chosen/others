#include <windows.h>

long long fadd_loop(volatile long long* pCnt, long long maxN)
{
  long long it = 0;
  while (InterlockedIncrement64(pCnt) < maxN)
    ++it;
  return it;
}
