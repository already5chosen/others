#include <windows.h>
#include "report_t.h"

report_t fadd_loop_sqrsum(volatile long long* pCnt, long long maxN)
{
  long long it = 0;
  long long sqrSum = 0;
  long long cnt, prevCnt = *pCnt;
  while ((cnt=InterlockedIncrement64(pCnt)) < maxN) {
    ++it;
    long long d = prevCnt - cnt;
    sqrSum += d * d;
    prevCnt = cnt;
  }
  report_t res = { {it, sqrSum } };
  return res;
}
