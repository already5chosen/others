#include "report_t.h"

report_t fadd_loop(volatile long long* pCnt, long long maxN)
{
  long long it = 0;
  long long nBounce = 0;
  long long cnt, nxtCnt = *pCnt;
  while ((cnt=__atomic_add_fetch(pCnt, 1, __ATOMIC_RELAXED)-1) < maxN) {
    ++it;
    nBounce += (cnt != nxtCnt);
    nxtCnt = cnt + 1;
  }
  report_t res = { {it, nBounce } };
  return res;
}
