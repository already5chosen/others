#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
#include <array>

int main()
{
  typedef std::unordered_map<uint64_t, uint64_t> map_t;
  map_t aSrc, aDst;

  unsigned long long tot = 0;
  for (int nDigits = 1; nDigits <= 19; ++nDigits) {
    aSrc.clear();
    for (int ch0 = 1; ch0 < 10; ++ch0) {
      int r  = ch0 % nDigits;
      unsigned key = (1u << r)*2 + (r==0);
      aSrc[key] += 1;
    }

    map_t* src = &aSrc;
    map_t* dst = &aDst;
    for (int prefixLen = 1; prefixLen < nDigits; ++prefixLen) {
      dst->clear();
      for (auto pSrc = src->cbegin(); pSrc != src->cend(); ++pSrc) {
        uint64_t key = (*pSrc).first;
        unsigned oneChild = key & 1;
        unsigned key0 = key >> 1;
        unsigned key1 = key >> (nDigits+1);
        for (int suffix = 0; suffix < 10; ++suffix) {
          unsigned nxtKey0 = 1u << (suffix % nDigits);
          unsigned nxtKey1 = 0;
          for (int prefixRem = 0; prefixRem < nDigits; ++prefixRem) {
            if ((key0 >> prefixRem) & 1) {
              int nxtRem = (prefixRem*10+suffix) % nDigits;
              unsigned nxtBit = 1u << nxtRem;
              nxtKey1 |= (nxtKey0 & nxtBit);
              nxtKey0 |= nxtBit;
              if ((key1 >> prefixRem) & 1)
                nxtKey1 |= nxtBit;
            }
          }

          int nxtOneChild = oneChild;
          if (nxtKey0 & 1) { // at least 1 child
            if (nxtKey1 & 1)
              continue; // more than one child
            if (oneChild)
              continue; // more than one child
            nxtOneChild = 1;
          }
          // at most 1 child
          uint64_t nxtKey = ((uint64_t(nxtKey1) << nDigits) | nxtKey0)*2+nxtOneChild;
          (*dst)[nxtKey] += (*pSrc).second;
        }
      }
      map_t* tmp = src; src = dst; dst = tmp;
    }

    unsigned long long cnt = 0;
    for (auto pSrc = src->cbegin(); pSrc != src->cend(); ++pSrc)
      if ((*pSrc).first & 1)
        cnt += (*pSrc).second;

    tot += cnt;
    printf("%2d %20llu %20llu\n", nDigits, cnt, tot);
  }

  return 0;
}
