#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <immintrin.h>

#ifdef _MSC_VER
__declspec(align(32))
#endif
static float sincos_tab[(512+128+32+8)*2*3+1]
#ifdef __GNUC__
__attribute__ ((aligned (32)))
#endif
;
static unsigned char nibrev_tab[64];

#ifdef __AVX2__
#define MULADD( x, y, acc) _mm256_fmadd_ps (x, y, acc)
#define MULnADD(x, y, acc) _mm256_fnmadd_ps(x, y, acc)
#else
#define MULADD( x, y, acc) _mm256_add_ps(acc, _mm256_mul_ps(x,y))
#define MULnADD(x, y, acc) _mm256_sub_ps(acc, _mm256_mul_ps(x,y))
#endif

void cfft2048_r4dif(float srcDst[2048*2])
{
#ifdef _MSC_VER
  __declspec(align(32))
#endif
  float wbuf[2048*2]
#ifdef __GNUC__
  __attribute__ ((aligned (32)))
#endif
  ;
  const float* pTab = sincos_tab;
  float* dst = wbuf;
  float* src = srcDst;
  for (int i = 0; i < 2048/4/8; ++i) {
    __m256 c00 = _mm256_loadu_ps(&src[8*0+512*2*0]); // re0 im0 re1 im1 re2 im2 re3 im3
    __m256 c01 = _mm256_loadu_ps(&src[8*1+512*2*0]); // re4 im4 re5 im5 re6 im6 re7 im7
    __m256 c10 = _mm256_loadu_ps(&src[8*0+512*2*1]);
    __m256 c11 = _mm256_loadu_ps(&src[8*1+512*2*1]);
    __m256 c20 = _mm256_loadu_ps(&src[8*0+512*2*2]);
    __m256 c21 = _mm256_loadu_ps(&src[8*1+512*2*2]);
    __m256 c30 = _mm256_loadu_ps(&src[8*0+512*2*3]);
    __m256 c31 = _mm256_loadu_ps(&src[8*1+512*2*3]);

    __m256 re0 = _mm256_shuffle_ps(c00, c01, 0x88); // re0 re1 re4 re5  re2 re3 re6 re7
    __m256 im0 = _mm256_shuffle_ps(c00, c01, 0xDD); // im0 im1 im4 im5  im2 im3 im6 im7
    __m256 re1 = _mm256_shuffle_ps(c10, c11, 0x88);
    __m256 im1 = _mm256_shuffle_ps(c10, c11, 0xDD);
    __m256 re2 = _mm256_shuffle_ps(c20, c21, 0x88);
    __m256 im2 = _mm256_shuffle_ps(c20, c21, 0xDD);
    __m256 re3 = _mm256_shuffle_ps(c30, c31, 0x88);
    __m256 im3 = _mm256_shuffle_ps(c30, c31, 0xDD);

    __m256 reH0 = _mm256_add_ps(re0, re2);
    __m256 reH2 = _mm256_sub_ps(re0, re2);
    __m256 imH0 = _mm256_add_ps(im0, im2);
    __m256 imH2 = _mm256_sub_ps(im0, im2);

    __m256 reH1 = _mm256_add_ps(re1, re3);
    __m256 reH3 = _mm256_sub_ps(re1, re3);
    __m256 imH1 = _mm256_add_ps(im1, im3);
    __m256 imH3 = _mm256_sub_ps(im1, im3);

    re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
    re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
    im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
    im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

    re1 =  _mm256_add_ps(reH2, imH3); // x1 = xh2 - j*xh3;
    im1 =  _mm256_sub_ps(imH2, reH3); // x1 = xh2 - j*xh3;
    re3 =  _mm256_sub_ps(reH2, imH3); // x3 = xh2 + j*xh3;
    im3 =  _mm256_add_ps(imH2, reH3); // x3 = xh2 + j*xh3;

    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;

    _mm256_store_ps(&dst[8*0+512*2*0], re0);
    _mm256_store_ps(&dst[8*1+512*2*0], im0);

    __m256 re1w = MULnADD(im1, sin1, _mm256_mul_ps(re1, cos1)); // rew = re*cos - im*sin
    __m256 im1w = MULADD (re1, sin1, _mm256_mul_ps(im1, cos1)); // imw = im*cos + re*sin

    __m256 re2w = MULnADD(im2, sin2, _mm256_mul_ps(re2, cos2));
    __m256 im2w = MULADD (re2, sin2, _mm256_mul_ps(im2, cos2));

    __m256 re3w = MULnADD(im3, sin3, _mm256_mul_ps(re3, cos3));
    __m256 im3w = MULADD (re3, sin3, _mm256_mul_ps(im3, cos3));

    _mm256_store_ps(&dst[8*0+512*2*1], re1w);
    _mm256_store_ps(&dst[8*1+512*2*1], im1w);

    _mm256_store_ps(&dst[8*0+512*2*2], re2w);
    _mm256_store_ps(&dst[8*1+512*2*2], im2w);

    _mm256_store_ps(&dst[8*0+512*2*3], re3w);
    _mm256_store_ps(&dst[8*1+512*2*3], im3w);

    src += 8*2;
    dst += 8*2;
  }

  dst = wbuf;
  for (int i = 0; i < 2048/16/8; ++i) {
    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;
    for (int k = 0; k < 4; ++k) {
      __m256 re0 = _mm256_loadu_ps(&dst[8*0+128*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
      __m256 im0 = _mm256_loadu_ps(&dst[8*1+128*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
      __m256 re1 = _mm256_loadu_ps(&dst[8*0+128*2*1]);
      __m256 im1 = _mm256_loadu_ps(&dst[8*1+128*2*1]);
      __m256 re2 = _mm256_loadu_ps(&dst[8*0+128*2*2]);
      __m256 im2 = _mm256_loadu_ps(&dst[8*1+128*2*2]);
      __m256 re3 = _mm256_loadu_ps(&dst[8*0+128*2*3]);
      __m256 im3 = _mm256_loadu_ps(&dst[8*1+128*2*3]);

      __m256 reH0 = _mm256_add_ps(re0, re2);
      __m256 reH2 = _mm256_sub_ps(re0, re2);
      __m256 imH0 = _mm256_add_ps(im0, im2);
      __m256 imH2 = _mm256_sub_ps(im0, im2);

      __m256 reH1 = _mm256_add_ps(re1, re3);
      __m256 reH3 = _mm256_sub_ps(re1, re3);
      __m256 imH1 = _mm256_add_ps(im1, im3);
      __m256 imH3 = _mm256_sub_ps(im1, im3);

      re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
      re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
      im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
      im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

      re1 =  _mm256_add_ps(reH2, imH3); // x1 = xh2 - j*xh3;
      im1 =  _mm256_sub_ps(imH2, reH3); // x1 = xh2 - j*xh3;
      re3 =  _mm256_sub_ps(reH2, imH3); // x3 = xh2 + j*xh3;
      im3 =  _mm256_add_ps(imH2, reH3); // x3 = xh2 + j*xh3;

      _mm256_store_ps(&dst[8*0+128*2*0], re0);
      _mm256_store_ps(&dst[8*1+128*2*0], im0);

      __m256 re1w = MULnADD(im1, sin1, _mm256_mul_ps(re1, cos1));
      __m256 im1w = MULADD (re1, sin1, _mm256_mul_ps(im1, cos1));

      __m256 re2w = MULnADD(im2, sin2, _mm256_mul_ps(re2, cos2));
      __m256 im2w = MULADD (re2, sin2, _mm256_mul_ps(im2, cos2));

      __m256 re3w = MULnADD(im3, sin3, _mm256_mul_ps(re3, cos3));
      __m256 im3w = MULADD (re3, sin3, _mm256_mul_ps(im3, cos3));

      _mm256_store_ps(&dst[8*0+128*2*1], re1w);
      _mm256_store_ps(&dst[8*1+128*2*1], im1w);

      _mm256_store_ps(&dst[8*0+128*2*2], re2w);
      _mm256_store_ps(&dst[8*1+128*2*2], im2w);

      _mm256_store_ps(&dst[8*0+128*2*3], re3w);
      _mm256_store_ps(&dst[8*1+128*2*3], im3w);

      dst += 2048/4*2;
    }
    dst -= (2048-8)*2;
  }

  dst = wbuf;
  for (int i = 0; i < 2048/64/8; ++i) {
    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;
    for (int k = 0; k < 16; ++k) {
      __m256 re0 = _mm256_loadu_ps(&dst[8*0+32*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
      __m256 im0 = _mm256_loadu_ps(&dst[8*1+32*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
      __m256 re1 = _mm256_loadu_ps(&dst[8*0+32*2*1]);
      __m256 im1 = _mm256_loadu_ps(&dst[8*1+32*2*1]);
      __m256 re2 = _mm256_loadu_ps(&dst[8*0+32*2*2]);
      __m256 im2 = _mm256_loadu_ps(&dst[8*1+32*2*2]);
      __m256 re3 = _mm256_loadu_ps(&dst[8*0+32*2*3]);
      __m256 im3 = _mm256_loadu_ps(&dst[8*1+32*2*3]);

      __m256 reH0 = _mm256_add_ps(re0, re2);
      __m256 reH2 = _mm256_sub_ps(re0, re2);
      __m256 imH0 = _mm256_add_ps(im0, im2);
      __m256 imH2 = _mm256_sub_ps(im0, im2);

      __m256 reH1 = _mm256_add_ps(re1, re3);
      __m256 reH3 = _mm256_sub_ps(re1, re3);
      __m256 imH1 = _mm256_add_ps(im1, im3);
      __m256 imH3 = _mm256_sub_ps(im1, im3);

      re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
      re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
      im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
      im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

      re1 =  _mm256_add_ps(reH2, imH3); // x1 = xh2 - j*xh3;
      im1 =  _mm256_sub_ps(imH2, reH3); // x1 = xh2 - j*xh3;
      re3 =  _mm256_sub_ps(reH2, imH3); // x3 = xh2 + j*xh3;
      im3 =  _mm256_add_ps(imH2, reH3); // x3 = xh2 + j*xh3;


      _mm256_store_ps(&dst[8*0+32*2*0], re0);
      _mm256_store_ps(&dst[8*1+32*2*0], im0);

      __m256 re1w = MULnADD(im1, sin1, _mm256_mul_ps(re1, cos1));
      __m256 im1w = MULADD (re1, sin1, _mm256_mul_ps(im1, cos1));

      __m256 re2w = MULnADD(im2, sin2, _mm256_mul_ps(re2, cos2));
      __m256 im2w = MULADD (re2, sin2, _mm256_mul_ps(im2, cos2));

      __m256 re3w = MULnADD(im3, sin3, _mm256_mul_ps(re3, cos3));
      __m256 im3w = MULADD (re3, sin3, _mm256_mul_ps(im3, cos3));

      _mm256_store_ps(&dst[8*0+32*2*1], re1w);
      _mm256_store_ps(&dst[8*1+32*2*1], im1w);

      _mm256_store_ps(&dst[8*0+32*2*2], re2w);
      _mm256_store_ps(&dst[8*1+32*2*2], im2w);

      _mm256_store_ps(&dst[8*0+32*2*3], re3w);
      _mm256_store_ps(&dst[8*1+32*2*3], im3w);

      dst += 2048/16*2;
    }
    dst -= (2048-8)*2;
  }

  dst = wbuf;
  {
    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;
    for (int k = 0; k < 64; ++k) {
      __m256 re0 = _mm256_loadu_ps(&dst[8*0+8*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
      __m256 im0 = _mm256_loadu_ps(&dst[8*1+8*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
      __m256 re1 = _mm256_loadu_ps(&dst[8*0+8*2*1]);
      __m256 im1 = _mm256_loadu_ps(&dst[8*1+8*2*1]);
      __m256 re2 = _mm256_loadu_ps(&dst[8*0+8*2*2]);
      __m256 im2 = _mm256_loadu_ps(&dst[8*1+8*2*2]);
      __m256 re3 = _mm256_loadu_ps(&dst[8*0+8*2*3]);
      __m256 im3 = _mm256_loadu_ps(&dst[8*1+8*2*3]);

      __m256 reH0 = _mm256_add_ps(re0, re2);
      __m256 reH2 = _mm256_sub_ps(re0, re2);
      __m256 imH0 = _mm256_add_ps(im0, im2);
      __m256 imH2 = _mm256_sub_ps(im0, im2);

      __m256 reH1 = _mm256_add_ps(re1, re3);
      __m256 reH3 = _mm256_sub_ps(re1, re3);
      __m256 imH1 = _mm256_add_ps(im1, im3);
      __m256 imH3 = _mm256_sub_ps(im1, im3);

      re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
      re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
      im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
      im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

      re1 =  _mm256_add_ps(reH2, imH3); // x1 = xh2 - j*xh3;
      im1 =  _mm256_sub_ps(imH2, reH3); // x1 = xh2 - j*xh3;
      re3 =  _mm256_sub_ps(reH2, imH3); // x3 = xh2 + j*xh3;
      im3 =  _mm256_add_ps(imH2, reH3); // x3 = xh2 + j*xh3;

      _mm256_store_ps(&dst[8*0+8*2*0], re0);
      _mm256_store_ps(&dst[8*1+8*2*0], im0);

      __m256 re1w = MULnADD(im1, sin1, _mm256_mul_ps(re1, cos1));
      __m256 im1w = MULADD (re1, sin1, _mm256_mul_ps(im1, cos1));

      __m256 re2w = MULnADD(im2, sin2, _mm256_mul_ps(re2, cos2));
      __m256 im2w = MULADD (re2, sin2, _mm256_mul_ps(im2, cos2));

      __m256 re3w = MULnADD(im3, sin3, _mm256_mul_ps(re3, cos3));
      __m256 im3w = MULADD (re3, sin3, _mm256_mul_ps(im3, cos3));

      _mm256_store_ps(&dst[8*0+8*2*1], re1w);
      _mm256_store_ps(&dst[8*1+8*2*1], im1w);

      _mm256_store_ps(&dst[8*0+8*2*2], re2w);
      _mm256_store_ps(&dst[8*1+8*2*2], im2w);

      _mm256_store_ps(&dst[8*0+8*2*3], re3w);
      _mm256_store_ps(&dst[8*1+8*2*3], im3w);

      dst += 2048/64*2;
    }
  }

  src = wbuf;
  __m256 cs  = _mm256_broadcast_ss(pTab);              // sqrt(0.5)
  __m256 csN = _mm256_sub_ps(_mm256_setzero_ps(), cs); // -sqrt(0.5)
  for (int i = 0; i < 64; ++i) {
    __m256 re0 = _mm256_loadu_ps(&src[8*0+512*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
    __m256 im0 = _mm256_loadu_ps(&src[8*1+512*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
    __m256 re1 = _mm256_loadu_ps(&src[8*0+512*2*1]);
    __m256 im1 = _mm256_loadu_ps(&src[8*1+512*2*1]);
    __m256 re2 = _mm256_loadu_ps(&src[8*0+512*2*2]);
    __m256 im2 = _mm256_loadu_ps(&src[8*1+512*2*2]);
    __m256 re3 = _mm256_loadu_ps(&src[8*0+512*2*3]);
    __m256 im3 = _mm256_loadu_ps(&src[8*1+512*2*3]);
    dst = &srcDst[(ptrdiff_t)nibrev_tab[i] * 4*2];

    __m256 re01a = _mm256_unpacklo_ps(re0, re1); // re0[0] re0[1] re1[0] re1[1] re2[0] re2[1] re3[0] re3[1]
    __m256 re01b = _mm256_unpackhi_ps(re0, re1); // re4[0] re4[1] re5[0] re5[1] re6[0] re6[1] re7[0] re7[1]
    __m256 re23a = _mm256_unpacklo_ps(re2, re3); // re0[2] re0[3] re1[2] re1[3] re2[2] re2[3] re3[2] re3[3]
    __m256 re23b = _mm256_unpackhi_ps(re2, re3); // re4[2] re4[3] re5[2] re5[3] re6[2] re6[3] re7[2] re7[3]
    __m256 im01a = _mm256_unpacklo_ps(im0, im1); // im0[0] im0[1] im1[0] im1[1] im2[0] im2[1] im3[0] im3[1]
    __m256 im01b = _mm256_unpackhi_ps(im0, im1); // im4[0] im4[1] im5[0] im5[1] im6[0] im6[1] im7[0] im7[1]
    __m256 im23a = _mm256_unpacklo_ps(im2, im3); // im0[2] im0[3] im1[2] im1[3] im2[2] im2[3] im3[2] im3[3]
    __m256 im23b = _mm256_unpackhi_ps(im2, im3); // im4[2] im4[3] im5[2] im5[3] im6[2] im6[3] im7[2] im7[3]

    __m256 re01 = _mm256_permute2f128_ps(re01a, re23a, 0x20); // re0[0] re0[1] re1[0] re1[1] re0[2] re0[3] re1[2] re1[3]
    __m256 re23 = _mm256_permute2f128_ps(re01a, re23a, 0x31); // re2[0] re2[1] re3[0] re3[1] re2[2] re2[3] re3[2] re3[3]
    __m256 re45 = _mm256_permute2f128_ps(re01b, re23b, 0x20); // re4[0] re4[1] re5[0] re5[1] re4[2] re4[3] re5[2] re5[3]
    __m256 re67 = _mm256_permute2f128_ps(re01b, re23b, 0x31); // re6[0] re6[1] re7[0] re7[1] re6[2] re6[3] re7[2] re7[3]
    __m256 im01 = _mm256_permute2f128_ps(im01a, im23a, 0x20); // im0[0] im0[1] im1[0] im1[1] im0[2] im0[3] im1[2] im1[3]
    __m256 im23 = _mm256_permute2f128_ps(im01a, im23a, 0x31); // im2[0] im2[1] im3[0] im3[1] im2[2] im2[3] im3[2] im3[3]
    __m256 im45 = _mm256_permute2f128_ps(im01b, im23b, 0x20); // im4[0] im4[1] im5[0] im5[1] im4[2] im4[3] im5[2] im5[3]
    __m256 im67 = _mm256_permute2f128_ps(im01b, im23b, 0x31); // im6[0] im6[1] im7[0] im7[1] im6[2] im6[3] im7[2] im7[3]

    __m256 reH0 = _mm256_add_ps(re01, re45);
    __m256 reH2 = _mm256_sub_ps(re01, re45);
    __m256 imH0 = _mm256_add_ps(im01, im45);
    __m256 imH2 = _mm256_sub_ps(im01, im45);

    __m256 reH1 = _mm256_add_ps(re23, re67);
    __m256 reH3 = _mm256_sub_ps(re23, re67);
    __m256 imH1 = _mm256_add_ps(im23, im67);
    __m256 imH3 = _mm256_sub_ps(im23, im67);

    re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
    re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
    im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
    im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

    re1 =  _mm256_add_ps(reH2, imH3); // x1 = xh2 - j*xh3;
    im1 =  _mm256_sub_ps(imH2, reH3); // x1 = xh2 - j*xh3;
    re3 =  _mm256_sub_ps(reH2, imH3); // x3 = xh2 + j*xh3;
    im3 =  _mm256_add_ps(imH2, reH3); // x3 = xh2 + j*xh3;

    __m256 c0 = _mm256_unpacklo_ps(re0, im0); // re0[0] im0[0] re0[1] im0[1] re0[2] im0[2] re0[3] im0[3]
    __m256 c1 = _mm256_unpackhi_ps(re0, im0); // re1[0] im1[0] re1[1] im1[1] re1[2] im1[2] re1[3] im1[3]
    __m256 res0 = _mm256_add_ps(c0, c1);
    __m256 res1 = _mm256_sub_ps(c0, c1);
    _mm256_store_ps(&dst[1024*2*0+256*2*0], res0);
    _mm256_store_ps(&dst[1024*2*1+256*2*0], res1);

    // w1^1 = sqrt(0.5) -j*sqrt(0.5)
    __m256 re1w = _mm256_mul_ps(_mm256_add_ps(re1, im1), cs);
    __m256 im1w = _mm256_mul_ps(_mm256_sub_ps(im1, re1), cs);
    __m256 c2 = _mm256_unpacklo_ps(re1,  im1);  // re2[0] im2[0] re2[1] im2[1] re2[2] im2[2] re2[3] im2[3]
    __m256 c3 = _mm256_unpackhi_ps(re1w, im1w); // re3[0] im3[0] re3[1] im3[1] re3[2] im3[2] re3[3] im3[3]
    __m256 res2 = _mm256_add_ps(c2, c3);
    __m256 res3 = _mm256_sub_ps(c2, c3);

    // w1^2 = -j
    __m256 re2w = im2;
    __m256 im2w = _mm256_sub_ps(_mm256_setzero_ps(), re2);
    __m256 c4 = _mm256_unpacklo_ps(re2,  im2);  // re4[0] im4[0] re4[1] im4[1] re4[2] im4[2] re4[3] im4[3]
    __m256 c5 = _mm256_unpackhi_ps(re2w, im2w); // re5[0] im5[0] re5[1] im5[1] re5[2] im5[2] re5[3] im5[3]
    __m256 res4 = _mm256_add_ps(c4, c5);
    __m256 res5 = _mm256_sub_ps(c4, c5);

    // w1^3 = -sqrt(0.5) -j*sqrt(0.5)
    __m256 re3w = _mm256_mul_ps(_mm256_sub_ps(im3, re3), cs);
    __m256 im3w = _mm256_mul_ps(_mm256_add_ps(re3, im3), csN);
    __m256 c6 = _mm256_unpacklo_ps(re3,  im3);  // re6[0] im6[0] re6[1] im6[1] re6[2] im6[2] re6[3] im6[3]
    __m256 c7 = _mm256_unpackhi_ps(re3w, im3w); // re7[0] im7[0] re7[1] im7[1] re7[2] im7[2] re7[3] im7[3]
    __m256 res6 = _mm256_add_ps(c6, c7);
    __m256 res7 = _mm256_sub_ps(c6, c7);

    _mm256_store_ps(&dst[1024*2*0+256*2*1], res2);
    _mm256_store_ps(&dst[1024*2*1+256*2*1], res3);

    _mm256_store_ps(&dst[1024*2*0+256*2*2], res4);
    _mm256_store_ps(&dst[1024*2*1+256*2*2], res5);

    _mm256_store_ps(&dst[1024*2*0+256*2*3], res6);
    _mm256_store_ps(&dst[1024*2*1+256*2*3], res7);

    src += 8*2;
  }

  // memcpy(srcDst, wbuf, sizeof(wbuf));
}

void cifft2048_r4dif(float srcDst[2048*2])
{
#ifdef _MSC_VER
  __declspec(align(32))
#endif
  float wbuf[2048*2]
#ifdef __GNUC__
  __attribute__ ((aligned (32)))
#endif
  ;
  const float* pTab = sincos_tab;
  float* dst = wbuf;
  float* src = srcDst;
  for (int i = 0; i < 2048/4/8; ++i) {
    __m256 c00 = _mm256_loadu_ps(&src[8*0+512*2*0]); // re0 im0 re1 im1 re2 im2 re3 im3
    __m256 c01 = _mm256_loadu_ps(&src[8*1+512*2*0]); // re4 im4 re5 im5 re6 im6 re7 im7
    __m256 c10 = _mm256_loadu_ps(&src[8*0+512*2*1]);
    __m256 c11 = _mm256_loadu_ps(&src[8*1+512*2*1]);
    __m256 c20 = _mm256_loadu_ps(&src[8*0+512*2*2]);
    __m256 c21 = _mm256_loadu_ps(&src[8*1+512*2*2]);
    __m256 c30 = _mm256_loadu_ps(&src[8*0+512*2*3]);
    __m256 c31 = _mm256_loadu_ps(&src[8*1+512*2*3]);

    __m256 re0 = _mm256_shuffle_ps(c00, c01, 0x88); // re0 re1 re4 re5  re2 re3 re6 re7
    __m256 im0 = _mm256_shuffle_ps(c00, c01, 0xDD); // im0 im1 im4 im5  im2 im3 im6 im7
    __m256 re1 = _mm256_shuffle_ps(c10, c11, 0x88);
    __m256 im1 = _mm256_shuffle_ps(c10, c11, 0xDD);
    __m256 re2 = _mm256_shuffle_ps(c20, c21, 0x88);
    __m256 im2 = _mm256_shuffle_ps(c20, c21, 0xDD);
    __m256 re3 = _mm256_shuffle_ps(c30, c31, 0x88);
    __m256 im3 = _mm256_shuffle_ps(c30, c31, 0xDD);

    __m256 reH0 = _mm256_add_ps(re0, re2);
    __m256 reH2 = _mm256_sub_ps(re0, re2);
    __m256 imH0 = _mm256_add_ps(im0, im2);
    __m256 imH2 = _mm256_sub_ps(im0, im2);

    __m256 reH1 = _mm256_add_ps(re1, re3);
    __m256 reH3 = _mm256_sub_ps(re1, re3);
    __m256 imH1 = _mm256_add_ps(im1, im3);
    __m256 imH3 = _mm256_sub_ps(im1, im3);

    re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
    re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
    im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
    im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

    re1 =  _mm256_sub_ps(reH2, imH3); // x1 = xh2 + j*xh3;
    im1 =  _mm256_add_ps(imH2, reH3); // x1 = xh2 + j*xh3;
    re3 =  _mm256_add_ps(reH2, imH3); // x3 = xh2 - j*xh3;
    im3 =  _mm256_sub_ps(imH2, reH3); // x3 = xh2 - j*xh3;

    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;

    _mm256_store_ps(&dst[8*0+512*2*0], re0);
    _mm256_store_ps(&dst[8*1+512*2*0], im0);

    __m256 re1w = MULADD (im1, sin1, _mm256_mul_ps(re1, cos1));
    __m256 im1w = MULnADD(re1, sin1, _mm256_mul_ps(im1, cos1));

    __m256 re2w = MULADD (im2, sin2, _mm256_mul_ps(re2, cos2));
    __m256 im2w = MULnADD(re2, sin2, _mm256_mul_ps(im2, cos2));

    __m256 re3w = MULADD (im3, sin3, _mm256_mul_ps(re3, cos3));
    __m256 im3w = MULnADD(re3, sin3, _mm256_mul_ps(im3, cos3));

    _mm256_store_ps(&dst[8*0+512*2*1], re1w);
    _mm256_store_ps(&dst[8*1+512*2*1], im1w);

    _mm256_store_ps(&dst[8*0+512*2*2], re2w);
    _mm256_store_ps(&dst[8*1+512*2*2], im2w);

    _mm256_store_ps(&dst[8*0+512*2*3], re3w);
    _mm256_store_ps(&dst[8*1+512*2*3], im3w);

    src += 8*2;
    dst += 8*2;
  }

  dst = wbuf;
  for (int i = 0; i < 2048/16/8; ++i) {
    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;
    for (int k = 0; k < 4; ++k) {
      __m256 re0 = _mm256_loadu_ps(&dst[8*0+128*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
      __m256 im0 = _mm256_loadu_ps(&dst[8*1+128*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
      __m256 re1 = _mm256_loadu_ps(&dst[8*0+128*2*1]);
      __m256 im1 = _mm256_loadu_ps(&dst[8*1+128*2*1]);
      __m256 re2 = _mm256_loadu_ps(&dst[8*0+128*2*2]);
      __m256 im2 = _mm256_loadu_ps(&dst[8*1+128*2*2]);
      __m256 re3 = _mm256_loadu_ps(&dst[8*0+128*2*3]);
      __m256 im3 = _mm256_loadu_ps(&dst[8*1+128*2*3]);

      __m256 reH0 = _mm256_add_ps(re0, re2);
      __m256 reH2 = _mm256_sub_ps(re0, re2);
      __m256 imH0 = _mm256_add_ps(im0, im2);
      __m256 imH2 = _mm256_sub_ps(im0, im2);

      __m256 reH1 = _mm256_add_ps(re1, re3);
      __m256 reH3 = _mm256_sub_ps(re1, re3);
      __m256 imH1 = _mm256_add_ps(im1, im3);
      __m256 imH3 = _mm256_sub_ps(im1, im3);

      re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
      re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
      im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
      im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

      re1 =  _mm256_sub_ps(reH2, imH3); // x1 = xh2 + j*xh3;
      im1 =  _mm256_add_ps(imH2, reH3); // x1 = xh2 + j*xh3;
      re3 =  _mm256_add_ps(reH2, imH3); // x3 = xh2 - j*xh3;
      im3 =  _mm256_sub_ps(imH2, reH3); // x3 = xh2 - j*xh3;

      _mm256_store_ps(&dst[8*0+128*2*0], re0);
      _mm256_store_ps(&dst[8*1+128*2*0], im0);

      __m256 re1w = MULADD (im1, sin1, _mm256_mul_ps(re1, cos1));
      __m256 im1w = MULnADD(re1, sin1, _mm256_mul_ps(im1, cos1));

      __m256 re2w = MULADD (im2, sin2, _mm256_mul_ps(re2, cos2));
      __m256 im2w = MULnADD(re2, sin2, _mm256_mul_ps(im2, cos2));

      __m256 re3w = MULADD (im3, sin3, _mm256_mul_ps(re3, cos3));
      __m256 im3w = MULnADD(re3, sin3, _mm256_mul_ps(im3, cos3));

      _mm256_store_ps(&dst[8*0+128*2*1], re1w);
      _mm256_store_ps(&dst[8*1+128*2*1], im1w);

      _mm256_store_ps(&dst[8*0+128*2*2], re2w);
      _mm256_store_ps(&dst[8*1+128*2*2], im2w);

      _mm256_store_ps(&dst[8*0+128*2*3], re3w);
      _mm256_store_ps(&dst[8*1+128*2*3], im3w);

      dst += 2048/4*2;
    }
    dst -= (2048-8)*2;
  }

  dst = wbuf;
  for (int i = 0; i < 2048/64/8; ++i) {
    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;
    for (int k = 0; k < 16; ++k) {
      __m256 re0 = _mm256_loadu_ps(&dst[8*0+32*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
      __m256 im0 = _mm256_loadu_ps(&dst[8*1+32*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
      __m256 re1 = _mm256_loadu_ps(&dst[8*0+32*2*1]);
      __m256 im1 = _mm256_loadu_ps(&dst[8*1+32*2*1]);
      __m256 re2 = _mm256_loadu_ps(&dst[8*0+32*2*2]);
      __m256 im2 = _mm256_loadu_ps(&dst[8*1+32*2*2]);
      __m256 re3 = _mm256_loadu_ps(&dst[8*0+32*2*3]);
      __m256 im3 = _mm256_loadu_ps(&dst[8*1+32*2*3]);

      __m256 reH0 = _mm256_add_ps(re0, re2);
      __m256 reH2 = _mm256_sub_ps(re0, re2);
      __m256 imH0 = _mm256_add_ps(im0, im2);
      __m256 imH2 = _mm256_sub_ps(im0, im2);

      __m256 reH1 = _mm256_add_ps(re1, re3);
      __m256 reH3 = _mm256_sub_ps(re1, re3);
      __m256 imH1 = _mm256_add_ps(im1, im3);
      __m256 imH3 = _mm256_sub_ps(im1, im3);

      re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
      re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
      im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
      im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

      re1 =  _mm256_sub_ps(reH2, imH3); // x1 = xh2 + j*xh3;
      im1 =  _mm256_add_ps(imH2, reH3); // x1 = xh2 + j*xh3;
      re3 =  _mm256_add_ps(reH2, imH3); // x3 = xh2 - j*xh3;
      im3 =  _mm256_sub_ps(imH2, reH3); // x3 = xh2 - j*xh3;

      _mm256_store_ps(&dst[8*0+32*2*0], re0);
      _mm256_store_ps(&dst[8*1+32*2*0], im0);

      __m256 re1w = MULADD (im1, sin1, _mm256_mul_ps(re1, cos1));
      __m256 im1w = MULnADD(re1, sin1, _mm256_mul_ps(im1, cos1));

      __m256 re2w = MULADD (im2, sin2, _mm256_mul_ps(re2, cos2));
      __m256 im2w = MULnADD(re2, sin2, _mm256_mul_ps(im2, cos2));

      __m256 re3w = MULADD (im3, sin3, _mm256_mul_ps(re3, cos3));
      __m256 im3w = MULnADD(re3, sin3, _mm256_mul_ps(im3, cos3));

      _mm256_store_ps(&dst[8*0+32*2*1], re1w);
      _mm256_store_ps(&dst[8*1+32*2*1], im1w);

      _mm256_store_ps(&dst[8*0+32*2*2], re2w);
      _mm256_store_ps(&dst[8*1+32*2*2], im2w);

      _mm256_store_ps(&dst[8*0+32*2*3], re3w);
      _mm256_store_ps(&dst[8*1+32*2*3], im3w);

      dst += 2048/16*2;
    }
    dst -= (2048-8)*2;
  }

  dst = wbuf;
  {
    __m256 cos1 = _mm256_load_ps(&pTab[0*8]);
    __m256 sin1 = _mm256_load_ps(&pTab[1*8]);
    __m256 cos2 = _mm256_load_ps(&pTab[2*8]);
    __m256 sin2 = _mm256_load_ps(&pTab[3*8]);
    __m256 cos3 = _mm256_load_ps(&pTab[4*8]);
    __m256 sin3 = _mm256_load_ps(&pTab[5*8]);
    pTab += 6*8;
    for (int k = 0; k < 64; ++k) {
      __m256 re0 = _mm256_loadu_ps(&dst[8*0+8*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
      __m256 im0 = _mm256_loadu_ps(&dst[8*1+8*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
      __m256 re1 = _mm256_loadu_ps(&dst[8*0+8*2*1]);
      __m256 im1 = _mm256_loadu_ps(&dst[8*1+8*2*1]);
      __m256 re2 = _mm256_loadu_ps(&dst[8*0+8*2*2]);
      __m256 im2 = _mm256_loadu_ps(&dst[8*1+8*2*2]);
      __m256 re3 = _mm256_loadu_ps(&dst[8*0+8*2*3]);
      __m256 im3 = _mm256_loadu_ps(&dst[8*1+8*2*3]);

      __m256 reH0 = _mm256_add_ps(re0, re2);
      __m256 reH2 = _mm256_sub_ps(re0, re2);
      __m256 imH0 = _mm256_add_ps(im0, im2);
      __m256 imH2 = _mm256_sub_ps(im0, im2);

      __m256 reH1 = _mm256_add_ps(re1, re3);
      __m256 reH3 = _mm256_sub_ps(re1, re3);
      __m256 imH1 = _mm256_add_ps(im1, im3);
      __m256 imH3 = _mm256_sub_ps(im1, im3);

      re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
      re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
      im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
      im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

      re1 =  _mm256_sub_ps(reH2, imH3); // x1 = xh2 + j*xh3;
      im1 =  _mm256_add_ps(imH2, reH3); // x1 = xh2 + j*xh3;
      re3 =  _mm256_add_ps(reH2, imH3); // x3 = xh2 - j*xh3;
      im3 =  _mm256_sub_ps(imH2, reH3); // x3 = xh2 - j*xh3;

      _mm256_store_ps(&dst[8*0+8*2*0], re0);
      _mm256_store_ps(&dst[8*1+8*2*0], im0);

      __m256 re1w = MULADD (im1, sin1, _mm256_mul_ps(re1, cos1));
      __m256 im1w = MULnADD(re1, sin1, _mm256_mul_ps(im1, cos1));

      __m256 re2w = MULADD (im2, sin2, _mm256_mul_ps(re2, cos2));
      __m256 im2w = MULnADD(re2, sin2, _mm256_mul_ps(im2, cos2));

      __m256 re3w = MULADD (im3, sin3, _mm256_mul_ps(re3, cos3));
      __m256 im3w = MULnADD(re3, sin3, _mm256_mul_ps(im3, cos3));

      _mm256_store_ps(&dst[8*0+8*2*1], re1w);
      _mm256_store_ps(&dst[8*1+8*2*1], im1w);

      _mm256_store_ps(&dst[8*0+8*2*2], re2w);
      _mm256_store_ps(&dst[8*1+8*2*2], im2w);

      _mm256_store_ps(&dst[8*0+8*2*3], re3w);
      _mm256_store_ps(&dst[8*1+8*2*3], im3w);

      dst += 2048/64*2;
    }
  }

  src = wbuf;
  __m256 cs  = _mm256_broadcast_ss(pTab);              // sqrt(0.5)
  __m256 csN = _mm256_sub_ps(_mm256_setzero_ps(), cs); // -sqrt(0.5)
  for (int i = 0; i < 64; ++i) {
    __m256 re0 = _mm256_loadu_ps(&src[8*0+512*2*0]); // re0 re1 re4 re5  re2 re3 re6 re7
    __m256 im0 = _mm256_loadu_ps(&src[8*1+512*2*0]); // im0 im1 im4 im5  im2 im3 im6 im7
    __m256 re1 = _mm256_loadu_ps(&src[8*0+512*2*1]);
    __m256 im1 = _mm256_loadu_ps(&src[8*1+512*2*1]);
    __m256 re2 = _mm256_loadu_ps(&src[8*0+512*2*2]);
    __m256 im2 = _mm256_loadu_ps(&src[8*1+512*2*2]);
    __m256 re3 = _mm256_loadu_ps(&src[8*0+512*2*3]);
    __m256 im3 = _mm256_loadu_ps(&src[8*1+512*2*3]);
    dst = &srcDst[(ptrdiff_t)nibrev_tab[i] * 4*2];

    __m256 re01a = _mm256_unpacklo_ps(re0, re1); // re0[0] re0[1] re1[0] re1[1] re2[0] re2[1] re3[0] re3[1]
    __m256 re01b = _mm256_unpackhi_ps(re0, re1); // re4[0] re4[1] re5[0] re5[1] re6[0] re6[1] re7[0] re7[1]
    __m256 re23a = _mm256_unpacklo_ps(re2, re3); // re0[2] re0[3] re1[2] re1[3] re2[2] re2[3] re3[2] re3[3]
    __m256 re23b = _mm256_unpackhi_ps(re2, re3); // re4[2] re4[3] re5[2] re5[3] re6[2] re6[3] re7[2] re7[3]
    __m256 im01a = _mm256_unpacklo_ps(im0, im1); // im0[0] im0[1] im1[0] im1[1] im2[0] im2[1] im3[0] im3[1]
    __m256 im01b = _mm256_unpackhi_ps(im0, im1); // im4[0] im4[1] im5[0] im5[1] im6[0] im6[1] im7[0] im7[1]
    __m256 im23a = _mm256_unpacklo_ps(im2, im3); // im0[2] im0[3] im1[2] im1[3] im2[2] im2[3] im3[2] im3[3]
    __m256 im23b = _mm256_unpackhi_ps(im2, im3); // im4[2] im4[3] im5[2] im5[3] im6[2] im6[3] im7[2] im7[3]

    __m256 re01 = _mm256_permute2f128_ps(re01a, re23a, 0x20); // re0[0] re0[1] re1[0] re1[1] re0[2] re0[3] re1[2] re1[3]
    __m256 re23 = _mm256_permute2f128_ps(re01a, re23a, 0x31); // re2[0] re2[1] re3[0] re3[1] re2[2] re2[3] re3[2] re3[3]
    __m256 re45 = _mm256_permute2f128_ps(re01b, re23b, 0x20); // re4[0] re4[1] re5[0] re5[1] re4[2] re4[3] re5[2] re5[3]
    __m256 re67 = _mm256_permute2f128_ps(re01b, re23b, 0x31); // re6[0] re6[1] re7[0] re7[1] re6[2] re6[3] re7[2] re7[3]
    __m256 im01 = _mm256_permute2f128_ps(im01a, im23a, 0x20); // im0[0] im0[1] im1[0] im1[1] im0[2] im0[3] im1[2] im1[3]
    __m256 im23 = _mm256_permute2f128_ps(im01a, im23a, 0x31); // im2[0] im2[1] im3[0] im3[1] im2[2] im2[3] im3[2] im3[3]
    __m256 im45 = _mm256_permute2f128_ps(im01b, im23b, 0x20); // im4[0] im4[1] im5[0] im5[1] im4[2] im4[3] im5[2] im5[3]
    __m256 im67 = _mm256_permute2f128_ps(im01b, im23b, 0x31); // im6[0] im6[1] im7[0] im7[1] im6[2] im6[3] im7[2] im7[3]

    __m256 reH0 = _mm256_add_ps(re01, re45);
    __m256 reH2 = _mm256_sub_ps(re01, re45);
    __m256 imH0 = _mm256_add_ps(im01, im45);
    __m256 imH2 = _mm256_sub_ps(im01, im45);

    __m256 reH1 = _mm256_add_ps(re23, re67);
    __m256 reH3 = _mm256_sub_ps(re23, re67);
    __m256 imH1 = _mm256_add_ps(im23, im67);
    __m256 imH3 = _mm256_sub_ps(im23, im67);

    re0 =  _mm256_add_ps(reH0, reH1); // x0 = xh0 +   xh1;
    re2 =  _mm256_sub_ps(reH0, reH1); // x2 = xh0 -   xh1;
    im0 =  _mm256_add_ps(imH0, imH1); // x0 = xh0 +   xh1;
    im2 =  _mm256_sub_ps(imH0, imH1); // x2 = xh0 -   xh1;

    re1 =  _mm256_sub_ps(reH2, imH3); // x1 = xh2 + j*xh3;
    im1 =  _mm256_add_ps(imH2, reH3); // x1 = xh2 + j*xh3;
    re3 =  _mm256_add_ps(reH2, imH3); // x3 = xh2 - j*xh3;
    im3 =  _mm256_sub_ps(imH2, reH3); // x3 = xh2 - j*xh3;

    __m256 c0 = _mm256_unpacklo_ps(re0, im0); // re0[0] im0[0] re0[1] im0[1] re0[2] im0[2] re0[3] im0[3]
    __m256 c1 = _mm256_unpackhi_ps(re0, im0); // re1[0] im1[0] re1[1] im1[1] re1[2] im1[2] re1[3] im1[3]
    __m256 res0 = _mm256_add_ps(c0, c1);
    __m256 res1 = _mm256_sub_ps(c0, c1);
    _mm256_store_ps(&dst[1024*2*0+256*2*0], res0);
    _mm256_store_ps(&dst[1024*2*1+256*2*0], res1);

    // w1^1 = sqrt(0.5) + j*sqrt(0.5)
    __m256 re1w = _mm256_mul_ps(_mm256_sub_ps(re1, im1), cs);
    __m256 im1w = _mm256_mul_ps(_mm256_add_ps(im1, re1), cs);
    __m256 c2 = _mm256_unpacklo_ps(re1,  im1);  // re2[0] im2[0] re2[1] im2[1] re2[2] im2[2] re2[3] im2[3]
    __m256 c3 = _mm256_unpackhi_ps(re1w, im1w); // re3[0] im3[0] re3[1] im3[1] re3[2] im3[2] re3[3] im3[3]
    __m256 res2 = _mm256_add_ps(c2, c3);
    __m256 res3 = _mm256_sub_ps(c2, c3);

    // w1^2 = j
    __m256 im2w = re2;
    __m256 re2w = _mm256_sub_ps(_mm256_setzero_ps(), im2);
    __m256 c4 = _mm256_unpacklo_ps(re2,  im2);  // re4[0] im4[0] re4[1] im4[1] re4[2] im4[2] re4[3] im4[3]
    __m256 c5 = _mm256_unpackhi_ps(re2w, im2w); // re5[0] im5[0] re5[1] im5[1] re5[2] im5[2] re5[3] im5[3]
    __m256 res4 = _mm256_add_ps(c4, c5);
    __m256 res5 = _mm256_sub_ps(c4, c5);

    // w1^3 = -sqrt(0.5) + j*sqrt(0.5)
    __m256 re3w = _mm256_mul_ps(_mm256_add_ps(re3, im3), csN);
    __m256 im3w = _mm256_mul_ps(_mm256_sub_ps(re3, im3), cs);
    __m256 c6 = _mm256_unpacklo_ps(re3,  im3);  // re6[0] im6[0] re6[1] im6[1] re6[2] im6[2] re6[3] im6[3]
    __m256 c7 = _mm256_unpackhi_ps(re3w, im3w); // re7[0] im7[0] re7[1] im7[1] re7[2] im7[2] re7[3] im7[3]
    __m256 res6 = _mm256_add_ps(c6, c7);
    __m256 res7 = _mm256_sub_ps(c6, c7);

    _mm256_store_ps(&dst[1024*2*0+256*2*1], res2);
    _mm256_store_ps(&dst[1024*2*1+256*2*1], res3);

    _mm256_store_ps(&dst[1024*2*0+256*2*2], res4);
    _mm256_store_ps(&dst[1024*2*1+256*2*2], res5);

    _mm256_store_ps(&dst[1024*2*0+256*2*3], res6);
    _mm256_store_ps(&dst[1024*2*1+256*2*3], res7);

    src += 8*2;
  }

  // memcpy(srcDst, wbuf, sizeof(wbuf));
}

static void setSinCos(float* dst, int x)
{ // quarter-aware sin/cos is a little more precise
  // ideally, we should use sinPi/cosPi, but it's not available in out standard library
  if (x < 512) {
    dst[0] = (float) cos(M_PI*x/1024);
    dst[8] = (float)-sin(M_PI*x/1024);
  } else if (x < 1024) {
    dst[0] = (float)-sin(M_PI*(x-512)/1024);
    dst[8] = (float)-cos(M_PI*(x-512)/1024);
  } else {
    dst[0] = (float)-cos(M_PI*(x-1024)/1024);
    dst[8] = (float) sin(M_PI*(x-1024)/1024);
  }
}

void cfft2048_r4dif_int(void)
{
  static const char il_tab[8] = { 0, 1, 4, 5, 2, 3, 6, 7};
  float* dst = sincos_tab;
  for (int pass = 0, tab_step = 1; pass < 4; ++pass) {
    int i_max = 512/tab_step;
    for (int i = 0; i < i_max; ++i) {
      int ih = i / 8;
      int il = il_tab[i % 8];
      for (int k = 0; k < 3; ++k) {
        int x = (i*(k+1)*tab_step) %  2048;
        int idx = (ih * 3 + k)*8*2 + il;
        setSinCos(&dst[idx], x);
      }
    }
    dst += i_max*2*3;
    tab_step *= 4;
  }
  dst[0] = (float)sqrt(0.5);

  for (unsigned i = 0; i < 64; ++i) {
    nibrev_tab[i] = (unsigned char)((i % 4) * 16 + ((i/4) % 4)*4 + i/16);
  }
}