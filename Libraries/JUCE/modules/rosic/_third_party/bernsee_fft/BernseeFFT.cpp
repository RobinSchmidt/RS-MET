#include "BernseeFFT.h"

void smbFft(float *fftBuffer, long fftFrameSize, long sign)
{
  float wr, wi, arg, *p1, *p2, temp;
  float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
  long i, bitm, j, le, le2, k, logN;
  logN = (long)(log((double)fftFrameSize)/log(2.)+.5);

  for (i = 2; i < 2*fftFrameSize-2; i += 2) {

    for (bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {

      if (i & bitm) j++;
      j <<= 1;

    }

    if (i < j) {

      p1 = fftBuffer+i; p2 = fftBuffer+j;
      temp = *p1; *(p1++) = *p2;
      *(p2++) = temp; temp = *p1;
      *p1 = *p2; *p2 = temp;

    }

  }

  for (k = 0, le = 2; k < logN; k++) {

    le <<= 1;
    le2 = le>>1;
    ur = 1.0;
    ui = 0.0;
    arg = (float) M_PI / (le2>>1);
    wr = cos(arg);
    wi = sign*sin(arg);

    for (j = 0; j < le2; j += 2) {

      p1r = fftBuffer+j; p1i = p1r+1;
      p2r = p1r+le2; p2i = p2r+1;

      for (i = j; i < 2*fftFrameSize; i += le) {

        tr = *p2r * ur - *p2i * ui;
        ti = *p2r * ui + *p2i * ur;
        *p2r = *p1r - tr; *p2i = *p1i - ti;
        *p1r += tr; *p1i += ti;
        p1r += le; p1i += le;
        p2r += le; p2i += le;

      }

      tr = ur*wr - ui*wi;
      ui = ur*wi + ui*wr;
      ur = tr;

    }

  }

}