//#include "rosic_Transformations.h"
//using namespace rosic;

// We need to include the kiss-fft code here
#include "../_third_party/kiss_fft_v1_2_6/kiss_fft.c"

void rosic::discreteFourierTransform(Complex *in, Complex *out, int length, bool isInverse)
{
  // memory allocation:
  kiss_fft_cfg cfg     = kiss_fft_alloc(length, isInverse, 0, 0);
  kiss_fft_cpx* cx_in  = new kiss_fft_cpx[length];
  kiss_fft_cpx* cx_out = new kiss_fft_cpx[length];

  // convert rosic::Complex into the kiss_fft format:
  int k;
  for(k=0; k<length; k++)
  {
    // put kth sample in cx_in[k].r and cx_in[k].i:
    cx_in[k].r = (kiss_fft_scalar) in[k].re;
    cx_in[k].i = (kiss_fft_scalar) in[k].im;
  }

  // do the actual fft:
  kiss_fft(cfg, cx_in, cx_out);

  // convert back to rosic::Complex:
  for(k=0; k<length; k++)
  {
    // transformed. DC is in cx_out[0].r and cx_out[0].i
    out[k].re = cx_out[k].r;
    out[k].im = cx_out[k].i;
  }

  // free memeory:
  delete[] cx_in;
  delete[] cx_out;
  free(cfg);
}

void rosic::radixTwoFastFourierTransform(Complex *in, Complex *out, int length, bool isInverse)
{
  // setup the internal variables for Stephan Bernsee's routine:
  //#define M_PI 3.14159265358979323846
  long    fftFrameSize = length;
  double* fftBuffer    = new double[2*fftFrameSize];
  double  sign;
  if( isInverse == true )
    sign = -1.0;
  else 
    sign = +1.0;

  // copy the complex input data into the fftBuffer:
  int n;
  for(n=0; n<length; n++)
  {
    fftBuffer[2*n]   = in[n].re;
    fftBuffer[2*n+1] = in[n].im;
  }

  //---------------------------------------------------------------------------
  // the FFT-routine as it is published in DFT a pied (with some minor tweaks):

  double wr, wi, arg, *p1, *p2, temp;
  double tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
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
    arg = PI / (le2>>1);
    wr = cos(arg);
    wi = sign*sin(arg);
    //rosic::sinCos(arg, &wi, &wr);
    //wi *= sign;

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

  //---------------------------------------------------------------------------
  // the routine has finished
  
  // copy the result into rosic::Complex output buffer
  for(n=0; n<length; n++)
  {
    out[n].re = fftBuffer[2*n];
    out[n].im = fftBuffer[2*n+1];
  }

  // free temporarily allocated memory:
  delete[] fftBuffer;
}


