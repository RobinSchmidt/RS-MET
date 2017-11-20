#ifndef BernseeFFT_h
#define BernseeFFT_h

#include <math.h>
#define M_PI 3.14159265358979323846

/**
FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)

Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the time domain data in 
fftBuffer[0...2*fftFrameSize-1]. The FFT array takes and returns the cosine and sine parts in an 
interleaved manner, ie. fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize 
must be a power of 2. It expects a complex input signal (see footnote 2), ie. when working with 
'common' audio signals our input signal has to be passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. 
In that case, the transform of the frequencies of interest is in fftBuffer[0...fftFrameSize]. */
void smbFft(float *fftBuffer, long fftFrameSize, long sign);

#endif