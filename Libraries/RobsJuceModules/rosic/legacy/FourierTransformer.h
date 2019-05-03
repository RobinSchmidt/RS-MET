#ifndef FourierTransformer_h
#define FourierTransformer_h

#include <stdio.h>       //for debugging
//#include "VstTools.h"
#include "Definitions.h"

/**

This class can be used to transform a block of data into the frequency domain
via the fast fourier transform (FFT). It is designed to be easy to use, in that
it shields the user from the underlying complex arithmetic (and the resulting
interleaved buffer-represantation of complex numbers), from the inherent
redundancies of the DFT of real signals and normalization issues.


This class is built around the free FFT-routine in fft4g.c
....blah blah blah still to come

*/

class FourierTransformer  
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          FourierTransformer();  //< Constructor.
 virtual ~FourierTransformer();  //< Destructor.

 //---------------------------------------------------------------------------
 //parameter settings:

 virtual void setBlockSize(int newBlockSize);

 //---------------------------------------------------------------------------
 // forward FFT's:

 virtual void getReAndImFromSig(double *Sig, double *ReAndIm);   
 /**< Calculates real and imaginary part of the spectrum as interleaved
      double buffer: buf[2]=re[1], buf[3]=im[1], buf[4]=re[2], 
      buf[5]=im[2],... in general: buf[2*k]=re[k], buf[2k+1]=im[k],
      k=1,..,(n/2)-1 where n is the FFT-size. The first two elements of the
      buffer have a special meaning: buf[0] is the (purely real) DC and buf[1]
      is the (purely real) coefficient for the Nyquist frequency. The other
      fields contain the real and imaginary parts of the positive frequencies
      only (interleaved) because the negative frequencies are redundant 
      (they are conjugate symmetric). To compensate for the energy loss while
      discarding the negative frequencies, the positive frequencies are
      multiplied by a factor of 2. */

 virtual void getMagAndPhsFromSig(double *Sig, double *Mag, double *Phs);
 /**< Calculates magnitude and phase from a signal, where *Sig should be of
      length n, where n is the block-size as chosen with setBlockSize() 
      *Mag and *Phs should be of length n/2. */

 virtual void getMagFromSig(double *Sig, double *Mag);
 /**< Calculates the magnitude only from a signal (for analyzer-stuff). */

 //---------------------------------------------------------------------------
 //inverse FFT's:

 virtual void getSigFromReAndIm(double *ReAndIm, double *Sig);
 /**< Calculates a time signal from and interleaved buffer containing the
      real and imaginary parts of the positive frequencies (the negative
      frequencies are assumed to be conjugate symmetric). */

 virtual void getSigFromMagAndPhs(double *Mag, double *Phs, double *Sig);
 /**< Calculates a time signal from its magnitudes and phases, *Mag and *Phs
      should be of length n/2, *Sig is of length n where n is the blcok-size
      as chosen with setBlockSize(). */

 //===========================================================================

protected:

 //static const int maxBlockSize = 8192;    // maximum FFT-Size
 static const int maxBlockSize = 1048576;

 static const int bitRevWorkSize = 2048;  
 //< size of the working area for bit reversal, 
 //  has to be >= 2+sqrt(maxBlockSize/4) !!!

 double complexBuf[2*maxBlockSize];  
  // internal buffer for the in-place calculation of the FFT, twice as long
  // as the (maximum) FFT-size, because it holds complex numbers (as 
  // interleaved re/im representation).

 int    blockSize; // the actual FFT-size, has to be a power of 2

 double blockSizeRec;          
  // reciprocal of the blocksize - has to be used as scale factor 
  // for the inverse FFT

 int    ip[bitRevWorkSize];    
  // working area for bit reversal, ip[0]=0 indicates, that the table of
  // twiddle factors has to be recalculated

 double w[maxBlockSize];       
  // table of the twiddle factors - will be calculated before the actual 
  // FFT if ip[0]=0 ...can be shrinked to >= (maxBlockSize*5)/8-1

};

#endif // FourierTransformer_h
