#include "FourierTransformer.h"
#include "fft4g.c"

FourierTransformer::FourierTransformer()
{
 blockSize    = 4;
 blockSizeRec = 1.0/4.0;
 ip[0]        = 0;

}

FourierTransformer::~FourierTransformer()
{

}

//----------------------------------------------------------------------------
//parameter settings:

void FourierTransformer::setBlockSize(int BlockSize)
{
 if(BlockSize!=blockSize)  
 {
  blockSize    = (int) BlockSize;
  blockSizeRec = 1.0/blockSize;
  ip[0]        = 0;
 }

 //if-branch, to reset the ip[0] only if the blockSize really changes - this
 //leads to a recalculation of the w[]-table in the next call of the fft-routine
}

//----------------------------------------------------------------------------
//signal processing:

void FourierTransformer::getReAndImFromSig(double *Sig, double *ReAndIm)
{
 //copy the signal into the complex processing buffer by interleaving 
 //it with zeros:
 int i;
 for(i=0; i<blockSize; i++)
 {
  complexBuf[2*i]   = Sig[i];
  complexBuf[2*i+1] = 0;
 }

 //perform the actual FFT on the complex buffer (in place):
 cdft(2*blockSize, -1, complexBuf, ip, w);

 //copy the first half of the complex spectrum (the positive
 //frequencies) into the output buffer, thereby multiply the coeffs
 //with a factor of two to compensate for the energy loss from
 //dicarding the negative frequencies. the first two values in the 
 //array are treated seperately because they have a special meaning:
 //ReAndIm[0] represents the (purely real) DC offset and
 //ReAndIm[1] represents the (purely real) coefficient for Nyquist-frequency:
 ReAndIm[0] = complexBuf[0];
 ReAndIm[1] = complexBuf[blockSize]; //=(2*N)/2
 for(i=2; i<blockSize; i++)
  ReAndIm[i] = 2*complexBuf[i];
}

void FourierTransformer::getMagAndPhsFromSig(double *Sig, double *Mag, double *Phs)
{
 double im, re;

 //copy the signal into the complex processing buffer by interleaving 
 //it with zeros:
 int i;
 for(i=0; i<blockSize; i++)
 {
  complexBuf[2*i]   = Sig[i];
  complexBuf[2*i+1] = 0;
 }

 //perform the actual FFT on the complex buffer (in place):
 cdft(2*blockSize, -1, complexBuf, ip, w);

 /*
 //for debug: print the complex buffer:
 printf("%s", "complex buffr after FFT: \n");
 for(i=0; i<2*blockSize; i++)
  printf("%.5f %s", complexBuf[i], "\n");
 */


 //from the interleaved re/im reprensentation in the complex buffer
 //calculate the magnitudes and the phases of the fourier-coeffs, the 
 //first output buffer will contain the blockSize/2 magnitude values and
 //the second output buffer will contain the blockSize/2 phase values
 //the first field in the phs-buffer is used for the (purely real)
 //fourier coefficient at the Nyquist frequency
 Mag[0]           = complexBuf[0];          //the purely real DC magnitude
 Phs[0]           = complexBuf[blockSize];  //the purely real Nyquist magnitude is
                                            //stored i the first field of the phs-array
 for(i=1; i<(blockSize/2); i++)
 {
  re = complexBuf[2*i];
  im = complexBuf[2*i+1];

  Mag[i] = 2*sqrt(re*re + im*im);

  if(re==0 && im==0)
   Phs[i] = 0;
  else
   Phs[i] = atan2(im, re);
 }

 /*
 printf("%s", "magnitudes: \n");
 for(i=0; i<(blockSize/2); i++)
  printf("%.5f %s", Mag[i], "\n");
 printf("%s", "phases: \n");
 for(i=0; i<(blockSize/2); i++)
  printf("%.5f %s", Phs[i], "\n");
  */
}


void FourierTransformer::getMagFromSig(double *Sig, double *Mag)
{
 double im, re;

 //copy the signal into the complex processing buffer by interleaving 
 //it with zeros:
 int i;
 for(i=0; i<blockSize; i++)
 {
  complexBuf[2*i]   = Sig[i];
  complexBuf[2*i+1] = 0;
 }

 //perform the actual FFT on the complex buffer (in place):
 cdft(2*blockSize, -1, complexBuf, ip, w);

 Mag[0]           = complexBuf[0];          //the purely real DC magnitude
 for(i=1; i<(blockSize/2); i++)
 {
  re = complexBuf[2*i];
  im = complexBuf[2*i+1];

  Mag[i] = 2*sqrt(re*re + im*im);
 }
}


void FourierTransformer::getSigFromReAndIm(double *ReAndIm, double *Sig)
{
 //copy the real and imaginary part into the complex buffer, thereby
 //duplicating it conjugate symmetric into the upper half of the
 //spectrum (the negative frequencies) and multiply by 0.5 to compensate
 //for the energy-gain due to theis duplication. ReAndIm[0] and ReAndIm[1]
 //play the special role:
 //ReAndIm[0] represents the (purely real) DC offset and
 //ReAndIm[1] represents the (purely real) coefficient for Nyquist-frequency:
 complexBuf[0]           = ReAndIm[0];
 complexBuf[1]           = 0.0;
 complexBuf[blockSize]   = ReAndIm[1];
 complexBuf[blockSize+1] = 0.0;
 int i;
 for(i=2; i<blockSize; i+=2)
 {
  //first half of the buffer is filled with the original spectrum
  //time 0.5:
  complexBuf[i]   = 0.5*ReAndIm[i];
  complexBuf[i+1] = 0.5*ReAndIm[i+1]; 
  //second half of the buffer is filled according to the
  //conjugate symmetry of the DFT:
  complexBuf[blockSize+i]   =  0.5*ReAndIm[blockSize-i];
  complexBuf[blockSize+i+1] = -0.5*ReAndIm[blockSize-i+1];
 }
 
 /*
 //for debug: print the complex buffer:
 printf("%s", "generated complex buffer from ReAndIm: \n");
 for(i=0; i<2*blockSize; i++)
  printf("%.5f %s", complexBuf[i], "\n");
 */
 
 //perform the actual IFFT on the complex buffer (in place):
 cdft(2*blockSize, 1, complexBuf, ip, w);

 //from the resulting (complex) signal, take the real part only
 //(the imaginary part should be zero anyway, because of the
 //conjugate symmetry of the spectrum), thereby dividing it by
 //the blockSize to normalize it according to the IDFT-formula:
 for(i=0; i<blockSize; i++)
  Sig[i] = blockSizeRec*complexBuf[2*i];
}

void FourierTransformer::getSigFromMagAndPhs(double *Mag, double *Phs, double *Sig)
{
 //generate the complex buffer from the magnitude and phase values:
 complexBuf[0]           = Mag[0];   //DC magnitude (DC is purely real)
 complexBuf[1]           = 0.0;
 complexBuf[blockSize]   = Phs[0];   //Nyquist magnitude (Nyquist is purely real)
 complexBuf[blockSize+1] = 0.0;

 int i;
 for(i=1; i<(blockSize/2); i++)
 {
  //fill the first half of the complex buffer:
 //apply a factor of 0.5 to all fourier coeffs in order to preserve energy
 //when duplicating the coefficients according to their conjugate symmetry:
  complexBuf[2*i]   = 0.5 * Mag[i] * cos(Phs[i]);   //real part
  complexBuf[2*i+1] = 0.5 * Mag[i] * sin(Phs[i]);   //imaginary part

  //second half of the buffer is filled according to the
  //conjugate symmetry of the DFT:
  complexBuf[2*blockSize-2*i]   =  complexBuf[2*i];
  complexBuf[2*blockSize-2*i+1] = -complexBuf[2*i+1];
 }

 /*
 //for debug: print the complex buffer:
 printf("%s", "generated complex buffer from MagAndPhs: \n");
 for(i=0; i<2*blockSize; i++)
  printf("%.5f %s", complexBuf[i], "\n");
  */

 //perform the actual IFFT on the complex buffer (in place):
 cdft(2*blockSize, 1, complexBuf, ip, w);

 //from the resulting (complex) signal, take the real part only
 //(the imaginary part should be zero anyway, because of the
 //conjugate symmetry of the spectrum), thereby dividing it by
 //the blockSize to normalize it according to the IDFT-formula:
 for(i=0; i<blockSize; i++)
  Sig[i] = blockSizeRec*complexBuf[2*i];

}