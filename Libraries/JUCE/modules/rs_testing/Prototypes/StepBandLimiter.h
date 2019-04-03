#pragma once


/** UNDER CONSTRUCTION

A generic implementation of the "BLEP" (= (B)and (L)imited st(EP)) technique for anti-aliasing
step discontinuities in a signal and/or its derivatives. Discontinuities in the signal itself show 
up as steps and produce spectra with a 6 dB/oct spectral roll-off, discontinuties in the first 
derivative show up as corners and produce spectra with a 12 dB/oct roll--off and so on. The basic 
idea is to take the difference between a bandlimited step (i.e. an integrated sinc function) and a 
naive step (this difference is called the residual) and add that residual to the output signal.

....

*/

template<class TSig, class TTim> // types for signal values and continuous time
class rsStepBandLimiter
{

public:




  inline TSig getSampleBlep(TSig in)
  {
    TSig y = delayBuf[bufIndex] + blep[bufIndex];
    blep[bufIndex] = 0;      // clear blep at this position - it has been consumed
    delayBuf[wrap(bufIndex+delayLength)] = in;       // write input into delayline
    bufIndex = wrap(bufIndex + 1);
    return y;
  }
  // blep only - suitable if you have discontinuities only in the signal itself

  //inline T getSampleBlamp(T in)     // 
  //inline T getSampleBlepBlamp(T in)


  inline int wrap(int i)
  {
    return i & mask;
  }


protected:


  int delayLength = 8;    // sincLength 
  int mask        = 7;    // == delayLength - 1
  int bufIndex    = 0;

  int stepSize = 1;  
  // Step-size to move from one sample to the next in the blit, blep, etc. buffers, i.e. the 
  // sample-rate at which the bandlimited functions are sampled/tabulated with respect to their 
  // zero corssings. A step-size of 1 means, they are sampled at the zero-crossings of the sinc, 
  // 2 means at zero-crossings and halfway in between (i.e. at the maxima of the underlying sine), 
  // etc.

  //int sincLength = 8;    // redundant with delayLength ...get rid of one of them
  //

  std::vector<TTim> blit, blitDrv, blep, blamp;
  // Buffers/tables for bandlimited impulse (windowed sinc), its derivative, first integral 
  // (bandlimited step) and second integral (banlimited ramp)

  std::vector<TSig> delayBuf; // buffer for delayed input signal values

};