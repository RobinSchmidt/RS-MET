#ifndef rosic_WaveShaper_h
#define rosic_WaveShaper_h

// rosic-indcludes:
#include "../filters/rosic_LowpassHighpass.h"
#include "../filters/rosic_EllipticSubBandFilterDirectForm.h" 

namespace rosic
{

  /**

  This is a waveshaping distortion unit with several predefined transfer functions.

  todo:
  -include compressor shape (maybe deprecate hardclip then (redundant))

  */

  class WaveShaper
  {

  public:

    /** Enumeratation of the transfer function shapes. */
    enum transferFunctions
    {
      LINEAR,
      TANH,
      HARDCLIP,
      //CUBIC,
      QUINTIC,
      //HEPTIC,
      //SIN,

      NUM_TRANSFER_FUNCTIONS
    };

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    WaveShaper(); 

    /** Destructor. */
    ~WaveShaper(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the transfer function to use @see transferFunctions. */
    void setTransferFunction(int newTransferFunction);

    /** Sets the drive (gain) for the input (expected in dB). */
    void setDrive(double newDrive) { driveFactor = dB2amp(newDrive); }

    /** Sets the DC offset to be added to the signal (expected as raw value). */
    void setDcOffset(double newDcOffset) { dcOffset = newDcOffset; }

    /** Sets the oversampling factor for the internal signal processing - a factor of one (or 
    lower) will result in no oversampling. */
    void setOversampling(int newOversamplingFactor);

    /** Sets the global output volume (expected in dB). */
    void setOutputLevel(double newOutputLevel) { outVolFactor  = dB2amp(newOutputLevel); }

    /** Sets the amount of the effect in percent (scales the difference between the original and 
    distorted signal by a factor). */
    void setAmount(double newAmount) { amount = newAmount; }

    /** Sets the value of the derivative (slope) of the pentic tranfer function at (0,intercept). 
    A value of 15/8 = 1.875 will yield a curve that has a vanishing second derivative at (1,1) in 
    addition to the vanishing first derivative (when intercept==0). A value of 1.5 will yield the 
    same curve as with the cubic transfer function. */
    void setPenticSlopeAtZero(double newSlope) 
    { slope = newSlope; calculateQuinticCoefficients(); }

    /** Sets the interception with the y-axis for the pentic curve. */
    void setPenticInterceptY(double newIntercept) 
    { intercept = newIntercept; calculateQuinticCoefficients(); }

    /** Sets the ratio between dry and wet between 0...1. */
    //void setDryWetRatio(double newDryWet) 
    //{ equalPowerGainFactors(newDryWet, &dry, &wet, 0.0, 1.0); }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns true when the currently selected transfer function supports a slope parameter, 
    false otherwise. */
    bool doesTransferFunctionSupportSlope() { return transferFunctionIndex == QUINTIC; }

    /** Returns true when the currently selected transfer function supports an y-intercept 
    parameter, false otherwise. */
    bool doesTransferFunctionSupportInterceptY() { return transferFunctionIndex == QUINTIC; }

    /** Retruns the value of the waveshaping transfer function at some value x. */
    INLINE double transferFunctionAt(double x);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Processes a block of samples at a time. */
    inline void processBlock(double *inOutL, double *inOutR, double numSampleFrames);

    /** Calculates one stereo output frame at a time. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the resampling filters. */
    void reset();

    //=============================================================================================

  protected:

    /** Calculates the coefficients ofr the quintic polynomial according to the specifications. */
    void calculateQuinticCoefficients();

    double driveFactor, dcOffset, outVolFactor, amount;

    double a0, a1, a2, a3, a4, a5; // polynomial coefficients

    double slope, intercept;

    int    transferFunctionIndex, oversampling;

    EllipticSubBandFilterDirectForm upsamplerL, upsamplerR;
    EllipticSubBandFilterDirectForm antiAliasFilterL, antiAliasFilterR;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double WaveShaper::transferFunctionAt(double x)
  {
    double y;
    switch( transferFunctionIndex )
    {
    case TANH:     y = tanhApprox(x);         break;
    case HARDCLIP: y = clip(x, -1.0, 1.0);    break;
      /*
    case CUBIC:
      {
        x = clip(x, -1.0, 1.0);
        y = x*(1.5-0.5*x*x);
      }
      */
    case QUINTIC:
      {
        if( x > 1.0 )
          return 1.0;
        else if( x < -1.0 )
          return -1.0;
        else
        {
          double x2 = x*x;
          y         = x * (a1 + a3*x2 + a5*x2*x2) + x2 * (a4*x2 + a2) + a0;
        }
      }
      break;

      /*
    case HEPTIC:
      {
        x  = clip(x, -1.0, 1.0);
        x2 = x*x;
        return (1.0/8.0)*x * (8 + 11*x2 - 18*x2*x2 + 7*x2*x2*x2);
      }
    case SIN:      return sin(x);
    */

    default: y = x;
    }

    return y;
  }

  INLINE void WaveShaper::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {    
    double tmpL = driveFactor * (*inOutL) + dcOffset;    
    double tmpR = driveFactor * (*inOutR) + dcOffset;
    if( oversampling > 1 )
    {
      // calculate n frames of the oversampled, distorted, and anti-alias filtered signal
      // (we do oversampling by a factor of n):
      tmpL = upsamplerL.getSample(oversampling*tmpL);
      tmpR = upsamplerR.getSample(oversampling*tmpR);
      tmpL = transferFunctionAt(tmpL); 
      tmpR = transferFunctionAt(tmpR);
      tmpL = antiAliasFilterL.getSample(tmpL);
      tmpR = antiAliasFilterR.getSample(tmpR);

      for(int i=1; i<oversampling; i++)
      {
        tmpL = upsamplerL.getSample(0.0);
        tmpR = upsamplerR.getSample(0.0);
        tmpL = transferFunctionAt(tmpL);
        tmpR = transferFunctionAt(tmpR);
        tmpL = antiAliasFilterL.getSample(tmpL);
        tmpR = antiAliasFilterR.getSample(tmpR);
      }
    }
    else
    {
      tmpL = transferFunctionAt(tmpL);
      tmpR = transferFunctionAt(tmpR);
    }

    // establish the difference between distorted and clean signal:
    tmpL -= *inOutL;
    tmpR -= *inOutR;  

    // obtain the output by adding the scaled difference to the input:
    tmpL = *inOutL + (0.01*amount)*tmpL;
    tmpR = *inOutR + (0.01*amount)*tmpR;

    // finally, apply volume scaling:
    *inOutL = tmpL*outVolFactor;
    *inOutR = tmpR*outVolFactor;
  }

  void WaveShaper::processBlock(double *inOutL, double *inOutR, double numSampleFrames)
  {
    for(int n=0; n<numSampleFrames; n++)
      getSampleFrameStereo(&inOutL[n], &inOutR[n]);
  }

} // end namespace rosic

#endif // rosic_WaveShaper_h