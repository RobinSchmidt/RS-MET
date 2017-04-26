#ifndef rosic_Compressor_h
#define rosic_Compressor_h

// rosic-indcludes:
#include "rosic_DynamicsProcessorBase.h"

namespace rosic
{

  /**

  This class implements a classic compressor with threshold, ratio, attack, release and some
  additional stuff.

  */

  class Compressor : public DynamicsProcessorBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a dynamics processor with a given maximum number of samples 
    lookahead. */
    Compressor(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~Compressor();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the threshold above which the characteristic curve snaps off (in decibels). */
    void setThreshold(double newThreshold) 
    { threshold = newThreshold; updateAutoGainFactor(); }

    /** Sets the ratio which determines by how many decibels the input must be above the 
    threshold to have the output signal to hit 1 dB above the threshold - its reciprocal is
    equal to the slope of the characteristic curve above the threshold. */
    void setRatio(double newRatio) 
    { ratio = newRatio; oneOverRatio = 1.0/ratio; updateAutoGainFactor(); }

    /** Toggles automatic gain compensation on/off. When switched on, it will apply a gain such
    that the characteristic curve goes through the point (0, 0) on a dBIn/dBOut plot. */
    void setAutoGain(bool shouldBeActive) 
    { autoGainActive = shouldBeActive; updateAutoGainFactor(); }

    /** Switches the compressor into limiter mode - levels above the threshold will be mapped to
    the threshold itself - this corresponds to an infinite ratio. */
    void setLimiterMode(bool shouldLimit) { limiterMode = shouldLimit; updateAutoGainFactor(); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

  protected:

    /** Returns the value of the transfer curve at a given x - x is assumed to be the current
    amplitude (envelope) of some signal, expressed in dB. */
    INLINE double transferCurveAt(double x);

    /** Updates the member 'autoGainFactor' according to current settings. */
    void updateAutoGainFactor();

    double threshold, ratio, oneOverRatio, autoGainFactor;
    bool   autoGainActive, limiterMode;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double Compressor::transferCurveAt(double x)
  {
    if( x <= threshold )
      return x;
    else
    {
      if( limiterMode )
        return threshold;
      else
        return threshold +  oneOverRatio * (x-threshold);
    }
  }

  INLINE void Compressor::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // establish the mono sum of the (non-delayed) input signal and obtain its envelope:
    double tmpL    = inputGainFactor * (*inOutL);
    double tmpR    = inputGainFactor * (*inOutR);  
    double inM     = ONE_OVER_SQRT2  * (tmpL + tmpR);
    double inLevel = levelDetector.getLevel(inM);

    // calculate the compressor gain:
    double outLevel       = transferCurveAt(inLevel);
    double compressorGain = dB2amp(outLevel-inLevel);

    // apply the attack/release smoothing:
    //compressorGain   = attackReleaseEnveloper.getSample(compressorGain);
    double totalGain = compressorGain*outputGainFactor*autoGainFactor;

    // apply delay, gain an mix dry/wet (this can be almalgameted into a single operation):
    applyLookAheadDelay(&tmpL, &tmpR, &tmpL, &tmpR);
    tmpL    *= (dry + wet*totalGain);
    tmpR    *= (dry + wet*totalGain);
    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_Compressor_h
