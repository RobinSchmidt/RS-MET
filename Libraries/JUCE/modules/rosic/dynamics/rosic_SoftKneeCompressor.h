#ifndef rosic_SoftKneeCompressor_h
#define rosic_SoftKneeCompressor_h

// rosic-indcludes:
#include "rosic_Compressor.h"
#include "rosic_BesselFilterForGainSignal.h"
#include "../math/rosic_Interpolation.h"

namespace rosic
{

  /**

  This class extends the Compressor class by a soft-knee characteristic. The knee is given as a
  transition width between the two slopes of the characteristic line (in the dBIn-dBout plot). The
  curve between threshold-0.5*knee and threshold+0.5*knee will follow a cubic polynomial that 
  matches the classic (hard-knee) curve at these points and also matches the first derivative 
  there.

  */

  class SoftKneeCompressor : public Compressor
  {

  public:

    /** Enumeratation of the transtion curves between the two slopes of the compressors characteristic 
    curve. */
    enum kneeShapes
    {
      CUBIC,
      QUARTIC
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a dynamics processor with a given maximum number of samples 
    lookahead. */
    SoftKneeCompressor(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~SoftKneeCompressor();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the threshold above which the characteristic curve snaps off (in decibels). */
    void setThreshold(double newThreshold)
    { Compressor::setThreshold(newThreshold); calculateCoefficients(); }

    /** Sets the ratio which determines by how many decibels the input must be above the 
    threshold to have the output signal to hit 1 dB above the threshold - its reciprocal is
    equal to the slope of the characteristic curve above the threshold. */
    void setRatio(double newRatio) 
    { Compressor::setRatio(newRatio); calculateCoefficients(); }

    /** Switches the compressor into limiter mode - levels above the threshold will be mapped to
    the threshold itself - this corresponds to an infinite ratio. */
    void setLimiterMode(bool shouldLimit) 
    { Compressor::setLimiterMode(shouldLimit); calculateCoefficients(); }

    /** Sets the transition width between the two slopes on the dBIn-dBout plot (in decibels). */
    void setKneeWidth(double newKneeWidth);

    /** Sets the shape of the transition between threshold-knee/2 and  threshold+knee/2. 
    @see kneeShapes */
    void setKneeShape(int newKneeShape);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

  protected:

    /** Returns the value of the transfer curve at a given x - x is assumed to be the current
    amplitude (envelope) of some signal, expressed in dB. */
    INLINE double transferCurveAt(double x);

    /** Computes the polynomial coefficients. */
    void calculateCoefficients();

    /** Calculates the polynomial coefficients for the quartic transition. */
    void calculateQuarticCoefficients();

    /** Calculates the polynomial coefficients for the cubic transition. */
    void calculateCubicCoefficients();

    double kneeWidth;
    double a0, a1, a2, a3, a4; // polynomial coefficients
    int    kneeShape;
    bool   antiAlias; 

    BesselFilterForGainSignal gainFilter;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double SoftKneeCompressor::transferCurveAt(double x)
  {
    if( x < threshold-0.5*kneeWidth )
      return x;
    else if( x > threshold+0.5*kneeWidth )
    {
      if( limiterMode )
        return threshold;
      else
        return threshold + oneOverRatio * (x-threshold); 
    }
    else
      return evaluateQuartic(x, a0, a1, a2, a3, a4);
  }

  INLINE void SoftKneeCompressor::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // establish the mono sum of the (non-delayed) input signal and obtain its envelope:
    double tmpL    = inputGainFactor * (*inOutL);
    double tmpR    = inputGainFactor * (*inOutR);  
    double inM     = ONE_OVER_SQRT2  * (tmpL + tmpR);
    //inM *= ONE_OVER_SQRT2; // for test only
    double inLevel = levelDetector.getLevel(inM);

    // calculate the compressor gain:
    double outLevel       = transferCurveAt(inLevel);
    double compressorGain = dB2amp(outLevel-inLevel);

    // apply the attack/release smoothing and the addition Bessel filter for anti-aliasing:
    //compressorGain = -attackReleaseEnveloper.getSample(-compressorGain);
    if( antiAlias )
      compressorGain = gainFilter.getSample(compressorGain); 
    double totalGain = compressorGain*outputGainFactor*autoGainFactor;
    
    // apply delay, gain an mix dry/wet (this can be almalgameted into a single operation):
    applyLookAheadDelay(&tmpL, &tmpR, &tmpL, &tmpR);
    tmpL    *= (dry + wet*totalGain);
    tmpR    *= (dry + wet*totalGain);
    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_SoftKneeCompressor_h
