#ifndef rosic_SoftKneeExpander_h
#define rosic_SoftKneeExpander_h

//// rosic-indcludes:
//#include "rosic_Expander.h"
//#include "rosic_BesselFilterForGainSignal.h"
//#include "../math/rosic_Interpolation.h"

namespace rosic
{

  /**

  This class extends the Expander class by a soft-knee characteristic. The knee is given as a
  transition width between the two slopes of the characteristic line (in the dBIn-dBout plot). The
  curve between threshold-0.5*knee and threshold+0.5*knee will follow a cubic polynomial that 
  matches the classic (hard-knee) curve at these points and also matches the first derivative 
  there.

  */

  class SoftKneeExpander : public Expander
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
    SoftKneeExpander(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~SoftKneeExpander();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the threshold below which the characteristic curve snaps off (in decibels). */
    void setThreshold(double newThreshold)
    { Expander::setThreshold(newThreshold); calculateCoefficients(); }

    /** Sets the ratio which determines by how many decibels the output will be below the 
    threshold when the input signal is 1 dB below the threshold - this is equal to the slope of the 
    characteristic curve below the threshold. */
    void setRatio(double newRatio) 
    { Expander::setRatio(newRatio); calculateCoefficients(); }

    /** Switches the expander into gate mode - levels below the threshold will be mapped to
    -inf (in dB) or 0 in linear amplitude units - this corresponds to an infinite ratio. */
    void setGateMode(bool shouldGate) 
    { Expander::setGateMode(shouldGate); calculateCoefficients(); }

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

    /** Returns the value of the transfer curve at a given x - x is assumed to be the current
    amplitude (envelope) of some signal, expressed in dB. */
    INLINE double transferCurveAt(double x);

  protected:

    /** Computes the polynomial coefficients. */
    void calculateCoefficients();

    /** Calculates the polynomial coefficients for the quintic transition. */
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

  INLINE double SoftKneeExpander::transferCurveAt(double x)
  {
    if( x > threshold+0.5*kneeWidth )
      return x;
    else if( x <= threshold-0.5*kneeWidth )
    {
      if( gateMode )
        return NEG_INF;
      else
        return threshold + ratio * (x-threshold);
    }
    else
      return RAPT::rsEvaluateQuartic(x, a0, a1, a2, a3, a4);
  }

  INLINE void SoftKneeExpander::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // establish the mono sum of the (non-delayed) input signal and obtain its envelope:
    double tmpL    = inputGainFactor * (*inOutL);
    double tmpR    = inputGainFactor * (*inOutR);  
    double inM     = SQRT2_INV  * (tmpL + tmpR);
    double inLevel = levelDetector.getLevel(inM);

    // calculate the compressor gain:
    double outLevel     = transferCurveAt(inLevel);
    double expanderGain = RAPT::rsDbToAmp(outLevel-inLevel);

    // apply the attack/release smoothing and the addition Bessel filter for anti-aliasing:
    //expanderGain = -attackReleaseEnveloper.getSample(-expanderGain);
    if( antiAlias )
      expanderGain = gainFilter.getSample(expanderGain); 
    double totalGain = expanderGain*outputGainFactor;
    
    // apply delay, gain an mix dry/wet (this can be almalgameted into a single operation):
    applyLookAheadDelay(&tmpL, &tmpR, &tmpL, &tmpR);
    tmpL    *= (dry + wet*totalGain);
    tmpR    *= (dry + wet*totalGain);
    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_SoftKneeExpander_h
