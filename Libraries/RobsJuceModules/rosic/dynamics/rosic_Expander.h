#ifndef rosic_Expander_h
#define rosic_Expander_h

//// rosic-indcludes:
//#include "rosic_DynamicsProcessorBase.h"

namespace rosic
{

  /**

  This class implements a classic expander with threshold, ratio, attack, release and some
  additional stuff.

  */

  class Expander : public DynamicsProcessorBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a dynamics processor with a given maximum number of samples 
    lookahead. */
    Expander(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~Expander();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the threshold below which the characteristic curve snaps off (in decibels). */
    void setThreshold(double newThreshold) { threshold = newThreshold; }

    /** Sets the ratio which determines by how many decibels the output will be below the 
    threshold when the input signal is 1 dB below the threshold - this is equal to the slope of the 
    characteristic curve below the threshold. */
    void setRatio(double newRatio) { ratio = newRatio; }

    /** Switches the expander into gate mode - levels below the threshold will be mapped to
    -inf (in dB) or 0 in linear amplitude units - this corresponds to an infinite ratio. */
    void setGateMode(bool shouldGate) { gateMode = shouldGate; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

  protected:

    /** Returns the value of the transfer curve at a given x - x is assumed to be the current
    amplitude (envelope) of some signal, expressed in dB. */
    INLINE double transferCurveAt(double x);

    double threshold, ratio;
    bool   gateMode;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double Expander::transferCurveAt(double x)
  {
    if( x >= threshold )
      return x;
    else
    {
      if( gateMode )
        return NEG_INF;
      else
        return threshold + ratio * (x-threshold) ;
    }
  }

  INLINE void Expander::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // establish the mono sum of the (non-delayed) input signal and obtain its envelope:
    double tmpL    = inputGainFactor * (*inOutL);
    double tmpR    = inputGainFactor * (*inOutR);  
    double inM     = SQRT2_INV  * (tmpL + tmpR);
    double inLevel = levelDetector.getLevel(inM);

    // calculate the expander gain:
    double outLevel    = transferCurveAt(inLevel);
    double expanderGain = RAPT::rsDbToAmp(outLevel-inLevel);

    // apply the attack/release smoothing:
    //expanderGain     = attackReleaseEnveloper.getSample(expanderGain);
    double totalGain = expanderGain*outputGainFactor;

    // apply delay, gain an mix dry/wet (this can be almalgameted into a single operation):
    applyLookAheadDelay(&tmpL, &tmpR, &tmpL, &tmpR);
    tmpL    *= (dry + wet*totalGain);
    tmpR    *= (dry + wet*totalGain);
    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_Expander_h
