#ifndef rosic_PhaseStereoizer_h
#define rosic_PhaseStereoizer_h

//// rosic-indcludes:
//#include "../filters/rosic_LowpassHighpass.h"
//#include "../filters/rosic_QuadratureNetwork.h"

namespace rosic
{

  /**

  This is pseudo stereo creator based on a qudrature network to create left and right channel
  signals that are 90 degrees out of phase with respect to one another.

  todo: maybe perform dry/wet mixing on M/S representation rather than L/R

  */

  class PhaseStereoizer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    PhaseStereoizer();

    /** Destructor. */
    ~PhaseStereoizer();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the phase-offset between left and right channel in degrees. */
    void setPhaseOffset(double newPhaseOffset);

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet)
    { dryWetRatio = newDryWet; RAPT::rsEqualPowerGainFactors(dryWetRatio, &dry, &wet, 0.0, 1.0);  }

    /** Adjusts the ratio between the mid and side signal where 0.0 means mid only and 1.0 means
    side only. */
    void setMidSideRatio(double newRatio)
    { RAPT::rsEqualPowerGainFactors(newRatio, &mid, &side, 0.0, 1.0); }

    /** Sets the cutoff frequency of a lowpass filter for the wet side signal. */
    void setLowpassCutoff(double newCutoff)
    { filterM.setLowpassCutoff(newCutoff); filterS.setLowpassCutoff(newCutoff); }

    /** Sets the cutoff frequency of a highpass filter for the wet side signal. */
    void setHighpassCutoff(double newCutoff)
    { filterM.setHighpassCutoff(newCutoff); filterS.setHighpassCutoff(newCutoff); }

    /** Switches the channels of the wet signal */
    void setSwitchWetLeftForRight(bool shouldSwitch) { channelSwitch = shouldSwitch; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the ratio between dry and wet signal in percent wet. */
    double getDryWetRatio() const { return dryWetRatio; }

    /** Informs whether the wet signal is channel-switched left for right. */
    bool isWetSignalChannelSwitched() const { return channelSwitch; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo-ouput frame. */
    INLINE void getSampleFrameStereo(double* inOutL,  double* inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

    //=============================================================================================

    // embedded modules:
    rsQuadratureNetwork quadratureNetwork;
    LowpassHighpass   filterM, filterS;

  protected:

    double dry, wet, mid, side;
    double dryWetRatio;
    double phaseOffset;
    double sinFactor, cosFactor;
    bool   channelSwitch;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void PhaseStereoizer::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {
    // obtain the mono-sum:
    double tmp = 0.5 * (*inOutL + *inOutR);

    // create two versions of the mono-sum with 90 degree phase shift:
    double wetL, wetR;
    quadratureNetwork.getOutputSamplePair(tmp, &wetL, &wetR);

    // re-adjust the phase of the right channel so as to have the desired phase relationship
    // with respect to the left channel:
    wetR = cosFactor*wetL + sinFactor*wetR;

    // apply mid/side, filtering and dry/wet:
    double wetM = mid  * (wetL + wetR);
    double wetS = side * (wetL - wetR);
    //wetM        = filterM.getSample(wetM);
    wetS        = filterS.getSample(wetS);
    wetL        = 0.5  * (wetM + wetS);
    wetR        = 0.5  * (wetM - wetS);

    //*inOutL = wetL;
    //*inOutR = wetR;

    //double inM  = ONE_OVER_SQRT2 * (*inOutL + *inOutR);
    double inS  = SQRT2_INV * (*inOutL - *inOutR);
    double outM = wetM;
    double outS = dry*inS + wet*wetS;
    *inOutL     = SQRT2_INV * (outM + outS);
    *inOutR     = SQRT2_INV * (outM - outS);

    //*inOutL     = dry*(*inOutL) + wet*wetL;
    //*inOutR     = dry*(*inOutR) + wet*wetR;
  }

} // end namespace rosic

#endif // rosic_PhaseStereoizer_h
