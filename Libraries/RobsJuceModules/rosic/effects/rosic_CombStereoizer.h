#ifndef rosic_CombStereoizer_h
#define rosic_CombStereoizer_h

//// rosic-indcludes:
//#include "../delaylines/rosic_IntegerDelayLine.h"
//#include "../filters/rosic_LpfHpfApf.h"

namespace rosic
{

  /**

  This is a delayline-based pseudo stereo creator.

  */

  class CombStereoizer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    CombStereoizer();   

    /** Destructor. */
    ~CombStereoizer(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { dryWetRatio = newDryWet; RAPT::rsEqualPowerGainFactors(dryWetRatio, &dry, &wet, 0.0, 1.0);  }

    /** Switches the channels of the wet signal */
    void setSwitchWetLeftForRight(bool shouldSwitch);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the ratio between dry and wet signal in percent wet. */
    double getDryWetRatio() const { return dryWetRatio; }

    /** Informs whether the wet signal is channel-switched left for right. */
    bool isWetSignalChannelSwitched();

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
    IntegerDelayLine delayLine;
    LpfHpfApf        wetFilter;

  protected:

    double dry, wet;
    double dryWetRatio;
    bool   channelSwitch;

  };

  //-----------------------------------------------------------------------------
  //from here: definitions of the functions to be inlined, i.e. all functions
  //which are supposed to be called at audio-rate (they can't be put into
  //the .cpp file):

  INLINE void CombStereoizer::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {
    // obtain the mono-sum:
    double tmp = SQRT2_INV * (*inOutL + *inOutR);

    // delay the mono-sum:
    tmp = delayLine.getSample(tmp);

    // filter the delayed signal:
    tmp = wetFilter.getSample(tmp);

    // weight the dry and wet signal and mix it with positve phase to the left channel and 
    // with negative phase to the right channel:
    if( channelSwitch == false )
    {
      *inOutL = dry * (*inOutL) + wet * tmp;
      *inOutR = dry * (*inOutR) - wet * tmp;
    }
    else
    {
      *inOutL = dry * (*inOutL) - wet * tmp;
      *inOutR = dry * (*inOutR) + wet * tmp;
    }
  }

} // end namespace rosic

#endif // rosic_CombStereoizer_h
