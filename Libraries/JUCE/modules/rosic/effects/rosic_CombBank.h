#ifndef rosic_CombBank_h
#define rosic_CombBank_h

//// rosic-indcludes:
//#include "../filters/rosic_CombResonator.h"

namespace rosic
{

  /**

  This class implements a bank of 12 parallel comb filters. The 12 comb filters can be set up in 
  terms of pitch-offsets with respect to some fundamental pitch (which is also one of the user 
  parameters).

  \todo: 
  -include stereo-pitch-offset as in CombResonatorStereo

  */

  class CombBank
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    CombBank();   

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { equalPowerGainFactors(newDryWet, &dry, &wet, 0.0, 1.0); }

    /** Sets the level for the wet signal. */
    void setLevel(double newLevel) { gain = dB2amp(newLevel) / sqrt((double)numActiveCombs); }

    /** Sets the reference pitch (as MIDI note number) - all the individual pitch-offsets of each 
    comb are with respect to this pitch. */
    void setReferencePitch(double newPitch);

    /** Sets the reference frequency (in Hz) - all the individual pitch-offsets of each comb are 
    with respect to this pitch. */
    void setReferenceFrequency(double newFrequency) 
    { setReferencePitch(RAPT::rsFreqToPitch(newFrequency)); }

    /** Sets a offset between the the pitches for each comb-pair (for left and right channel)
    in semitones. */
    void setDetune(double newDetune);

    /** Sets the panorama position for the 1st comb (between -1...+1). */
    void setPan1(double newPan1) { pan1 = 0.5*(newPan1+1.0); }

    /** Sets the panorama position for the 2nd comb (between -1...+1). */
    void setPan2(double newPan2) { pan2 = 0.5*(newPan2+1.0); }

    /** Initializes all the pitch offsets so as to realize the 12 tone equal temperament. */
    void initPitchOffsets();

    /** Sets the pitch-offset for the comb with the given index. */
    //void setPitchOffset(int index, double newOffset);

    /** Sets an pitch-offset between the two stereo channels. */
    //void setStereoPitchOffset();

    /** Sets the stereo spread (between -100% and +100%) where -100% will pan all the odd indexed
    combs hard left (and the even indexed hard right), for +100% it's vice versa. */
    void setStereoSpread(double newStereoSpread);

    /** Sets the time for the combs to decay to -60 dB (in seconds). */
    void setDecayTime(double newDecayTime);

    /** Scales the decay-time for the high-frequency band. */
    void setHighDecayScale(double newScale);

    /** Scales the decay-time for the low-frequency band. */
    void setLowDecayScale(double newScale);

    /** Sets the crossover-frequency between mid and high frequencies. */
    void setHighCrossoverFreq(double newFreq);

    /** Sets the crossover-frequency between low and mid frequencies. */
    void setLowCrossoverFreq(double newFreq);

    /** Switches the resonator into a mode where it produces only odd harmonics. */
    void setOddOnlyMode(bool shouldCreateOnlyOddHarmonics);

    /** Sets the polarity of the wet signal negative (or positive, if false). */
    void setNegativePolarity(bool shouldBeNegative);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

    //=============================================================================================

  protected:

    static const int maxNumCombs = 12;
    CombResonator combsL[maxNumCombs], combsR[maxNumCombs];
    double referencePitch, detune;
    double pitchOffsets[maxNumCombs];
    double stereoSpread;
    double pan1, pan2, gainOddL, gainEvenL, gainOddR, gainEvenR, gain;
    double dry, wet;
    int    numActiveCombs;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void CombBank::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double accuOddL  = 0.0;
    double accuOddR  = 0.0;
    double accuEvenL = 0.0;
    double accuEvenR = 0.0;

    for(int i=0; i<maxNumCombs; i+=2)
    {
      accuEvenL += combsL[i]  .getSample(*inOutL);
      accuOddL  += combsL[i+1].getSample(*inOutL);
      accuEvenR += combsR[i]  .getSample(*inOutR);
      accuOddR  += combsR[i+1].getSample(*inOutR);
    }

    accuEvenL *= gainEvenL;
    accuOddL  *= gainOddL;
    accuEvenR *= gainEvenR;
    accuOddR  *= gainOddR;

    double wetL  = (accuEvenL+accuOddL);
    double wetR  = (accuEvenR+accuOddR);
    wetL         = (1.0-pan1)*wetL + (1.0-pan2)*wetR;
    wetR         =      pan1 *wetL +      pan2 *wetR;
    wetL        *= gain;
    wetR        *= gain;

    *inOutL = dry*(*inOutL) + wet*wetL;
    *inOutR = dry*(*inOutR) + wet*wetR;


    //*inOutL = dry * (*inOutL) + gain*wet * (accuEvenL+accuOddL);
    //*inOutR = dry * (*inOutR) + gain*wet * (accuEvenR+accuOddR);
  }

} // end namespace rosic

#endif // #ifndef rosic_CombBank_h
