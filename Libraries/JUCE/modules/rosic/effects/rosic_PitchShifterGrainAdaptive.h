#ifndef rosic_PitchShifterGrainAdaptive_h
#define rosic_PitchShifterGrainAdaptive_h

//// rosic-indcludes:
//#include "rosic_PitchShifter.h"
//#include "../analysis/rosic_PitchDetector.h"
//#include "../analysis/rosic_FormantPreserver.h"

namespace rosic
{

  /**

  This class augments the PitchShifter class by a PitchDetector and (optionally) adjusts the
  grainLength on the fly according to the detected pitch period.

  */

  class PitchShifterGrainAdaptive : public PitchShifter
  {

  public:

    enum grainLengthUnits
    {
      MILLISECONDS = 0,
      PITCH_CYCLES,
      BEATS
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    PitchShifterGrainAdaptive();

    /** Destructor. */
    ~PitchShifterGrainAdaptive();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Overrides setSampleRate in order to update the additional embedded modules. */
    void setSampleRate(double newSampleRate);

    /** Sets the current tempo in beats per minute - required when syncing the grain length to beats. */
    void setBeatsPerMinute(double newBpm);

    /** Sets the length of the grains for the non-adaptive mode (i.e. the grain length which is used when adaption is switched off). */
    void setGrainLengthInMilliseconds(double newGrainLength);

    /** Sets the length of the grains for the adaptive mode as a number of cycles per grain. */
    void setGrainLengthInPitchCycles(double newNumCycles);

    /** Sets the length of the grains for the adaptive mode as a number of beats per grain. */
    void setGrainLengthInBeats(double newNumBeats);

    /** Switches the unit in which the grain length is adjusted. */
    void setGrainLengthUnit(int newUnit);

    /** Switches forman preservation on/off. */
    void setFormantPreserve(bool shouldPreserveFormants);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the length of the grains for the non-adaptive mode (i.e. the grain length which is
    used when adaption is switched off). */
    double getGrainLengthInMilliseconds() const { return grainLengthInMilliseconds; }

    /** Returns the length of the grains for the adaptive mode as a number of cycles per grain. */
    double getGrainLengthInPitchCycles() const { return grainLengthInPitchCycles; }

    /** Returns the length of the grains for the beat-sync mode as a number of betas per grain. */
    double getGrainLengthInBeats() const { return grainLengthInBeats; }

    /** Informs about the current unit for the grain length. */
    int getGrainLengthUnit() const { return grainLengthUnit; }

    /** Informs, whether formant preservation is on (true) or off (false). */
    bool getFormantPreserve() const { return formantPreserve; }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo-ouput frame. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    //=====================================================================================================================================

  protected:

    double grainLengthInMilliseconds, grainLengthInPitchCycles, grainLengthInBeats;
    double tempoInBpm;
    int    grainLengthUnit;
    bool   formantPreserve;

    /** This embedded object is used to detect the pitch of the incoming signal. */
    PitchDetector pitchDetector;

    /** A smoother for the grain length adaption. */
    FourPoleFilter smoother;

    /** This object removes and later re-applies the formants. */
    FormantPreserver formantPreserver;

  private:

    /** This inherited method has been made private because it should not be used anymore - use setNonAdaptiveGrainLength() and 
    setNumPeriodsPerGrain instead. */
    void setGrainLength(double newGrainLength);

    /** This inherited method has been made private because it should not be used anymore - use getNonAdaptiveGrainLength() and 
    getNumPeriodsPerGrain() instead. */
    void getGrainLength();

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions: 

  INLINE void PitchShifterGrainAdaptive::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    //PitchShifter::getSampleFrameStereo(inOutL, inOutR);

    double length, period;

    if( grainLengthUnit == PITCH_CYCLES )
    {
      // estimate  the pitch-period (in seconds)
      period = pitchDetector.estimatePeriod(*inOutL+*inOutR);
      length = 1000.0 * grainLengthInPitchCycles * period;
    }
    else if( grainLengthUnit == BEATS )
    {
      length = 1000.0 * beatsToSeconds(grainLengthInBeats, tempoInBpm);
        // shouldn't be necessarry, but leaving it out leads to harsh noise sometimes - apparently
        // when the unit changes (load preset Shimmer, hit preset-minus button while feeding
        // sound into it to reprodcue it)
    }
    else
      length = grainLengthInMilliseconds;

    if( length != PitchShifter::grainLength )
    {
      // smooth the instantaneous grain-length to avoid glitches (the exp(-smooth(-log(blah))) can
      // be seen as conversion from length to pitch -> smoothing -> back-conversion:
      length = smoother.getSample(length);
      //length = exp(-smoother.getSample(-log(length)));
      PitchShifter::setGrainLength(length);
    }

    // do the pitch shift (with formant preservation if this is desired):
    if( formantPreserve == true )
    {
      //formantPreserver.removeFormants(inL, inR, outL, outR);
      //int dummy = 0;
      double whiteL, whiteR;
      formantPreserver.removeFormants(inOutL, inOutR, &whiteL, &whiteR);
      PitchShifter::getSampleFrameStereo(&whiteL, &whiteR);
      formantPreserver.reApplyFormants(&whiteL, &whiteR, inOutL, inOutR);
    }
    else
      PitchShifter::getSampleFrameStereo(inOutL, inOutR);
  }

} // end namespace rosic

#endif // rosic_PitchShifterGrainAdaptive_h
