#ifndef rosic_FrequencyShifter_h
#define rosic_FrequencyShifter_h

namespace rosic
{

  //===============================================================================================
  // class FreqShifterHalfbandFilter:

  /**

  This is a halfband filter specifically designed to the very special and demanding requirements 
  for use inside a frequency shifter algorithm such as flat passband up to sampleRate/4 (and no 
  rolloff before because that would result in highpassed freq-shifted outputs) and a very high 
  slope beyond sampleRate/4 (transition width should be on the order of 20 Hz). A 24th order 
  elliptic filter was found to be adequate. The filter is realized by a cascade of 3 8th order 
  direct-form filters.

  ...mmh first performance measurements seem to indicate that using EngineersFilter directly
  is more efficient (but these were done with a debug build - we need some more testing here...)

  -with feedback sometimes a parasitic oscillation at the nyquist frequency builds up (for 
   example: shift = -26.68, feedback >= 71.6 - maybe include a nyquist blocker into the feedback 
   path ...done

  */

  class FreqShifterHalfbandFilter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    FreqShifterHalfbandFilter(); 

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates output sample at a time. */
    INLINE double getSample(double in)
    { return stage3.getSample(stage2.getSample(stage1.getSample(in))); }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the 3 8th order sections. */
    void reset() { stage1.reset(); stage2.reset(); stage3.reset(); }

    //=============================================================================================

  protected:

    rsDirectFormFilterDD  stage1, stage2, stage3;

  };


  //===============================================================================================
  // class FrequencyShifter:

  /**

  This is a frequency shifter based on the article 'An Efficient Precise Frequency Shifter' by 
  Jens Groh in the CSound magazine.

  ToDo:
  -Sometimes there is some kind of sizzle after tweaking the shift. Maybe we need a mutex to ensure 
   that all oscillators are triggered at the exact same sample instant?
  -Use an optimized elliptic halfband filter (convert biquad-cascade structure to 24th order direct 
   form or series of two 12th order direct forms, if the former will be unstable)
  -Maybe introduce feedback
  -Obtain lower sideband (subtract instead of add?). Maybe make a function 
   processFrame(double in, double* outLowerSide, double* outUpperSide).

  */

  class FrequencyShifter
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    FrequencyShifter(); 

    /** Destructor. */
    ~FrequencyShifter(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the frequency shift in Hz. */
    void setFrequencyShift(double newShift) { shiftInHz = newShift; setupOscillators(); }

    /** Sets a feedback factor by which the frequency shifted output sample from the previous call
    will be fed back to the input. */
    void setFeedbackFactor(double newFactor) { feedbackFactor = newFactor; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates output sample at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the subband filters. */
    void reset();

    //=============================================================================================

  protected:

    /** Sets up the oscillators according to the desired frequency shift. */
    void setupOscillators();

    double yOld;
    double shiftInHz, sampleRate, feedbackFactor;

    rsEngineersFilterMono halfbandFilter1, halfbandFilter2;
    //FreqShifterHalfbandFilter halfbandFilter1, halfbandFilter2;  // Old - but ma be used again
    SineOscillator   cosOsc1, sinOsc1, cosOsc2, sinOsc2;
    NyquistBlocker   nyquistBlocker;

    MutexLock mutex; // Try to get rid! Thread stuff should not be done on the DSP code level.
  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE double FrequencyShifter::getSample(double in)
  { 
    mutex.lock();

    in += feedbackFactor * yOld;

    double tmp1 = in * cosOsc1.getSample();
    tmp1        = halfbandFilter1.getSample(tmp1);
    tmp1       *= cosOsc2.getSample();

    double tmp2 = in * sinOsc1.getSample();
    tmp2        = halfbandFilter2.getSample(tmp2);
    tmp2       *= sinOsc2.getSample();

    yOld = nyquistBlocker.getSample( 2.0*(tmp1+tmp2) );

    mutex.unlock();

    //if( fabs(yOld) > 10.0 )
    //  DEBUG_BREAK;

    return yOld;
  }


  //===============================================================================================
  // class FrequencyShifterStereo:

  /**

  This is a stereo-version of the FrequencyShifter class.

  */

  class FrequencyShifterStereo
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:
 
    /** Constructor. */
    FrequencyShifterStereo(); 

    /** Destructor. */
    ~FrequencyShifterStereo(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate) 
    { shifterL.setSampleRate(newSampleRate); shifterR.setSampleRate(newSampleRate); }

    /** Sets the frequency shift in Hz. */
    void setFrequencyShift(double newShift) 
    { 
      shift = newShift; 
      shifterL.setFrequencyShift(shift+0.5f*offset); 
      shifterR.setFrequencyShift(shift-0.5f*offset); 
    }

    /** Sets the offset of the frequency shift between left and right channel in Hz. */
    void setStereoOffset(double newOffset) 
    { 
      offset = newOffset; 
      shifterL.setFrequencyShift(shift+0.5f*offset); 
      shifterR.setFrequencyShift(shift-0.5f*offset); 
    }

    /** Sets a feedback factor by which the frequency shifted output sample from the previous call
    will be fed back to the input. */
    void setFeedbackFactor(double newFactor) 
    { shifterL.setFeedbackFactor(newFactor); shifterR.setFeedbackFactor(newFactor); }

    /** Sets the feedback in percent. */
    void setFeedbackInPercent(double newFeedback) { setFeedbackFactor(0.01*newFeedback); }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { RAPT::rsEqualPowerGainFactors(newDryWet, &dry, &wet, 0.0, 1.0); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one stereo output frame at a time. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states. */
    void reset() { shifterL.reset(); shifterR.reset(); }

    //=============================================================================================

  protected:

    double dry, wet;
    double shift, offset;

    FrequencyShifter shifterL, shifterR;

  };

  //-----------------------------------------------------------------------------------------------
  // definitions of inlined functions:

  INLINE void FrequencyShifterStereo::getSampleFrameStereo(double* inOutL, double* inOutR)
  { 
    double tmpL = shifterL.getSample(*inOutL);
    double tmpR = shifterR.getSample(*inOutR);

    *inOutL     = dry*(*inOutL) + wet*tmpL;
    *inOutR     = dry*(*inOutR) + wet*tmpR;
  }

} // end namespace rosic

#endif // rosic_FrequencyShifter_h