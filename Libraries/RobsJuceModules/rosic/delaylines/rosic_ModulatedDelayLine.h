#ifndef rosic_ModulatedDelayLine_h
#define rosic_ModulatedDelayLine_h

//// rosic-indcludes:
//#include "rosic_FractionalDelayLine.h"
//#include "../generators/rosic_SineOscillator.h"
////#include "../filters/rosic_LowpassHighpass.h"
//#include "../filters/rosic_Equalizer.h"

namespace rosic
{

  /**

  This class ...

  */

  class ModulatedDelayLine : public FractionalDelayLine
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
    ModulatedDelayLine(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~ModulatedDelayLine();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets up the tempo in  beats per minute - overriden from FractionalDelayLine to update the 
    LFO frequencies, if necessarry. */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets a global gain value (as raw amplitude factor). */
    void setGlobalGainFactor(double newFactor) { g = newFactor; }

    /** Sets the feedforward factor (as raw amplitude factor). */
    void setFeedforwardFactor(double newFactor) { feedforward = newFactor; }

    /** Sets the feedback factor (as raw amplitude factor). */
    void setFeedbackFactor(double newFactor) { feedback = newFactor; }

    /** Sets the feedback amount in percent. */
    void setFeedbackInPercent(double newPercentage) { feedback = 0.01*newPercentage; }

    /** Sets the blend factor (as raw amplitude factor). */
    void setBlendFactor(double newFactor) { blend = newFactor; }

    /** Sets the panorama position between -1...+1. */
    void setPan(double newPan) { pan = newPan; RAPT::rsEqualPowerGainFactors(pan, &gL, &gR, -1.0, 1.0); }


    /** Sets the cycle-length/period by which the delay-times are modulated (in seconds or 
    beats). */
    void setDelayModulationCycleLength(double newCycleLength);

    /** Sets the depth by which the delay-times are modulated in normlized unit between 0...1 where
    1 is the maximum permissible modulation amount (the delay modulates down to zero in this
    case). 

    ToDo: write a function that expresses the depth in terms of pitch-deviation in semitones
    ...has to take frequency into account which means recalculation also for cycle-length and bpm 
    changes
    
    */
    void setDelayModulationDepth(double newDepth);



    /** Sets the start phase for the delay-modulation of the left channel (in degrees). */
    void setDelayModulationPhaseLeft(double newPhase);

    /** Sets the start phase for the delay-modulation of the right channel (in degrees). */
    void setDelayModulationPhaseRight(double newPhase);

    /** Switches tempo-sync for delay-modulation on or off. */
    void setDelayModulationSyncMode(bool shouldBeSynced);


    /** Sets the cycle-length/period by which the amplitudes are modulated (in seconds or 
    beats). */
    void setAmplitudeModulationCycleLength(double newCycleLength);

    /** Sets the depth by which the amplitudes are modulated in normlized unit between 0...1 where
    1 is the maximum permissible modulation amount (the amplitude modulates down to zero in this
    case). */
    void setAmplitudeModulationDepth(double newDepth);

    /** Sets the start phase for the amp-modulation of the left channel (in degrees). */
    void setAmplitudeModulationPhaseLeft(double newPhase);

    /** Sets the start phase for the amp-modulation of the right channel (in degrees). */
    void setAmplitudeModulationPhaseRight(double newPhase);

    /** Switches tempo-sync for delay-modulation on or off. */
    void setAmplitudeModulationSyncMode(bool shouldBeSynced);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the global gain value (in dB). */
    //double getGlobalGain() const { return amp2dB(g); }

    /** Returns the global gain as value raw amplitude factor. */
    double getGlobalGainFactor() const { return g; }

    /** Returns the feedforward factor (as raw amplitude factor). */
    double getFeedforwardFactor() const { return feedforward; }

    /** Returns the feedback factor (as raw amplitude factor). */
    double getFeedbackFactor() const { return feedback; }

    /** Returns the amount of feedback factor in percent. */
    double getFeedbackInPercent() const { return 100.0*feedback; }

    /** Returns the blend factor (as raw amplitude factor). */
    double getBlendFactor() const { return blend; }

    /** Returns the panorama position between -1...+1.. */
    double getPan() const { return pan; }


    /** Returns the cycle-length/period by which the delay-times are modulated (in seconds or 
    beats). */
    double getDelayModulationCycleLength() const { return delayModCycle; }

    /**Returns the depth by which the delay-times are modulated in normlized unit between 0...1 
    where 1 is the maximum permissible modulation amount (the delay modulates down to zero in this
    case). */
    double getDelayModulationDepth() const { return delayModDepth; }

    /** Returns the start phase for the delay-modulation of the left channel (in degrees). */
    double getDelayModulationPhaseLeft() const 
    { return RAPT::rsRadiantToDegree(delayOscL.getStartPhase()); }

    /** Returns the start phase for the delay-modulation of the right channel (in degrees). */
    double getDelayModulationPhaseRight() const 
    { return RAPT::rsRadiantToDegree(delayOscR.getStartPhase()); }

    /** Returns the cycle-length/period by which the amplitudes are modulated (in seconds or 
    beats). */
    double getAmplitudeModulationCycleLength() const { return ampModCycle; }

    /**Returns the depth by which the amplitudes are modulated. */
    double getAmplitudeModulationDepth() const { return ampModDepth; }

    /** Returns the start phase for the amp-modulation of the left channel (in degrees). */
    double getAmplitudeModulationPhaseLeft() const 
    { return RAPT::rsRadiantToDegree(ampOscL.getStartPhase()); }

    /** Returns the start phase for the amp-modulation of the right channel (in degrees). */
    double getAmplitudeModulationPhaseRight() const 
    { return RAPT::rsRadiantToDegree(ampOscR.getStartPhase()); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. If the accumulate flag is true, it
    will add it's output to the input samples. */
    INLINE void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR,
      bool accumulate);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the content of the delayline to zeros and resets the filter's state. */
    void resetBuffers();

    /** Re-triggers the 2 oscillators for the delay-modulation to their respective start-phase. */
    void triggerDelayModOscillators();

    /** Re-triggers the 2 oscillators for the amp-modulation to their respective start-phase. */
    void triggerAmpModOscillators();

    /** Re-triggers all 4 oscillators to their respective start-phase. */
    void triggerOscillators();

    //---------------------------------------------------------------------------------------------
    // embedded public modules:

    Equalizer equalizer;

    //=============================================================================================

  protected:

    /** Sets up the frequency of the LFO for the delaytime modulation according to the values of 
    sync, bpm and delayModCycle. */
    void setupDelayModulationFrequency();

    /** Sets up the frequency of the LFO for the amplitude modulation according to the values of 
    sync, bpm and delayModCycle. */
    void setupAmplitudeModulationFrequency();


    double delayModDepth; // modulation depth for the delaytime
    double ampModDepth;   // modulation amount for the amplitude
    double feedforward;   // feedforward factor
    double feedback;      // feedback factor
    double blend;         // blend factor
    double g;             // global gain factor
    double pan;           // pan user parameter (-1...+1)
    double gL, gR;        // gain factors for left and right for pan
    double delayModCycle; // cycle-length (in seconds or beats) for delaytime modulation
    double ampModCycle;   // cycle-length (in seconds or beats) for amplitude modulation
    bool   delayModSync;  // flag to indicate tempo-sync of delay-LFO
    bool   ampModSync;    // flag to indicate tempo-sync of amp-LFO

    //LowpassHighpass filter;
    SineOscillator  delayOscL, delayOscR, ampOscL, ampOscR;
    Interpolator    interpolatorL, interpolatorR;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  void ModulatedDelayLine::getSampleFrameStereo(double *inL, double *inR,
    double *outL, double *outR, bool accumulate)
  {
    // compute the 3 instantaneous delays:
    double d  = delayInSamples;                              // instantaneous delay for the unmodulated tap
    double dL = d + d*delayModDepth*delayOscL.getSample();   // instantaneous delay for the left tap
    double dR = d + d*delayModDepth*delayOscR.getSample();   // instantaneous delay for the right tap

    // read out the unmodulated signal from the delayline:
    int    iPart    = floorInt(d);                // integer part of the delay
    double fPart    = d - iPart;                  // fractional part of the delay
    fPart           = 1.0 - fPart;                // fractional part of the read-position
    if( fPart >= 1.0 )                            // wrap to open interval [0.0...1.0[
    {
      fPart  = 0.0;
      iPart += 1;
    }
    int    readPos  = wrapAround(tapIn-iPart-1);  // integer part of the read-position
    double yM       = interpolator.getSample(fPart, &(delayBuffer[readPos]));

    // same procedure for the 2 modulated taps:
    iPart     = floorInt(dL);
    fPart     = dL - iPart;
    fPart     = 1.0 - fPart;
    if( fPart >= 1.0 )
    {
      fPart  = 0.0;
      iPart += 1;
    }
    readPos   = wrapAround(tapIn-iPart-1);
    double yL = interpolatorL.getSample(fPart, &(delayBuffer[readPos]));

    iPart     = floorInt(dR);
    fPart     = dR - iPart;
    fPart     = 1.0 - fPart;
    if( fPart >= 1.0 )
    {
      fPart  = 0.0;
      iPart += 1;
    }
    readPos   = wrapAround(tapIn-iPart-1);
    double yR = interpolatorR.getSample(fPart, &(delayBuffer[readPos]));

    // establish the mono sum of the input signal:
    double inM = SQRT2_INV * (*inL + *inR);
    //double inM = 0.5 * (*inL + *inR);    // debug

    // add the feedback of the unmodulated output:
    double tmp = inM + feedback*yM;

    // filter that signal with the embedded equalizer and feed it into the delayline:
    tmp                = equalizer.getSample(tmp);
    delayBuffer[tapIn] = tmp;

    // !!!! sample for the interpolator ???!!!!


    // apply blend and feedforward factors to the left and right output signals:
    yL = feedforward*yL + blend*inM;
    yR = feedforward*yR + blend*inM;

    // apply the amplitude modulation, gain- and pan-factors:
    yL *= g * gL * ( 1.0 + 0.5*ampModDepth*(ampOscL.getSample()+1.0) );
    yR *= g * gR * ( 1.0 + 0.5*ampModDepth*(ampOscR.getSample()+1.0) );

    // increment tap-pointers:
    tapIn  = wrapAround(tapIn+1);
    tapOut = wrapAround(tapOut+1);  // is not used here, however...

    if( accumulate == true )
    {
      *outL = *inL + yL;
      *outR = *inR + yR;
    }
    else
    {
      *outL = yL;
      *outR = yR;
    }

  }

} // end namespace rosic

#endif // #ifndef rosic_ModulatedDelayLine_h
