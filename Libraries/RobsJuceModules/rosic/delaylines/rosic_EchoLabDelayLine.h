#ifndef rosic_EchoLabDelayLine_h
#define rosic_EchoLabDelayLine_h

//// rosic-indcludes:
////#include "rosic_FractionalDelayLine.h"
////#include "../filters/rosic_Equalizer.h"
//#include "../filters/rosic_EqualizerStereo.h"

namespace rosic
{

  /**

  This class ...

  */

  class EchoLabDelayLine
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
    EchoLabDelayLine(int maximumDelayInSamples = 65536);

    /** Destructor */
    ~EchoLabDelayLine();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the delay-time in seconds or beats (depending on whether sync is active). */
    void setDelayTime(double newDelayTime);

    /** Switches the tempo-sync on or off. */
    void setSyncMode(bool shouldTempoSync);

    /** Sets up the tempo in  beats per minute. */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets a global gain value (as raw amplitude factor). */
    void setGlobalGainFactor(double newFactor) { g = newFactor; }

    /** Sets the feedback factor (as raw amplitude factor). */
    void setFeedbackFactor(double newFactor) { feedback = newFactor; }

    /** Sets the feedback amount in percent. */
    void setFeedbackInPercent(double newPercentage) { feedback = 0.01*newPercentage; }

    /** Sets the panorama position between -1...+1. */
    void setPan(double newPan) { pan = newPan; RAPT::rsEqualPowerGainFactors(pan, &gL, &gR, -1.0, 1.0); }

    /** Switches the delayline into/out of ping-pong mode. */
    void setPingPongMode(bool shouldBePingPong) { pingPongMode = shouldBePingPong; }

    /** Mutes the output of the delayline. */
    void setMute(bool shouldBeMuted) { mute = shouldBeMuted; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the delay-time in seconds or beats (depending on whether sync is active). */
    double getDelayTime() const { return delayTime; }

    /** Returns true when tempo-sync is active, false otherwise. */
    int isInSyncMode() const { return tempoSync; }

    /** Returns the global gain as value raw amplitude factor. */
    double getGlobalGainFactor() const { return g; }

    /** Returns the feedback factor (as raw amplitude factor). */
    double getFeedbackFactor() const { return feedback; }

    /** Returns the amount of feedback factor in percent. */
    double getFeedbackInPercent() const { return 100.0*feedback; }

    /** Returns the panorama position between -1...+1.. */
    double getPan() const { return pan; }

    /** Returns true when the delayline is in ping-pong mode. */
    bool isInPingPongMode() const { return pingPongMode; }

    /** Returns true when the delayline is muted. */
    bool isMuted() const { return mute; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. If the accumulate flag is true, it
    will add it's output to the input samples. */
    INLINE void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR,
      bool accumulate);

    /** Calculates one output stereo sample-frame at a time without acquiring the locks for the 
    embedded filters. Always wrap calls to this functions into 
    acquireFilterLocks() / releaseFilterLocks(). */
    INLINE void getSampleFrameStereoWithoutLocks(double *inL, double *inR, double *outL, 
      double *outR,  bool accumulate);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Aquires the mutex-lock for the embedded filters. to be deprecated.... */
    void acquireFilterLocks() 
    { 
      //inputEqualizer.acquireLock(); 
      //feedbackEqualizer.acquireLock(); 
    }

    /** Releases the mutex-lock for the embedded filters. to be deprecated.... */
    void releaseFilterLocks() 
    { 
      //inputEqualizer.releaseLock(); 
      //feedbackEqualizer.releaseLock(); 
    }

    /** Sets the content of the delayline to zeros and resets the filter's state. */
    void clearBuffers();

    //---------------------------------------------------------------------------------------------
    // embedded public modules:

    EqualizerStereo inputEqualizer, feedbackEqualizer;

    //=============================================================================================

  protected:

    /** Wraps an integer (read/write) position into the permitted range (0...length-1). */
    INLINE int wrapAround(int position);

    /** Sets up the delay-time in samples according to the chosen delayTime, sync-mode and 
    sample-rate user parameters. */
    void setupDelayInSamples();

    //static const int interpolatorMargin = 1;

    float *delayBuffer1, *delayBuffer2;

    int    tapIn, tapOut;
    int    length;         // nominal length (without interpolator margin, maximum delay will be length-1)

    double delayInSamples;
    double delayTime;      // in seconds or beats
    double sampleRate;
    double bpm;
    double feedback;        // feedback factor
    double g;               // global gain factor
    double pan;             // pan user parameter (-1...+1)
    double gL, gR;          // gain factors for left and right for pan

    bool   tempoSync;
    bool   pingPongMode;    // flag to indicate ping-pong mode
    bool   mute;            // flag to indicate that this delayline is muted

    bool   pointersInvalid; // flag to indicate that delayBuffer1 and/or delayBuffer2 could not be 
                            // successfully allocated in the constructor

    //WarpedAllpassInterpolator interpolatorL, interpolatorR;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE int EchoLabDelayLine::wrapAround(int position)
  {
    while( position >= length )
      position -= length;  
    while( position < 0 )
      position += length;   
    return position;
  }

  INLINE void EchoLabDelayLine::getSampleFrameStereo(double *inL, double *inR,
    double *outL, double *outR, bool accumulate)
  {
    acquireFilterLocks();
    getSampleFrameStereoWithoutLocks(inL, inR, outL, outR, accumulate);
    releaseFilterLocks();
  }

  INLINE void EchoLabDelayLine::getSampleFrameStereoWithoutLocks(double *inL, double *inR,
    double *outL, double *outR, bool accumulate)
  {
    if( mute || pointersInvalid )
    {
      if( accumulate == false )
      {
        *outL = 0.0;
        *outR = 0.0;
      }
      // else: add nothing -> do nothing
      return;
    }
      

    double inM, y1, y2, yL, yR, tmp;

    // establish the mono sum of the input signal:
    inM = SQRT2_INV * (*inL + *inR);
    if( pingPongMode == false )
    {
      // read the ouput form the delayline:
      y1  = delayBuffer1[tapOut];

      // add the feedback to the input:
      tmp = inM + RAPT::rsClip(feedback*y1, -2.0, 2.0);

      // filter that signal with the embedded equalizer and feed it into the delayline:
      tmp                 = feedbackEqualizer.getSample(tmp);
      delayBuffer1[tapIn] = (float) tmp;

      // apply the output-equalizer, gain- and pan-factors:
      y1 = inputEqualizer.getSample(y1);
      yL = g*gL*y1;
      yR = g*gR*y1;
    }
    else
    {
      // apply output-equalizer to the input - due to linearity of all involved blocks, this is 
      // equivalent to applying it to the output (but this way, we need only one instance of the 
      // EQ):
      inM = inputEqualizer.getSample(inM);

      // read the ouput from the 1st and 2nd delayline (could be optimized - calculation of coeff 
      // in Interpolator is done twice):
      y1 = delayBuffer1[tapOut];
      y2 = delayBuffer2[tapOut];

      // apply cross-feedback (2nd -> 1st) and EQ and feed the result into the first delayline:
      tmp                 = inM + RAPT::rsClip(feedback*y2, -2.0, 2.0);
      tmp                 = feedbackEqualizer.getSample(tmp);
      delayBuffer1[tapIn] = (float) tmp;

      // apply cross-feedback (1st -> 2nd) and feed the result into the 2nd delayline:
      delayBuffer2[tapIn] = (float) (feedback*y1);

      // establish the left and right output signals from the outputs of the two delaylines:
      yL = g * (gL*y1 + gR*y2);
      yR = g * (gL*y2 + gR*y1);
    }

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

#endif // #ifndef rosic_EchoLabDelayLine_h
