#ifndef rosic_LowFrequencyOscillator_h
#define rosic_LowFrequencyOscillator_h

// rosic-includes:
#include "../generators/rosic_WaveTable.h"
#include "../others/rosic_SlewRateLimiter.h"

namespace rosic
{

  // todo: include a waveform generator

  /**

  This class defines the user parameters for the LowFrequencyOscillator class.

  */

  class LowFrequencyOscillatorParameters
  {

  public:

    // user parameters:
    double cycleLength;    // length of one cycle (in beats or seconds)
    double startPhase;     // start phase in the waveform in degrees
    double stereoPhase;    // phase-offset between left and right channel
    double attack;         // attack time of the slewrate limiter in milliseconds
    double release;        // release time of the slewrate limiter in milliseconds
    double rise;           // rise time of the fade-in ramp in milliseconds
    double bpm;            // tempo in beats per minute
    double sampleRate;     // samplerate
    bool   mute;           // deactivates the LFO
    bool   sync;           // switches tempo sync on/off

    LowFrequencyOscillatorParameters()
    {
      cycleLength = 0.25;
      startPhase  = 0.0;
      stereoPhase = 0.0;
      attack      = 1.0;
      release     = 1.0;
      rise        = 1.0;
      bpm         = 120.0;
      sampleRate  = 44100.0;
      mute        = false;
      sync        = true;
    }

  };

  /**

  This class implements a low-frequency oscillator which can produce standard and custom waveforms
  and apply a couple of waveform manipulations. It also contains a slewrate limiter to smooth out
  the discontinuities in the waveform.

  \todo: maybe include a (subsonic) highpass
  pass a flag to the constructor to indicate whether or not we need to allocate memory for the
  table - avoids memory fragmentation on construction of polyphonic instruments (in which otherwise
  each voice-LFO would allocate its own memory and later de-allocate it and switch to shared
  memory)

  */

  class LowFrequencyOscillator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    LowFrequencyOscillator();

    /** Destructor. */
    ~LowFrequencyOscillator();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the wavetable object which should be used by this oscillator. */
    void setWaveTableToUse(WaveTable *newTableToUse);

    /** Selects one of the standard waveforms. @see WaveTable waveforms. */
    //void setWaveform(int newWaveform);

    /** Passes a custom waveform to be used for the LFO. */
    //void setWaveform(double *newWaveformBuffer, int newLength, char* name);

    /** Passes a custom waveform to be used for the LFO. */
    //void setWaveform(float *newWaveformBuffer, int newLength, char* name);

    /** Sets the cycle-length in seconds or beats (depending on whether sync is active). */
    void setCycleLength(double newCycleLength);

    /** Sets up the tempo in beats per minute (for sync mode). */
    void setBeatsPerMinute(double newBPM);

    /** Turns tempo-synchronization on/off - when on, the cycle-length is interpreted as value in
    beats, when off, it's interpreted as value in seconds. */
    void setTempoSync(bool shouldSync);

    /** Sets the start-phase in degrees (range: 0...360). */
    void setStartPhase(double newStartPhase);

    /** Sets the phase offset between left and right channel in degrees (range: 0...360). */
    void setStereoPhase(double newStereoPhase);

    /** Sets the attack-time (in milliseconds) for the slewrate limiter - this is time which it
    takes to rise 63% for upward jumps in the waveform. */
    void setUpwardSlewRate(double newSlewRate)
    { slewRateLimiterL.setAttackTime(newSlewRate);  slewRateLimiterR.setAttackTime(newSlewRate); }

    /** Sets the release-time (in milliseconds) for the slewrate limiter - this is time which it
    takes to fall to 37% for downward jumps in the waveform. */
    void setDownwardSlewRate(double newSlewRate)
    { slewRateLimiterL.setReleaseTime(newSlewRate); slewRateLimiterR.setReleaseTime(newSlewRate); }

    /** Set this to "true" when this instance should be muted. */
    void setMute(bool shouldBeMuted);


    // still to do (currently implemented as dummy functions for the rosof-editors):
    void setTimeReverse(bool shouldReverse)
    {
      if( shouldReverse )
      {
        // do stuff
      }
    }
    void setPolarityInversion(bool shouldInvert)
    {
      if( shouldInvert )
      {
        // do stuff
      }
    }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the cycle-length in seconds or beats (depending on whether sync is active). */
    double getCycleLength() const { return parameters->cycleLength; }

    /** Returns the frequency in Hz. */
    double getFrequency() const;

    /** Returns true when tempo-sync is active, false otherwise. */
    int isInSyncMode() const { return parameters->sync; }

    /** Returns the start-phase in degrees (range: 0...360). */
    double getStartPhase() const { return parameters->startPhase; }

    /** Fills the targetBuffer with values suitable for displaying the current waveform in a
    display. The targetBuffer is assumed to be of size [numSamplesToShow] where the index
    represents the sample-number on the display (which is the x-coordinate in pixels). */
    //void getWaveformForDisplay(double* targetBuffer, int numSamplesToShow);

    /** Returns the name of the waveform as zero-terminated c-string - if an audiofile is used,
    this 'name' will represent the relative path of the file. */
    //char* getWaveformName() const;

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output sample at a time. */
    INLINE double getSample();

    /** Calculates a stereo output sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Calculates the phase increment. */
    INLINE void calculateIncrement();

    /** Resets the table-pointer to its start-position. */
    void trigger();

    /** Resets the table pointer to a position given as phase between 0...360 degrees. */
    void triggerWithPhase(double phase);

    /** Resets the slewratelimiter's state and retriggers the LFO. */
    void reset();

    //---------------------------------------------------------------------------------------------
    // embedded publically accessible objects:

    /** This is a pointer to the wavetable object to be used. */
    WaveTable *waveTable;
    bool      tableIsOwned; // flag to indicate that we use our own memory for the table (as opposed
                            // to using shared memory

    //=============================================================================================

  protected:

    double positionFrac;
    double incrementFrac;
    int    positionInt;
    int    incrementInt;

    SlewRateLimiter slewRateLimiterL, slewRateLimiterR;
    // \todo: scrap the embedded object and do the slewrate-limiting directly here - this optimizes
    // memory occupation for multiple voices as we only the y1 member from the SlewRateLimiter
    // per voice - everything else can be shared

    /** A pointer to the parameters which are potentially shared by among instances. */
    LowFrequencyOscillatorParameters* parameters;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void LowFrequencyOscillator::calculateIncrement()
  {
    double cycleLengthInSeconds;
    if( parameters->sync )
      cycleLengthInSeconds = RAPT::rsBeatsToSeconds(parameters->cycleLength, parameters->bpm);
    else
      cycleLengthInSeconds = parameters->cycleLength;

    double frequency      = 1.0 / cycleLengthInSeconds;
    double phaseIncrement = frequency * (double) waveTable->tableLength / parameters->sampleRate;

    incrementInt  = (int) phaseIncrement;
    incrementFrac = phaseIncrement - (double) incrementInt;
  }

  INLINE double LowFrequencyOscillator::getSample()
  {
    if( (parameters->mute == true) || (waveTable == NULL) )
      return 0.0;

    // wraparound the integer part of the position-pointer if necesarry:
    /*
    while( positionInt >= waveTable->tableLength )
      positionInt -= waveTable->tableLength;
    while( positionInt < 0 )
      positionInt += waveTable->tableLength;
    */
    positionInt = RAPT::rsWrapAround(positionInt, waveTable->tableLength);
    double tmp  = waveTable->getValueAt(positionInt, positionFrac);

    // apply rise-envelope before slewratelimiter to smooth discontinuities due to resetting the
    // envelope as well:
    // tmpL *= riseEnvelope.getSample();

    // apply slewrate limiter:
    tmp = slewRateLimiterL.getSample(tmp);

    // increment position-pointer:
    positionInt  += incrementInt;
    positionFrac += incrementFrac;
    if( positionFrac >= 1.0 )
    {
      positionFrac -= 1.0;
      positionInt  += 1;
    }

    return tmp;
  }


  INLINE void LowFrequencyOscillator::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    if( (parameters->mute == true) || (waveTable == NULL) )
    {
      *inOutL = *inOutR = 0.0;
      return;
    }

    // wraparound the integer part of the position-pointer if necesarry:
    positionInt = RAPT::rsWrapAround(positionInt, waveTable->tableLength);
    /*
    while( positionInt >= waveTable->tableLength )
      positionInt -= waveTable->tableLength;
    while( positionInt < 0 )
      positionInt += waveTable->tableLength;
      */

    int positionInt2 = positionInt + roundToInt(waveTable->tableLength*parameters->stereoPhase/360.0);
       // to be optimized - integer offset can be precalulated
    positionInt2 = RAPT::rsWrapAround(positionInt2, waveTable->tableLength);

    double tmpL = waveTable->getValueAt(positionInt,  positionFrac);
    double tmpR = waveTable->getValueAt(positionInt2, positionFrac);

    // apply rise-envelope before slewratelimiter to smooth discontinuities due to resetting the
    // envelope as well:
    // tmpL *= riseEnvelope.getSample();
    // tmpR *= riseEnvelope.getSample();

    // apply slewrate limiter:
    tmpL = slewRateLimiterL.getSample(tmpL);
    tmpR = slewRateLimiterR.getSample(tmpR);

    // increment position-pointer:
    positionInt  += incrementInt;
    positionFrac += incrementFrac;
    if( positionFrac >= 1.0 )
    {
      positionFrac -= 1.0;
      positionInt  += 1;
    }

    *inOutL = tmpL;
    *inOutR = tmpR;
  }


} // end namespace rosic

#endif // rosic_LowFrequencyOscillator_h
