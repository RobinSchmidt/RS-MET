#ifndef rosic_SamplePlayer_h
#define rosic_SamplePlayer_h

//// rosic-indcludes:
//#include <vector>
//#include "../basics/rosic_HelperFunctions.h"
//#include "../basics/rosic_SampleBuffer.h"
//#include "../basics/rosic_SamplePlaybackParameters.h"
//#include "../filters/rosic_LowpassHighpass.h"

namespace rosic
{

  /**

  This class implements a sample player....


  */

  class SamplePlayer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SamplePlayer();

    /** Destructor */
    ~SamplePlayer();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the sample-buffer object which will be used for playback. */
    void setSampleBufferToUse(SampleBuffer *newSampleBufferToUse); 

    /** Sets the waveform data. */
    bool setSampleData(float** newSampleData, int newNumSamples, int newNumChannels); 

    /** Mutes or un-mutes this sample player. */
    void setMute(bool shouldBeMuted) { parameters->setMute(shouldBeMuted); }

    /** Sets the flag that this SamplePlayer should play solo - the flag can be retrieved by an 
    outlying class via isSolo to switch other tone generators off. */
    void setSolo(bool shouldBeSolo) { parameters->setSolo(shouldBeSolo); }

    /** Sets the root-key of this sample which is supposed to be the note which was played when
    the sample was recorded. */
    void setRootKey(double newRootKey);

    /** Sets the playback tuning in semitones. */
    void setTune(double newTune) 
    { parameters->setTune(newTune); calculateKeyAndVelocityDependentParameters(); }

    /** Sets the keytracking for the playback pitch in cents/key - a value of 100.0 means normal 
    chromatic playback. */
    void setTuneByKey(double newTuneByKey) 
    { parameters->setTuneByKey(newTuneByKey); calculateKeyAndVelocityDependentParameters(); }

    /** Sets the velocity tracking for the playback pitch in cents/step. */
    void setTuneByVel(double newTuneByVel) 
    { parameters->setTuneByVel(newTuneByVel); calculateKeyAndVelocityDependentParameters(); }

    /** Sets the playback level in decibels. */
    void setLevel(double newLevel) 
    { parameters->setLevel(newLevel); calculateKeyAndVelocityDependentParameters(); }

    /** Sets the key dependency of the playback level. */
    void setLevelByKey(double newLevelByKey) 
    { parameters->setLevelByKey(newLevelByKey); calculateKeyAndVelocityDependentParameters(); }

    /** Sets the velocity dependency of the playback level. */
    void setLevelByVel(double newLevelByVel) 
    { parameters->setLevelByVel(newLevelByVel); calculateKeyAndVelocityDependentParameters(); }

    /** Sets the detuning of the actual fundamental frequency from the root-key in cents. */
    //void setRootDetune(double newRootDetune);

    /** Activates or deactivates the loop. */
    void setLoopMode(bool shouldBeLooped) { parameters->setLoopMode(shouldBeLooped); }

    /** Sets the start-point (in samples) of the loop (if any). */
    void setLoopStart(double newLoopStart);

    /** Sets the end-point (in samples) of the loop (if any). */
    void setLoopEnd(double newLoopEnd);

    /** Sets the length of the loop (by adjusting the loopEnd). */
    void setLoopLength(double newLength);

    /** Switch snapping to zero crossings of the loop-locators on/off. */
    void setSnapLoopToZeros(bool shouldSnap);

    /** Sets the nominal frequency to be played. The actual playback-frequency will additionally 
    take into account the detune-factor as set by setDetune(). */
    INLINE void setPlaybackFrequencyNominal(double newFrequency);

    /** Sets the name (i.e. the relative path) of the current sample) - this functions is only
    for conviniently handling preset management in a plugIn-context */
    void setSampleName(char* newSampleName);

    /** Sets the current key and velocity, so the player can set up it's key and velocity
    dependent parameters. */
    void setKeyAndVel(int newKey, int newVel)
    { key = newKey; vel = newVel; calculateKeyAndVelocityDependentParameters(); }

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    /** Returns the current sample-rate. */
    double getSampleRate();

    /** Returns a pointer to a pointer to the stored data-buffers. When this double-pointer is 
    de-referenced as a two-dimensional array, the first index indicates the channel and the
    second index indicates the sample-numer. 
    WARNING: if you retrieve this pointer and do something with the pointer afterwards 
    (derefence it) make sure, that no new data-array is passed inside this time-slice (by 
    another thread). To be on the safe side, you can call lockSampleData() before you 
    obtain the pointer and unlockSampleData() after you have finished your work with the 
    data. */
    //float** getSampleData();

    /** Returns the number of channels. */
    int getNumChannels();

    /** Returns the number of samples. */
    int getNumSamples();

    /** @see setSnapLoopToZeros() */
    bool getLoopSnapToZeroMode();

    /** Returns the name (i.e. the relative path) of the sample as a zero-terminated string. */
    char* getSampleName();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double* outL, double* outR);

    //---------------------------------------------------------------------------------------------
    // matser/slave configuration:

    /** Adds a  slave instance to the vector of slaves - this will also cause the 'isMaster' flag 
    of the slave to be set to false, redirect the slaves parameters-pointer to the one of this 
    instance and delete the old (now unused) parameters pointer of the slave. */
    void addSlave(SamplePlayer* newSlave);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Calculates the phase increment from the nominal frequency and the detune-factor. */
    INLINE void calculateIncrement();

    /** Resets the phase pointer to the start sample. */
    void reset();

    /** Snaps the loop start and end locatros to the nearest upward zero crossings. */
    void snapLoopToZeros();

    /** Aquires a mutex-lock for access to the sample-data. */
    void lockSampleData();

    /** Releases th mutex-lock for access to the sample-data. */
    void unlockSampleData();

    //---------------------------------------------------------------------------------------------
    // embedded public objects:

    SamplePlaybackParameters *parameters;
    SampleBuffer             *theBuffer;

    //=============================================================================================

  protected:

    /** Calculates the gain-factor from the desired playback level, taking also into account key-
    and velocity dependencies. */
    //void calculateGainFactor();

    /** Triggers a calculation of the key and/or velocity dependent parameters. */
    void calculateKeyAndVelocityDependentParameters();

    double sampleRate;
    double position;
    double increment;
    double frequencyNominal;  // nominal playback frequency (without detune factor)
    double frequencyDetuned;  // final frequency including detuning
    double gainFactor;        // gain factor, including key- and velocity dependence
    int    key;               // currently played key/note
    int    vel;               // current velocity
    bool   loopSnap;


    LowpassHighpass filterL, filterR;





    /** A vector of pointers to other instances of this class which shall be kept in sync to this 
    instance with regard to their parameters. */
    std::vector<SamplePlayer*> slaves;

    /** A flag which indicates whether or not this instance is a master which controls other 
    instances of this class - this will also determine whether or not this objects will delete the
    pointer to the parameter set on destruction. By default, instances will be constructed as 
    master, later they can be re-configured as slaves by adding them as slaves via addSlave to 
    another instance. */
    bool isMaster;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void SamplePlayer::setPlaybackFrequencyNominal(double newFrequency)
  {
    if( newFrequency != frequencyNominal && newFrequency > 0.0 )
    {
      frequencyNominal = newFrequency;
      calculateIncrement();
    }
  }

  INLINE void SamplePlayer::calculateIncrement()
  {
    // take transposition and detuning into account for the calcluation of the
    // final frequency:
    //frequencyDetuned = parameters->detuneFactor * frequencyNominal + parameters->detuneHz;
    //frequencyDetuned = frequencyNominal;
    frequencyDetuned = frequencyNominal + parameters->getDetuneHz();

    // calculate the new phase increment:
    increment  = frequencyDetuned / parameters->getFundamentalFrequency();
    increment *= parameters->getRecordingSampleRate() / sampleRate;
  }


  INLINE void SamplePlayer::getSampleFrameStereo(double* outL, double* outR)
  {
    // \todo: optimize SamplePlayer::getSampleFrameStereo


    //
    if( theBuffer == NULL || parameters->isMuted() )
    {
      *outL = 0.0;
      *outR = 0.0;
      return;
    }

    if( position >= parameters->getLoopEnd() )
    {
      // wraparound the position-pointer:
      if( parameters->getLoopMode() == SamplePlaybackParameters::FORWARD_LOOP )
        position -= ( parameters->getLoopEnd() - parameters->getLoopStart() );
    }

    if( position >= theBuffer->getNumSamples()-1 ) // getPlaybackEnd instead of getNumSamples
    {
      *outL = 0.0;
      *outR = 0.0;
      return;
    }

    *outL = gainFactor * theBuffer->getSampleLinearAt(position, 0);
    *outR = gainFactor * theBuffer->getSampleLinearAt(position, 1);

    position += increment;
  }

} // end namespace rosic

#endif // #ifndef rosic_SamplePlayer_h
