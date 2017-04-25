#ifndef rosic_SampleModulator_h
#define rosic_SampleModulator_h

// rosic-indcludes:
//#include "../basics/rosic_BinarySemaphore.h"
#include "../basics/rosic_Interpolator.h"
#include "../infrastructure/rosic_MutexLock.h"
//#include "../_third_party/kvr/aciddose/admutexw32.h"
//#include "../_third_party/stk/Mutex.h"

namespace rosic
{

  /**

  This class implements a modulation generator based on a sampled waveform. This waveform may or 
  may not be a single cycle - in the case where it is not, a loop start and end-point may be 
  specified and it can be told from outside, how many cycles are inside this loop (from that value, 
  the readout-speed will be calculated). Also, it can be specified from outside, from which 
  position the readout should start, when the modulator is triggered - and this start-poisition can
  be different for left and right channel.

  */

  class SampleModulator  
  {

  public:


    /** This enumeration lists the available loop-modes. */
    enum loopModes
    {
      NO_LOOP = 0,
      FORWARD_LOOP
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    SampleModulator();
    /**< Constructor. */

    ~SampleModulator();
    /**< Destructor */

    void copyDataFrom(const SampleModulator& source);
    /**< Copies the data (the content of all member variables) from the passed source into this
         instance of SampleModulator. */

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setSampleName(char* newSampleName, int newLength);
    /**< Sets the name (i.e. the relative path) of the current sample) - this function is only
         for conveniently handling preset management in a plugIn-context */

    void setSampleData(double* newSampleData, int newNumSamples); 
    /**< Sets the waveform data. It will be copied into internal buffers. */

    void setSampleData(float* newSampleData, int newNumSamples); 
    /**< Sets the waveform data from single precision float. It will be casted to double and
         stored in internal buffers. */

    void setStartSample(int newStartSample);
    /**< Sets the sample at which the readout will start when the modulator is triggered for the 
         both channels (left and right). */

    //void setStartSampleLeft(int newLeftStartSampleLeft);
    /**< Sets the sample at which the readout will start when the modulator is triggered for the 
         left channel modulation signal. */

    //void setStartSampleRight(int newStartSampleRight);
    /**< Sets the sample at which the readout will start when the modulator is triggered for the 
         right channel modulation signal. */

    void setLoopMode(int newLoopMode);
    ///< Turns looping on or off.

    void setLoopStartSample(int newLoopStartSample);
    /**< Sets the sample at which the loop starts (the sample itself is inside the loop). */

    void setLoopEndSample(int newLoopEndSample);
    /**< Sets the sample at which the loop end (the sample itself is inside the loop). */

    void setNumCyclesInLoop(int newNumberOfCyclesInLoop);
    /**< Sets the number number of cycles which the loop represents. This is needed for the 
         calculation of the readout-frequency. */

    void setSyncMode(bool shouldBeSynced);
    /**< Switches into sync mode. This will have the effect that ... */

    void setBeatsPerMinute(double newBpm);
    /** Sets the beats per minute value in order to allow beat syncronous modulations (call 
    setSyncMode() to toggle sync-mode on or off). */

    void setNumCyclesPerTimeUnit(double newNumCycles);
    /**< Sets the number of cycles per time unit where the time unit is either seconds (in 
         which case this is the frequency in Hz) or whole notes. Which time unit of the two 
         actually will be used, depends on the value of the sync-flag (which can be set via 
         setSyncMode()). */

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    double getSampleRate();
    /**< Returns the current sample-rate. */

    char* getSampleName();
    /**< Returns the name (i.e. the relative path) of the sample as a zero-terminated string. */

    int getNumSamples();
    /**< Returns the number of samples. */

    double* getSampleData();
    /**< Returns a pointer to the currently loaded sample-data. 
         WARNING: if you retrieve this pointer and do something with the pointer afterwards 
         (derefence it) make sure, that no new data-array is passed inside this time-slice (by 
         another thread). To be on the safe side, you can call lockSampleData() before you 
         obtain the pointer and unlockSampleData() after you have finished your work with the 
         data. */

    int getStartSample();
    /**< Returns the start sample for the left channel. */

    //int getStartSampleRight();
    /**< Returns the start sample for the right channel. */

    int getLoopMode();
    /**< Returns the loop-mode. */

    int getLoopStartSample();
    /**< Returns the loop-start sample. */

    int getLoopEndSample();
    /**< Returns the loop-start sample. */

    int getNumCyclesInLoop();
    /**< Returns the number of cycles contained in the loop. */

    bool isInSyncMode();
    /**< Informs, whether the modulator is in tempo-sync mode or not. */

    double getNumCyclesPerTimeUnit();
    /**< Returns the number of cycles per time unit (seconds or whole notes). */

    int getInterpolationMethod();
    /**< Returns the interpolation method. */


    ///....more to come.....

    //---------------------------------------------------------------------------------------------
    // event-handling:

    void trigger();
    /**< Triggers the modulator, i.e. sets the read-pointers to their start-positions. */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE double getSample();
    /**< Calculates one output-sample at a time. */

    INLINE void getSampleFrameStereo(double *outL, double *outR);
    /**< Calculates one stereo output sample-frame at a time. */

    //---------------------------------------------------------------------------------------------
    // others:

    void suspendAudioProcessing();
    /**< Switches audio processing in the getSample-functions off auch that they will always 
         return a zero signal and - more importantly - will not periodically lock the mutex. */

    void resumeAudioProcessing();
    /**< Switches audio processing in the getSample-functions on - the audio processing functions 
         will periodically lock the mutex when switched on. */

    void lockSampleData();
    /**< Aquires a mutex-lock for access to the sample-data. */

    void unlockSampleData();
    /**< Releases th mutex-lock for access to the sample-data. */

  //===============================================================================================

  protected:

    static const int interpolatorMargin = 8;
    /**< The allocated memory will be a bit larger than the required delayline-length in order to 
         make life easier for the interpolator (such that the interpolator is concerned with 
         buffer-wraparounds). This is the number of samples which the buffer is longer. */

    double   sampleRate;
    double   position;
    double   increment;
    double   bpm;
    //double   periodInWholeNotes;
    double   numCyclesPerTimeUnit;
    int      numCyclesInLoop;
    int      numSamples;
    int      start;
    //int      startR;
    int      loopStart;
    int      loopEnd;
    int      loopLength;
    double*  sampleData;
    bool     loopIsOn;
    bool     syncMode;

    char*    sampleName;
    int      sampleNameLength;

    //BinarySemaphore mutex;
    MutexLock mutex;
    //ADMUTEXW32 mutex;
    //Mutex mutex;

    void updateIncrement();
    /**< Calculates the phase-increment per sample. */

    bool audioProcessingIsSuspended;

    // for debugging:
    double* sampleDataEnd;

    SampleModulator(const SampleModulator& modulatorToCopy);
    /**< Make a copy-constructor unavailable because this class needs deep copying (because of the 
         pointers in the MutexLocks). If you need to create a new object based on an existing one,
         simply create a new object with the default constructor and then use copyDataFrom(). */

  public:

    bool     endIsReached;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE double SampleModulator::getSample()
  {
    double out = 0.0;

    if( audioProcessingIsSuspended || sampleData == NULL || numSamples == 0 )
      return 0.0;


    // dummy-code for debugging:
    /*
    int dummy;
    if( mutex.isLocked == 1 ) // this should happen on wavefile loading
    {
      dummy = 1;
    }
    */


    // when accessing the arrays, we need to aquire a mutex-lock, because the memory is allocated 
    // dynamically.....
    mutex.lock();

    ///< \todo: optimize all these typecasts away (see SampleOscillator)....

    // wraparound:
    //while( position >= (double) numSamples )
    //  position -= (double) numSamples;
    while( position >= (double) (loopEnd+1) )
      position -= (double) loopLength;

    if( position < 0.0 ) // should actually never happen - check for debug
      position = 0.0;

    int    positionInt  = (int) floor(position);
    double positionFrac = position - positionInt;

    /*
    // dummy-code for debugging:
    if( mutex.isLocked == 0 ) // this should never happen here
    {
      dummy = 1; 
    }
    */

    out = Interpolator::getSampleLinear(positionFrac, &(sampleData[positionInt]) );

    position += increment;

    mutex.unlock();

    return out;
  }

  INLINE void SampleModulator::getSampleFrameStereo(double *outL, double *outR)
  {

    if( audioProcessingIsSuspended || sampleData == NULL || numSamples == 0 )
    {
      *outL = 0.0;
      *outR = 0.0;
      return;
    }

    // when accessing the arrays, we need to aquire a mutex-lock, because the memory is allocated 
    // dynamically.....
    mutex.lock();

    // access the array....

    mutex.unlock();
  }


} // end namespace rosic

#endif // #ifndef rosic_SampleModulator_h
