#ifndef romos_ProcessingStatus_h
#define romos_ProcessingStatus_h

// some definitions for shorthand notation - get rid of these macros:
#define NUM_PLAYING_VOICES processingStatus.getNumPlayingVoices()
#define PLAYING_VOICE_INDICES processingStatus.getPlayingVoiceIndices()
//#define NUM_ALLOCATED_VOICES processingStatus.getNumAllocatedVoices()
//#define ALLOCATED_BLOCKSIZE processingStatus.getBufferSize()

//namespace romos
//{

/** This class represents the ProcessingStatus in which the modular system is used. Modules can 
inquire information about their processing ProcessingStatus by accessing the global object 
processingStatus. This information comprises the system sample-rate, the tempo in BPM, the notes 
that are currently active (for each voice), etc. Some modules access these informations in order to
make them available inside the modular system. */

class ProcessingStatus
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ProcessingStatus();

  /** Destructor. */
  ~ProcessingStatus();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the samplerate which is reported (i.e. can be inquired from) all modules in the 
  system - this function should be called in response to changes in the outside system samplerate 
  (i.e. the plugin-host sets a new samplerate or similar). */
  void setSystemSampleRate(double newSampleRate);

  /** Sets the tempo in BPM - should be called whenever the host changes the song-tempo. */
  void setTempo(double newTempo);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the sample rate under which the modular system is currently run. */
  INLINE double getSystemSampleRate() const { return systemSampleRate; }

  /** Returns the reciprocal of the sample rate under which the modular system is currently run. */
  INLINE double getSystemSamplePeriod() const { return systemSamplePeriod; }

  /** Returns 2*PI/systemSampleRate - this is the factor to convert from a frequency given in Hz 
  to a normalized radian frequency which is commonly denoted as "omega" in the DSP field. */
  INLINE double getFreqToOmegaFactor() const { return freqToOmegaFactor; }

  /** Returns the maximum blocksize that modules can take in their processBlock function. This is 
  determined by the amount of memory, they have allocated. */
  INLINE int getBufferSize() const { return bufferSize; }
   // rename to getMaxBlockSize

  /** Returns the number of voices for which memory is allocated for in/out signals. */
  //INLINE int getNumAllocatedVoices() const { return voiceAllocator.getNumVoices(); }

  /** Returns the number of currently playing voices. */
  INLINE int getNumPlayingVoices() const { return voiceAllocator.getNumPlayingVoices(); }

  /** Returns a pointer to an array with the indices of the currently playing voices. The array's 
  contents are valid for indices from 0...getNumPlayingVoices()-1. Above index 
  getNumPlayingVoices()-1, there may be garbage from previously playing voices. */
  INLINE const int* const getPlayingVoiceIndices() const 
  { return voiceAllocator.getPlayingVoiceIndices(); }

  /** Returns the number of memory slots that are required for each pin. For example if you have a 
  module with 2 inputs and 3 outputs and N is the number returned by this function, you would 
  allocate 2*N doubles for the input signals and 3*N doubles for the output signals. We must 
  allocate enough memory to hold a full buffer for all voices. */
  int getRequiredMemorySlotsPerPin() const { return voiceAllocator.getNumVoices() * bufferSize; }
  // comment seems out of date - there are no input buffers anymore. inputs are just pointers to 
  // another module's output 


  // maximum value for the buffersize to allocate. this affects also the size of the WorkArea
  static const int maxBufferSize = 512;


protected:

  /** Sets the maximum blocksize that modules can take in their processBlock function. This 
  determines the amount of memory, that modules need to allocate for their output signals (and 
  perhaps internal state). Larger blockSizes result in less function call overhead but more memory 
  consumption. The optimal value might depened on the machine, so we make this value 
  user-adjustable. This function should be called only from the top-level module, which 
  subsequently recursively triggers memory re-allocation for all child modules. Thus the 
  TopLevelModule class is responsible for keeping the value of our "allocatedBlockSize" member in 
  sync with the actually allocated blocksize inside the modules. */
  void setBufferSize(int newBufferSize);


private:

  double systemSampleRate, systemSamplePeriod, tempo;
  double freqToOmegaFactor;          // == (2*PI)/sampleRate
  // maybe have more conversion factors here like beatsToSecondsFactor, etc.
  // systemSampleRate should also include any global oversampling factor - maybe rename it to 
  // internal samplerate or just samplerate maybe later we let modules have local oversampling 
  // factors, too

  // stuff for block processing
  int bufferSize;
  // defines the maximum blocksize a module can take - larger blockSizes result in more memory 
  // being used, smaller blocksizes result in more function-call overhead -> find the sweet spot 
  // (might depend on the machine -> let the user choose)

  //bool processingSuspended;

};

extern ProcessingStatus processingStatus;

//}

#endif
