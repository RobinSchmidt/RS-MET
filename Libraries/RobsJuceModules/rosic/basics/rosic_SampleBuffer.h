#ifndef rosic_SampleBuffer_h
#define rosic_SampleBuffer_h

//// rosic-indcludes:
//#include "rosic_HelperFunctions.h"
//#include "rosic_Interpolator.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /**

  This class can be used to store and retrieve audio-data and to manipulate it in various ways.

  todo: estimateFundamental, setPhaseRandomize

  */

  class SampleBuffer  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SampleBuffer();

    /** Destructor */
    virtual ~SampleBuffer();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the waveform data. It will be copied into internal buffers. The return value informs 
    whether this operation was succesful, i.e. if the required memory could be allocated. */
    bool setSampleData(float** newSampleData, int newNumSamples, int newNumChannels, 
      float newSampleRate = 44100.0f); 

    /** Copies the data (the content of all member variables) from the passed source into this
    instance of SampleBuffer. */
    void copyDataFrom(const SampleBuffer& source);

    /** Sets the waveform data from single precision float. It will be casted to double and
    stored in internal buffers. @see setSampleData(double**, int, int) */
    //bool setSampleData(float** newSampleData, int newNumSamples, int newNumChannels); 

    //---------------------------------------------------------------------------------------------
    // inquiry and analysis:

    /** Returns the number of channels. */
    int getNumChannels() const { return numChannels; }

    /** Returns the number of samples. */
    int getNumSamples() const { return numSamples; }

    /** Finds and returns the position of an upward zero-crossing inside the signal that is closest
    to the passed position. This position is generally not an integer as the exact zero-crossing 
    will in general fall in between two samples. For multichannel data, the first channel is 
    used. */
    double findNearestUpwardZeroCrossing(double position);

    //---------------------------------------------------------------------------------------------
    // audio data manipulation:

    /** Fills the whole buffer with all zeros. */
    void fillWithZeros();

    /** Sets the value of a sample at a given index and channel (both indexes start from 0). */
    //INLINE void setSampleAt(float newSampleValue, int sampleNumber, int channel = 0);

    //---------------------------------------------------------------------------------------------
    // audio data retrieval:

    /** Returns a pointer to a pointer to the stored data-buffers. When this double-pointer is 
    de-referenced as a two-dimensional array, the first index indicates the channel and the
    second index indicates the sample-numer. 
    WARNING: if you retrieve this pointer and do something with the pointer afterwards 
    (derefence it) make sure, that no new data-array is passed inside this time-slice (by 
    another thread). To be on the safe side, you can call lockSampleData() before you 
    obtain the pointer and unlockSampleData() after you have finished your work with the 
    data. */
    //float** getSampleData() { return channelPointers; }

    /** Returns a pointer to the data in one channel - for more convenient handling of one channel
    signals (no double de-referencing nececsarry). 
    WARNING: if you retrieve this pointer and do something with the pointer afterwards 
    (derefence it) make sure, that no new data-array is passed inside this time-slice (by 
    another thread). To be on the safe side, you can call lockSampleData() before you 
    obtain the pointer and unlockSampleData() after you have finished your work with the 
    data.  */
    //float* getChannelSampleData(int channel = 0);

    /** Returns the sample at a given index and channel (both indexes start from 0). */
    INLINE double getSampleAt(int sampleNumber, int channel = 0);

    /** Returns the sample at a given fractional index using linear interpolation. */
    INLINE double getSampleLinearAt(double sampleNumber, int channel = 0);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Aquires a mutex-lock for access to the sample-data. */
    void lockSampleData() { mutex.lock(); }

    /** Releases th mutex-lock for access to the sample-data. */
    void unlockSampleData() { mutex.unlock(); }

    //=============================================================================================

  protected:

    /** Re-allocates the memory for the actual sample data and the channel pointers if the new 
    number of samples and/or channels differs from the current settings. It will also set up the 
    channel pointers to the beginnings of the channel-data sub-arrays. It returns true when memory 
    could be allocated as required, false otherwise. In the latter case, it will free all resources, 
    set the pointers to NULL and numChannels and numSamples to 0. */
    bool reAllocateMemoryIfNecessary(int newNumSamples, int newNumChannels);

    /** Frees the allocated memory and sets the pointers to NULL. */
    void freeMemory();

    static const int interpolatorMargin = 1;

    int     numChannels;
    int     numSamples;
    float   sampleRate;
    float*  sampleData;
    float** channelPointers;

    MutexLock mutex;  // get rid of this! threading should be handled in jura_framework, and/or
                      // make a subclass SampleBufferSafe that wraps all accesses into mutex-locks

    /** Make a copy-constructor unavailable because this class needs deep copying (because of the 
    pointers in the MutexLocks). If you need to create a new object based on an existing one,
    simply create a new object with the default constructor and then use copyDataFrom(). */
    SampleBuffer(const SampleBuffer& sampleBufferToCopy);

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double SampleBuffer::getSampleAt(int sampleNumber, int channel)
  {
    if( sampleData == NULL || numSamples <= 0 || numChannels <= 0)
      return 0.0;

    mutex.lock();
    float out = channelPointers[channel][sampleNumber];
    mutex.unlock();

    return out;
  }

  INLINE double SampleBuffer::getSampleLinearAt(double sampleNumber, int channel)
  {
    if( sampleData == NULL || numSamples <= 0 )
      return 0.0;

    if( channel >= numChannels )
      channel = numChannels-1;

    int    iPart;
    double fPart;
    double out; 

    mutex.lock();
    iPart = (int) sampleNumber;
    fPart = sampleNumber - (double) iPart;
    //float  out   = Interpolator::getSampleLinear(fPart, &(channelPointers[channel][iPart]) );
    out = (1.0-fPart)*channelPointers[channel][iPart] + fPart*channelPointers[channel][iPart+1];
    mutex.unlock();

    return out;
  }

} // end namespace rosic

#endif // #ifndef rosic_SampleBuffer_h
