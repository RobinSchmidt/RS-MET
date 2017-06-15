#ifndef rosic_NoiseGenerator_h
#define rosic_NoiseGenerator_h

//// rosic-indcludes:
////#include "../math/rosic_ElementaryFunctionsReal.h"
//#include "../infrastructure/rosic_MutexLock.h"
////#include "../others/rosic_RandomNumberGenerator01.h"
//#include "../transforms/rosic_FourierTransformerBluestein.h"
//#include <new>


namespace rosic
{

  /**

  This class implements a colored noise generator based on specifying a magnitude spectrum in terms
  of a slope, randomizing the phases and generating a sample from that.

  \todo:
  -generate stereo outputs by using another read-pointer with distance of half the buffersize
  -seems to be slow on construction (Quadrifex)
  -random number generator seems to be not so random
  -parasitic clicks (at buffer wraparounds?) even with the standard rand function
  -introduce brickwall high- and lowpass filtering in the spectral domain

  */

  class NoiseGenerator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - the parameter determines how much memory will be allocated for the sample to
    be played. If you pass 0 here, you may later make the gerator use shared memory by means of
    setSharedMemoryAreaToUse(). */
    NoiseGenerator(int bufferLengthToAllocate = 262144);

    /** Destructor */
    ~NoiseGenerator();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Passes a memory area that is to be used by this object. */
    void setSharedMemoryAreaToUse(void *startAddress, int sizeInBytes);

    /** Sets the sample-rate - relevant for the high- and lowpass filters. */
    void setSampleRate(double newSampleRate);

    /** Sets the slope of the magnitude spectrum in dB/oct. */
    void setSpectralSlope(double newSlope);

    /** Sets the lowest frequency to be generated (in Hz) - the actual cutoff frequency will be the
    frequency of the FFT-bin associated with a frequency strictly lower than the value passed
    here. */
    void setLowestFrequency(double newLowestFrequency);

    /** Sets the highest frequency to be generated (in Hz) - the actual cutoff frequency will be the
    frequency of the FFT-bin associated with a frequency strictly higher than the value passed
    here. */
    void setHighestFrequency(double newHighestFrequency);

    /** Sets the seed for the random number generator that creates the phase spectrum. */
    void setRandomPhaseSeed(int newSeed);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Aquires the mutex-lock for accessing the (possibly shared) memory buffer. */
    void acquireLock() { mutex.lock(); }

    /** Releases the mutex-lock for accessing the (possibly shared) memory buffer. */
    void releaseLock() { mutex.unlock(); }

    /** Resets the noise generator's  sample index. */
    void trigger();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time with aquiring the mutex lock for the
    possibly shared memory. Should be used when the app potentially (re)allocates memory during the
    lifetime of this object. */
    INLINE double getSampleThreadSafe();

    /** Calculates one output stereo sample-frame at a time without aquiring the mutex lock for the
    possibly shared memory. If you are potentially doing memory (re)allocations during the lifetime
    of this object (via setSharedMemoryAreaToUse), either wrap calls to this function into
    acquireLock()/releaseLock() or use getSampleFrameStereoThreadSafe(). */
    INLINE double getSample();

    //=============================================================================================

  protected:

    /** Creates the noise sequence for subsequent playback. */
    void createNoiseSequence();

    /** Frees the allocated memory, in case it is not shared memory. */
    void freeMemoryIfNotShared();

    double *buffer;          // buffer containing the noise sequence for playback
    int    length;           // length of the buffer
    int    readIndex;        // current index in the buffer to read from
    double sampleRate;       // the samplerate at which this generator runs
    double slope;            // slope of the magnitude spectrum
    double lowestFreq;       // lowest frequency to be generated
    double highestFreq;      // highest frequency to be generated
    int    seed;             // seed for the random number generator for the phase spectrum
    bool   pointerInvalid;   // flag to indicate that the 'buffer' poinet is invalid
    bool   memoryIsShared;   // indicates whether we use shared memory

    MutexLock mutex;         // mutex-lock for accessing the possibly shared memory

    //=============================================================================================

  private:

    // make assignment operator and copy constructor unavailable because this class contains
    // pointer members:
    NoiseGenerator& operator=(const NoiseGenerator& /*other*/) { return *this; }
    NoiseGenerator(const NoiseGenerator& /*other*/) { }

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  double NoiseGenerator::getSampleThreadSafe()
  {
    mutex.lock();
    double result = getSample();
    mutex.unlock();
    return result;
  }

  INLINE double NoiseGenerator::getSample()
  {
    if( pointerInvalid )
      return 0.0;

    readIndex++;
    if( readIndex >= length )
      readIndex = 0;
    return buffer[readIndex];
  }

} // end namespace rosic

#endif // #ifndef rosic_NoiseGenerator_h
