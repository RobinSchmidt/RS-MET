#ifndef rosic_WaveformBuffer_h
#define rosic_WaveformBuffer_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"
//#include "../datastructures/rosic_String.h"

namespace rosic
{

  /**

  This is a class for buffering a waveform and copiying that waveform into another buffer of 
  possibly different length by interpolating. The buffer can also have a name, which could be the 
  name of the file, in case the waveform originates from an audiofile.

  \todo: 
  -support stereo wavefiles (or even multichannel)
  -enable the class itself to read in the buffer from a file

  */

  class WaveformBuffer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WaveformBuffer();   

    /** Destructor. */
    ~WaveformBuffer();   

    //---------------------------------------------------------------------------------------------
    // parameter-settings:

    /** Sets the waveform from an arbitrary waveform buffer and copies the content into the 
    internal waveform buffer. The 'name' should represent something that facilitates a 
    storage/recall of the waveform - for example the relative path of an audio-file. It will be 
    returned by getWaveformName when a custom waveform is used. This function is mainly for 
    conveniently handling preset management. */
    void setWaveform(double* newWaveform, int newLength, String newName);

    /** Similar to setWaveform(float*, int), just fo single precision buffers. */
    void setWaveform(float* newWaveform, int newLength, String newName);

    /** Initializes the internal buffer with all zeros. */
    void initWaveform(int newLength = 2048, String newName = String("Empty"));

    //---------------------------------------------------------------------------------------------
    // waveform retrieval:

    /** Copies the waveform into the passed buffer, possibly interpolating if necessary (when the 
    target buffer has different length from our internal buffer). */
    void getWaveform(double *targetBuffer, int length);

    /** Returns the name of the waveform as zero-terminated c-string. This will be a pointer to 
    the original char-array which is maintained here, so DON'T attempt to DELETE it */
    const char* getWaveformName() const { return name.getRawString(); }

    //=============================================================================================

  protected:

    /** Allocates memory for the internal buffer of the given size - the function will check if 
    re-allocation is actually necessary and re-allocate only when the new size is different from 
    the current size. */
    void allocateMemory(int newNumSamples);

    double *buffer;
    int    numSamples;
    int    numChannels;
    String name;

  };

} // end namespace rosic

#endif 
