#include "Utilities.h"
#include "InfiniteDataStream.h"
using namespace RSLib;


/**

Abstract baseclass for time-stretching of audio data via different underlying algorithmic engines,
like, for example Olli Parviainen's SoundTouch or zPlane's Elastique. Subclasses wrap the actual
timestretching engines, such that from the client side, we just have to pass a pointer to the 
complete input data (for random access), set up the user parameters and then regularly call 
getBuffer(float *buffer, int numFrames).  Occassionaly, in between such calls to getBuffer, there 
may be calls to setStretchAndPitchFactor and/or jumpToInputFrame(int frameIndex) and the subclasses
do all the required reset/warm-up/crossfade/etc. mumbo-jumbo for the timestretching-enginess such 
that we do not see any latency from the client side. 

*/

class TimeStretcher
{

public:

  /** Constructor. You must pass the desired length for the output buffer. 512 seems to be the 
  sweet spot for the Elastique engine. */
  TimeStretcher(int outputBufferSize = 512);

  /** Destructor. */
  ~TimeStretcher();

  /** Sets the pointer to the full input stream (to which we need random access). */
  virtual void setInputStream(float *newInputStream, int newLength);

   /** Sets the stretch factor. */
  void setStretchFactor(float newStretchFactor);

  /** Sets the pitch factor. Note that the actual pitch-factor will be slightly off, due to 
  Elastique's implementation. */
  void setPitchFactor(float newPitchFactor);

  /** Sets the stretch and pitch factor at the same time. */
  virtual void setStretchAndPitchFactor(float newStretchFactor, float newPitchFactor) = 0;
  
  /** Jumps to a given sample-frame in the input material. */
  virtual void jumpToInputFrame(int frameIndex) = 0;

  /** Produces a buffer of output samples of the given length. */
  void getBuffer(float *buffer, int numFrames);

  /** Returns the current position in the input material as it appears from the outside, 
  compensating the initial (warm-up) latency of the time-stretching engine and divergence of the 
  numbers of produced outputs and consumed inputs due to stretching */
  int getApparentInputFrameIndex() const { return (int) readIndexApparent; }

  /** Returns the length of the underlying input audio stream (in sample-frames). */
  int getInputLength() const { return inStream.getDataLength(); }

  /** A complete reset. */
  virtual void reset();

protected:

  /** Must be overriden by subclasses to fill our outBuffer member with fresh data. */
  virtual void fillOutputBuffer() = 0;

  /** Must be overriden by subclasses to reset additional variables defined in subclasses, in 
  particular, the algorithmic engine. */
  virtual void resetSubclassVariables() = 0;

  /** Retrieves an output buffer with the current settings and stores it in our tmpBuffer member.
  This can subsequently be used for crossfading when the strethcing engine has been set up (and 
  warmed up) with new settings and/or we jump to another position in the input material.  */
  void prepareForCrossfade();

  /** Performs a crossfade between the contents of our outBuffer and tmpBuffer members and stores
  the result in outBuffer. Assuming that our member outBuffer is fully filled up with valid output 
  samples with new settings and/or position in the input material (as is the case after a call to 
  warmUpElastique), and the tmpBuffer is filled up with with valid output samples that have been 
  generated with the old settings, this provides a smooth transition between the old and new
  settings/postion. */
  void applyCrossfade();

  // user parameters:
  float stretchFactor, pitchFactor;

  // internal buffering:
  int outBufferLength;
  float *outBuffer, *tmpBuffer;

  // bookkeeping:
  rsInfiniteDataStream<float> inStream;  // holds the input material
  long double readIndexApparent;         // as seen from outside -> latency compensated
  int         readIndexActual;           // position, where we actually have to read data
  int         outFrameIndex;             // position in outBuffer, where we read (in getBuffer)

};

//-------------------------------------------------------------------------------------------------

/**

TimeStretcher subclass using zPlane's Elastique as algorithmic engine.

*/

class TimeStretcherElastique : public TimeStretcher
{

public:

  TimeStretcherElastique(int outputBufferSize = 512);
  ~TimeStretcherElastique();
  virtual void setStretchAndPitchFactor(float newStretchFactor, float newPitchFactor);
  virtual void jumpToInputFrame(int frameIndex);

protected:

  virtual void fillOutputBuffer();
  virtual void resetSubclassVariables();

  /** Feeds an appropriate number of input samples into the Elastique engine, such that the 1st 
  buffer produced by getBuffer() immediately starts with nonzero samples. */
  void warmUpElastique();

  // internal buffering:
  int maxInBufferLength;  // inquired from Elastique
  float *inBuffer;
  float **ppIn, **ppOut;

  // the wrapped Elastique object:
  CElastiqueProIf *pcElastiqueHandle;

};

//-------------------------------------------------------------------------------------------------

/**

TimeStretcher subclass using Olli Parviainen's SoundTouch as algorithmic engine.

*/

class TimeStretcherSoundTouch : public TimeStretcher
{

public:

  TimeStretcherSoundTouch(int outputBufferSize = 512);
  ~TimeStretcherSoundTouch();
  virtual void setStretchAndPitchFactor(float newStretchFactor, float newPitchFactor);
  virtual void jumpToInputFrame(int frameIndex);

protected:

  virtual void fillOutputBuffer();
  virtual void resetSubclassVariables();

  /** Feeds an appropriate number of input samples into the SoundTouch engine, such that the 1st 
  buffer produced by getBuffer() immediately starts with nonzero samples. */
  void warmUpSoundTouch();

  // internal buffering:
  int inBufferLength;
  float *inBuffer;

  // the wrapped SoundTouch object:
  SoundTouch soundTouch;

};

//-------------------------------------------------------------------------------------------------

// idea: implement a time-stretcher that is not based repeating and overlapping grains but on
// using two different block/FFT-sizes for input and output: the input buffer is of length N, the 
// output buffer is of length M, where M = stretchFactor*N, or N = M/stretchFactor. if this is 
// non-integer, use jittering of the buffer-positions (i.e. occasionaly skip or insert a single 
// sample) for compensation. -> requires arbitrary length FFT (Bluestein algorithm)
// -> we need a decent method for interpolating/decimating the input spectrum

// or: use fixed analysis/synthesis FFT-sizes (synthesis blocksize should be twice the analysis
// blocksize to avoid circular convolution artifacts), do the pitch-shift in the frequency domain and 
// after iFFT, do the timestretch in the time domain (via polynomial interpolation). aliasing can be 
// avoided by just zeroing out high-frequency bins before iFFT resynthesis.
