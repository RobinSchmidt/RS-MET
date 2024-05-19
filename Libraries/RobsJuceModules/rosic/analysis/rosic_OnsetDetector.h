#ifndef rosic_OnsetDetector_h
#define rosic_OnsetDetector_h

namespace rosic
{

/** A simple class for representing an onset. */

class Onset
{

public:

  /** Standard constructor - initializes members to zeros. */
  Onset()
  {
    timeInSamples   = 0;
    strength        = 0.f;
    bpm             = 0.f;
    beatProbability = 0.f;
    isBeat          = false;
  }

  int   timeInSamples;   // time of the onset in samples 
  float strength;        // strength of the onset
  float bpm;             // estimated instantaneous tempo in bpm
  float beatProbability; // probability that this onset is a beat
  bool  isBeat;          // flag to mark this onset as a beat

};

/**

This class can be used to detect onsets in audio data. It is based on computing the spectral
flux and finding local maxima of the time-varying flux function (subject to some constraints).

References:
-Simon Dixon: Onset Detection Revisited (DAFX 2006 paper)

*/

class OnsetDetector
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Standard constructor. */
  OnsetDetector();

  /** Standard destructor. */
  ~OnsetDetector();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the sample-rate of the signal to be analyzed - required for setting up blocksizes and
  interpreting the time correctly. */
  void setSampleRate(int newSampleRate);

  /** Resets the internal state. */
  void reset();

  //-----------------------------------------------------------------------------------------------
  // processing:

  /** Accepts the entire input signal at once for processing. The signal should be a mono-signal
  in single precision floating point format. When the function returns (which may take a while
  for long signals), the detected onsets can be found in the member 'onsets' and can be retrieved 
  via getOnsets(). */
  void processSignalAtOnce(float *sampleData, int numSamples, int sampleRate);

  /** Prepares the object for block-wise processing (allocating the required buffers, etc.). You
  should pass an expected signal length here (in samples) which will determine the amount of
  memory that will be reserved - if the signal will be longer (i.e. more blocks than expected are
  passed to feedSignalBlock), we will need to re-allocate memory during the analysis process
  which will degrade efficiency. */
  void prepareForBlockProcessing(int expectedSignalLength, int sampleRate);

  /** Feeds one block of input samples into the detector. */
  void feedSignalBlock(float* sampleData, int numSamples);

  /** When you have passed all the blocks of the signal via feedSignalBlock, this function should
  be called to trigger the post-processing of the spectral flux and rms data that was generated
  during these block-wise calls. When the function returns, the detected onsets can be found in
  the member 'onsets'.  */
  void finishBlockProcessing();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** A call to findOnsets triggers a potentially long calculation. You may want from time to
  time request the progress from another thread (presumably a user-interface thread). That's what
  this function is made for. */
  //float getProgressInPercent();

  /** Returns a vector containing the detected onsets. */
  std::vector<Onset> getOnsets() { return onsets; }
  // ToDo: Maybe return a const reference

  //===============================================================================================

protected:

  /** Creates an array of type 'float' (via new), initializes it with zeros and returns the
  pointer to it - the caller is responsible to evetually delete it. */
  float* createAllZeroFloatArray(int numValues);

  /** Clears all the buffers (initializes them to all zeros. */
  void clearBuffers();

  /** Creates the window function for windowing the blocks. */
  void createWindow();

  /** Computes the frequency dependent weights for the spectral flux. */
  void computeSpectralWeights();

  /** Creates a complex signal-block (as required by the FFT routine) from a real signal-block by
  filling the real parts with the windowed 'realSignal' and the imaginary parts with zeros.
  The realSignal is assumed to have a length of our 'blockSize' member. 'destination' must be a
  buffer that has twice this length. */
  void createComplexBlockForTransform(float *realSignal, float *destination);

  /** Computes the FFT magnitudes from the passed 'complexSpectrum' and stores them in
  'magnitudes'. We assume, that the input complex input spectrum has conjugate symmetry (because
  it is the spectrum of a real signal), so the array of the magnitudes will only have to be of
  size numBins/2 in order to store all the non-redundant magnitudes. 'complexSpectrum' is assumed
  to be of size 2*numBins (because it contains all the redundant bins and is complex). */
  void computeMagnitudes(float *complexSpectrum, float *magnitudes, int numBins);

  /** Computes the RMS-value of the passed block. */
  float computeBlockRms(float *block);

  /** Computes one value of the weighted spectral flux from the two passed magnitude spectra. */
  float computeSpectralFluxValue(float *magnitudes, float *oldMagnitudes,
    float *spectralWeights);

  /** Computes the time-varying spectral flux for the whole signal at once. */
  void computeSpectralFlux();

  /** Finds the onsets from the maxima in the spectral flux (subject to some constraints) and
  stores them in the member 'onsets'. */
  void findOnsetsFromFluxMaxima();

  /** Fits a quadratic parabola defined by y(x) = a*x^2 + b*x + c
  to the 3 points given in the vectors x and y (x and y should be of length 3) */
  void fitQuadratic(float *x, float *y, float &a, float &b, float &c);


  static const int maxBlockSize = 4096;

  int    blockSize;      // 1024 @ 44.1/48, 2048 @ 88.2/96, 4096 @ 176.4/192
  int    hopSize;        // blockSize/8
  int    sampleRate;     // the sample-rate of the signal to be analyzed

  int    length;         // length of the signal to be analyzed
  float* window;         // the pre-computed window-function
  float* signal;         // the (whole) signal to analyze
  int    numBlocks;      // number of blocks (equals number of values in the spectral flux array)

  float* circularBuffer;     // a circular buffer used for analyzing a signal block-wise
  float* linearBuffer;       // a linear buffer for extracting chunks from the circular buffer
  float* complexSpectrum;    // complex spectrum of current block
  float* magnitudes;         // FFT magnitudes of current block
  float* magnitudesOld;      // FFT magnitudes of previous block
  float* weights;            // frequency-dependent weights for the spectral flux
  int    bufferIndex;        // current index in the circular buffer
  int    nextBlockEnd;       // when we hit this index in the circular buffer, a new block 
                             // has been accumulated and should be processed

  std::vector<float> rms;       // RMS values of the blocks
  std::vector<float> flux;      // values of the spectral flux
  std::vector<Onset> onsets;    // the detected onsets

};

}

#endif
