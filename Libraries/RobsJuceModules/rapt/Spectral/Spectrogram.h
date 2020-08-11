#ifndef RAPT_SPECTROGRAM_H
#define RAPT_SPECTROGRAM_H

//=================================================================================================

/** A class for representing complex spectrograms. These are 2-dimensional arrays of complex values
where each value represents the short time Fourier transform at a particular time-instant and
frequency. The data is accessed via the regular matrix(i, j) syntax of matrices. The first index
i is the frame-index and the second index j is the bin index within FFT frame i.

class is not yet used - shall be used for the return type of
rsSpectrogramPocessor::getComplexSpectrogram */

template<class T>
class rsComplexSpectrogram
{

public:


  /** Constructor. Allocates memory for the given number of frames and bins */
  rsComplexSpectrogram(int numFrames, int numBins) : stftData(numFrames, numBins) {}


  rsComplexSpectrogram(rsMatrix<std::complex<T>>&& newData) : stftData(std::move(newData)) {}
  // should be invoked in arithmetic operators +,-,etc. when client code writes:
  // rsComplexSpectrogram<double> specSum = spec1 + spec1;
  // -> test this

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  inline void setNumFramesAndBins(int numFrames, int numBins)
  { stftData->setSize(numFrames, numBins); }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of STFT frames in this spectrogram. */
  inline int getNumFrames() const { return stftData->getNumRows(); }

  /** Retruns the number of FFT bins in this spectrogram. */
  inline int getNumBins() const { return stftData->getNumColumns(); }

  /** Read and write access to STFT value at the given frame- and bin-index. */
  inline T& operator()(const int frameIndex, const int binIndex)
  { return stftData(frameIndex, binIndex); }


  //-----------------------------------------------------------------------------------------------
  // \name Manipulation

  /** Scales the data at given binIndex by the given scaleFactor in all frames. The type for the
  scale factor should either be T or std::complex<T>. */
  template<class TScaler>
  void scaleBin(int binIndex, TScaler scaleFactor)
  {
    rsAssert(binIndex >= 0 && binIndex < getNumBins());
    for(int i = 0; i < getNumFrames(); i++)
      stftData(i, binIndex) *= scaleFactor;
  }


  // maybe make a scaleFrame function as well - could be useful for applying amplitude envelopes



  //-----------------------------------------------------------------------------------------------
  // \name Operators

  bool operator==(const rsComplexSpectrogram<T>& rhs) const
  { return this->stftData == rhs.stftData; }

  bool operator!=(const rsComplexSpectrogram<T>& rhs) const
  { return this->stftData != rhs.stftData; }

  rsComplexSpectrogram<T> operator+(const rsComplexSpectrogram<T>& rightOperand) const
  { return rsComplexSpectrogram<T>(stftData + rightOperand.stftData); }

  rsComplexSpectrogram<T> operator-(const rsComplexSpectrogram<T>& rightOperand) const
  { return rsComplexSpectrogram<T>(stftData - rightOperand.stftData); }


  /** Multiplies two spectrograms element-wise. */
  rsComplexSpectrogram<T> operator*(const rsComplexSpectrogram<T>& rightOperand) const
  {  return rsComplexSpectrogram<T>(stftData.getElementwiseProduct(rightOperand.stftData)); }
  // having the multiplication operator do element-wise multiplication is one of the reasons why
  // we don't just inherit from rsMatrix and instead use a matrix as member variable - matrix
  // multiplication would make no sense for spectrograms

  //rsComplexSpectrogram<T> operator/(const rsComplexSpectrogram<T>& rightOperand) const
  //{  return rsComplexSpectrogram<T>(stftData.getElementwiseQuotient(rightOperand.stftData)); }


  // todo: implement +=,-=,*=,/= operators



protected:

  rsMatrix<std::complex<T>> stftData;

  // maybe this should also have additional data-members for hopSize, blockSize, sampleRate - these
  // values are needed in order to correctly interpret the data - but maybe not - it would not be
  // clear from which of the two object the result of +,-,.. etc. would get their values and we
  // don't want to be restricted to add,subtract,... spectrograms with matching block/hopSizes
  // ...or do we? it may most likely be an error trying to combine spectrograms with non-matched
  // values - maybe be a bit too restrictive - we may always relax restrictions later but not
  // tighten them when client code already relies on certain things
  // hmm - i think hopSize and blockSize should be members but sampleRate probably not

};

//=================================================================================================

/** This is .... under construction ...

-maybe rename to rsSpectrogramProcessor - a spectrogram itself is actually just a matrix of
 complex values...or magnitude/phase values....but what about re-assignment? there may be
 re-assigned spectrograms
-maybe let rsComplexSpectrogram be a private subclass of rsMatrix<double<T>>, implement the ()
 operator as simple delegation and replace getNumRows/getNumColumns with getNumFrames/getNumBins
-i think, i have implemented analysis and (re)synthesis in a single class rather than two separate
 classes (like with the sinusoidal model) in order to allow for identity resynthesis with arbitrary
 windows (the demodulation procedure needs to know both windows - analysis and synthesis)
-implement move contructors for rsArrayTools and rsMatrixOld  */

template<class T>
class rsSpectrogramProcessor
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsSpectrogramProcessor();

  /** Destructor. */
  ~rsSpectrogramProcessor();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the size (in samples) of the blocks/frames for analysis and synthesis (todo: allow
  different sizes for analysis and synthesis) */
  void setBlockSize(int newBlockSize);

  /** Sets the FFT size to be used. It should be >= the block size. If it's greater, the blocks
  will be zero padded to the desired FFT size. */
  void setTrafoSize(int newSize);

  /** Sets FFT size and block size at the same time (setting them simultaneously avoids potential
  temporary violations of trafoSize >= blockSize (which triggers an assertion) during setup). */
  void setBlockAndTrafoSize(int newBlockSize, int newTrafoSize);



  /** Sets the time-delta between successive blocks/frames in samples (value used for analysis and
  synthesis - todo: allow different values for both) */
  void setHopSize(int newHopSize) { hopSize = newHopSize; }

  /** Should be one of the type in RAPT::rsWindowFunction::windowTypes */
  void setAnalysisWindowType(rsWindowFunction::WindowType newType)
  {
    analysisWindowType = newType;
    updateAnalysisWindow();
  }

  void setSynthesisWindowType(rsWindowFunction::WindowType newType)
  {
    synthesisWindowType = newType;
    updateSynthesisWindow();
  }

  // todo: maybe have distinct setAnalysisBlockSize, setSynthesisBlockSize, setAnalysisHopSize, ..

  //updateAnalysisWindow();

  /** Sets the time origin for each analysis window to the center of the respective window. This is
  relevant for the phase spectrum. The natural time origin for each DFT buffer is the sample zero
  of the buffer - but the center often seems more meaningful. By default, this option is set to
  true. */
  void setTimeOriginAtWindowCenter(bool shouldBeAtCenter)
  {
    timeOriginAtWindowCenter = shouldBeAtCenter;
  }

  /** Selects whether or not the output signal should be demodulated by the product of the
  analysis- and synthesis window. When this product does not add up to unity when summed over all
  blocks, the raw resynthesized output will show an amplitude modulation given by these summed
  window products. The demodulation step is for undoing this modulation such that an identity
  analysis/resynthesis roundtrip is possible even with windows that do not satsify these nice
  "add up to unity" properties. With, for example, the (the "zn" version of) the Hanning window
  and a hop-size equal to blockSize/4, such demodulation should not be necessarry because in this
  case, the windows indeed do add up to unity, so you may save some computation time by not doing
  this step. That's why it can be switched on and off by client code - and for testing purposes.

  ToDo: even with the zn-hanning window it doesn't yet work right because it currently seems to
  add up not to unity but to some other constant (1.5 in case of H = B/4 - we seem to miss a
  scaling somewhere..... also, at the moment, there are transient artifacts at the start and
  end (fade in/out due to not overlapping enough windows): workaround: prepend and append one
  block-length of zero samples (which can be removed after resynthesis). */
  void setOutputDemodulation(bool shouldDemodulate)
  {
    demodulateOutput = shouldDemodulate;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  //size_t getFftSize() const { return blockSize * zeroPaddingFactor; }

  //size_t getFftSize() const { return trafoSize; }
  int getFftSize() const { return trafoSize; }
  // rename to getTrafoSize


  //size_t getNumNonRedundantBins() const { return getFftSize()/2 + 1; }
  int getNumNonRedundantBins() const { return getFftSize()/2 + 1; }

  int getHopSize() const { return hopSize; }

  /** Converts a physical frequency in Hz to an FFT bin index for the given sample-rate. This
  depends on our trafoSize setting and will return a real number. Client code should round the
  returned itself as appropriate. */
  T frequencyToBinIndex(T freq, T sampleRate) const { return trafoSize*freq/sampleRate; }

  /** Returns the amplitude scale factor by which the outputs of shortTimeSpectrum() has to be
  scaled to obtain the actual amplitude of a real sinusoid. This is used internally but made
  available mostly for testing purposes. */
  T getAnalysisScaler() const { return 2. / rsArrayTools::sum(&analysisWindow[0], blockSize); }
    // the factor two comes from the fact that two complex sinusoids sum up to one real sinusoid
    // and we count only one of the two symmetric magnitude values

  /** Computes the required number of frames, given the length of the signal and the hopsize.
  We assume that the peak value of the first window w[B/2] coincides with the 1st sample x[0] and
  the peak of the last window either coincides with the last sample x[N-1] or is behind the last
  sample. This ensures that each sample is covered within the range where overlapping windows
  add up to a constant (assuming a properly chosen window shape).  */
  static int getNumFrames(int numSamples, int hopSize);

  /** Assuming that the analysis-window wa and synthesis-windows ws are chosen in a way to sum up
  to a constant for given blocksize B and hopsize H, this function computes that constant and
  returns it. */
  static T getWindowSum(T *wa, T *ws, int B, int H);


  /** Computes the number of required zero valued samples to prepend and append for given block
  size and transform size. We use a zero-padding such that the actual audio data sits in the middle
  of a transform buffer with zero padding to the left and right */
  static void getLeftRightPaddingAmount(int blockSize, int trafoSize, int* left, int* right);



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  //rsMatrixOld<std::complex<T>> getComplexSpectrogram(const T *signal, int numSamples);

  /** Computes the complex spectrogram of the given signal x of the given length in samples. */
  rsMatrix<std::complex<T>> getComplexSpectrogram(const T* signal, int numSamples);
  // rename to signalToComplexSpectrogram


  /** Computes a short-time FFT spectrum ... */
  void shortTimeSpectrum(const T* signal, int numSamples, int blockCenter,
    std::complex<T>* spectrum);

  /** Writes a segment from the signal x (assumend to be of length N) centered at n into the
  complex buffer X - completely with windowing, zero-padding and shifting (if applicable), ready
  to pass to the FFT routine. Used internally by shortTimeSpectrum. */
  void prepareTrafoBuffer(const T* x, int N, int n, std::complex<T> *X);


  /** Computes the amplitude modulation signal that would be imposed on a signal when doing an
  analysis/synthesis roundtrip with the given analysis- and synthesis-windows and block- and
  hop-size. This can be used for demodulating a resynthesized signal, making the
  analysis/resynthesis/demodulation roundtrip a identity-operation (up to roundoff).
  Note: In typical phase-vocoder implementations, these window-shapes and block- and hopsizes are
  chosen such that this amplitude modulation signal comes out as a constant value (or, even better,
  unity). To make the implementation more flexible with regard to the choice of these parameters,
  this function can be can be used to get the modulation signal and divide by it for
  demodulation. */
  std::vector<T> getRoundTripModulation(int numFrames);

  // is this formula also correct for odd fft sizes? verify?


  //-----------------------------------------------------------------------------------------------
  /** \name Static Processing Functions */

  /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero
  value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again.
  That means, the window has a period length of N. Such a window is suitable for applications
  where it is important that suitably overlapped windows sum up to a constant, like when identity
  resynthesis is required. */
  //static void hanningWindowZN(T *w, int N);
  // maybe have a version NZ, ZZ, NN
  // moved to rsWindowFunction


  /** Zeroes out all bins above "highestBinToKeep" (in all frames) */
  static void lowpass(rsMatrix<std::complex<T>>& spectrogram, int highestBinToKeep);
  // todo: take as argument rsComplexSpectrogram instead of rsMatrix

  //static void lowpass(rsMatrix<std::complex<T>>& spectrogram, T cutoffBin);
  // not yet tested


  /** Zeroes out all bins below "lowestBinToKeep" (in all frames) */
  static void highpass(rsMatrix<std::complex<T>>& spectrogram, int lowestBinToKeep);

  //static void highpass(rsMatrix<std::complex<T>>& spectrogram, T cutoffBin);

  /** Zeroes out all bin-values values except those in the range lowBin <= bin <= highBin. If
  lowBin == 0 or highBin == numBins-1, you can also get lowpass- and highpass-filters,
  respectively. */
  static void bandpass(rsMatrix<std::complex<T>>& spectrogram, int lowBin, int highBin);





  /** Given a complex spectrogram, this function synthesizes a signal using given synthesis
  window, blockSize and hopSize. You must also pass the analysis window that was used - this is
  needed internally for amplitude demodulation.
  a synthesis
  window w of length B (which is the blocksize) with hopsize H. The number of equals the number
  of columns in the matrix s - each row is one short-time spectrum (of positive frequencies only
  due to symmetry). */
  std::vector<T> synthesize(const rsMatrix<std::complex<T>>& spectrogram);
  // rename to spectrogramToSignal

  /** Given a signal and its complex spectrogram, this function computes the matrix of time
  reassignments for each time/frequency value. The rampedWindow array should be a time-ramped
  version of the window that was used in the sepctrogram analysis. */
  static rsMatrix<T> timeReassignment(T *signal, int numSamples,
    const rsMatrix<std::complex<T>> &spectrogram, T *rampedWindow, int blockSize,
    int hopSize);
  // not yet implemented - probably needs to be made non-static

  /** Given a signal and its complex spectrogram, this function computes the matrix of frequency
  reassignments for each time/frequency value. The derivativeWindow array should be the derivative
  of the window that was used in the sepctrogram analysis. */
  static rsMatrix<T> frequencyReassignment(T *signal, int numSamples,
    const rsMatrix<std::complex<T>> &spectrogram, T *derivativeWindow, int blockSize,
    int hopSize);
  // not yet implemented - probably needs to be made non-static

  /** Used inside synthesize() - this is the raw resynthesis in the sense that any amplitude
  modulation artifacts that result from overlapping grains that don't sum up to unity are not
  yet compensated for. If you want that compensation, use synthesize(). For properly chosen
  analysis and synthesis windows and hop-size and block-size, the demodulation step can be
  legitimately skipped because the overlapped windows add up to unity - but this is not the case
  in general. */
  std::vector<T> synthesizeRaw(const rsMatrix<std::complex<T>> &spectrogram);



  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Sets up the analysis/synthesis parameters to default values. */
  //void init();


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Internal */

  void updateAnalysisWindow();

  void updateSynthesisWindow();

  static void fillWindowArray(T* w, int length, rsWindowFunction::WindowType type);

  /** Given a buffer of complex values, this functions swaps the first and second half which is the
  same as a circular shift by half of the buffer size. This shifts the time-origin of buffer from
  sample zero to the center of the buffer. That means that the center of the buffer becomes the
  reference for measuring the phase - a sinusoid that starts at the center sample will have zero
  phase. This function is applied before transforming a windowed grain to the frequency domain and
  after syntehsizing a grain from a short time spectrum. */
  void swapForZeroPhase(std::complex<T>* X, int length);
  // todo: comment about the even/odd window size issue, maybe rename to shiftOriginToCenter or
  // something
  // make a unit test that ensures that applying this function twice gives an identity operation
  // -if this is not the case for odd windows, we may have to use two functions - one for forward,
  // the other for backward shift ..maybe shiftZeroToCenter, shiftCenterToZero
  // or applyAnalysisPhaseAdjustment, applySynthesisPhaseAdjustment ..or TimeAdjustment
  // of shiftForAnalysis, shiftForSynthesis

  /** Performs the forward FFT on given complex buffer X of length M */
  void fft(std::complex<T>* X, int M);

  /** Performs the inverse FFT on given complex buffer X of length M */
  void ifft(std::complex<T>* X, int M);


  /** \name Data */

  int trafoSize = 0;    // FFT size - initialized in constructor
  int blockSize = 0;    // initialized in constructor which also generates the window functions
  int hopSize   = 128;
  // maybe we should also distiguish between analysis and synthesis hop- and block-size

  rsWindowFunction::WindowType analysisWindowType  = rsWindowFunction::WindowType::hanningZN;
  rsWindowFunction::WindowType synthesisWindowType = rsWindowFunction::WindowType::hanningZN;

  std::vector<T> analysisWindow, synthesisWindow;

  bool timeOriginAtWindowCenter = true;
  bool demodulateOutput = true;

  /** The Fourier transformer object. */
  rsFourierTransformerBluestein<T> transformer;

  // buffers used for swapping first and second half of window for zero-phase windowing
  //std::vector<std::complex<T>> swapBufferAna, swapBufferSyn;
  //std::vector<std::complex<T>> swapBuffer;
};

#endif
