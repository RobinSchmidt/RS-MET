#ifndef RAPT_PHASEVOCODER_H
#define RAPT_PHASEVOCODER_H

/** This is .... under construction ...


-implement move contructors for rsArray and rsMatrix  */

template<class T>
class rsSpectrogram
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsSpectrogram();

  /** Destructor. */
  ~rsSpectrogram();


  /** \name Setup */

  /** Sets the size (in samples) of the blocks/frames for analysis and synthesis (todo: allow 
  different sizes for analysis and synthesis) */
  void setBlockSize(int newBlockSize);

  /** Sets the time-delta between successive blocks/frames in samples (value used for analysis and
  synthesis - todo: allow different values for both) */
  void setHopSize(int newHopSize) { hopSize = newHopSize; }

  /** Sets the amout of zero-padding as an integer factor. 
  ...padding should be applied symmetrically at front and back such that the phase stays measured
  at the center of the frame */
  void setZeroPaddingFactor(int newFactor) { zeroPaddingFactor = newFactor; }
  // get rid - replace by setTransformSize ot setFftSize

  void setTrafoSize(int newSize) { trafoSize = newSize; }


  /** Should be one of the type in RAPt::rsWindowFunction::windowTypes */
  void setAnalysisWindowType(int newType) 
  { 
    analysisWindowType = newType; 
    updateAnalysisWindow();
  }

  void setSynthesisWindowType(int newType) 
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



  /** \name Inquiry */



  size_t getFftSize() const { return blockSize * zeroPaddingFactor; }

  size_t getNumNonRedundantBins() const { return getFftSize()/2 + 1; }

  int getHopSize() const { return hopSize; }






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



  /** \name Audio Processing */

  //rsMatrix<std::complex<T>> complexSpectrogram(const T *signal, int numSamples);


  /** Computes a short-time FFT spectrum ... */
  void shortTimeSpectrum(const T* signal, int numSamples, int blockCenter, std::complex<T>* spectrum);

  /** Computes the complex spectrogram of the given signal x of the given length in samples. */
  rsMatrix<std::complex<T>> complexSpectrogram(const T* signal, int numSamples);




  /** \name Static Processing Functions */

  /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero
  value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again.
  That means, the window has a period length of N. Such a window is suitable for applications
  where it is important that suitably overlapped windows sum up to a constant, like when identity
  resynthesis is required. */
  //static void hanningWindowZN(T *w, int N);
  // maybe have a version NZ, ZZ, NN
  // moved to rsWindowFunction




  /** Given a complex spectrogram, this function synthesizes a signal using given synthesis
  window, blockSize and hopSize. You must also pass the analysis window that was used - this is
  needed intrenally for amplitude demodulation.
  a synthesis
  window w of length B (which is the blocksize) with hopsize H. The number of equals the number
  of columns in the matrix s - each row is one short-time spectrum (of positive frequencies only
  due to symmetry). */
  std::vector<T> synthesize(const rsMatrix<std::complex<T>> &spectrogram);

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
  legitimately skipped because the overallped windows add up to unity - but this is not the case
  in general. */
  std::vector<T> synthesizeRaw(const rsMatrix<std::complex<T>> &spectrogram);

  /** Computes the amplitude modulation signal that would be imposed on a signal when doing an
  analysis/synthesis roundtrip with the given analysis- and synthesis-windows and block- and
  hop-size. This can be used for demodulating a resynthesized signal, making the
  analysis/resynthesis/demodulation roundtrip a identity-operation (up to roundoff).
  Note: In typical phase-vocoder implementations, these window-shapes and block- and hopssizes are
  chosen such that this amplitude modulation signal comes out as a constant value (or, even better,
  unity). To make the implementation more flexible with regard to the choice of these parameters,
  this function can be can be used to get the modulation signal and divide by it for
  demodulation. */
  std::vector<T> getRoundTripModulation(int numFrames);


    // is this formula also correct for odd fft sizes? verify?


  /** \name Misc */

  /** Sets up the analysis/synthesis parameters to default values. */
  //void init();


protected:

  /** \name Internal */

  void updateAnalysisWindow();

  void updateSynthesisWindow();

  static void fillWindowArray(T* w, int length, int type);

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


  int zeroPaddingFactor = 1;
  int blockSize = 0;    // initialized in constructor which also generates the window functions   
  int hopSize   = 128;
  // int trafoSize = 0;  // FFT size - use later
  // maybe we should also distiguish between analysis and synthesis hop-and block-size


  int analysisWindowType  = rsWindowFunction::HANNING_WINDOW_ZN;
  int synthesisWindowType = rsWindowFunction::HANNING_WINDOW_ZN;

  std::vector<T> analysisWindow, synthesisWindow;


  bool timeOriginAtWindowCenter = true;

  /** The Fourier transformer object. */
  rsFourierTransformerBluestein<T> transformer;

  // buffers used for swapping first and second half of window for zero-phase windowing
  //std::vector<std::complex<T>> swapBufferAna, swapBufferSyn;
  //std::vector<std::complex<T>> swapBuffer;
};

#endif
