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
  void setBlockSize(int newBlockSize) 
  { 
    if(newBlockSize != blockSize) {
      blockSize = newBlockSize;
      updateAnalysisWindow();
      updateSynthesisWindow();
      // update fourier transformer(s)
    }
  }

  /** Sets the time-delta between successive blocks/frames in samples (value used for analysis and
  synthesis - todo: allow different values for both) */
  void setHopSize(int newHopSize) { hopSize = newHopSize; }

  /** Sets the amout of zero-padding as an integer factor. 
  ...padding should be applied symmetrically at front and back such that the phase stays measured
  at the center of the frame */
  void setZeroPaddingFactor(int newFactor) { zeroPaddingFactor = newFactor; }

  /** Should be one of the type in RAPt::rsWindowFunction::windowTypes */
  void setAnalysisWindowType(int newType) 
  { 
    analysisWindowType = newType; 
    updateAnalysisWindow();
  }

  //updateAnalysisWindow();



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

  rsMatrix<std::complex<T>> complexSpectrogram(const T *signal, int numSamples) const;



  /** \name Static Processing Functions */

  /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero
  value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again.
  That means, the window has a period length of N. Such a window is suitable for applications
  where it is important that suitably overlapped windows sum up to a constant, like when identity
  resynthesis is required. */
  //static void hanningWindowZN(T *w, int N);
  // maybe have a version NZ, ZZ, NN

  /** Computes a short-time FFT spectrum ... */
  static void shortTimeSpectrum(const T* signal, int numSamples, int blockCenter, const T* window,
    int blockSize, int fftSize, std::complex<T> *spectrum);

  /** Computes the complex spectrogram of the given signal using the given analysis window (the
  length of which should be equal to the given blockSize), using the given hopSize. The optional
  padFactor parameter determines the amount of zero-padding for the FFT (it should be chosen such
  that fftSiize = padFactor*blockSize is a power of two). */
  static rsMatrix<std::complex<T>> complexSpectrogram(const T* signal, int numSamples,
    const T* window, int blockSize, int hopSize, int padFactor = 1);

  /** Given a complex spectrogram, this function synthesizes a signal using given synthesis
  window, blockSize and hopSize. You must also pass the analysis window that was used - this is
  needed intrenally for amplitude demodulation.
  a synthesis
  window w of length B (which is the blocksize) with hopsize H. The number of equals the number
  of columns in the matrix s - each row is one short-time spectrum (of positive frequencies only
  due to symmetry). */
  static std::vector<T> synthesize(const rsMatrix<std::complex<T>> &spectrogram,
    T *synthesisWindow, int blockSize, int hopSize, T *analysisWindow);

  /** Given a signal and its complex spectrogram, this function computes the matrix of time
  reassignments for each time/frequency value. The rampedWindow array should be a time-ramped
  version of the window that was used in the sepctrogram analysis. */
  static rsMatrix<T> timeReassignment(T *signal, int numSamples,
    const rsMatrix<std::complex<T>> &spectrogram, T *rampedWindow, int blockSize,
    int hopSize);

  /** Given a signal and its complex spectrogram, this function computes the matrix of frequency
  reassignments for each time/frequency value. The derivativeWindow array should be the derivative
  of the window that was used in the sepctrogram analysis. */
  static rsMatrix<T> frequencyReassignment(T *signal, int numSamples,
    const rsMatrix<std::complex<T>> &spectrogram, T *derivativeWindow, int blockSize,
    int hopSize);

  /** Used inside synthesize() - this is the raw resynthesis in the sense that any amplitude
  modulation artifacts that result from overlapping grains that don't sum up to unity are not
  yet compensated for. If you want that compensation, use synthesize().  */
  static std::vector<T> synthesizeRaw(const rsMatrix<std::complex<T>> &spectrogram,
    T *window, int blockSize, int hopSize);

  /** Computes the amplitude modulation signal that would be imposed on a signal when doing an
  analysis/synthesis roundtrip with the given analysis- and synthesis-windows and block- and
  hop-size. This can be used for demodulating a resynthesized signal, making the
  analysis/resynthesis/demodulation roundtrip a identity-operation (up to roundoff).
  Note: In typical phase-vocoder implementations, these window-shapes and block- and hopssizes are
  chosen such that this amplzide modulation signal comes out as a constant value (or, even better,
  unity). To make the implementation more flexible with regard to the choice of these parameters,
  this function can be can be used to get the modulation signal and divide by it for
  demodulation. */
  static std::vector<T> getModulation(T *analysisWindow, T *synthesisWindow,
    int blockSize, int hopSize, int numFrames);


    // is this formula also correct for odd fft sizes? verify?


  /** \name Misc */

  /** Sets up the analysis/synthesis parameters to default values. */
  //void init();


protected:

  /** \name Internal */

  void updateAnalysisWindow();

  void updateSynthesisWindow();

  static void fillWindowArray(T* w, int length, int type);

  /** \name Data */

  int blockSize = 0;    // initialized in constructor which also generates the window functions   
  int hopSize   = 128;
  // maybe we should also distiguish between analysis and synthesis hop-and block-size

  int zeroPaddingFactor = 1;
  int analysisWindowType  = rsWindowFunction::HANNING_WINDOW_ZN;
  int synthesisWindowType = rsWindowFunction::HANNING_WINDOW_ZN;

  std::vector<T> analysisWindow, synthesisWindow;




  //T fs;   // samplerate

  //int Nb;      // block size
  //int Nh;      // hop size
  //int Nf;      // FFT size

  // Ba, Ha, Ka: // analysis blocksize, hopsize, number of frequencies

  //rsFourierTransformerBluestein<T> analysisTransformer, synthesisTransformer;


};

#endif
