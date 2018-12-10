#ifndef RAPT_PHASEVOCODER_H
#define RAPT_PHASEVOCODER_H

/** This is .... under construction ...

todo:
-maybe factor out a class rsSpectrogramProcessor - the phase vocoder is actually a higher level
 of analysis on top of the spectrogram
 -actually, the original phase vocoder relies on multiplication of the signal with various complex
  sinusoids ("heterodyning") and a filter bank - maybe rename the class here to rsSpectrogram and 
  additionally implement a true phase vocoder (see DAFX, Ch. 8 - Time-frequency processing)
  -maybe that true phase vocoder may be extended by allowing the complex sinusoids to be not 
  necessarily harmonic and vary in frequency (and maybe also allowing the filters to have 
  time-variying bandwidths) - this could be used to track the frequency of the partial being 
  analyzed (but requires a preliminary knowlegde of the frequency trajectory of the partial=
  ...maybe we can firsat obtain preliminary (rough) frequency tracks by spectrogram processing and
  the refine them by a time-varying phase vocoder approach...and then use that data to refine 
  further etc. until it converges
-implement move contructors for rsArray and rsMatrix  */

template<class T>
class rsPhaseVocoder
{

public:

  enum windowTypes
  {
    RECTANGULAR_WINDOW = 0,
    HANNING_WINDOW
    // HAMMING_WINDOW,
    // BLACKMAN_HARRIS_WINDOW,
  };

  /** \name Construction/Destruction */

  /** Constructor. */
  rsPhaseVocoder();

  /** Destructor. */
  ~rsPhaseVocoder();


  /** \name Setup */

  void setBlockSize(int newBlockSize) { blockSize = newBlockSize; }

  void setHopSize(int newHopSize) { hopSize = newHopSize; }

  void setZeroPaddingFactor(int newFactor) { zeroPaddingFactor = newFactor; }



  /** \name Inquiry */



  size_t getFftSize() const { return blockSize * zeroPaddingFactor; }

  size_t getNumNonRedundantBins() const { return getFftSize()/2 + 1; }






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

  /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero
  value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again.
  That means, the window has a period length of N. Such a window is suitable for applications
  where it is important that suitably overlapped windows sum up to a constant, like in the
  phase-vocoder. */
  static void hanningWindowZN(T *w, int N);

  /** Computes a short-time FFT spectrum ... */
  static void shortTimeSpectrum(T *signal, int numSamples, int blockCenter, T *window,
    int blockSize, int fftSize, std::complex<T> *spectrum);

  /** Computes the complex spectrogram of the given signal using the given analysis window (the
  length of which should be equal to the given blockSize), using the given hopSize. The optional
  padFactor parameter determines the amount of zero-padding for the FFT (it should be chosen such
  that fftSiize = padFactor*blockSize is a power of two). */
  static rsMatrix<std::complex<T>> complexSpectrogram(T *signal, int numSamples,
    T *window, int blockSize, int hopSize, int padFactor = 1);

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
  analysis/resynthesis/demodulation roundtrip a unity-operation.
  Note: In typical phase-vocoder implementations, these window-shapes and block- and hopssize are
  chose such that this amplzide modulation signal comes out as a constant value (or, even better,
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


  /** \name Data */

  int blockSize = 512;
  int hopSize   = 128;
  int zeroPaddingFactor = 1;
  int windowType = HANNING_WINDOW;


  //T fs;   // samplerate

  //int Nb;      // block size
  //int Nh;      // hop size
  //int Nf;      // FFT size

  // Ba, Ha, Ka: // analysis blocksize, hopsize, number of frequencies

  // FourierTransformerBluestein 


};

#endif
