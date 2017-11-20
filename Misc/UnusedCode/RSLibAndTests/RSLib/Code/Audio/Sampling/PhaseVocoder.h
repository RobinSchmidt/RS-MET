#ifndef RS_PHASEVOCODER_H
#define RS_PHASEVOCODER_H

namespace RSLib
{

  /**

  This is ....


  // todo:
  // maybe factor out a class rsSpectrogramProcessor - the phase vocoder is actually a higher level
  // of analysis on top of the spectrogram
  // implement move contructors for rsArray and rsMatrix

  */

  class RSLib_API rsPhaseVocoder
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsPhaseVocoder();  

    /** Destructor. */
    ~rsPhaseVocoder();  


    /** \name Setup */



    /** \name Inquiry */


    /** Computes the required number of frames, given the length of the signal and the hopsize.
    We assume that the peak value of the first window w[B/2] coincides with the 1st sample x[0] and
    the peak of the last window either coincides with the last sample x[N-1] or is behind the last 
    sample. This ensures that each sample is covered within the range where overlapping windows 
    add up to a constant (assuming a properly chosen window shape).  */
    static int getNumFrames(int numSamples, int hopSize);

    /** Assuming that the analysis-window wa and synthesis-windows ws are chosen in a way to sum up
    to a constant for given blocksize B and hopsize H, this function computes that constant and 
    returns it. */
    static double getWindowSum(double *wa, double *ws, int B, int H);



    /** \name Audio Processing */

    /** Creates a Hanning window that starts with a zero value in w[0] and ends with a nonzero 
    value in w[N-1] = w[1], such that the nominal and nonexistent value w[N] would be zero again. 
    That means, the window has a period length of N. Such a window is suitable for applications 
    where it is important that suitably overlapped windows sum up to a constant, like in the 
    phase-vocoder. */
    static void hanningWindowZN(double *w, int N);

    /** Computes a short-time FFT spectrum ... */
    static void shortTimeSpectrum(double *signal, int numSamples, int blockCenter, double *window, 
      int blockSize, int fftSize, rsComplexDbl *spectrum);

    /** Computes the complex spectrogram of the given signal using the given analysis window (the 
    length of which should be equal to the given blockSize), using the given hopSize. The optional
    padFactor parameter determines the amount of zero-padding for the FFT (it should be chosen such 
    that fftSiize = padFactor*blockSize is a power of two). */
    static rsMatrix<rsComplexDbl> complexSpectrogram(double *signal, int numSamples, 
      double *window, int blockSize, int hopSize, int padFactor = 1);

    /** Given a complex spectrogram, this function synthesizes a signal using given synthesis 
    window, blockSize and hopSize. You must also pass the analysis window that was used - this is
    needed intrenally for amplitude demodulation.
    a synthesis 
    window w of length B (which is the blocksize) with hopsize H. The number of equals the number 
    of columns in the matrix s - each row is one short-time spectrum (of positive frequencies only 
    due to symmetry). */
    static std::vector<double> synthesize(const rsMatrix<rsComplex<double>> &spectrogram,
      double *synthesisWindow, int blockSize, int hopSize, double *analysisWindow);

    /** Given a signal and its complex spectrogram, this function computes the matrix of time 
    reassignments for each time/frequency value. The rampedWindow array should be a time-ramped
    version of the window that was used in the sepctrogram analysis. */
    static rsMatrix<double> timeReassignment(double *signal, int numSamples,
      const rsMatrix<rsComplex<double>> &spectrogram, double *rampedWindow, int blockSize, 
      int hopSize);

    /** Given a signal and its complex spectrogram, this function computes the matrix of frequency
    reassignments for each time/frequency value. The derivativeWindow array should be the derivative 
    of the window that was used in the sepctrogram analysis. */
    static rsMatrix<double> frequencyReassignment(double *signal, int numSamples,
      const rsMatrix<rsComplex<double>> &spectrogram, double *derivativeWindow, int blockSize, 
      int hopSize);

    /** Used inside synthesize() - this is the raw resynthesis in the sense that any amplitude 
    modulation artifacts that result from overlapping grains that don't sum up to unity are not 
    yet compensated for. If you want that compensation, use synthesize().  */
    static std::vector<double> synthesizeRaw(const rsMatrix<rsComplex<double>> &spectrogram, 
      double *window, int blockSize, int hopSize);

    /** Computes the amplitude modulation signal that would be imposed on a signal when doing an
    analysis/synthesis roundtrip with the given analysis- and synthesis-windows and block- and 
    hop-size. This can be used for demodulating a resynthesized signal, making the 
    analysis/resynthesis/demodulation roundtrip a unity-operation.
    Note: In typical phase-vocoder implementations, these window-shapes and block- and hopssize are
    chose such that this amplzide modulation signal comes out as a constant value (or, even better,
    unity). To make the implementation more flexible with regard to the choice of these parameters,
    this function can be can be used to get the modulation signal and divide by it for 
    demodulation. */
    static std::vector<double> getModulation(double *analysisWindow, double *synthesisWindow,
      int blockSize, int hopSize, int numFrames);


    /** \name Misc */

    /** Sets up the analysis/synthesis parameters to default values. */
    void init();


  protected:

    /** \name Internal */


    /** \name Data */

    //double fs;   // samplerate

    //int Nb;      // block size
    //int Nh;      // hop size
    //int Nf;      // FFT size

    // Ba, Ha, Ka: // analysis blocksize, hopsize, number of frequencies

    // FourierTransformerBluestein 


  };

}

#endif
