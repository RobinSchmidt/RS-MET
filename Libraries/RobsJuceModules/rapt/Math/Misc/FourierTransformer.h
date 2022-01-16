#ifndef RAPT_FOURIERTRANSFORMERRADIX2_H
#define RAPT_FOURIERTRANSFORMERRADIX2_H

/** This class performs a fast Fourier Transform on a block of complex numbers. The length of the 
block has to be a power of 2. */

template<class T>
class rsFourierTransformerRadix2
{

public:

  /** The direction of the transform. */
  enum directions
  {
    FORWARD,
    INVERSE
  };

  /** These are the possible normalization modes. */
  enum normalizationModes
  {
    NORMALIZE_ON_FORWARD_TRAFO, // divide by blockSize on forward FFT
    NORMALIZE_ON_INVERSE_TRAFO, // divide by blockSize on inverse FFT (default)
    ORTHONORMAL_TRAFO,          // divide by sqrt(blockSize) on both transforms
    NEVER_NORMALIZE             // no normalization at all
  };


  /** \name Construction/Destruction */

  /** Constructor. */
  rsFourierTransformerRadix2();

  /** Destructor. */
  ~rsFourierTransformerRadix2();


  /** \name Setup */

  /** FFT-size, has to be a power of 2 and >= 2. */
  void setBlockSize(int newBlockSize);

  /** Sets the direction of the transform (@see: directions). This will affect the sign of the
  exponent (or equivalently: theimaginary part) in the twiddling factors and the normalization
  constant. */
  void setDirection(int newDirection);

  /** When you switch bewteen usage of this object for real or complex signals, you will need to
  call this switch-function in between which just triggers a re-computation of the twiddle
  factors (which must be different for the two cases). */
  //void setRealSignalMode(bool willBeUsedForRealSignals);
    // check, if this is obsolete

  /** Sets the mode for normalization of the output (@see: normalizationModes). */
  void setNormalizationMode(int newNormalizationMode);


  /** \name Inquiry */

  /** Returns the size of FFT blocks. */
  int getBlockSize() const { return N; }


  /** \name Transforms */

  /** Transforms a buffer of complex numbers into its (forward or inverse) fourier transform.
  The inBuffer will remain intact. Both, inBuffer and outBuffer must be of the size which was
  specified when setting up the blockSize with setBlockSize(). */
  void transformComplexBuffer(const std::complex<T> *inBuffer, std::complex<T> *outBuffer);

  /** Does the same thing as transformComplexBuffer but performes in-place calculation
  (overwrites the input buffer). */
  void transformComplexBufferInPlace(std::complex<T> *buffer);

  /** Transforms a real signal into its corresponding (conjugate symmetric) complex spectrum
  using an algorithm which exploits the symmetries for more efficiency. When the input array is
  an array of Ts of length N, the output array will be an array of complex numbers (class
  Complex) of length N/2 with the (purely real) transform value of bin N/2 stored in the
  imaginary part of the first array entry (outSpectrum[0].im = real(N/2)). */
  void transformRealSignal(const T *inSignal, std::complex<T> *outSpectrum);

  /** Calculates real and imaginary part of the spectrum as interleaved T buffer:
  buf[2]=re[1], buf[3]=im[1], buf[4]=re[2], buf[5]=im[2],... in general: buf[2*k]=re[k],
  buf[2k+1]=im[k], k=1,..,(N/2)-1 where N is the FFT-size. The first two elements of the buffer
  have a special meaning: buf[0] is the (purely real) DC and buf[1] is the (purely real)
  coefficient for the Nyquist frequency. The other fields contain the real and imaginary parts of
  the positive frequencies only (interleaved) because the negative frequencies are redundant
  (they are conjugate symmetric). */
  void transformRealSignal(const T *signal, T *reAndIm);

  /** Calculates spectral magnitudes and phases from a signal, where *signal should be of
  length N, where N is the block-size as chosen with setBlockSize() *magnitudes and *phases
  should be of length N/2. The first values of the output arrays have a special meaning:
  magnitudes[0] is the (purely real) DC and phases[0] is the (purely real) coefficient for the
  Nyquist frequency. */
  void getRealSignalMagnitudesAndPhases(const T *signal, T *magnitudes, T *phases);
    // hmm...this is a somewhat odd convention

  /** Calculates the magnitudes only from a signal (useful for analyzer-stuff). */
  void getRealSignalMagnitudes(const T *signal, T *magnitudes);

  /** Transforms a complex conjugate symmetric spectrum (i.e. a spectrum of a real signal) into
  the corresponding real signal. */
  void transformSymmetricSpectrum(const std::complex<T> *inSpectrum, T *outSignal);

  /** Calculates a time signal from and interleaved buffer containing the real and imaginary
  parts of the positive frequencies (the negative frequencies are assumed to be conjugate
  symmetric). */
  void transformSymmetricSpectrum(const T *reAndIm, T *signal);

  /** Calculates a real time signal from its magnitudes and phases, *magnitudes and *phases
  should be of length N/2, *signal is of length N where N is the block-size as chosen with
  setBlockSize(). The first values of the input arrays have a special meaning: magnitudes[0] is
  the (purely real) DC and phases[0] is the (purely real) coefficient for the Nyquist
  frequency. */
  void getRealSignalFromMagnitudesAndPhases(const T *magnitudes, const T *phases, T *signal);


  /** \name Static Member Functions */

  /** Returns the physical frequency in Hz that corresponds to the given 'binIndex' for a given
  'fftSize' and 'sampleRate'. */
  static T binIndexToFrequency(int binIndex, int fftSize, T sampleRate)
  {
    return binIndex*sampleRate/fftSize;
  }

  static T binIndexToOmega(int binIndex, int fftSize)
  {
    return T(2*PI*binIndex)/fftSize;
  }

  /** Returns the normalization factor to applied after the forward and/or inverse transform, 
  according to the block-size, direction of the transform and normalization mode. */
  static T getNormalizationFactor(int blockSize, int direction, int normalizationMode);


protected:

  /** \name Internal Functions */

  /** Updates the normalizationFactor member variable acording to a new blockSize, direction or
  normalizationMode. */
  void updateNormalizationFactor();




  /** \name Data */

  int N;                  /**< blocksize of the FFT. */
  int logN;               /**< Base 2 logarithm of the blocksize. */
  int direction;          /**< The direction of the transform (@see: directions). */
  int normalizationMode;  /**< The normalization mode (@see: normalizationModes. */
  T normalizationFactor;  /**< The normalization factor (can be 1, 1/N or 1/sqrt(N)). */

  // work-area stuff for Ooura's fft-routines:
  T *w;        /**< Table of the twiddle-factors. */
  int *ip;     /**< Work area for bit-reversal (index pointer?). */

  // our own temporary storage area:
  std::complex<T>* tmpBuffer;

};

//===============================================================================================

/** This class works similar to the FourierTransformerRadix2 class, except that it accepts 
arbitrary blockSizes. It uses Bluestein's algorithm to evaluate the DFT of a given block which 
works for arbitrary block-sizes and preserves the desirable O(N*log(N)) complexity of FFT 
algorithms (although the implicit constants are larger than for radix-2 algorithms).

\todo implement normalization */

template<class T>
class rsFourierTransformerBluestein
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  rsFourierTransformerBluestein();

  /** Destructor. */
  ~rsFourierTransformerBluestein();

  //---------------------------------------------------------------------------------------------
  // \name Setup

  /** FFT-size, can be an arbitrary integer > 1. Powers of 2 will be most efficient, other
  numbers increase the CPU-load significantly, but O(N*log(N))-scaling is preserved. */
  void setBlockSize(int newBlockSize);

  /** Sets the direction of the transform (@see: FourierTransformerRadix2::directions). This
  will affect the sign of the exponent (or equivalently: theimaginary part) in the twiddling
  factors and the normalization constant. */
  void setDirection(int newDirection);

  /** Sets the mode for normalization of the output
  (@see: FourierTransformerRadix2::normalizationModes). */
  void setNormalizationMode(int newNormalizationMode);

  //---------------------------------------------------------------------------------------------
  // \name Processing

  /** Transforms a buffer of complex numbers into its (forward or inverse) fourier transform.
  The inBuffer will remain intact. Both, inBuffer and outBuffer must be of the size which was
  specified when setting up the blockSize with setBlockSize(). */
  void transformComplexBuffer(std::complex<T> *inBuffer, std::complex<T> *outBuffer);

  /** Does the same thing as transformComplexBuffer but performes in-place calculation
  (overwrites the input buffer). */
  void transformComplexBufferInPlace(std::complex<T> *buffer);

  /** Convenience function to be used if you need just a single FFT. If you need a full series
  of FFTs (i.e. multiple blocks), this is not recommended for efficiency reasons. */
  static void fft(std::complex<T>* buffer, int length, bool normalize)
  {
    rsFourierTransformerBluestein<T> trafo;
    trafo.setBlockSize(length);
    trafo.setDirection(rsFourierTransformerRadix2<T>::FORWARD);
    if(normalize)
      trafo.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_FORWARD_TRAFO);
    trafo.transformComplexBufferInPlace(&buffer[0]);
  }
  // needs test

  /** Convenience function for a single inverse FFT. */
  static void ifft(std::complex<T>* buffer, int length, bool normalize)
  {
    rsFourierTransformerBluestein<T> trafo;
    trafo.setBlockSize(length);
    trafo.setDirection(rsFourierTransformerRadix2<T>::INVERSE);
    if(normalize)
      trafo.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_INVERSE_TRAFO);
    trafo.transformComplexBufferInPlace(&buffer[0]);
  }
  // needs test

  //=============================================================================================

protected:

  /** Generates the chirp signal which is needed modulate the input signal (and to retrieve the
  output-spectrum), as well as the spectrum of the conjugate chirp signal which is needed to
  multiply the spectrum of the chirp-modulated signal. */
  void generateChirp();

  /** Updates the normalizationFactor member variable acording to a new blockSize, direction or
  normalizationMode. */
  void updateNormalizationFactor();

  std::complex<T> *h;        /**< The spectrum of the of the h-values. */
  std::complex<T> *c;        /**< The modulating chirp-signal. */
  std::complex<T> *y;        /**< Internal buffer of size M. */

  int N;                  /**< The blocksize of the Bluestein-FFT. */
  int M;                  /**< The enlarged blocksize for the embedded radix-2 FFT. M ist the
                          smallest power of two, such that M >= 2*N-1. */

  int direction;          /**< The direction of the transform (@see: directions). */
  int normalizationMode;  /**< The normalization mode (@see: normalizationModes. */

  T normalizationFactor;

  bool blockSizeIsPowerOfTwo;
  /**< Indicates if we have the special case of a power of two blocksize - in this case, the
  whole transformation is significantly simplified. */

  rsFourierTransformerRadix2<T> transformerRadix2;
  /**< This embedded object is used to perform the radix-2 forward and inverse transformations
  which occur as part of the Bluestein algorithm. */

};

#endif
