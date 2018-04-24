#ifndef RS_MISCAUDIO_H
#define RS_MISCAUDIO_H

namespace RSLib
{

  /**

  This file defines functions for processing audio in non-realtime.

  */

  //-----------------------------------------------------------------------------------------------
  // signal generators:

  /** The classical standard waveforms that are found in synthesizers. \todo move inot some file 
  containing enumerations and definitions for use throughout the rsAudio module */
  enum standardWaveforms
  {
    SINE = 0,
    SAW,
    SQUARE,
    TRIANGLE,

    NUM_STANDARD_WAVEFORMS
  };

  /** Synthesizes a standard waveform at the desired frequency and samplerate. The phase is 
  expected in radians and there is an option to do anti-aliased synthesis (by means of adding
  the sinusoidal components up to the Nyquist frequency). @see standardWaveshapes */
  RSLib_API void synthesizeWaveform(double *x, int N, int shape, double frequency, 
    double sampleRate, double phase = 0.0, bool antiAlias = false);

  RSLib_API void synthesizePulseWave(double *x, int N, double frequency, double dutyCycle, 
    double sampleRate, double phase = 0.0, bool antiAlias = false);

  RSLib_API void synthesizeDecayingSine(double *x, int N, double frequency, double amplitude, 
    double decayTime, double startPhase, double sampleRate);

  RSLib_API void synthesizeModal(double *x, int N, rsVectorDbl frequencies, rsVectorDbl amplitudes, 
    rsVectorDbl decayTimes, rsVectorDbl startPhases, double sampleRate);

  RSLib_API void synthesizeModalPluckedString(double *x, int N, double frequency, 
    double sampleRate, double decayTime, double decayExponent, double amplitudeExponent, 
    double inharmonicity, double phase, double evenAmplitudeScaler);

  RSLib_API void synthesizeModalRodFreeFree(double *x, int N, double frequency, double sampleRate, 
    double decayTime, double decayExponent, double amplitudeExponent, double phase);

  //-----------------------------------------------------------------------------------------------
  // interpolation:

  // comment these, move into class rsResampler:

  RSLib_API void upsampleLinear(double *in, int inLength, double *out, int upsamplingFactor);

  RSLib_API void upsampleHermiteAsymmetric1(double *in, int inLength, double *out, 
    int upsamplingFactor, double shape);

  RSLib_API void upsampleHermiteAsymmetricM(double *in, int inLength, double *out, 
    int upsamplingFactor, int M, double shape);

  //-----------------------------------------------------------------------------------------------
  // IIR filters:

  /** Applies a Butterworth filter to the given signal 'x' of length 'N' and stores the result in
  'y' (which may be equal to x). 'mode' ...*/
  RSLib_API void filterButterworth(double *x, double *y, int N, double frequency, 
    double sampleRate, int mode, int prototypeOrder, double gain = 0.0, 
    bool forwardBackward = false);



  //-----------------------------------------------------------------------------------------------
  // others:

  /** Applies an attack/release envelope follower to the signal 'x' of length 'N' and stores the 
  result in 'y' (which may be equal to x). */
  RSLib_API void estimateEnvelope(double *x, double *y, int N, double sampleRate, 
    double attackTime, double releaseTime, int mode, bool forwardBackward = false);


  //-----------------------------------------------------------------------------------------------
  // spectral processing:

  //void fft(Complex *signalBlock, int blockSize, Complex *spectrum, int fftSize);

  /** Transforms the 'signalBlock' into the frequency domain via the FFT algorithm. The 'fftSize' 
  should be greater than or equal to the 'blockSize' where in the former case zero padding will be 
  used. When the 'fftSize' is a power of two, the efficient radix-2 FFT will be used, otherwise the
  Bluestein algorithm will be used. */
  RSLib_API void fft(double *signalBlock, int blockSize, rsComplexDbl *spectrum, int fftSize);

  /** Inverse Foruier transforms the 'spectrum' to yield a complex signal block of the same 
  length. ...to be tested */
  RSLib_API void ifft(rsComplexDbl *spectrum, int fftSize, rsComplexDbl *signalBlock);

  /** Inverse Foruier transforms the 'spectrum' (which is assumed to be conjugate symmetric) to 
  yield a real signal block (of the same length). ...to be tested */
  RSLib_API void ifftReal(rsComplexDbl *spectrum, int fftSize, double *signalBlock);

  /** Transforms the 'signalBlock' into the frequency domain via a call to fft and then extracts the 
  magnitudes and phases. The arrays 'magnitudes' and 'phases' should be of length fftSize/2+1 for 
  even 'fftSize' and (fftSize-1)/2+1 odd 'fftSize'. You may pass a NULL pointer for the phases, if 
  you are not interested in them. */
  RSLib_API void fftMagnitudesAndPhases(double *signalBlock, int blockSize, double *magnitudes, 
    double *phases, int fftSize);

  /** Computes the real cepstrum of a real-valued signal. The real cepstrum is defined as the 
  inverse dicrete Fourier transform of the logarithm of the absolute value of the discrete 
  Fourier transfrom of the signal: output = IDFT(log(abs(DFT(input)))). Note that the real 
  cepstrum of a finite length sequence is generally an infinite-length sequence. This function 
  returns a time-aliased version of this infinite-length sequence. You can reduce this 
  time-aliasing by zero-padding the input sequence. */
  RSLib_API void signalToRealCepstrum(double *signal, int numSamples, double *cepstrum);

  /** Transforms a real cepstrum back into the signal domain. @see signalToRealCepstrum  */
  RSLib_API void realCepstrumToSignal(double *cepstrum, int numSamples, double *signal);

  /** Computes a minimum-phase reconstruction of some input sequence. The resulting sequence will 
  have the same spectral magnitudes as the original signal, but phases will be adjusted such that a
  greater portion of the signal energy will be concentrated near the beginning of the sequence. */
  RSLib_API void minimumPhaseReconstruction(double *input, int numSamples, double *output);


  /** Computes cross-correlation between the signals x and y and stores the result in 'result' 
  which should be of length xLength+yLength-1. If x and y are are equal (point to the same buffer,
  the result will be the autocorrelation of x. 
  
  \todo: TEST THIS!!!
  
  */
  RSLib_API void crossCorrelation(double *x, int xLength, double *y, int yLength, double *result);

  /** Estimates the fundamental frequency of some signal using the autocorrelation function. */
  //RSLib_API double estimateFundamental(double *x, int N, double sampleRate, 
  //  double minExpected = 20.0, double maxExpected = 10000.0);

  RSLib_API void estimateModalParameters(double *x, int N, rsVectorDbl *frequencies, 
    rsVectorDbl *amplitudes, rsVectorDbl *decayTimes, double sampleRate);


  /** Given a correlation sequence r of length N, this function returns the subsample-precision 
  position of the maximum value of that sequence. The subsample position is obtained by fitting
  a quadratic parabola to the maximum value and its left and right neighbours and solving for the
  maximum of the parabola. */
  RSLib_API double rsMaxCorrelationLag(double *r, int N);

  /** Given two signals x1, x2 of length N, this function applies a Hamming window and computes the
  cross-correlation of the windowed signals and then locates the position of the maximum of that 
  sequence with subsample precision. The deBias parameter determines whether a biased or unbiased 
  sample cross-correlation should be used. I'm not yet sure, which is better... */
  RSLib_API double rsMaxCorrelationLag(double *x1, double *x2, int N, bool deBias = false);

  /** Given two signals x1, x2 of length N, this function computes the amount (in samples) by which
  x2 has to be shifted to obtain a best match to x1. If you have signals of different length, say 
  N1 and N2, you may still use this function by passing rsMin(N1, N2) for the N parameter. The 
  shift amount can be positive or negative, corresponding to righward and leftward shifts, 
  respectively. For the deBias parameter, @see rsMaxCorrelationLag() */
  RSLib_API double rsGetShiftForBestMatch(double *x1, double *x2, int N, bool deBias = false);

}

#endif
