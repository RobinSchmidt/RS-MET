#ifndef rosic_NonRealtimeProcesses_h
#define rosic_NonRealtimeProcesses_h

//// standard library includes
//#include <vector>
//
//// rosic-indcludes:
//#include "../basics/rosic_FunctionTemplates.h"
//#include "../math/rosic_ElementaryFunctionsReal.h"
//#include "../filters/rosic_EngineersFilter.h"
//#include "../generators/rosic_ModalSynthesizer.h"
//#include "../transforms/rosic_FourierTransformerBluestein.h"

namespace rosic
{

  /**

  This file defines functions for processing audio in non-realtime.

  */

  //-----------------------------------------------------------------------------------------------
  // signal generators:

  /** Synthesizes a standard waveform at the desired frequency and samplerate. The phase is 
  expected in radians and there is an option to do anti-aliased synthesis (by means of adding
  the sinusoidal components up to the Nyquist frequency). @see standardWaveshapes */
  void synthesizeWaveform(double *x, int N, int shape, double frequency, double sampleRate, 
    double phase = 0.0, bool antiAlias = false);

  void synthesizePulseWave(double *x, int N, double frequency, double dutyCycle, double sampleRate, 
    double phase = 0.0, bool antiAlias = false);

  void synthesizeDecayingSine(double *x, int N, double frequency, double amplitude, 
    double decayTime, double startPhase, double sampleRate);

  void synthesizeModal(double *x, int N, Vector frequencies, Vector amplitudes, Vector decayTimes, 
    Vector startPhases, double sampleRate);

  void synthesizeModalPluckedString(double *x, int N, double frequency, double sampleRate, 
    double decayTime, double decayExponent, double amplitudeExponent, double inharmonicity, 
    double phase, double evenAmplitudeScaler);

  void synthesizeModalRodFreeFree(double *x, int N, double frequency, double sampleRate, 
    double decayTime, double decayExponent, double amplitudeExponent, double phase);

  //-----------------------------------------------------------------------------------------------
  // interpolation:

  // comment these:

  void upsampleLinear(double *in, int inLength, double *out, int upsamplingFactor);

  void upsampleHermiteAsymmetric1(double *in, int inLength, double *out, int upsamplingFactor, 
    double shape);

  void upsampleHermiteAsymmetricM(double *in, int inLength, double *out, int upsamplingFactor, 
    int M, double shape);

  //-----------------------------------------------------------------------------------------------
  // IIR filters:

  /** Applies a Butterworth filter to the given signal 'x' of length 'N' and stores the result in
  'y' (which may be equal to x). 'mode' ...*/
  void filterButterworth(double *x, double *y, int N, double frequency, double sampleRate, 
    int mode, int prototypeOrder, double gain = 0.0, bool forwardBackward = false);



  //-----------------------------------------------------------------------------------------------
  // others:

  /** Applies an attack/release envelope follower to the signal 'x' of length 'N' and stores the 
  result in 'y' (which may be equal to x). */
  void estimateEnvelope(double *x, double *y, int N, double sampleRate, double attackTime, 
    double releaseTime, int mode, bool forwardBackward = false);


  //-----------------------------------------------------------------------------------------------
  // spectral processing:

  //void fft(Complex *signalBlock, int blockSize, Complex *spectrum, int fftSize);

  /** Transforms the 'signalBlock' into the frequency domain via the FFT algorithm. The 'fftSize' 
  should be greater than or equal to the 'blockSize' where in the former case zero padding will be 
  used. When the 'fftSize' is a power of two, the efficient radix-2 FFT will be used, otherwise the
  Bluestein algorithm will be used. */
  void fft(double *signalBlock, int blockSize, Complex *spectrum, int fftSize);

  /** Inverse Foruier transforms the 'spectrum' to yield a complex signal block of the same 
  length. ...to be tested */
  void ifft(Complex *spectrum, int fftSize, Complex *signalBlock);

  /** Inverse Foruier transforms the 'spectrum' (which is assumed to be conjugate symmetric) to 
  yield a real signal block (of the same length). ...to be tested */
  void ifftReal(Complex *spectrum, int fftSize, double *signalBlock);

  /** Transforms the 'signalBlock' into the frequency domain via a call to fft and then extracts the 
  magnitudes and phases. The arrays 'magnitudes' and 'phases' should be of length fftSize/2+1 for 
  even 'fftSize' and (fftSize-1)/2+1 odd 'fftSize'. You may pass a NULL pointer for the phases, if 
  you are not interested in them. */
  void fftMagnitudesAndPhases(double *signalBlock, int blockSize, double *magnitudes, 
    double *phases, int fftSize, bool scale = false);

  /** Computes the real cepstrum of a real-valued signal. The real cepstrum is defined as the 
  inverse dicrete Fourier transform of the logarithm of the absolute value of the discrete 
  Fourier transfrom of the signal: output = IDFT(log(abs(DFT(input)))). Note that the real 
  cepstrum of a finite length sequence is generally an infinite-length sequence. This function 
  returns a time-aliased version of this infinite-length sequence. You can reduce this 
  time-aliasing by zero-padding the input sequence. */
  void signalToRealCepstrum(double *signal, int numSamples, double *cepstrum);

  /** Transforms a real cepstrum back into the signal domain. @see signalToRealCepstrum  */
  void realCepstrumToSignal(double *cepstrum, int numSamples, double *signal);

  /** Computes a minimum-phase reconstruction of some input sequence. The resulting sequence will 
  have the same spectral magnitudes as the original signal, but phases will be adjusted such that a
  greater portion of the signal energy will be concentrated near the beginning of the sequence. */
  void minimumPhaseReconstruction(double *input, int numSamples, double *output);


  /** Computes cross-correlation between the signals x and y and stores the result in 'result' 
  which should be of length xLength+yLength-1. If x and y are are equal (point to the same buffer,
  the result will be the autocorrelation of x. 
  
  \todo: TEST THIS!!! not yet ready for production
  
  */
  void crossCorrelation(double *x, int xLength, double *y, int yLength, double *result);

  /** Estimates the fundamental frequency of some signal using the autocorrelation function. */
  double estimateFundamental(double *x, int N, double sampleRate, double minExpected = 20.0, 
    double maxExpected = 10000.0);



  void estimateModalParameters(double *x, int N, Vector *frequencies, Vector *amplitudes, 
    Vector *decayTimes, double sampleRate);


}

#endif
