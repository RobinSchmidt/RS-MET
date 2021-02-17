#ifndef RAPT_MISCAUDIO_H
#define RAPT_MISCAUDIO_H

/** This file defines functions for processing audio in non-realtime. 

\todo maybe wrap in class

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
template<class T>
void synthesizeWaveform(T *x, int N, int shape, T frequency, T sampleRate, T phase = 0.0, 
  bool antiAlias = false);

template<class T>
void synthesizePulseWave(T *x, int N, T frequency, T dutyCycle, T sampleRate, T phase = 0.0, 
  bool antiAlias = false);

// what about TriSaw

template<class T>
void synthesizeDecayingSine(T *x, int N, T frequency, T amplitude, T decayTime, T startPhase, 
  T sampleRate);

template<class T>
void synthesizeModal(T *x, int N, std::vector<T> frequencies, std::vector<T> amplitudes,
  std::vector<T> decayTimes, std::vector<T> startPhases, T sampleRate);
// todo: pass vectors as references

template<class T>
void synthesizeModalPluckedString(T *x, int N, T frequency, T sampleRate, T decayTime, 
  T decayExponent, T amplitudeExponent, T inharmonicity, T phase, T evenAmplitudeScaler);

template<class T>
void synthesizeModalRodFreeFree(T *x, int N, T frequency, T sampleRate, T decayTime, 
  T decayExponent, T amplitudeExponent, T phase);

//-----------------------------------------------------------------------------------------------
// interpolation:

// comment these, move into class rsResampler:

template<class T>
void upsampleLinear(T *in, int inLength, T *out, int upsamplingFactor);

template<class T>
void upsampleHermiteAsymmetric1(T *in, int inLength, T *out, int upsamplingFactor, T shape);

template<class T>
void upsampleHermiteAsymmetricM(T *in, int inLength, T *out, int upsamplingFactor, int M, 
  T shape);

//-----------------------------------------------------------------------------------------------
// IIR filters:

/** Applies a Butterworth filter to the given signal 'x' of length 'N' and stores the result in
'y' (which may be equal to x). 'mode' ...*/
template<class T>
void filterButterworth(T *x, T *y, int N, T frequency, T sampleRate, int mode, int prototypeOrder, 
  T gain = 0.0, bool forwardBackward = false);



//-----------------------------------------------------------------------------------------------
// others:

/** Applies an attack/release envelope follower to the signal 'x' of length 'N' and stores the
result in 'y' (which may be equal to x). */
template<class T>
void estimateEnvelope(T *x, T *y, int N, T sampleRate, T attackTime, T releaseTime, int mode, 
  bool forwardBackward = false);


//-----------------------------------------------------------------------------------------------
// spectral processing:

//void fft(Complex *signalBlock, int blockSize, Complex *spectrum, int fftSize);

/** Transforms the 'signalBlock' into the frequency domain via the FFT algorithm. The 'fftSize'
should be greater than or equal to the 'blockSize' where in the former case zero padding will be
used. When the 'fftSize' is a power of two, the efficient radix-2 FFT will be used, otherwise the
Bluestein algorithm will be used. */
template<class T>
void fft(T *signalBlock, int blockSize, std::complex<T> *spectrum, int fftSize);

/** Inverse Foruier transforms the 'spectrum' to yield a complex signal block of the same
length. ...to be tested */
template<class T>
void ifft(std::complex<T> *spectrum, int fftSize, std::complex<T> *signalBlock);

/** Inverse Foruier transforms the 'spectrum' (which is assumed to be conjugate symmetric) to
yield a real signal block (of the same length). ...to be tested */
template<class T>
void ifftReal(std::complex<T> *spectrum, int fftSize, double *signalBlock);

/** Transforms the 'signalBlock' into the frequency domain via a call to fft and then extracts the
magnitudes and phases. The arrays 'magnitudes' and 'phases' should be of length fftSize/2+1 for
even 'fftSize' and (fftSize-1)/2+1 odd 'fftSize'. You may pass a NULL pointer for the phases, if
you are not interested in them. */
template<class T>
void fftMagnitudesAndPhases(T *signalBlock, int blockSize, T *magnitudes, T *phases, int fftSize);

/** Computes the real cepstrum of a real-valued signal. The real cepstrum is defined as the
inverse dicrete Fourier transform of the logarithm of the absolute value of the discrete
Fourier transfrom of the signal: output = IDFT(log(abs(DFT(input)))). Note that the real
cepstrum of a finite length sequence is generally an infinite-length sequence. This function
returns a time-aliased version of this infinite-length sequence. You can reduce this
time-aliasing by zero-padding the input sequence. */
template<class T>
void signalToRealCepstrum(T *signal, int numSamples, T *cepstrum);

/** Transforms a real cepstrum back into the signal domain. @see signalToRealCepstrum  */
template<class T>
void realCepstrumToSignal(T *cepstrum, int numSamples, T *signal);

/** Computes a minimum-phase reconstruction of some input sequence. The resulting sequence will
have the same spectral magnitudes as the original signal, but phases will be adjusted such that a
greater portion of the signal energy will be concentrated near the beginning of the sequence. */
template<class T>
void minimumPhaseReconstruction(T *input, int numSamples, T *output);


/** Computes cross-correlation between the signals x and y and stores the result in 'result'
which should be of length xLength+yLength-1. If x and y are are equal (point to the same buffer,
the result will be the autocorrelation of x.

\todo: TEST THIS!!!*/
template<class T>
void crossCorrelation(T *x, int xLength, T *y, int yLength, T *result);
// get rid - is redundant

/** Estimates the fundamental frequency of some signal using the autocorrelation function. */
//RSLib_API double estimateFundamental(double *x, int N, double sampleRate, 
//  double minExpected = 20.0, double maxExpected = 10000.0);

template<class T>
void estimateModalParameters(T *x, int N, std::vector<T> *frequencies,
  std::vector<T> *amplitudes, std::vector<T> *decayTimes, T sampleRate);


/** Given a correlation sequence r of length N, this function returns the subsample-precision
position of the maximum value of that sequence. The subsample position is obtained by fitting
a quadratic parabola to the maximum value and its left and right neighbours and solving for the
maximum of the parabola. */
template<class T>
T rsMaxCorrelationLag(T *r, int N);

/** Given two signals x1, x2 of length N, this function applies a Hamming window and computes the
cross-correlation of the windowed signals and then locates the position of the maximum of that
sequence with subsample precision. The deBias parameter determines whether a biased or unbiased
sample cross-correlation should be used. I'm not yet sure, which is better... */
template<class T>
T rsMaxCorrelationLag(T *x1, T *x2, int N, bool deBias = false);

/** Given two signals x1, x2 of length N, this function computes the amount (in samples) by which
x2 has to be shifted to obtain a best match to x1. If you have signals of different length, say
N1 and N2, you may still use this function by passing rsMin(N1, N2) for the N parameter. The
returned shift amount can be positive or negative, corresponding to rightward and leftward shifts,
respectively. For the deBias parameter, @see rsMaxCorrelationLag() */
template<class T>
T rsGetShiftForBestMatch(T *x1, T *x2, int N, bool deBias = false);




#endif
