#ifndef RAPT_TRANSFORMS_H
#define RAPT_TRANSFORMS_H

// \todo wrap into class rsTransforms and name the particular transforms only fourier, kronecker,
// hadamard, walshHadamerd, etc. - the DFT can be named fourierNaive or fourierSlow, the arbitrary
// size version fourierBluestein, also have fourierRadix2_DIT, _DIF, fourierRadix4, fourierRadix8,
// fourierRadix3, wavelets: gabor, daubechies, ... -> wrap the old functions into an 
//   RS_DEPRECATED( T oldFunc(T p) { return newFunc(p); } )
// macro. See, how juce does that and imitate that. Maybe the templatizeation introduces another
// complication, so maybe a special macro RS_DEPRECATE_T for function templates is needed (but
// it may use the regular RS_DEPRECATED macro

/** FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)

Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the time domain data in
fftBuffer[0...2*fftFrameSize-1]. The FFT array takes and returns the cosine and sine parts in an
interleaved manner, ie. fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
must be a power of 2. It expects a complex input signal (see footnote 2), ie. when working with
'common' audio signals our input signal has to be passed as {in[0],0.,in[1],0.,in[2],0.,...} asf.
In that case, the transform of the frequencies of interest is in fftBuffer[0...fftFrameSize]. */
template<class T>
void smbFft(T *fftBuffer, long fftFrameSize, long sign);
// fourierBernsee

/** A direct implementation of the DFT with the slow O(N^2) complexitiy scaling. Its main utility
is to check the outputs of the FFT algorithms for debugging.
\todo: check this function.
*/
template<class T>
void rsDFT(std::complex<T> *buffer, int N);

/** A general purpose FFT-routine for complex inputs. */
template<class T>
void rsFFT(std::complex<T> *buffer, int N);

/** Computes FFT-magnitudes and (optionally) also phases of a given real signal of length N. The
resulting magnitude and phase arrays are both of length N/2 due to the symmetry of the Fourier
transform of a real signal. */
template<class T>
void rsMagnitudeAndPhase(T *signal, int N, T *magnitudes, T *phases = NULL);
// rename to fourierMagPhs or fourierPolar

/** Inverse transformation of rsFFT. */
template<class T>
void rsIFFT(std::complex<T> *buffer, int N);

/** A radix-2 FFT-routine. */
template<class T>
void rsRadix2FFT(std::complex<T> *buffer, int N);

/** An FFT-routine for arbitrary input sizes that uses the Bluestein algorithm. */
//template<class T>
//void rsBluesteinFFT(rsComplex<T> *buffer, int N);

/** Fast generalized Hadamard transform of length N array A with seed matrix coeffs a,b,c,d. N must
be a power of 2.  */
template<class T>
void rsFGHT(T* A, int N, T a, T b, T c, T d);
// ToDo: maybe provide a version that takes N at compile time (as template parameter) to make it 
// easier for the compiler to unroll the loops.
// maybe call this rsFastKroneckerTrafo2x2

/** Inverse of rsFGHT. */
template<class T>
void rsIFGHT(T* A, int N, T a, T b, T c, T d) 
{ T s = T(1) / (a*d - b*c); rsFGHT(A, N, s*d, -s*b, -s*c, s*a); }



#endif
