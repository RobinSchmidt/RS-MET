#ifndef RAPT_TRANSFORMS_H
#define RAPT_TRANSFORMS_H

/** A collection of linear transforms for signal vectors where by "linear transform" we mean here 
in general something like a matrix-vector product of the form:
  
  y = A * x

where x and y are input- and output vectors respectively and A is the matrix that defines the 
transformation. Many of the important linear transforms used in signal processing do not require 
the explicit computation of a matrix-vector product, which is an O(N^2) process, but can use more 
efficient algorithms due to the special structure of the matrix. The most notable example of 
these is certainly the fast Fourier transform (FFT) which can be computed by an O(N*log(N)) 
algorithm. But there are many others, each with its own unique features and application domain. 

Note that the Fourier transforms here are mostly convenience functions with no attempt to be 
particularly fast or precise. They may be good enough for prototyping, but production code that 
needs a lot of FFTs of the same size should use rsFourierTransformerRadix2 because that class 
precomputes the twiddle factors for the given size, avoiding some computations in the actual FFT 
and also avoiding error accumulation in the twiddle factors due to the recursion in the algos 
used here. */

class rsLinearTransforms
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Fourier Transforms

  /** A radix-2 FFT using a decimation in frequency (DIF) algorithm. The third parameter WN is a
  primitive N-th root of unity that the caller must pass. When T is a complex type, WN should be 
  e^(-2*i*pi/N) for a forward FFT or e^(+2*i*pi/N) for an inverse FFT. I have opted to require 
  client code to pass this primitive N-th root of unity to make it possible to instantiate the 
  routine also unchanged for the case when T is a modular integer type in which case the FFT is 
  also known as NTT for "number theoretic transform". See:
    https://ccrma.stanford.edu/~jos/st/Number_Theoretic_Transform.html
  The algorithm was adapted from algorithm 4.2 in the book "Inside the FFT black box" and then
  simplified. All twiddle factors are computed on the fly via recursion. For a complex FFT, the 
  caller just needs one single complex exponential (i.e. a sin/cos pair) evaluation and then call 
  this routine with that basic twiddle factor. The recursion is perhaps not very advisable from a 
  numeric point of view due to floating point error accumulation but for number theoretic 
  transform (i.e. FFT using the roots of unity of modular arithmetic), that's no issue, so that 
  seems desirable when the algo is instantiated for that. For quick and dirty prototype code, the 
  numeric precision is probably still good enough and production code should use optimized routines
  with (precisely) precomputed twiddle factors anyway, because that's more accurate and presumably
  also more efficient (maybe, but it will need to access two arrays instead of just one). This is 
  the way, it's done in class rsFourierTransformerRadix2, for example (which uses Ouura FFT). The 
  advantage of this implementation is that it requires no additional storage space for the twiddle
  factors. */
  template<class T>
  static void fourierRadix2DIF(T *x, int N, T WN);

  /** Invokes fourierRadix2DIF with WN = e^(-2*i*pi/N). This result in the regular FFT for complex
  number types. */
  template<class T>
  static void fourierRadix2DIF(std::complex<T> *x, int N)
  { fourierRadix2DIF(x, N, exp(std::complex<T>(T(0), T(-2.0*PI/N)))); }

  /** Inverse transform of fourierRadix2DIF() for complex types. Invokes fourierRadix2DIF with 
  WN = e^(+2*i*pi/N) and then scales the output by 1/N. */
  template<class T>
  static void fourierInvRadix2DIF(std::complex<T> *x, int N)
  { fourierInvUnscaledRadix2DIF(x, N); rsArrayTools::scale(x, N, T(1)/T(N)); }

  /** Invokes fourierRadix2DIF with WN = e^(+2*i*pi/N). This result in the inverse FFT for complex
  number types, except for the scaling by 1/N. It may be occasionally useful to bypass the scaling,
  for example, because you have already done the scaling after a forward FFT (which is actually 
  a very typical thing do in FFT analysis, because there we want the actual amplitudes of the 
  sines which do not depend on N, but DSP textbooks typically have the convention the other way
  around). That's why the unscaled transform has been factored out.  */
  template<class T>
  static void fourierInvUnscaledRadix2DIF(std::complex<T> *x, int N)
  { fourierRadix2DIF(x, N, exp(std::complex<T>(T(0), T(+2.0*PI/N)))); }


  //-----------------------------------------------------------------------------------------------
  // \name Fourier-alike Transforms (using non-sinusoidal waveshapes?)

  /** Fast Hadamard transform without scaling or sequency based ordering. The algorithm works in 
  place without allocating heap memory and was ported from here:
    https://en.wikipedia.org/wiki/Fast_Walsh%E2%80%93Hadamard_transform
  Calling hadamard(x, N) gives the same result as calling kronecker2x2(x, N, 1, 1, 1, -1), but the
  algorithm here is a bit simpler. In particular, the multiplications could be thrown away. */
  template<class T>
  static void hadamard(T* x, int N);
  // needs test
  // https://en.wikipedia.org/wiki/Hadamard_transform



  //-----------------------------------------------------------------------------------------------
  // \name Kronecker Transforms

  /** Fast Kronecker transform of a length N array x using a 2x2 seed matrix with coeffs
  a,b,c,d. N must be a power of 2. It is basically the transform described here:
    http://www.rs-met.com/documents/dsp/GeneralizedHadamardTransform.pdf
  The result will again end up in x. The algorithm works in place without allocating heap memory
  and is an adaption of the algorithm used in hadamard(). I don't really know, if "Kronecker 
  transform" is a common name for this kind of transform. I realized, that this "Generalized 
  Hadamard Transform" can be further generalized, using an NxN (or even MxN) matrix as seed and
  I found it appropriate to call the resulting family of transforms (fast) Kronecker transforms, 
  because they are all based on the Kronecker product. Their complexity is generally O(N*log(N)).*/
  template<class T>
  static void kronecker2x2(T* x, int N, T a, T b, T c, T d);
  // ToDo: 
  // -Maybe provide a version that takes N at compile time (as template parameter) to make it 
  //  easier for the compiler to unroll the loops. It's used in FDN reverb, so it may be worthwhile
  //  to aggressively optimize for that.
  // -Implement kronkecker3x3, kroneckerNxN, kroneckerMxN (see rosic_EffectsTests.cpp for 
  //  prototypes)

  /** Inverse transform of kronecker2x2(). It is obtained by simply using the very same fast 
  Kronecker transform algorithm but with an inverted seed matrix. */
  template<class T>
  static void kroneckerInv2x2(T* x, int N, T a, T b, T c, T d)
  { T s = T(1) / (a*d - b*c); kronecker2x2(x, N, s*d, -s*b, -s*c, s*a); }

  /** Fast Kronecker transform of a length N array x using the given 3x3 seed matrix. N must be a
  power of 3. @see kronecker2x2()  */
  template<class T>
  static void kronecker3x3(T* v, int N, const rsMatrix3x3<T>& A);

  // ToDo:
  //template<class T>
  //static void kroneckerInv3x3(T* v, int N, const rsMatrix3x3<T>& A)
  //{ kronecker3x3(v, N, A.getInverse()); }

};

//=================================================================================================
// The old, soon to be deprecated functions:

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
// rename to fourierBernsee, move to prototypes

/** A direct implementation of the DFT with the slow O(N^2) complexitiy scaling. Its main utility
is to check the outputs of the FFT algorithms for debugging.
\todo: check this function.
*/
template<class T>
void rsDFT(std::complex<T> *buffer, int N);
// maybe move to prototypes

/** A general purpose FFT-routine for complex inputs. */
template<class T>
void rsFFT(std::complex<T> *buffer, int N);

/** Computes FFT-magnitudes and (optionally) also phases of a given real signal of length N. The
resulting magnitude and phase arrays are both of length N/2 due to the symmetry of the Fourier
transform of a real signal. */
template<class T>
void rsMagnitudeAndPhase(T *signal, int N, T *magnitudes, T *phases = NULL);
// this is a convenience function mainly for prototyping that allocates heap memory
// rename to fourierMagPhs or fourierPolar - maybe move to prototypes


template<class T> // replacement: rsLinearTransforms::fourierInvRadix2DIF
RS_DEPRECATED(void rsIFFT(std::complex<T> *buffer, int N));

template<class T> RS_DEPRECATED_WITH_BODY(
void rsRadix2FFT(std::complex<T>* x, int N),
{ rsLinearTransforms::fourierRadix2DIF(x, N); })  // use that directly instead!

template<class T> RS_DEPRECATED_WITH_BODY(
  void rsFGHT(T* A, int N, T a, T b, T c, T d),
  { rsLinearTransforms::kronecker2x2(A, N, a, b, c, d); })    // use that directly instead!

template<class T> RS_DEPRECATED_WITH_BODY(
  void rsIFGHT(T* A, int N, T a, T b, T c, T d),
  { rsLinearTransforms::kroneckerInv2x2(A, N, a, b, c, d); }) // use that directly instead!

/** An FFT-routine for arbitrary input sizes that uses the Bluestein algorithm. */
//template<class T>
//void rsBluesteinFFT(rsComplex<T> *buffer, int N);


#endif
