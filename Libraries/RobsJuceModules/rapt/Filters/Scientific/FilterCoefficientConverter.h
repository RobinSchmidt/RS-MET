#ifndef RAPT_FILTERCOEFFICIENTCONVERTER_H_INCLUDED
#define RAPT_FILTERCOEFFICIENTCONVERTER_H_INCLUDED

/** This class contains static functions to convert between the coefficients for various filter 
representations and realization structures. */

template<class T>
class rsFilterCoefficientConverter
{

  typedef std::complex<T> Complex; // preliminary

public:

  /** Converts direct form FIR coefficients to FIR lattice filter reflection coefficients. */
  static void directFormToLatticeFir(T* directFormCoeffs, int order, T* reflectionCoeffs);
  // allocates heap memory

  /** Converts FIR lattice filter reflection coefficients to direct form FIR coefficients. */
  static void latticeToDirectFormFir(T* reflectionCoeffs, int order, T* directFormCoeffs);
  // allocates heap memory

  /** Converts complex poles and zeros into coefficients for a biquad cascade in which each stage 
  implements the difference equation:
  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2] \f]
  The arrays of poles and zeros are understood to contain only one pole (zero) for each complex 
  conjugate pair of poles and zeros of the actual filter. If there is a first order stage present 
  (a real pole/zero) than this real pole/zero should be the last entry in the array and the flag 
  lastStageIsFirstOrder should be set to true. ...to be deprecated */
  static void polesAndZerosToBiquadCascade(Complex* poles, int numPoles, Complex* zeros, 
    int numZeros, T* b0, T* b1, T* b2, T* a1, T* a2, bool lastStageIsFirstOrder);

  static void polesAndZerosToBiquadCascade(Complex* poles, Complex* zeros, int order,
    T* b0, T* b1, T* b2, T* a1, T* a2);
  // maybe define an alias pz2sos

  /**  !!! NOT YET FULLY IMPLEMETED !!!
  Converts the coefficients of a cascade biquad of biquad filters each of which having the form:
  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2] \f]
  into coefficients for a direct form filter of the form:
  \f[ y[n] = sum_{k=0}^K b_k x[n-k] - sum_{m=1}^M a_m y[n-m] \f]
  The b0, b1, b2, a1, a2 arrays should contain the coefficients for the individual biquad stages. 
  Then, when B is the number of biquad stages, the b array will have to be of size 
  K=2*B+1 (b_0,..., b_K) and the a array will have to be of size M=2*B (a_1,...,a_M) */
  static void biquadCascadeToDirectForm(int numBiquads, T* b0, T* b1, T* b2,
    T* a1, T* a2, T* b, T* a);
  // allocates heap memory

  /** Calculates the magnitude-response of a digital biquad filter with coefficients b0, b1, b2, 
  a0, a1, a1 at the normalized radian frequency 'omega'.  */
  static T getBiquadMagnitudeAt(T b0, T b1, T b2, T a1, T a2, T omega);
  // \todo remove - function is redundant  with the function in FilterAnalyzer

  /** Normalizes the biquad stages described by the given coefficients in such a way that each 
  stage has unit magnitude at the normalized radian frequency 'omega'. If a gainFactor is passed, 
  it will normalize to this gain factor.  */
  static void normalizeBiquadStages(T* b0, T* b1, T* b2, T* a1, T* a2, T omega, int numStages, 
    T gainFactor = 1.0);


  // todo: biquadToPartialFractions - converts a biquad to it partial fraction expansion

  //static void biquadToStateVariable(T* b0, T* b1, T* b2, T* a1, T* a2, ...);

  //static void biquadToPhasor(T* b0, T* b1, T* b2, T* a1, T* a2, T* rc, T* rs, T* wx, T* wy, T* wi);
  // update equation:
  // |x| = r * |c -s| * |x| + |in|
  // |y|       |s  c|   |y|   |0 |
  // output equation:
  // out = wx*x + wy*y + wi*in
  // where r is a decay factor <= 1; c,s are rotation matrix coeffs (sin/cos of angle w); wx,wy,wi 
  // are weights for the state-variables x,y and input -> 5 independent coeffs (r can be absorbed 
  // into the rotation matrix, turning it into a spiraling matrix)
  // the filter's state is a 2D vector - could be called state-vector filter
};

#endif
