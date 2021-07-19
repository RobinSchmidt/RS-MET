#ifndef RAPT_FUNCTIONITERATORS_H_INCLUDED
#define RAPT_FUNCTIONITERATORS_H_INCLUDED

/** This file contains a couple of classes that can be used to efficiently evaluate multiple
values of the same function at equidistant abscissa values, like: y[n] = f(x0 + n*dx), starting
at some x-axis value x0, giving y0 and then computing the successive values of y[n] at the other
x[n] by means of some recursive formula that doesn't require the full re-evaluation of the
function at each x-value. Such recursion formulas are available for the complex exponential and
thus also for sine/cosine and the real exponential, etc. */


/** Computes successive values of w[n] = a * z^n, such that w[0] = a, w[1] = a*z, w[2] = a*z^2,
w[3] = a*z^3, etc.. The parameter "a" is the multiplier and also equals the intial value. */

template<class T>
class rsComplexExponentialIterator
{

public:

  /** Constructor. You have to pass the two parameters "a" and "z" for the formula
  w[n] = a * z^n */
  rsComplexExponentialIterator(std::complex<T> a, std::complex<T> z);

  /** Resets the object's next output value to the given initialValue. The initial value equals
  the multiplier "a" in w[n] = a * z^n. */
  void resetValue(std::complex<T> initialValue);

  /** Sets the value of z in the computation formula  w[n] = a * z^n. */
  void setZ(std::complex<T> newZ);

  /** Returns one output value at a time (and updates the internal state variables for
  the next call). After intialization (or reset), the first call will return w[0] = a = a*z^0,
  the next call gives w[1] = a*z^1, then w[2] = a*z^2, and so on. */
  inline std::complex<T> getValue()
  {
    std::complex<T> ww = w;
    w *= z;
    return ww;
  };

  // \todo in "Matters Computational" (section about FFT), Jörg Arndt says that this recursion 
  // leads to exponential buildup of error and proposes another recursion that is numerically 
  // better behaved - try it, compare, maybe make two classes one with the simple but unstable 
  // and one with the more stable iteration ...oh - but this statement is concerned only with
  // sinusoids (i.e. abs(z) = 1)


protected:

  std::complex<T> w, z;

};

//=================================================================================================

/** Computes successive values of y[n] = a * sin(w*n + p). The parameter "a" is an overall
amplitude, "w" can be seen as a normalized radian frequency and "p" is a startphase. */

template<class T>
class rsSineIterator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Standard constructor. Directly after creation, calling geValue will just produce silence, 
  i.e. y[n] = 0 * sin(0*n + 0). If you want to produce an actual sinusoid, you need to call setup 
  after creation or alternatively, the parametrized constructor below. */
  rsSineIterator()
  {
    //a1 = T( 1.0806046117362795);
    //s1 = T(-0.84147098480789650);
    //s2 = T(-0.90929742682568171);
    // calling setup(1, 0, 1) would compute these values, but that would be more costly.
    // ...maybe do a direct in-class initialization - measure if that makes a difference in 
    // creation cost
    // maybe initialize all to zero

    //setup(0.05, 0.0, 1.0);  // ??! for debug?
  }

  /** Constructor. You have to pass the normalized radian frequency "w" and may optionally pass
  a start-phase "p" and an amplitude "a" and for the formula y[n] = a * sin(w*n + p). If you
  don't pass phase and/or amplitude, zero and unity will be used as default values respectively.
  After construction (or a call to setup), the first call to getValue() will return
  y[0] = a*sin(p), the next call gives y[1] = a*sin(w+p), then y[2] = a*sin(2*w+p), and so on. */
  rsSineIterator(T w, T p = 0.0, T a = 1.0) { setup(w, p, a); }

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up the sine oscillator such that the next call to getValue will return
  y[0] = a*sin(p), then the next call gives y[1] = a*sin(p + w), then y[2] = a*sin(p + 2*w), 
  then y[3] = a*sin(p + 3*w), and so on. So, the object will behave exactly in the same way as
  after creating it with the constructor using the same arguments */
  void setup(T w = 1.0, T p = 0.0, T a = 1.0);

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the current normalized radian frequency. Note that the value is computed/retrieved 
  from the recursion coefficient, so it may be subject to roundoff errors with respect to the value
  that was passed to setup (we don't store the user parameters here to keep the object size 
  small). */
  T getOmega() const { return acos(T(0.5)*a1); }  //  a1 = 2.0*cos(w);

  /*
  T getPhase() const 
  {
    T p = asin(s1);
    if(s1 < s2)
      p += T(PI);    // needs test
    return p;
  }
  */
  // ...nope! this does not yet work!
  // -is this the phase of the last or next call to getValue? maybe the function should be named 
  //  accordingly: getNextPhase, getPreviousPhase or getLastPhase - getNextPhase is 
  //  getLastPhase + getOmega
  // -what is the supposed output range? [-pi, +pi)? or [0, 2pi)? document this!

  // what about getAmplitude?

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns one output value at a time (and updates internal state variables for next call). */
  inline T getValue()
  {
    T tmp = a1*s1 - s2;
    s2 = s1;
    s1 = tmp;
    return tmp;
  }
  // todo: maybe make a version that takes an input sample - it should just be added to tmp:
  // T tmp = in + a1*s1 - s2;



protected:

  T a1 = T(0);             // recursion coefficient
  T s1 = T(0), s2 = T(0);  // past sine outputs

};

//=================================================================================================

/** under construction...works already but setup is still in an inefficient prototype stage... 

A class for iteratively evaluating a polynomial p(x). To use it, you need to instantiate it for the 
desired type T for coefficients and output and degree N. To initialize the iterative evaluation, 
you call setup with the polynomial coefficents, stepsize h and initial value x0. After that, the 
1st call to getValue will return p(x0), the 2nd p(x0+h), the 3rd p(x0+2*h) and so on. */

template<class T, int N>     // T: type for coeffs and output, N: degree of the polynomial
class rsPolynomialIterator
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setup(const T* newCoeffs, T newStepSize, T initialValue);


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns one output value at a time (and updates internal state variables for next call). */
  inline T getValue()
  {
    T r = y[N];                  // result
    for(int i = N; i > 0; i--) 
      y[i] += y[i-1];            // state update
    return r;
  }


protected:

  T y[N+1];

};

//=================================================================================================

/** Similar to rsPolynomialIterator, but instead of iteratively evaluating a polynomial p(x) itself, 
it evaluates the exponential function of that polynomial, i.e. exp(p(x)).

Warning: i think, this may be numerically unstable, especially for larger N...more tests needed.  */

template<class T, int N>
class rsExpPolyIterator : public rsPolynomialIterator<T, N>
{

public:

  void setup(const T* newCoeffs, T newStepSize, T initialValue)
  {
    rsPolynomialIterator<T, N>::setup(newCoeffs, newStepSize, initialValue);
    for(int i = 0; i <= N; i++) 
      this->y[i] = rsExp(this->y[i]);  
  }
  // todo: maybe make a version that lets the user specify a basis b..i think, we just need to 
  // multiply all y[i] by log(b) (inside the exp call)

  inline T getValue()
  {
    T r = y[N];                    // result
    for(int i = N; i > 0; i--) 
      this->y[i] *= this->y[i-1];  // state update
    return r;
  }

};

//=================================================================================================

/** Iteratively computes a sinusoidal frequency sweep using a cubic rsExpPolyIterator. That means,
the instantaneous phase and the logarithm of the instantaneous amplitude are given by cubic 
polynomials. The polynomials are set up in terms of user parameters that are compatible with the 
sinusoidal modeling framework. The intention is to use this class to synthesize a sinusoidal
partial efficiently in a realtime context. Computation of one output sample requires 3 complex 
multiplications (i.e. 12 real muls, 3 adds, 3 subs) and the class lends itself well to 
vectorization (i.e. usage with T = rsSimdVector<float, 16> or something). The driver class has to 
take care of not letting the accumulated roundoff error grow out of control. In a typical 
situation, one should probably re-initialize the state with a proper direct calculation every 
couple of hundreds or thousands of samples ...tbc...  */

template<class T>
class rsSineSweepIterator
{

public:

  /** Structure for the user parameters that determine the cubic polynomials for instantaneous 
  phase and logarithm of instantaneous amplitude. There are always two values indexed by 0 or 1, 
  standing for the start and the end of the synthesized sweep. */
  struct Parameters
  {
    T t0, t1; /**< time stamps in samples */
    T p0, p1; /**< unwrapped(!) phases in radians */
    T w0, w1; /**< omega = 2*pi*frequency/sampleRate, derivative of the phase */
    T l0, l1; /**< log(amplitude) */
    T r0, r1; /**< "raise", derivative of log of amplitude with respect to t in samples */
  };

  /** Sets up the initial state according to the user parameters. */
  void setup(const Parameters& params);

  // todo:
  //inline T getPhase()     const { return std::arg(core.y[N]); }
  //inline T getAmplitude() const { return std::abs(core.y[N]); }

  /** Computes a complex output value within which the imaginary part is the actual sine wave as 
  required by the sinusoidal modeling framework and the real part is a corresponding cosine 
  quadrature component that you get for free due to the way the algorithm works. And it updates 
  the state for the next call. */
  inline std::complex<T> getComplexValue() { return core.getValue(); }

  /** Convenience function to compute the sine only (and update the state). */
  inline T getSine()   { return getComplexValue().imag(); }

  /** Convenience function to compute the cosine only (and update the state). Note that you should 
  not call getSine and after that getCosine if you need both values because that would trigger two 
  state updates where you probably intend just one. If you need both values, use either 
  getComplexValue or getSineAndCosine. */
  inline T getCosine() { return getComplexValue().real(); }

  /** Computes both, the original sine and its quadrature component. This has no extra cost 
  compared to calling either getSine or getCosine to produce either of the outputs because
  both must be computed internally anyway. It's mostly meant as convenience function though because 
  it may be more efficient, if the driver just calls getComplexValue instead (that may save two 
  assignments). */
  inline void getSineAndCosine(T* sine, T* cosine)
  {
    std::complex<T> w = getComplexValue();
    *cosine = w.real();
    *sine   = w.imag();
  }


protected:

  rsExpPolyIterator<std::complex<T>, 3> core;
  // todo: implement and use a special optimized rsExpCubicIterator

};



#endif
