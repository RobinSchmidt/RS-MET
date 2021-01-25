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
  RS_INLINE std::complex<T> getValue()
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

  /** Standard constructor. Sets the iterator up for producing y[n] = a * sin(w*n + p) with
  a = 1, w = 1, p = 0. It initializes the recursion coefficient and state variables using
  precomuted constants, so no costly call to setup() is invoked. */
  rsSineIterator()
  {
    a1 = T( 1.0806046117362795);
    s1 = T(-0.84147098480789650);
    s2 = T(-0.90929742682568171);
    // calling setup(1, 0, 1) would compute these values, but that would be more costly.
    // maybe initialize all to zero

    setup(0.05, 0.0, 1.0);
  }

  /** Constructor. You have to pass the normalized radian frequency "w" and may optionally pass
  a start-phase "p" and an amplitude "a" and for the formula y[n] = a * sin(w*n + p). If you
  don't pass phase and/or amplitude, zero and unity will be used as default values respectively.
  After construction (or a call to setup), the first call to getValue() will return
  y[0] = a*sin(p), the next call gives y[1] = a*sin(w+p), then y[2] = a*sin(2*w+p), and so on. */
  rsSineIterator(T w, T p = 0.0, T a = 1.0) { setup(w, p, a); }



  /** Sets up the sine oscillator such that the next call to getValue will return
  y[0] = a*sin(p), then the next call gives y[1] = a*sin(w + p), then
  y[2] = a*sin(2*w + p), and so on. So, the object will behave exactly in the same way as
  after creating it with the constructor using the same arguments */
  void setup(T w = 1.0, T p = 0.0, T a = 1.0);


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  T getOmega() const { return acos(T(0.5)*a1); }  //  a1 = 2.0*cos(w);

  T getPhase() const 
  {
    T p = asin(s1);
    if(s1 < s2)
      p += T(PI);    // needs test
    return p;
  }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns one output value at a time (and updates internal state variables for next call). */
  RS_INLINE T getValue()
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

/*
todo:
rsSineCosineIterator : public rsComplexExponentialIterator
rsExponentialIterator, rsLinearIterator, rsQuadraticIterator, rsCubicIterator, 
rsPolynomialIterator,
for the polynomial iterators, see Salomon - Computer Graphics, page 275ff ("Fast Calculation of the 
Curve" and page 698 ("Forward Differences")
other functions: sinh, cosh, tanh, 1/cosh, 1/cosh^2 (gaussian-like?)

here's an interesting thread about a recursive sine oscillator:
https://dsp.stackexchange.com/questions/124/how-to-implement-a-digital-oscillator
especially the amplitude drift compensation approach with a taylor expansion of
1 / (sqrt(re^2 + im^2)) ~= (1/2) * (3 - (re^2 + im^2))
every 1000 (ot something) samples

*/
#endif
