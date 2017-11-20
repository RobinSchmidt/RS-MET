#ifndef RS_FUNCTIONITERATORS_H
#define RS_FUNCTIONITERATORS_H

namespace RSLib
{

  /** This file contains a couple of classes that can be used to efficiently evaluate multiple
  values of the same function at equidistant abscissa values, like: y[n] = f(x0 + n*dx), starting
  at some x-axis value x0, giving y0 and then computing the successive values of y[n] at the other
  x[n] by means of some recursive formula that doesn't require the full re-evaluation of the 
  function at each x-value. Such recrusion formulas are available for the complex exponential and 
  thus also for sine/cosine and the real exponential, etc. */


  /** Computes successive values of w[n] = a * z^n, such that w[0] = a, w[1] = a*z, w[2] = a*z^2,
  w[3] = a*z^3, etc.. The parameter "a" is the multiplier and also equals the intial value. */
  class RSLib_API rsComplexExponentialIterator
  {

  public:

    /** Constructor. You have to pass the two parameters "a" and "z" for the formula
    w[n] = a * z^n */
    rsComplexExponentialIterator(rsComplexDbl a, rsComplexDbl z);

    /** Resets the object's next output value to the given initialValue. The initial value equals
    the multiplier "a" in w[n] = a * z^n. */
    void resetValue(rsComplexDbl initialValue);

    /** Sets the value of z in the computation formula  w[n] = a * z^n. */
    void setZ(rsComplexDbl newZ);

    /** Returns one output value at a time (and updates the internal state variables for
    the next call). After intialization (or reset), the first call will return w[0] = a = a*z^0,
    the next call gives w[1] = a*z^1, then w[2] = a*z^2, and so on. */
    RS_INLINE rsComplexDbl getValue()
    {
      rsComplexDbl ww = w;
      w *= z;
      return ww;
    };

    // \todo in "Matters Computational" (section about FFT), Jörg Arndt says that this recursion 
    // leads to exponential buildup of error and proposes another recursion that is numerically 
    // better behaved - try it, compare, maybe make two classes one with the simple but unstable 
    // and one with the more stable iteration ...oh - but this statement is concerned only with
    // sinusoids (i.e. abs(z) = 1)


  protected:

    rsComplexDbl w, z;

  };


  /** Computes successive values of y[n] = a * sin(w*n + p). The parameter "a" is an overall
  amplitude, "w" can be seen as a normalized radian frequency and "p" is a startphase. */
  class RSLib_API rsSineIterator
  {

  public:

    /** Standard constructor. Sets the iterator up for producing y[n] = a * sin(w*n + p) with
    a = 1, w = 1, p = 0. It initializes the recursion coefficient and state variables using 
    precomuted constants, so no costly call to setup() is invoked. */
    rsSineIterator();

    /** Constructor. You have to pass the normalized radian frequency "w" and may optionally pass
    a start-phase "p" and an amplitude "a" and for the formula y[n] = a * sin(w*n + p). If you 
    don't pass phase and/or amplitude, zero and unity will be used as default values respectively.
    After construction (or a call to setup), the first call to getValue() will return 
    y[0] = a*sin(p), the next call gives y[1] = a*sin(w+p), then y[2] = a*sin(2*w+p), and so on. */
    rsSineIterator(double w, double p = 0.0, double a = 1.0);

    /** Sets up the sine oscillator such that the next call to getValue will return
    y[0] = a*sin(p), then the next call gives y[1] = a*sin(w + p), then 
    y[2] = a*sin(2*w + p), and so on. So, the object will behave exactly in the same way as
    after creating it with the constructor using the same arguments */
    void setup(double w = 1.0, double p = 0.0, double a = 1.0);

    /** Returns one output value at a time (and updates internal state variables for next call). */
    RS_INLINE double getValue()
    {
      double tmp = a1*s1 - s2;
      s2         = s1;
      s1         = tmp;
      return tmp;
    }

  protected:

    double a1;      // recursion coefficient
    double s1, s2;  // past sine outputs

  };


  /*
  class RSLib_API rsSineCosineIterator : public rsComplexExponentialIterator
  {

  };
  class RSLib_API rsExponentialIterator
  {

  };
  */


}

#endif
