#ifndef RAPT_NUMERICCALCULUS_H
#define RAPT_NUMERICCALCULUS_H

// todo: wrap into class
// -maybe rsNumCalc or rsNumericalCalculus - get rid of the then redundant "Numeric" in the 
//  function names - and also of the rs-prefixes

/** Cannot be used in-place yet: y and yd have to be distinct!
Given an array of strictly monotonically increasing but not necessarily equidistant abscissa
values in x and corresponding function values in y, this function fills the array yd with a
numeric approximation of the derivative for each x value. All arrays are of length N. To
compute the numeric derivative, we use a weighted average of the difference quotients left and
right to the data point:

               y[n] - y[n-1]          y[n+1] - y[n]
 yd[n] = wL * --------------- + wR * ---------------
               x[n] - x[n-1]          x[n+1] - x[n]

where the weights for the left and right difference quotients are determined by the distances
dxL = x[n]-x[n-1] and dxR = x[n+1]-x[n] as wL = dxR/(dxL+dxR), wR = dxL/(dxL+dxR), such that
the closer the x-axis value x[n-1] is to x[n], the more weight is given for the left quotient
and vice versa. The weights add up to unity. If extrapolateEnds is true, the function will use
linear extrapolation of the inner derivative values for the endpoints yd[0] and yd[N-1],
otherwise it will use the (divided) forward difference at 0 and the backward difference at
N-1. In a test with a sine function, the extrapolation gave more accurate results at the
endpoints compared to simple differences, so it's probably better to use extrapolation. */
template<class Tx, class Ty>
void rsNumericDerivative(const Tx *x, const Ty *y, Ty *yd, int N, bool extrapolateEnds = true);
// move into class rsNumericDifferentiatior and rename to derivative

// todo: make a numeric derivative routine that is the inverse of the trapezoidal integrator
// rsDifferentiateTrapezoidal, rename this one to rsWeightedCentralDifference


/** Computes the numerical integral of a function defined by data points, i.e. the function:
\f[ F(x) = \int_c^x f(t) dt \f] where the lower integration limit c can be passed as a parameter 
into the function. Usage is similar to rsNumericDerivative. The parameter c can also be seen as an 
integration constant and determines yi[0] that shifts the overall resulting function up or down 
along the y-axis. The algorithm uses a trapezoidal rule, i.e. it sums up the areas under the 
trapezoids defined by a piecewise linear interpolant that passes through the datapoints. */
template<class Tx, class Ty>
void rsNumericIntegral(const Tx *x, const Ty *y, Ty *yi, int N, Ty c = Ty(0));
// maybe rename to rsNumericIntegralTrapezoidal, or rsIntegrateTrapezoidal, use Tx, Ty for 
// datatypes maybe rename parameters to x, f, F
// move to rsNumericIntegrator and rename to integrateTrapezoidal or just trapezoidal, have a
// simpler version riemannSum which uses the midpoint of each interval as evaluation point


// Maybe rename to NumericAnalysis and include the interpolation stuff into this file as well 
// because some interpolation stuff depends on numeric derivatives but some numeric derivatives/
// integration stuff may depend on interpolation and if we templatize the functions, we need to 
// take care that everything is defined before it gets used.



//=================================================================================================

/** A class for computing numerical approximations to derivatives of functions. In some cases, the
functions are assumed to be given a function object (aka functor), in other cases as an array of 
datapoints. The meaning of the 3 template parameters is:
  Tx: type of the input to the function (abscissa type)
  Ty: type of the output of the function (ordinate type)
  F:  type of the function object (if applicable)
 For example, Tx could be "double", Ty could be "rsVector2D<double>" and F could be 
 "std::function<rsVector2D<double>(double)>". This would apply to functions that take real scalars 
 (double) as input and produce 2D vectors (of real numbers) as output. Such functions define the 
 parametric equation of 2D curves. They depend on a scalar parameter t (which is interpreted as 
 time) and produce a 2D vector as output for each t. In this case, the meaning of the first 
 derivative would be the velocity vector and the second derivative would be the acceleration 
 vector. 
 
 References:
   https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences
   http://web.media.mit.edu/~crtaylor/calculator.html  */

template<class Tx, class Ty, class F>
class rsNumericDifferentiator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Functor derivatives

  /** Numeric approximation of the first derivative of function f at the value x with approximation
  step-size h. Uses a central difference and is 2nd order accurate in h.  */
  static Ty derivative(const F& f, const Tx& x, const Tx& h)
  {
    return (f(x+h) - f(x-h)) / (Tx(2)*h);
  }

  /** Numeric approximation of the second derivative. 2nd order accurate in h. */
  static Ty secondDerivative(const F& f, const Tx& x, const Tx& h)
  {
    return (f(x-h) - Tx(2)*f(x) + f(x+h)) / (h*h);
    // the Tx(2) may seem weird because it multiplies the result f(x) which is of type Ty, so one 
    // might expect Ty(2) here. However, if Tx is a scalar type and Ty is a vector type, this is
    // totally appropriate.
  }

  /** Numeric approximation of the third derivative. 3rd order accurate in h (i think - verify). */
  static Ty thirdDerivative(const F& f, const Tx& x, const Tx& h)
  {
    return (-f(x-Tx(2)*h) + Tx(2)*f(x-h) - Tx(2)*f(x+h) + f(x+Tx(2)*h)) / (Tx(2)*h*h*h);
  }
  // needs test
  // coeffs found by:
  // http://web.media.mit.edu/~crtaylor/calculator.html
  // f_xxx = (-1*f[i-2]+2*f[i-1]+0*f[i+0]-2*f[i+1]+1*f[i+2])/(2*1.0*h**3)



  /** Computes 0th, 1st and 2nd derivative of f at x. This is more efficient than using the 
  separate functions. It uses 3 function evaluations whereas you would need 6, if you would call
  the function itself (1 evaluation) and derivative (2 evaluations) and secondDerivative (3 
  evaluations). The result is exactly the same - it just avoids to compute all the values twice. */
  static void derivativesUpTo2(const F& f, const Tx& x, const Tx& h, Ty* f0, Ty* f1, Ty* f2)
  {
    Ty fp = f(x+h);  // "plus"
    Ty fm = f(x-h);  // "minus"
    *f0 = f(x);
    *f1 = (fp - fm) / (Tx(2)*h);
    *f2 = (fm - Tx(2)*(*f0) + fp) / (h*h);
  }
  // needs test

  /** Computes 0th, 1st and 2nd derivative of f at x. Uses a 5-point stencil: -2,-1,0,1,2 */
  static void derivativesUpTo3(const F& f, const Tx& x, const Tx& h, 
    Ty* f0, Ty* f1, Ty* f2, Ty* f3)
  {
    // evaluate function at stencil points:
    Ty fm2 = f(x-2*h);  // "minus 2h", etc...
    Ty fm1 = f(x - h);
    Ty fc  = f(x    );  // centered ..get rid - assign to *f0 directly
    Ty fp1 = f(x + h);
    Ty fp2 = f(x+2*h);

    // form linear combinations to approximate derivatives:
    *f0 =                     fc;
    *f1 = ( fm2 -  8*fm1         +  8*fp1 - fp2) / (12*h);
    *f2 = (-fm2 + 16*fm1 - 30*fc + 16*fp1 - fp2) / (12*h*h);
    *f3 = (-fm2 +  2*fm1         -  2*fp1 + fp2) / (2*h*h*h);
  }
  // needs test
  // stencil: -2,-1,0,1,2
  // f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*h**1)
  // f_xx = (-1*f[i-2]+16*f[i-1]-30*f[i+0]+16*f[i+1]-1*f[i+2])/(12*1.0*h**2)
  // f_xxx = (-1*f[i-2]+2*f[i-1]+0*f[i+0]-2*f[i+1]+1*f[i+2])/(2*1.0*h**3)

  // todo:
  // -maybe it's better to not use x += h, x +- 2h, x +- 3h but instead 
  //  x +- h/3, x +- 2h/3, x +- h such that the total width of the stencil stays the same, 
  //  regardless of how many stencil-points we use. the goal is that the optimal choice of h 
  //  depends only on the problem, not on the number of stencil-points
  // -try stencils where the points are not distributed equidistantly but maybe exponentially, for
  //  example: x +- h/4, x +- h/2, x +- h for a 5-point stencil
  // -to avoid numerical error, it is desirable that the offsets (2, 2h, 3h, ..) are exactly 
  //  representable - also the final divisors should be exactly representable - maybe (inverse)
  //  powers of two are a good choice - at least, for the basic 3-point stencil x +- h, x +- 2h
  //  where the divisor is 1/(2h)
  // -test it with some standard functions like exp, log, sin, cos, tan, 1/x, 1/(1+x^2) - maybe 
  //  polynomials (in which case we should get exact results, if the number of stencil points 
  //  matches the degree)


  // make a similar function to compute up to 3rd derivative - uses 5 function evaluations
  // (use stencil: -2,-1,0,+1,+2)

  // maybe have more accurate formulas - formulas can be produced by
  // http://web.media.mit.edu/~crtaylor/calculator.html
  // ..but i have also implemented this stencil coeff computation myself somewhere in the 
  // experiments - clean that code up and move it over here, too

  // todo: figure out the accuracy experimentally - maybe this can be done by testing, how high 
  // degree a polynomial can be such that we still get perfect results - i think, a 5-point stencil
  // should/ be perfect for polynomials up to 5th degree (it's based on an interpolating polynomial
  // of degree 5)



  //-----------------------------------------------------------------------------------------------
  // \name Data derivatives

  // move the function rsNumericDerivative here

};



//=================================================================================================

/** just a stub, at the moment */

template<class Tx, class Ty>
class rsNumericIntegrator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setNumberOfSamplePoints(int newNumSamples) { numSamples = newNumSamples; }


  //-----------------------------------------------------------------------------------------------
  // \name Integration

  /** Computes the definite integral of f in the interval from a to b. */
  //Ty integrate(const std::function<Ty(Tx)>& f, Tx a, Tx b);


protected:

  int numSamples = 10;

};


#endif
