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
// move into class rsNumericDifferentiatior and rename to firstDerivative or just derivative

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
