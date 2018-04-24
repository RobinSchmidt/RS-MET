#ifndef RAPT_NUMERICCALCULUS_H
#define RAPT_NUMERICCALCULUS_H

// todo: wrap into class

/** Given an array of strictly monotonically increasing but not necessarily equidistant abscissa
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
N-1. In a test with a sine function, the extrapolation gave more accurate reults at the
endpoints compared to simple differences, so it's probably better to use extrapolation. */
template<class T>
void rsNumericDerivative(T *x, T *y, T *yd, int N, bool extrapolateEnds = true);

// maybe make the function a template, possibly even with different template types for 
// x, y, yd so it can be used for complex or multivariate data as well (in the latter case
// the derivative would be a matrix namely the Jacobian)

// maybe write a function that computes the numeric derivative a one particular datapoint and
// ideally also higher order derivatives at that point - that's more convenient to use because
// client code does not need to have buffers for all derivatives

// try another approach: fit a polynomial of arbitrary order to a number of datapoints around
// the n and return the derivative of the poynomial at that point (may this be equivalent to the
// approach above when using 3 points for a quadratic polynomial?)

// write numerical integration functions. these should possibly also generalize to multivariate 
// functions when intergration along the x-axis is replaced by path-integration - so maybe pass
// an array of "x"-values for the path

// write N-dimensional integration functions that return the amount of N+1 space contained in
// some hyperblock between x1, x2 (both of dimensionality N)

// Maybe rename to NumericAnalysis and include the interpolation stuff into this file as well 
// because some interpolation stuff depends on numeric derivatives but some numeric derivatives/
// integration stuff may depend on intrepolation and if we templatize the functions, we need to 
// take care that everything is defined before it gets used.

#endif
