#ifndef RAPT_NUMERICCALCULUS_H
#define RAPT_NUMERICCALCULUS_H

// todo: wrap into class
// -maybe rsNumCalc or rsNumericalCalculus - get rid of the then redundant "Numeric" in the 
//  function names

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

/** Computes the numerical integral of a function defined by data points. Usage is similar to
rsNumericDerivative... 
The parameter c is the integration constant and determines yi[0]. This shifts the overall 
resulting function up or down along the y-axis. */
template<class T>
void rsNumericIntegral(T *x, T *y, T *yi, int N, T c = T(0));
// maybe rename to rsNumericIntegralTrapezoidal


// Maybe rename to NumericAnalysis and include the interpolation stuff into this file as well 
// because some interpolation stuff depends on numeric derivatives but some numeric derivatives/
// integration stuff may depend on intrepolation and if we templatize the functions, we need to 
// take care that everything is defined before it gets used.



#endif
