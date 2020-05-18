#ifndef RAPT_CURVEFITTING_H
#define RAPT_CURVEFITTING_H

class rsCurveFitter
{

public:

/** Fits a weigthed sum of exponentials (with weights and exponents to be determined) to the data
points, such that y[n] = A[0]*exp(a[0]*n) + A[1]*exp(a[1]*n) + ... where n = 0, 1, 2, ...
This function may fail - if it does, it will return false and leave the A and a arrays untouched.
If it succeeds, it returns true and the A, a arrays are filled with the weights and exponents
respectively.
The numExponentials parameter is currently assumed to be numValues/2 - but for later
generalizations to least-squares fitting (where we will allow numExponentials <= numValues/2), it
was already included as parameter. */
template<class T>
static bool fitExponentialSum(T* y, int numValues, T* A, T* a, int numExponentials);


// add fitLine, fitPolynomial, fitRational, 
// fitPower: a * x^p + b ...p can be found by fitting a line to log-of-y-values

// fitFunctionSet -> should take a set of basis functions that are linearly combined
// maybe we can have a linear combination of basis functions in the numerator and in the 
// denominator, i.e. f(x) = (a0*f0(x) + ... + aN*fN(x)) / (b0*g0(x) + ... + bM*gM(x))
// this is a quite general set of functions but the coefficients can still be found by solving
// a linear system of equations - ahh..wait - no - this may work only when interpolating with
// a set of arbitrary basis-functions but nor least-squares fitting...or maybe it will?
// if not, maybe we can use a smoothed version of the data to find an interpolant and then
// use the resulting coeffs of the interpolant as start-values for gradient descent?

};

#endif
