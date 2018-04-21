#ifndef RAPT_CURVEFITTING_H
#define RAPT_CURVEFITTING_H

/** Fits a weigthed sum of exponentials (with weights and exponents to be determined) to the data
points, such that y[n] = A[0]*exp(a[0]*n) + A[1]*exp(a[1]*n) + ... where n = 0, 1, 2, ...
This function may fail - if it does, it will return false and leave the A and a arrays untouched.
If it succeeds, it returns true and the A, a arrays are filled with the weights and exponents
respectively.
The numExponentials parameter is currently assumed to be numValues/2 - but for later
generalizations to least-squares fitting (where we will allow numExponentials <= numValues/2), it
was already included as parameter. */
template<class T>
bool rsFitSumOfExponentials(T* y, int numValues, T* A, T* a, int numExponentials);

// wrap into a class rsCurveFitting

#endif
