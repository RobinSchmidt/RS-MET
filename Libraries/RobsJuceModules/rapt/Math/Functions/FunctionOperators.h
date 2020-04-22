#pragma once

/** This file contains operators in the mathematical sense of the term. In mathematics, an operator
is an object that takes a function as input and returns another function as output. An example for
an operator could be "take the derivative" which turns f(x) into f'(x) or "multiply by input 
argument" which turns f(x) into x*f(x). The input- and output functions are represented as 
std::function objects. */

/** Takes the segment of the function func from xMin to xMax and periodically repeats it. */
template<class T>
inline std::function<T(T)> rsMakePeriodic(const std::function<T(T)>& func, T xMin, T xMax)
{
  T P = xMax - xMin;  // period length
  return [=](T x) { T d = ceil((x-xMax) / P); return func(x-P*d); };
}
// maybe make a bivariate version of that

/** Returns the numerical derivative of f computed via a central difference using step-size h. */
template<class T>
inline std::function<T(T)> rsDerivative(const std::function<T(T)>& f, T h)
{
  return [=](T x) { return (f(x+h) - f(x-h)) / (2*h); };
}
// todo: 
// -maybe make versions that take the forward and backward difference
// -maybe let h be automatically selected - it sould be selected such that we get the best 
//  numerical accuracy - using smaller h reduces the approximation error but increases the the 
//  cancellation error (subtracting similar numbers) - so there should be an optimal choice 
//  somewhere - maybe there's a formula for that optimal choice based on the machine epsilon and
//  x and f(x)?
// -advantage of forward/backward difference: only one evaluation of f, may be used at points where
//  the function is undefined in either of the two directions from the point
// -maybe fit a parabola through x-h, x, x+h and take the derivative of that -> more accurate (?) 
 // ..or maybe the contribution of the (x,f(x)) point cancels and it actually gives the same result?

/** Extracts the even part of the function, i.e. the part with even symmetry (axial symmetry with 
respect to the y-axis). */
template<class T>
inline std::function<T(T)> rsEvenPart(const std::function<T(T)>& f)
{
  return [=](T x) { return (f(x) + f(-x)) / 2; };
}

/** Extracts the odd part of the function, i.e. the part with odd symmetry (point symmetry with 
respect to the origin). */
template<class T>
inline std::function<T(T)> rsOddPart(const std::function<T(T)>& f)
{
  return [=](T x) { return (f(x) - f(-x)) / 2; };
}

// -make scaled/shifted version(see class rsScaledAndShiftedSigmoid)

// todo: invert, numerical integral (definite and indefinite, the latter 
// needs a lower limit as parameter, the former both limits), 
// maybe implement a convolution (by some specific function), i.e.
// rsConvolve(const std::function<T(T)>& f, const std::function<T(T)>& g) - maybe the convolution
// integral should be limited to a finite range (i.e. finite instead of infinite integration 
// limits)
// maybe de-inline