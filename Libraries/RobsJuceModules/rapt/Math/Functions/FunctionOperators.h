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

/** Returns the numerical derivative computed via a central difference using step-size h. */
template<class T>
inline std::function<T(T)> rsDerivative(const std::function<T(T)>& f, T h)
{
  return [=](T x) { return (f(x+h) - f(x-h)) / (2*h); };
}
// todo: 
// -maybe make versions that take the forward and backward difference
// -maybe let h be automatically selected - it sould be selected such that we get the best 
//  numerical accuracy - using smaller h makes reduces the approximation error but increases the
//  the cancellation error (subtracting similar numbers) - so there should be an optimal choice 
//  somewhere - maybe there's a formula for that optimal choice based on the machine epsilon and
//  x and f(x)?

// todo: invert, numerical integral (definite and indefinite, the latter 
// needs a lower limit as parameter, the former both limits), 
// maybe de-inline