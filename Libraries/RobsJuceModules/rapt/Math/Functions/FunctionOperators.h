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
  T P = xMax - xMin;           // period length
  std::function<T(T)> f;
  f = [=] (T x) { 
    T d = ceil((x-xMax) / P);
    return func(x-P*d); };
  return f;
}

// todo: invert, numerical derivative, numerical integral (definite and indefinite, the latter 
// needs a lower limit as parameter, the former both limits), 