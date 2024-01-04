#pragma once

/** Minimizes the sum of squared differences between adjacent array elements under the constraint 
that the sums of adjacent array elements must be equal to given values. The input array s is the 
length N-1 array of the desired sums, the output array v is the length N value array, such that
v[i] + v[i+1] = s[i] and sum_i (v[i+1] - v[i])^2 = min. You may optionally pass an array of 
weights for the squared differences in the cost function - if you do, the w array must have the 
same length as s, if you don't, unit weights will be used for each squared difference. With 
weights, we will minimize sum_i w[i] * (v[i+1] - v[i])^2 subject to the (same) constraints that
v[i] + v[i+1] = s[i] for all i = 0,..,N-2 */
template<class T>
void rsMinSqrDifFixSum(T* v, int N, T* s, T* w = nullptr);
// allocates heap-memory - todo: make a version that takes a workspace array as parameter - may be
// nullptr too in which case the function allocates the workspace itself and cleans it up 
// afterwards (for more convenient use)
// maybe rename the function with suffix "Simple" or "Fast" or something to indicate that it uses 
// the simple pentadiagonal algorithm - and make a version of that function that uses class
// rsBandDiagonalSolver for better numeric precision



/** Stub.

A class with functions for general 1D function minimization. */

template<class T>
class rsMinimizer1D
{

public:


  /** UNDER CONSTRUCTION
  A minimization method that works analoguous to the bisection method of root-finding. You need to 
  pass the left and right limits of the interval between which the minimum shall be found. 
  ...TBC...
  See comments in implementation for algorithmic details. */
  static T bisection(const std::function<T(T)>& func, T xLeft, T xRight);



};

// This is only here to make the compile happy for the time being. There's a unit test under 
// construction that calls this fucntion but the implementation should actually live in the .cpp
// file but when we put it there, we get a linker error. The method should actually be instantiated
// in the same way as the corresponding method of class rsRootFinder. For this class, we also have 
// a unit test and it links just fine. We need to find out where and why exactly *that* gets 
// instantiated and then make sure, that it works the same way for rsMinimizer1D:
template<class T>
T rsMinimizer1D<T>::bisection(const std::function<T(T)>& f, T xL, T xR)
{
  rsError("Not yet implemented correctly!");
  return T(0);
}

