#ifndef RAPT_ROOTFINDER_H_INCLUDED
#define RAPT_ROOTFINDER_H_INCLUDED

/** This class implements various (one-dimensional) root finding algorithms, i.e. it finds 
solutions to the equation f(x) = 0 for an arbitrary given function f(x). For more generality, the
right hand side does not actually need to be equal to zero but can be any target value, so it 
actually finds solutions to f(x) = y for given f(x) and given y (which defaults to zero). This 
little addition to standard textbook root-finding algorithms doesn't change much algorithmically 
(we just need to subtract the target value in each function evaluation), yet adds a lot to the 
flexibility. Some algorithms require the user to pass an initial interval that is assumed to 
bracket the root, others require an initial estimate of the root.

References
(1) Numerical Recipies in C (2nd Edition), Chapter 9

*/

template<class T>
class rsRootFinder
{
public:

  /** Bisection takes an initial interval xLeft, xRight (assumed to bracket the root) and evaluates 
  the function at the midpoint of the interval. Depending on the function value, the midpoint 
  becomes either the new left or the new right border of the bracketing interval. So, in each 
  iteration, the size of the interval is halved. The order of convergence is linear (in each 
  iteration, we get one more correct binary digit in the root estimate) and convergence is 
  guaranteed. That means: the method is slow but safe. */
  static T bisection(const std::function<T(T)>& func, T xLeft, T xRight, T y = 0);
  // -Maybe instead of std::function use a second template parameter F
  // -Maybe declare the template parameters in front of the functions, not the class (like in 
  //  rsArrayTools)
  // -Let the function take a tolerance parameter (maybe defaulting 
  //  std::numeric_limits<T>::epsilon)
  // -Let the function take a maxNumIterations paremeter.
  // -When adding these additional parameters, make sure that their API and semantics matches that
  //  in rsMinimizer1D (in Optimization.h/cpp). 

  /** Similar to bisection but doesn't use the midpoint of the current bracketing interval, but the
  point where a line between (xLeft,yLeft), (xRight,yRight) crosses the x-axis. Convergence is
  typically (for smooth functions, that are well approximated by a line near the root) faster than
  for bisection and convergence is guranteed. But for some functions, the convergence may also be 
  slower than bisection. The method is also known as "regula falsi". */
  static T falsePosition(const std::function<T(T)>& func, T xLeft, T xRight, T y = 0);

  // todo: let caller pass a tolerance (maybe separate for x and y), state pathological conditions
  // when no convergence can be epxected (like when the function has poles like 1/x)


  /** For a given function f = f(x) and target vaule y, finds a reasonable left bracket xL such 
  that f(xL) <= y. You can pass an initial guess for xL as well as an initial distance d by which
  xL will be decremented, in case, the condition f(xL) doesn't hold. */
  //static inline T findLeftBracket(const std::function<T(T)>& f, T y, T xL = T(0), T d = T(1))
  //{
  //  while(f(xL) > y) { xL -= d; d *= 2; }
  //  return xL;
  //}
  //// hmm - but this function works only for increasing functions
  //// maybe the growth multiplier for d should be an optional user parameter defaulting to 2. 
  //// When doing so, we should assert that it is >= 1.

  //static inline T findRightBracket(const std::function<T(T)>& f, T y, T xR = T(0), T d = T(1))
  //{
  //  while(f(xR) < y) { xR += d; d *= 2; }
  //  return xR;
  //}



  //static T modifiedFalsePosition(std::function<T(T)>& func, T xLeft, T xRight, T y = 0);


  // secant, newton, ridders, brent, ...

  // see RSLib::UnivariateScalarFunction in file FunctionObjects.h/cpp
  

};


#endif
