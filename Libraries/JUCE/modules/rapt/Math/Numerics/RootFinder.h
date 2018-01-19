#ifndef RAPT_ROOTFINDER_H_INCLUDED
#define RAPT_ROOTFINDER_H_INCLUDED

/** This class implements various (one-dimensional) root finding algorithms, i.e. it finds 
solutions to the equation f(x) = 0 for an arbitrary given function f(x). For more generality, the
right hand side does not actually need to be equal to zero but can be any target value, so it 
actually find solutions to f(x) = y for given f(x) and given y. This little addition to standard 
textbook root-finding algorithms doesn't change much algorithmically (we just need to subtract the 
target value in each function evaluation), yet adds a lot to the flexibility. Some algorithms 
require the user to pass an initial interval that is assumed to bracket the root, others require an
initial estimate of the root.

References
(1) Numerical Recipies in C (2nd Edition), Chapter 9

*/

template<class Tx, class Ty>
class rsRootFinder
{
public:

  /** Bisection takes an initial interval xLeft, xRight (assumed to bracket the root) and evaluates 
  the function at the midpoint of the interval. Depending on the function value, the midpoint 
  becomes either the new left or the new right border of the bracketing interval. So, in each 
  iteration, the size of the interval is halved. The order of convergence is linear (in each 
  iteration, we get one more correct binary digit in the root estimate) and convergence is 
  guaranteed. */
  static Tx bisection(std::function<Ty(Tx)>& func, Tx xLeft, Tx xRight, Ty y = 0);

  /** Similar to bisection but doesn't use the midpoint of the current bracketing interval, but the
  point where a line between (xLeft,yLeft), (xRight,yRight) crosses the x-axis. Convergence is
  typically (for smooth functions, that are well approxiamted by a line near the root) faster than
  for bisection and convergence is guranteed. */
  static Tx falsePosition(std::function<Ty(Tx)>& func, Tx xLeft, Tx xRight, Ty y = 0);


  
  //Tx secant(Tx xLeft, Tx xRight, Ty y = 0);
  
  // newton, ridders, brent, ...

  // see RSLib::UnivariateScalarFunction in file FunctionObjects.h/cpp
  

};


#endif
