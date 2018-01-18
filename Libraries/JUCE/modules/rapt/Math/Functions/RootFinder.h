#ifndef RAPT_ROOTFINDER_H_INCLUDED
#define RAPT_ROOTFINDER_H_INCLUDED

/** This class implements various (one-dimensional) root finding algorithms, i.e. it finds 
solutions to the equation f(x) = 0 for an arbitrary given function f(x). For more generality, the
right hand side does not actually need to be equal to zero but can be any target value, so it 
actually find solutions to f(x) = y for given f(x) and given y. This little addition to standard 
textbook root-finding algorithms doesn't change much algorithmically (we just need to subtract the 
target value in each function evaluation), yet adds a lot to the flexibility. */

template<class Tx, class Ty>
class rsRootFinder
{
public:

  /** Given a function and two x-values xLeft, xRight that bracket the root, this function finds
  the root via the bisection algorithm. For more generality, the target-value for y does not need
  to be zero but can be passed as additional parameter. */
  static Tx bisection(std::function<Ty(Tx)>& func, Tx xLeft, Tx xRight, Ty y = 0);
  
  //Tx secant(Tx xLeft, Tx xRight, Ty y = 0);
  
  // newton, ridders, brent, ...

  // see RSLib::UnivariateScalarFunction in file FunctionObjects.h/cpp
  

};


#endif
