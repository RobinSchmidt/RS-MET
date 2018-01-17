#ifndef RAPT_ROOTFINDER_H_INCLUDED
#define RAPT_ROOTFINDER_H_INCLUDED

template<class Tx, class Ty>
class rsRootFinder
{
public:

  /**  */
  static Tx bisection(std::function<Ty(Tx)>& func, Tx xLeft, Tx xRight, Ty y = 0);
  
  //Tx secant(Tx xLeft, Tx xRight, Ty y = 0);
  
  // newton, ridders, brent, ...

  // see RSLib::UnivariateScalarFunction in file FunctionObjects.h/cpp
  

};


#endif
