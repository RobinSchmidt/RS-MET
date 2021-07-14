#pragma once

// just a stub at the moment

// maybe rename file to OdeSolvers (plural!)

//=================================================================================================

/** A class that tabulates and/or implements algorithms to compute the various coefficients that
occur in the different numerical methods to solve ordinary differential equations (ODEs). */

template<class T>
class rsOdeCoeffs
{

public:

  /** Computes the coefficients for an Adams-Bashforth multistep method. Directly implements the 
  formula given here:
    https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods
  without any attempt to optimize the efficiency. I recommend to use rsFraction<int> for the template
  parameter T. Adams-Bashforth methods are explicit multistep metods. The returned array will be of
  length "order". ...tbc... */
  static std::vector<T> adamsBashforth(int order);
  // Maybe rename "order" to numSteps. it's the number of steps, we look into the past

  /** Similar to adamsBashforth. See:
    https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Moulton_methods  
  Adams-Moulton methods are implicit methods. The returned array will be of length "order"+1. 
  ...tbc...   */
  static std::vector<T> adamsMoulton(int order);
  // Maybe rename "order" to numSteps. 


protected:

  // ToDo: 
  // -store the arrays of coeffs for the BDF formulas here (as static const), also the Butcher 
  //  tableaus of various Runge-Kutta type methods
  // -have functions with an API similar to the Adams-... coeff computation methods but instead of
  //  returning the vector as value, return it as const reference to our static members
  // -Or: let all functions have the API: adamsBashforth(int order, T* coeffs, T* work). Then 
  //  they are called all with the same API. The ones using tabulated values would just copy them
  //  into the output arrays. Client code does need to know or care, if values are computed or
  //  tabulated.
   
};

//=================================================================================================

/** Baseclass for all ODE solvers. Defines a common interface and implements the simplemost 
forward Euler method. Subclasses are supposed to extend that basic functionality by more 
sophisticated solver methods, error estimation, stepsize control, etc. */

template<class Tx, class Ty>
class rsOdeSolverBase   // renam to rsInitialValueSolver, we may also want rsBoundaryValueSolver
{

public:

  using Func = std::function<void(Tx x, const Ty& y, Ty& yd)>;


  /** Sets the function F to be used to compute y' = F(x,y). It takes the y' as output parameter 
  in order to avoid having to return in by value which would be costly, if Ty is a type that 
  allocates heap memory such as std::vector. Moreover, it allows us to pass an initial guess in 
  case F implements an implicit root finding algorithm, i.e. F is given implicitly as 
  F(x,y,y') = 0. The passed value for y' will be the y' value from the previous call to doStep. */
  void setFunction(const Func& newFunction) { F = newFunction; }

  void setStepSize(Tx newStepSize) { h = newStepSize; }
  // maybe it should have a 2nd boolean parameter "fixed" that can be used later by subclasses to 
  // switch between fixed and variable step sizes

  /** Performs one integration step at a time. The baseclass implementation just does the 
  simplemost thing: a forward Euler step. Subclasses can override this to implement more 
  sophisticated methods (such as Runge-Kutta, Adams-Bashforth, etc.). */
  virtual void doStep()
  {
    F(x, y, yd);  // compute derivative y' = F(x,y)
    x += h;       // advance x by stepsize h
    y += h*yd;    // advance y by a forward Euler step
  }
  // maybe it should take a bool adaptStepSize defaulting to false

  /** Resets x,y,y' to 0. Does not reset the setp-size (should it? maybe not). */
  virtual void reset() { x = 0; y = 0; yd = 0; }


  const Tx& getAbscissa() const { return x; }

  const Ty& getValue() const { return y; }

  const Ty& getDerivative() const { return yd; }



protected:

  // Current values:
  Tx h  = 1;  // step size
  Tx x  = 0;  // independent variable (maybe use t)
  Ty y  = 0;  // dependent variable
  Ty yd = 0;  // derivative of y with respect to x, aka y'
  Func F;     // function to compute y' = F(x,y) or F(x,y,y') = 0

};

//=================================================================================================

/** Implements explicit or implicit single step methods such as Runge-Kutta, ... */

template<class Tx, class Ty>
class rsOdeSolverSingleStep : public rsOdeSolverBase<Tx, Ty>
{

public:

  enum class StepMethod
  {                   // implicit
    forwardEuler,     //
    backwardEuler,    //   yes
    trapezoidal,      //   yes
    midpoint,         //
    heun,             //
    rungeKutta3,
    rungeKutta4,
    cashKarp
    // ...
  };



protected:

  // We need a std::function for solving the implicit equation that occurs in implicit methods. The
  // explicity methods can ignore that...or wait...do we? or do we have to solve a linear system in
  // case of an implicit method?

  // Butcher tableau:
  rsMatrix<Tx> A;
  std::vector<Tx> c, b;  // nodes and weights

  StepMethod stepMethod = StepMethod::forwardEuler;

  // maybe we should have a predictorMethod and a correctorMethod where use of the corrector is 
  // optional

};

//=================================================================================================

/** A class that implements multistep ODE solvers such as Adams-Bashforth, Adams-Moulton, Nyström,
backward difference formula (BDF), etc. */


template<class Tx, class Ty>
class rsOdeSolverMultiStep : public rsOdeSolverSingleStep<Tx, Ty>
{

public:

  using Vec  = std::vector<Ty>;
  //using Func = std::function<void(Tx x, const Ty& y, Ty& yd)>;
  // Type for the function to compute y' from y should take both as reference parameters, such that 
  // we can use std::vector as type T without having to re-allocate.

  
protected:



  Vec Y, Yd; // arrays of stored past values of y, y': Y[0] = y[n], Y[1] = y[n-1], ...


  /*
  Tx   x;     // the independent variable...maybe use t

  Func F;     // function to compute y' = F(x,y)
  //Func fp;    // 
  */

};

