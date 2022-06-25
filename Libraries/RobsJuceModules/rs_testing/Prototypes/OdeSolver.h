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

/** A class for solving initial value problems for ordinary differential equations (ODEs). The 
problem has to be given in the general form:

  y' = f(y)

where y is a state vector of the system and y' is its derivative with respect to the independent 
variable which is typically time, so we'll denote it by t. The function f shall be passed as an 
object of type std::function that takes the array representing y as first parameter and produces
the corresponding y' in its second (output) parameter. This class is intended for ODE systems with
arbitrary dimensionality. The stepper methods use a loop over the dimensions.

...tbc...  */

template<class T>
class rsInitialValueSolver2
{

public:

  using Func = std::function<void(const T* y, T* dy)>;

  static void stepForwardEuler(const Func& f, int N, T* y, T* v, T h);


  void init(int dimensionality, T* initialPosition, T* initialVelocity);


protected:

  T*   y = nullptr; // current state vector, i.e. position in phase-space y = y(t)
  T*   v = nullptr; // current derivative, i.e. velocity in phase space, v = y' = dy/dt
  T    h = 1;       // step size
  int  N = 0;       // dimensionality of the system, length of y and v
  Func f;           // function to compute y' = f(y)

};

// ToDo:
// -Define, what should be passed into f as "initial guess" for y' - namely, the computed y' from
//  the previous iteration. If the ODE system is given explicitly as y' = f(y), the f can just 
//  ignore the initial content of dy, but if the ODE given implicitly as f(y,y') = 0, f may use 
//  the initial content of dy as initial guess for e.g. Newton iteration
// -Maybe the y and v arrays should be owned by client code. For multistep methods, we may require
//  that they must be of length N*a, N*b where a and b are the numbers of past states for y and v
//  that the repective method uses, e.g. a=1, b=2 for 2nd order Adams-Bashforth or a=2, b=1 for 
//  BDF-2 (I think)







//=================================================================================================
// Another idea with a different API

/** Baseclass for all ODE solvers. Defines a common interface and implements the simplemost 
forward Euler method for reference. Subclasses are supposed to extend that very basic functionality 
by more sophisticated solver methods (higher order, implicit, multistep, predictor-corrector, etc.) 
and/or have more features like error estimation, stepsize adaption, etc. 

This implementation is intended for ODE systems of low dimensionality such that states can be held
in variables like rsVector2D or rsVector3D which do not require memory allocations or loops on 
assigment.

*/

template<class Tx, class Ty>
class rsOdeSolverBase   // rename to rsInitialValueSolver, we may also want an rsBoundaryValueSolver
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
  // switch between fixed and variable step sizes. In the case of variable stepsizes, this function
  // here will serve for initializing the stepsize

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

  /** Resets x,y,y' to all zeros. Does not reset the step-size (should it? maybe not). */
  virtual void reset() { x = 0; y = 0; yd = 0; }
  // maybe allow the user to pass values - but maybe that should be left to another function that
  // we may call setState. Maybe that function should avoid the assignment operator because we
  // anticipate that the class Ty may be something like std::vector<float> and we want to avoid
  // allocations in setState. ...not yet sure how to best solve that...


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
    forwardEuler,     //   no
    backwardEuler,    //   yes
    trapezoidal,      //   yes
    midpoint,         //    ?
    heun,             //    ?
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

