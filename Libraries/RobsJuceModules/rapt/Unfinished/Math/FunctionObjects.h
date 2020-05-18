#ifndef RAPT_FUNCTIONOBJECTS_H
#define RAPT_FUNCTIONOBJECTS_H

/** This class serves as baseclass for univariate scalar functions - that is, functions that take 
a scalar as input and produce a scalar at the output.

\todo the root-finding functions need some more proper stopping criteria. they are currently set
a bit ad hoc and do not guarantee that the function will always converge - do some research....

\todo generalize the root-finding (i.e. finding x for f(x) = 0) into value-finding (find x such
that f(x) = a). I think, we just need to modify the root-finding algorithms to subtract "a" after
each function evaluation. Make "a" an optional parameter that defaults to 0.

\todo implement methods for finding the definite integral between x=a and x=b, like
getIntegralRiemann(double a, double b);
getIntegralGaussian(double a, double b);
etc.
make a virtual function getIntegral(a, b); that defaults to one of the above numerical procedures
but can be overriden in subclasses, for functions where analytical expressions for the
antiderivative are known. see numerical integration methods in "Numerical Recipies in C".
....maybe it's possible to make a general function
getAntiDerivative(double x, double c); with integration constant c, that just returns
getIntegral(-inf, x) + c; ? ...don't know, if that makes sense

\todo include references for the root-finding methods (there was a paper called something like
"beyond Newton" - or something

// todo (maybe): implement functions like getDerivative(),
// getAntiDerivative(double integrationConstant = 0.0) that return another object of class
// UnivariateScalarFunction ...mmm...maybe not

*/

template<class T>
class UnivariateScalarFunction
{

public:

  //---------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  UnivariateScalarFunction() {}

  /** Destructor. */
  virtual ~UnivariateScalarFunction() {}

  //---------------------------------------------------------------------------------------------
  /** \name Evaluation */

  /** You must override this in your subclass to produce an output value for a given input value
  x. */
  virtual T getValueAt(T x) = 0;
    // maybe rename to getValue() or value()

  /** You may override this in your subclass to produce an output value of the first derivative
  of the function for a given input value x. If you don't override it, the baseclass
  implementation will invoke approximateDerivative. */
  virtual T getFirstDerivativeAt(T x);
    // rename to getDerivativeAt ...or maybe just getDerivative

  /** Similar to getFirstDerivativeAt, but for the second derivative. */
  virtual T getSecondDerivativeAt(T x);



  //virtual double getRootNear(double x0);


  // later put these into a protected or private section:

  /** Approximates the derivative of given order by using a finite difference approximation of
  the type: f'(x) ~ (f(x+eps) - f(x-eps)) / (2*eps). For higher derivatives, it uses the same
  approach recursively such that: f''(x) ~ (f'(x+eps) - f'(x-eps)) / (2*eps) where f' is
  approximated by an recursive call to this very function, likewise for higher derivatives.
  Note that the approximation becomes progressively worse for higher derivatives. */
  T approximateDerivativeAt(T x, int order, T eps);



  /** Finds a root between xL and xU.
  WARNING: This function is still preliminary - it doesn't always converge. */
  T findRootViaRidders(T xMin, T xMax);

  /** Finds a root of this function object, that is, a value x* for which (x) = 0 via a Newton
  iteration method. */
  T findRootViaNewtonNonRobust(T x0);

  /** Finds a root of this function object, that is, a value x* for which (x) = 0 via a
  Chebychev iteration method. */
  T findRootViaChebychevNonRobust(T x0);

};

template<class T>
class UnivariateScalarFunctionViaPointer : public UnivariateScalarFunction<T>
{

public:

  //---------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  UnivariateScalarFunctionViaPointer(T (*functionToUse) (T),
    T (*derivativeToUse)(T) = NULL);

//---------------------------------------------------------------------------------------------
/** \name Evaluation */

  virtual T getValueAt(T x);
  virtual T getFirstDerivativeAt(T x);

  //=============================================================================================

protected:

  // pointers to facilitate the use of simple c-style functions by wrapping them into an object:
  T (*functionPointer)   (T x);
  T (*derivativePointer) (T x);

};

/**

A class for representing polynomials of the form:
\f[ f(x) = a_0 + a_1 x + a_2 x^2 + \ldots + a_N x^N \f]
\todo: maybe rename to PolynomialReal - we may also want to have a class for complex polynomials

*/
/*
template<class T>
class Polynomial : public UnivariateScalarFunction
{

public:

  Polynomial(const unsigned int order = 0, const double* const coeffs = NULL);



  virtual ~Polynomial();




  //---------------------------------------------------------------------------------------------
  // setup:

  void zeroCoefficients();

  void setCoefficients(const double* const newCoeffs);

  //---------------------------------------------------------------------------------------------
  // evaluation:

  virtual double getValueAt(double x);
  virtual double getFirstDerivativeAt(double x);
  virtual double getSecondDerivativeAt(double x);

  //\todo: define operators for: multiplication, division, addition, subtraction, 
  //       (in)equality, etc.

  //=============================================================================================

protected:

  unsigned int order;
  double       *coeffs;

};
*/

/** This class serves as baseclass for multivariate scalar functions - that is, functions that take a
vector as input and produce a scalar at the output. */

template<class T>
class MultivariateScalarFunction
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultivariateScalarFunction(int numInputs);

  /** Destructor. */
  ~MultivariateScalarFunction();

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Informs, whether or not the function supports the evaluation of the gradient. The baseclass
  implmentation will return false. If you return true here in your subclass, you must also
  override getGradient. */
  virtual bool supportsGradient() { return false; }

  /** Returns the number of inputs of the function (the dimensionality of the input vector). */
  virtual int getNumInputs() { return numInputs; }

  //---------------------------------------------------------------------------------------------
  // evaluation:

  /** Sets the input vector - all subsequent evaluations of the function value and gradient are
  with respect to this input vector. */
  virtual void setInputVector(rsVectorDbl x) = 0;

  /** Computes an output number from a given input Vector. */
  virtual T getValue() = 0;

  /** Returns the local gradient at the given input vector x, that is, the vector of partial
  derivatives with respect to the elements of x. You should override this in your subclass if you
  return true in supportsGradient - the baseclass implementation will only return an empty dummy
  vector. */
  virtual rsVectorDbl getGradient() { return rsVectorDbl(); }

  //=============================================================================================

protected:

  int numInputs;

};


/** This class serves as baseclass for multivariate error functions - that is, functions that take 
a vector of parameters as input and produce a scalar error at the output. */

template<class T>
class MultivariateErrorFunction
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultivariateErrorFunction();

  /** Destructor. */
  ~MultivariateErrorFunction();

  //---------------------------------------------------------------------------------------------
  // evaluation:

  /** Sets the parameter vector - all subsequent evaluations of the function value and gradient
  are with respect to this input vector. */
  //virtual void setParameterVector(rsVectorDbl newParameters) = 0; 

  /** Computes an output number from a given input Vector. */
  virtual T getValue(rsVectorDbl p) = 0;

  /** Returns the local gradient at the given input vector p, that is, the vector of partial
  derivatives with respect to the elements of p. The baseclass implementation will approximate
  the gradient by a central difference which perturbs each element individually and evaluates
  the function itself via getValue - this requires 2*N evaluations of the function with N being
  the dimensionality of the parameter vector. See (1), page 147 for more details. This is
  expensive, so you really should override this function in your subclass if you have some better
  algorithm to compute the gradient for the problem at hand.  */
  virtual rsVectorDbl getGradient(rsVectorDbl p);
   // use const reference to vector of template type T

  /** Approximates the product v^T * H at the point p in parameter space where H denotes the
  local Hessian matrix. The approximation is based on a central difference of two local gradients
  at p + eps*v and p - eps*v, where eps is some small constant. Thus, the function will call
  getGradient two times. Note that in optimization algorithms like conjugate gradient, one of
  these two calls can be avoided when the central difference is replaced by a one-sided
  difference. See (1), page 158 for more details. */
  virtual rsVectorDbl getVectorTimesHessianApproximate(rsVectorDbl p, rsVectorDbl v);

};


/**

This class realizes a simple bivariate quadratic error function of the form
E = (1/2) * p^T * A * p - b^T * p + c  with p denoting the parameter vector to be optimized and
the constants c (scalar), b (vector) and A (matrix). This function object is intended to be used
as a test error function for testing optimization algorithms.

...move this into the test-suite - or generalize it to make it useful in other contexts (let A,
b, c be user variables) - rename to QuadraticForm or Quadric or something (look up, what it is called)

-get rid of the 1/2 factor

*/

template<class T>
class QuadraticTestErrorFunction : public MultivariateErrorFunction<T>
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  QuadraticTestErrorFunction();

  /** Destructor. */
  ~QuadraticTestErrorFunction();

  //---------------------------------------------------------------------------------------------
  // evaluation:

  /** Computes an output number from a given input Vector p. */
  virtual T getValue(rsVectorDbl p) { return  0.5*p*(A*p) - b*p + c; }

  /** Returns the local gradient at the given input vector p. */
  virtual rsVectorDbl getGradient(rsVectorDbl p) { return 0.5*(trans(A)+A)*p - b; }

  //=============================================================================================

protected:

  rsMatrixDbl A;
  rsVectorDbl b;
  T c;

};



/** This class serves as baseclass for multivariate vector functions - that is, functions that 
take a vector as input and produce a vector at the output. */

class MultivariateVectorFunction
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultivariateVectorFunction(int numInputs, int numOutputs);

  /** Destructor. */
  ~MultivariateVectorFunction();

  //---------------------------------------------------------------------------------------------
  // inquiry:

  // supportsJacobian

  /** Returns the number of inputs of the function (the dimensionality of the input vector). */
  virtual int getNumInputs() { return numInputs; }

  /** Returns the number of outputs of the function (the dimensionality of the output vector). */
  virtual int getNumOutputs() { return numOutputs; }

  //---------------------------------------------------------------------------------------------
  // evaluation:

  /** Computes an output vector from a given input Vector. */
  virtual rsVectorDbl getOutputVector(rsVectorDbl x) = 0;

  // virtual rsMatrixDbl getJacobian(rsVectorDbl x);

  //=============================================================================================

protected:

  int numInputs, numOutputs;

};

#endif 
