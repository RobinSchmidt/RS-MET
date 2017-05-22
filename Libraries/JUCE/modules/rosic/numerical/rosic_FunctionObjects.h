#ifndef rosic_FunctionObjects_h
#define rosic_FunctionObjects_h

//// rosic-indcludes:
//#include "../math/rosic_MatrixVectorFunctions.h"

namespace rosic
{

  /**

  This class serves as baseclass for univariate scalar functions - that is, functions that take a 
  scalar as input and produce a scalar at the output.

  */

  class UnivariateScalarFunction  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    UnivariateScalarFunction() {}

    /** Destructor. */
    virtual ~UnivariateScalarFunction() {}

    //---------------------------------------------------------------------------------------------
    // evaluation:

    /** You must override this in your subclass to produce an output value for a given input value 
    x. */
    virtual double getValueAt(double x) = 0;

    /** You may override this in your subclass to produce an output value of the first derivative 
    of the function for a given input value x. If you don't override it, the baseclass 
    implementation will invoke approximateDerivative. */
    virtual double getFirstDerivativeAt(double x);

    /** Similar to getFirstDerivativeAt, but for the second derivative. */
    virtual double getSecondDerivativeAt(double x);



    //virtual double getRootNear(double x0);


    // later put these into a protected or private section:

    /** Approximates the derivative of given order by using a finite difference approximation of 
    the type: f'(x) ~ (f(x+eps) - f(x-eps)) / (2*eps). For higher derivatives, it uses the same 
    approach recursively such that: f''(x) ~ (f'(x+eps) - f'(x-eps)) / (2*eps) where f' is 
    approximated by an recursive call to this very function, likewise for higher derivatives. 
    Note that the approximation becomes progressively worse for higher derivatives. */
    double approximateDerivativeAt(double x, int order, double eps);




    double findRootViaRidders(double xMin, double xMax);

    /** Finds a root of this function object, that is, a value x* for which (x) = 0 via a Newton 
    iteration method. */
    double findRootViaNewtonNonRobust(double x0);

    /** Finds a root of this function object, that is, a value x* for which (x) = 0 via a 
    Chebychev iteration method. */
    double findRootViaChebychevNonRobust(double x0);



    // todo (maybe): implement functions like getDerivative(), 
    // getAntiDerivative(double integrationConstant = 0.0) that return another object of class 
    // UnivariateScalarFunction ...mmm...maybe not

  };






  class UnivariateScalarFunctionViaPointer : public UnivariateScalarFunction  
  {

    //---------------------------------------------------------------------------------------------
    // construction/destruction: 

    /** Constructor. */
    UnivariateScalarFunctionViaPointer(double (*functionToUse)   (double), 
                                       double (*derivativeToUse) (double) = NULL);  

    //---------------------------------------------------------------------------------------------
    // evaluation:

    virtual double getValueAt(double x);
    virtual double getFirstDerivativeAt(double x);

    //=============================================================================================

  protected:

    // pointers to facilitate the use of simple c-style functions by wrapping them into an object:
    double (*functionPointer)   (double x);  
    double (*derivativePointer) (double x);

  };



  /**

  A class for representing polynomials of the form:
  \f[ f(x) = a_0 + a_1 x + a_2 x^2 + \ldots + a_N x^N \f]
  \todo: maybe rename to PolynomialReal - we may also want to have a class for complex polynomials

  */

  class Polynomial : public UnivariateScalarFunction  
  {

  public:

    Polynomial(const unsigned int order = 0, const double* const coeffs = NULL);  



    virtual ~Polynomial();  




    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets all coefficients to zero. */
    void zeroCoefficients();

    /** Sets all coefficients to the values passed in the array pointed to by "newCoeffs" - this 
    array should have a length of (at least) order+1, where order is the currently set order of 
    this object. */
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







  /**

  This class serves as baseclass for multivariate scalar functions - that is, functions that take a
  vector as input and produce a scalar at the output.

  */

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
    virtual void setInputVector(Vector x) = 0;

    /** Computes an output number from a given input Vector. */
    virtual double getValue() = 0;

    /** Returns the local gradient at the given input vector x, that is, the vector of partial 
    derivatives with respect to the elements of x. You should override this in your subclass if you 
    return true in supportsGradient - the baseclass implementation will only return an empty dummy
    vector. */
    virtual Vector getGradient() { return Vector(); }

    //=============================================================================================

  protected:

    int numInputs;

  };


  /**

  This class serves as baseclass for multivariate error functions - that is, functions that take a
  vector of parameters as input and produce a scalar error at the output.

  */

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
    //virtual void setParameterVector(Vector newParameters) = 0; 

    /** Computes an output number from a given input Vector. */
    virtual double getValue(Vector p) = 0;

    /** Returns the local gradient at the given input vector p, that is, the vector of partial 
    derivatives with respect to the elements of p. The baseclass implementation will approximate
    the gradient by a central difference which perturbs each element individually and evaluates 
    the function itself via getValue - this requires 2*N evaluations of the function with N being 
    the dimensionality of the parameter vector. See (1), page 147 for more details. This is 
    expensive, so you really should override this function in your subclass if you have some better 
    algorithm to compute the gradient for the problem at hand.  */
    virtual Vector getGradient(Vector p);

    /** Approximates the product v^T * H at the point p in parameter space where H denotes the 
    local Hessian matrix. The approximation is based on a central difference of two local gradients
    at p + eps*v and p - eps*v, where eps is some small constant. Thus, the function will call
    getGradient two times. Note that in optimization algorithms like conjugate gradient, one of 
    these two calls can be avoided when the central difference is replaced by a one-sided 
    difference. See (1), page 158 for more details. */
    virtual Vector getVectorTimesHessianApproximate(Vector p, Vector v);

  };


  /**

  This class realizes a simple bivariate quadratic error function of the form 
  E = (1/2) * p^T * A * p - b^T * p + c  with p denoting the parameter vector to be optimized and
  the constants c (scalar), b (vector) and A (matrix). This function object is intended to be used 
  as a test error function for testing optimization algorithms.

  */

  class QuadraticTestErrorFunction : public MultivariateErrorFunction
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
    virtual double getValue(Vector p) { return  0.5*p*(A*p) - b*p + c; }

    /** Returns the local gradient at the given input vector p. */
    virtual Vector getGradient(Vector p) { return 0.5*(trans(A)+A)*p - b; }

    //=============================================================================================

  protected:

    Matrix A;
    Vector b;
    double c;

  };



  /**

  This class serves as baseclass for multivariate vector functions - that is, functions that take a
  vector as input and produce a vector at the output.

  */

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
    virtual Vector getOutputVector(Vector x) = 0;

    // virtual Matrix getJacobian(Vector x);

    //=============================================================================================

  protected:

    int numInputs, numOutputs;

  };

} // end namespace rosic

#endif 
