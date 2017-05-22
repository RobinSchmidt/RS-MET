//#include <limits>
//
//#include "rosic_FunctionObjects.h"
//#include "../math/rosic_PolynomialAlgorithms.h"
//using namespace rosic;

//=================================================================================================
// class UnivariateScalarFunction  

//-------------------------------------------------------------------------------------------------    
// evaluation:

double UnivariateScalarFunction::getFirstDerivativeAt(double x)
{
  return approximateDerivativeAt(x, 1, rmax(0.001*x, 1.e-14));
}

double UnivariateScalarFunction::getSecondDerivativeAt(double x)
{
  double eps = rmax(0.001*x, 1.e-14);
  return (getFirstDerivativeAt(x+eps) - getFirstDerivativeAt(x-eps)) / (2*eps);
}

double UnivariateScalarFunction::approximateDerivativeAt(double x, int order, double eps)
{
  if( order == 0 )
    return getValueAt(x);
  else
  {
    return (  approximateDerivativeAt(x+eps, order-1, eps) 
            - approximateDerivativeAt(x-eps, order-1, eps)) / (2*eps);
  }
}

double UnivariateScalarFunction::findRootViaRidders(double xL, double xU)
{
  const double eps  = std::numeric_limits<double>::epsilon(); 
  const int    iMax = 1000;    // maximum number of iterations
  double fL = getValueAt(xL);  // function value at lower bound
  double fU = getValueAt(xU);  // function value at upper bound
  double fM;                   // function value at midpoint
  double xM;                   // midpoint of interval  (initialization unnecessary)
  double xN;                   // new x-value (to replace xL or xH)
  double fN;                   // new function value (to replace fL or fU)
  double sq;                   // square root in the denominator
  int    i  = 0;               // iteration counter

  bool upward = false;         // flag to distinguish between upward and downward zero crossings
  if( fU > fL )
    upward = true;

  if( fL*fU > 0.0 )
  {
    DEBUG_BREAK;  // brackets are not guaranteed to bracket a root
    return 0.0;
  }
  else if( fL*fU < 0.0 )
  {
    while( xU-xL > eps && i <= iMax )
    {
      xM = 0.5*(xL+xU); 
      fM = getValueAt(xM);                     // 1st evaluation
      if( fM == 0.0 )
        return xM;

      sq = sqrt(fM*fM-fL*fU);                  // square root in Eq. 9.2.4
      if( sq == 0.0 )
        return xM;

      xN = xM + (xM-xL) * sign(fL-fU)*fM / sq; // Eq. 9.2.4
      fN = getValueAt(xN);                     // 2nd evaluation
      if( fN == 0.0 )
        return xN;

      if( upward )
      {
        if( fN > 0.0 )
        {
          xU = xN; 
          fU = fN;
          if( fM < 0.0 ) 
          { 
            xL = xM; 
            fL = fM; 
          }
        }
        else if( fN < 0.0 )
        {
          xL = xN;
          fL = fN;
          if( fM > 0.0 )
          {
            xU = xM;
            fU = fM;
          }
        }
      }
      else  // downward zero crossing
      {
        if( fN < 0.0 )
        {
          xU = xN;
          fU = fN;
          if( fM > 0.0 )
          {
            xL = xM;
            fL = fM;
          }
        }
        else if( fN > 0.0 )
        {
          xL = xN;
          fL = fN;
          if( fM < 0.0 )
          {
            xU = xM;
            fU = fM;
          }
        }
      }

      i++;
    }
  }
  else // fL*fU == 0.0
  {
    if( fL == 0.0 )
      return xL;
    if( fU == 0.0 )
      return xU;
  }

  DEBUG_BREAK;  // did not converge - should not happen
  return 0.0; 
}

double UnivariateScalarFunction::findRootViaNewtonNonRobust(double x)
{
  const double eps  = std::numeric_limits<double>::epsilon(); 
  const int    iMax = 1000;    // maximum number of iterations
  double f, f1;                // f(x), f'(x)
  double d = 2*eps;            // delta for next iteration, initializer only to enter the loop
  int    i = 0;                // iteration counter
  while( fabs(d) > eps && i <= iMax )
  {
    f   = getValueAt(x);
    f1  = getFirstDerivativeAt(x);
    d   = -f/f1;  // what if f1 == 0?
    x  += d;
    i++;
  }
  return x;
}

double UnivariateScalarFunction::findRootViaChebychevNonRobust(double x)
{
  const double eps = std::numeric_limits<double>::epsilon(); 
  double f, f1, f2;  // f(x), f'(x), f''(x)
  double d = 2*eps;  
  int i    = 0;    
  int iMax = 1000; 
  while( fabs(d) > eps && i <= iMax )
  {
    f   = getValueAt(x);
    f1  = getFirstDerivativeAt(x);
    f2  = getSecondDerivativeAt(x);
    d   = -(f/f1 + (f2*f*f)/(2*f1*f1*f1));
    x  += d;
    i++;
  }
  return x;
}



//=================================================================================================
// class UnivariateScalarFunctionViaPointer

//-------------------------------------------------------------------------------------------------    
// construction/destruction:


UnivariateScalarFunctionViaPointer::UnivariateScalarFunctionViaPointer(double (*functionToUse)(double), 
                                                                       double (*derivativeToUse)(double))
{
  rassert( functionToUse != NULL );
  functionPointer   = functionToUse;
  derivativePointer = derivativeToUse;
}

//-------------------------------------------------------------------------------------------------    
// evaluation:

double UnivariateScalarFunctionViaPointer::getValueAt(double x)
{
  return functionPointer(x);
}

double UnivariateScalarFunctionViaPointer::getFirstDerivativeAt(double x)
{
  if( derivativePointer != NULL )
    return derivativePointer(x);
  else
    return UnivariateScalarFunction::getFirstDerivativeAt(x);
}


//=================================================================================================
// class Polynomial  

//-------------------------------------------------------------------------------------------------    
// construction/destruction:

Polynomial::Polynomial(const unsigned int initOrder, const double* const initCoeffs)
{
  order  = initOrder;
  coeffs = new double[order+1];
  if( initCoeffs != NULL )
    setCoefficients(initCoeffs);
  else
    zeroCoefficients();
}

Polynomial::~Polynomial()
{
  delete[] coeffs;
}

//-------------------------------------------------------------------------------------------------    
// setup:

void Polynomial::zeroCoefficients()
{
  for(unsigned int k = 0; k <= order; k++)
    coeffs[k] = 0.0;
}

void Polynomial::setCoefficients(const double* const newCoeffs)
{
  for(unsigned int k = 0; k <= order; k++)
    coeffs[k] = newCoeffs[k];
}

//-------------------------------------------------------------------------------------------------    
// inquiry:

double Polynomial::getValueAt(double x)
{
  return rosic::evaluatePolynomialAt(x, coeffs, order);
}

double Polynomial::getFirstDerivativeAt(double x)
{
  double f[2];  // f(x), f'(x)
  rosic::evaluatePolynomialAndDerivativesAt(x, coeffs, order, f, 1);
  return f[1];
}

double Polynomial::getSecondDerivativeAt(double x)
{
  double f[3];  // f(x), f'(x), f''(x)
  rosic::evaluatePolynomialAndDerivativesAt(x, coeffs, order, f, 2);
  return f[2];
}



//=================================================================================================
// class MultivariateScalarFunction  

//-------------------------------------------------------------------------------------------------    
// construction/destruction:

MultivariateScalarFunction::MultivariateScalarFunction(int numInputs)
{
  this->numInputs = numInputs;
}

MultivariateScalarFunction::~MultivariateScalarFunction()
{

}


//=================================================================================================
// class MultivariateErrorFunction  

//-------------------------------------------------------------------------------------------------    
// construction/destruction:

MultivariateErrorFunction::MultivariateErrorFunction()
{

}

MultivariateErrorFunction::~MultivariateErrorFunction()
{

}

rosic::Vector MultivariateErrorFunction::getGradient(rosic::Vector p)
{
  rosic::Vector g(p.dim);   // gradient vector to be computed
  double eps = 0.00001;     // epsilon for the approximation
  double pTmp;              // for temporary storage
  double ep;                // error at p[i] + eps
  double em;                // error at p[i] - eps
  for(int i=0; i<p.dim; i++)
  {
    pTmp   = p.v[i];
    p.v[i] = pTmp + eps;
    ep     = getValue(p);
    p.v[i] = pTmp - eps;
    em     = getValue(p);
    g.v[i] = (ep-em) / (2.0*eps);
    p.v[i] = pTmp;
  }
  return g;
}

rosic::Vector MultivariateErrorFunction::getVectorTimesHessianApproximate(rosic::Vector p, 
  rosic::Vector v)
{
  double eps  = 0.00001;                // epsilon for the approximation...
  eps        /= v.getEuclideanNorm();   // ...should be normalized by ||v||
  rosic::Vector gp   = getGradient(p + eps*v); // gradient at p + eps*v
  rosic::Vector gm   = getGradient(p - eps*v); // gradient at p - eps*v
  return (gp-gm) / (2.0*eps);
}

//=================================================================================================
// class QuadraticTestErrorFunction  

//-------------------------------------------------------------------------------------------------    
// construction/destruction:

QuadraticTestErrorFunction::QuadraticTestErrorFunction()
{
  b.setDimensionality(2);
  b.v[0] = -5;
  b.v[1] = +1;

  A.setSize(2, 2);
  A.m[0][0] =  1;
  A.m[0][1] = -5;
  A.m[1][0] = 10;
  A.m[1][1] = 15;

  c = 3.0;
}

QuadraticTestErrorFunction::~QuadraticTestErrorFunction()
{

}

//=================================================================================================
// class MultivariateVectorFunction  

//-------------------------------------------------------------------------------------------------    
// construction/destruction:

MultivariateVectorFunction::MultivariateVectorFunction(int numInputs, int numOutputs)
{
  this->numInputs  = numInputs;
  this->numOutputs = numOutputs;
}