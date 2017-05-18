#ifndef rosic_ComplexFunctions_h
#define rosic_ComplexFunctions_h

// this file can be removed...

// rosic includes:
#include "rosic_Complex.h" 
#include "rosic_RealFunctions.h"

namespace rosic
{

  /**

  A collection of functions for complex arguments.

  */

  Complex acos_c(Complex z);   
  /**< Complex inverse cosine (arccosine) function. */

  Complex acosh_c(Complex z);   
  /**< Complex inverse hyperbolic cosine function. */

  Complex asin_c(Complex z);   
  /**< Complex inverse sine (arcsine) function. */

  Complex asinh_c(Complex z);   
  /**< Complex inverse hyperbolic sine function. */

  Complex atan_c(Complex z);   
  /**< Complex inverse tangent (arctangent) function. */

  Complex atanh_c(Complex z);   
  /**< Complex inverse hyperbolic tangent function. */

  Complex cos_c(Complex z);   
  /**< Complex cosine function. */

  Complex cosh_c(Complex z);   
  /**< Complex hyperbolic cosine function. */

  Complex exp_c(Complex z);   
  /**< Calculates the complex exponential of a complex number. */

  Complex log_c(Complex z);   
  /**< Calculates the natural (base e) logarithm of a complex number - as the complex 
  logarithm is a multifunction, it returns the main value. */

  Complex pow_c(Complex basis, Complex exponent);
  /**< Raises a complex number to a complex power. */

  Complex sin_c(Complex z);   
  /**< Complex sine function. */

  Complex sinh_c(Complex z);   
  /**< Complex hyperbolic sine function. */

  Complex sqrt_c(Complex z);  
  /**< Calculates the (primitive) square root of a complex number. The second square root 
  is obtained by using the negative value. */

  Complex tan_c(Complex z);   
  /**< Complex tangent function. */

  Complex tanh_c(Complex z);   
  /**< Complex hyperbolic tangent function. */

} // end namespace rosic

#endif // #ifndef rosic_ComplexFunctions_h