#ifndef rosic_ComplexFunctions_h
#define rosic_ComplexFunctions_h

// rosic includes:
#include "../basics/rosic_HelperFunctions.h"
#include "rosic_Complex.h" 
#include "rosic_SpecialFunctionsReal.h" // todo - move this to a new file SpecialFunctionsComplex


namespace rosic
{

  /**

  A collection of functions for complex arguments.

  \todo: a few of the elliptic functions are still missing - implement them

  */

  /** Inverse complex Jacobian elliptic function cd with elliptic modulus k. */
  Complex acdC(Complex w, double k);   

  /** Complex inverse cosine (arccosine) function. */
  Complex acosC(Complex z);   

  /** Complex inverse hyperbolic cosine function. */
  Complex acoshC(Complex z);   

  /** Complex inverse sine (arcsine) function. */
  Complex asinC(Complex z);   

  /** Complex inverse hyperbolic sine function. */
  Complex asinhC(Complex z);   

  /** Inverse complex Jacobian elliptic function sn with elliptic modulus k. */
  Complex asnC(Complex w, double k);   

  /** Complex inverse tangent (arctangent) function. */
  Complex atanC(Complex z);   

  /** Complex inverse hyperbolic tangent function. */
  Complex atanhC(Complex z);   

  /** Complex Jacobian elliptic function cd with elliptic modulus k. */
  Complex cdC(Complex u, double k);   

  /** Complex cosine function. */
  Complex cosC(Complex z);   

  /** Complex hyperbolic cosine function. */
  Complex coshC(Complex z);   

  /** Calculates the complex exponential of a complex number. */
  Complex expC(Complex z);    

  /** Calculates the natural (base e) logarithm of a complex number - as the complex logarithm 
  is a multifunction, it returns the principal value. */
  Complex logC(Complex z);   

  /** Raises a complex number to a complex power. */
  Complex powC(Complex basis, Complex exponent);

  /** Complex sine function. */
  Complex sinC(Complex z);   

  /** Complex hyperbolic sine function. */
  Complex sinhC(Complex z);   

  /** Complex Jacobian elliptic function sn with elliptic modulus k. */
  Complex snC(Complex u, double k);   

  /** Calculates the (primitive) square root of a complex number. The second square root is 
  obtained by using the negative value. */
  Complex sqrtC(Complex z);  

  /** Complex tangent function. */
  Complex tanC(Complex z);   

  /** Complex hyperbolic tangent function. */
  Complex tanhC(Complex z);   

  /** Returns the number of finite values in the complex array "a" of length "N". */
  int getNumFiniteValues(Complex *a, int N);

  /** Given an array "z" of "N" complex numbers, this function copies those values from "z" to "zF" 
  that are finite. The return value returns the number of copied values, i.e. the effective (used) 
  length of the returned array "zF". The remaining values in "zF" are left as is. */
  int copyFiniteValues(Complex *z, Complex *zF, int N);

  /** Returns the product of finite values in the complex array "a" of length "N". */
  Complex productOfFiniteFactors(Complex *a, int N);

  /** Given an array "z" of "N" complex numbers, this function copies those values from "z" to "zL" 
  that are in the left half-plane (including the imaginary axis). The return value returns the 
  number of copied values, i.e. the effective (used) length of the returned array "zL". The 
  remaining values in "zL" are left as is. */
  int onlyLeftHalfPlane(Complex *z, Complex *zL, int N);

  /** Similar to onlyLeftHalfPlane - copies upper halfplane values (including real axis). */
  int onlyUpperHalfPlane(Complex *z, Complex *zU, int N);

  /** Zeros the imaginary parts in the passed complex array "z" when the absolute value of the 
  imaginary part is below the given threshold. */
  void zeroNegligibleImaginaryParts(Complex *z, int length, double threshold);

  /** Compares two complex numbers for a less-than condition by first comparing real parts and if 
  they are equal, comparing imaginary parts. */
  bool complexLessByReIm(const Complex left, const Complex right);

  /** Sorts an array of complex numbers according to the less-than criterion defined by the 
  function complexLessByReIm. */
  void sortComplexArrayByReIm(Complex *z, int length);

}

#endif