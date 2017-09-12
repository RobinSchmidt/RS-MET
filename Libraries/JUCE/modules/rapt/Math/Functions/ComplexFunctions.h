#ifndef RS_COMPLEXFUNCTIONS_H
#define RS_COMPLEXFUNCTIONS_H

namespace RSLib
{

  /**

  A collection of functions for complex arguments.

  \todo: a few of the elliptic functions are still missing - implement them

  */

  /** Inverse complex Jacobian elliptic function cd with elliptic modulus k. */
  template<class T>
  rsComplex<T> rsAcdC(rsComplex<T> w, double k);   

  /** Complex inverse cosine (arccosine) function. */
  template<class T>
  rsComplex<T> rsAcosC(rsComplex<T> z);   

  /** Complex inverse hyperbolic cosine function. */
  template<class T>
  rsComplex<T> rsAcoshC(rsComplex<T> z);   

  /** Complex inverse sine (arcsine) function. */
  template<class T>
  rsComplex<T> rsAsinC(rsComplex<T> z);   

  /** Complex inverse hyperbolic sine function. */
  template<class T>
  rsComplex<T> rsAsinhC(rsComplex<T> z);   

  /** Inverse complex Jacobian elliptic function sn with elliptic modulus k. */
  template<class T>
  rsComplex<T> rsAsnC(rsComplex<T> w, double k);   

  /** Complex inverse tangent (arctangent) function. */
  template<class T>
  rsComplex<T> rsAtanC(rsComplex<T> z);   

  /** Complex inverse hyperbolic tangent function. */
  template<class T>
  rsComplex<T> rsAtanhC(rsComplex<T> z);   

  /** Complex Jacobian elliptic function cd with elliptic modulus k. */
  template<class T>
  rsComplex<T> rsCdC(rsComplex<T> u, double k);   

  /** Complex cosine function. */
  template<class T>
  rsComplex<T> rsCosC(rsComplex<T> z);   

  /** Complex hyperbolic cosine function. */
  template<class T>
  rsComplex<T> rsCoshC(rsComplex<T> z);   

  /** Calculates the complex exponential of a complex number. */
  template<class T>
  rsComplex<T> rsExpC(rsComplex<T> z);    

  /** Calculates the natural (base e) logarithm of a complex number - as the complex logarithm 
  is a multifunction, it returns the principal value. */
  template<class T>
  rsComplex<T> rsLogC(rsComplex<T> z);   

  /** Raises a complex number to a complex power. */
  template<class T>
  rsComplex<T> rsPowC(rsComplex<T> basis, rsComplex<T> exponent);

  /** Complex sine function. */
  template<class T>
  rsComplex<T> rsSinC(rsComplex<T> z);   

  /** Complex hyperbolic sine function. */
  template<class T>
  rsComplex<T> rsSinhC(rsComplex<T> z);   

  /** rsComplex<T> Jacobian elliptic function sn with elliptic modulus k. */
  template<class T>
  rsComplex<T> rsSnC(rsComplex<T> u, double k);   

  /** Calculates the (primitive) square root of a complex number. The second square root is 
  obtained by using the negative value. */
  //template<class T>
  //rsComplex<T> rsSqrtC(rsComplex<T> z);  

  /** Complex tangent function. */
  template<class T>
  rsComplex<T> rsTanC(rsComplex<T> z);   

  /** Complex hyperbolic tangent function. */
  template<class T>
  rsComplex<T> rsTanhC(rsComplex<T> z);   

  /** Returns the number of finite values in the complex array "a" of length "N". */
  template<class T>
  int rsGetNumFiniteValues(rsComplex<T> *a, int N);

  /** Given an array "z" of "N" complex numbers, this function copies those values from "z" to "zF" 
  that are finite. The return value returns the number of copied values, i.e. the effective (used) 
  length of the returned array "zF". The remaining values in "zF" are left as is. */
  template<class T>
  int rsCopyFiniteValues(rsComplex<T> *z, rsComplex<T> *zF, int N);

  /** Returns the product of finite values in the complex array "a" of length "N". */
  template<class T>
  rsComplex<T> rsProductOfFiniteFactors(rsComplex<T> *a, int N);

  /** Given an array "z" of "N" complex numbers, this function copies those values from "z" to "zL" 
  that are in the left half-plane (including the imaginary axis). The return value returns the 
  number of copied values, i.e. the effective (used) length of the returned array "zL". The 
  remaining values in "zL" are left as is. */
  template<class T>
  int rsOnlyLeftHalfPlane(rsComplex<T> *z, rsComplex<T> *zL, int N);

  /** Similar to onlyLeftHalfPlane - copies upper halfplane values (including real axis). */
  template<class T>
  int rsOnlyUpperHalfPlane(rsComplex<T> *z, rsComplex<T> *zU, int N);

  /** Zeros the imaginary parts in the passed complex array "z" when the absolute value of the 
  imaginary part is below the given threshold. */
  template<class T>
  void rsZeroNegligibleImaginaryParts(rsComplex<T> *z, int length, double threshold);

  /** Applies complex conjugation to all values in the buffer. */
  template<class T>
  void rsConjugate(rsComplex<T> *z, int length);

  /** Compares two complex numbers for a less-than condition by first comparing real parts and if 
  they are equal, comparing imaginary parts. */
  template<class T>
  bool rsComplexLessByReIm(const rsComplex<T> left, const rsComplex<T> right);

  /** Sorts an array of complex numbers according to the less-than criterion defined by the 
  function complexLessByReIm. */
  template<class T>
  void rsSortComplexArrayByReIm(rsComplex<T> *z, int length);

  /** Returns true, when for all values in the passed z-array: 
  |z[i].im| <= |relativeTolerance*z[i].re|. That means that the absolute value of the imaginary part
  must be a factor "relativeTolerance" smaller than the real part for a number to be considered 
  purely real. */
  template<class T>
  bool rsAreAllValuesReal(rsComplex<T> *z, int length, T relativeTolerance = T(0));

}

#endif