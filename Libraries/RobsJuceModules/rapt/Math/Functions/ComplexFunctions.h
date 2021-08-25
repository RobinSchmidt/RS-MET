#ifndef RAPT_COMPLEXFUNCTIONS_H_INCLUDED
#define RAPT_COMPLEXFUNCTIONS_H_INCLUDED

/** A collection of functions for complex arguments.

\todo: a few of the elliptic functions are still missing - implement them 
\todo: remove the functions that are available in the standard library (epx,log,pow,sin,etc.)
  
*/


/** Inverse complex Jacobian elliptic function cd with elliptic modulus k. */
template<class T>
std::complex<T> rsAcdC(std::complex<T> w, T k);

/** Complex inverse cosine (arccosine) function. */
template<class T>
std::complex<T> rsAcosC(std::complex<T> z);

/** Complex inverse hyperbolic cosine function. */
template<class T>
std::complex<T> rsAcoshC(std::complex<T> z);

/** Complex inverse sine (arcsine) function. */
template<class T>
std::complex<T> rsAsinC(std::complex<T> z);

/** Complex inverse hyperbolic sine function. */
template<class T>
std::complex<T> rsAsinhC(std::complex<T> z);

/** Inverse complex Jacobian elliptic function sn with elliptic modulus k. */
template<class T>
std::complex<T> rsAsnC(std::complex<T> w, T k);

/** Complex inverse tangent (arctangent) function. */
template<class T>
std::complex<T> rsAtanC(std::complex<T> z);

/** Complex inverse hyperbolic tangent function. */
template<class T>
std::complex<T> rsAtanhC(std::complex<T> z);

/** Complex Jacobian elliptic function cd with elliptic modulus k. */
template<class T>
std::complex<T> rsCdC(std::complex<T> u, T k);

/** Complex cosine function. */
template<class T>
std::complex<T> rsCosC(std::complex<T> z);

/** Complex hyperbolic cosine function. */
template<class T>
std::complex<T> rsCoshC(std::complex<T> z);

/** Calculates the complex exponential of a complex number. */
//template<class T>
//std::complex<T> rsExpC(std::complex<T> z);

/** Calculates the natural (base e) logarithm of a complex number - as the complex logarithm
is a multifunction, it returns the principal value. */
//template<class T>
//std::complex<T> rsLogC(std::complex<T> z);
// not needed anymore (is part of the standard library)

/** Raises a complex number to a complex power. */
template<class T>
std::complex<T> rsPowC(std::complex<T> basis, std::complex<T> exponent);

/** Complex sine function. */
template<class T>
std::complex<T> rsSinC(std::complex<T> z);

/** Complex hyperbolic sine function. */
template<class T>
std::complex<T> rsSinhC(std::complex<T> z);

/** std::complex<T> Jacobian elliptic function sn with elliptic modulus k. */
template<class T>
std::complex<T> rsSnC(std::complex<T> u, T k);

/** Calculates the (primitive) square root of a complex number. The second square root is
obtained by using the negative value. */
//template<class T>
//std::complex<T> rsSqrtC(std::complex<T> z);  

/** Complex tangent function. */
template<class T>
std::complex<T> rsTanC(std::complex<T> z);

/** Complex hyperbolic tangent function. */
template<class T>
std::complex<T> rsTanhC(std::complex<T> z);

/** Returns true if real or imaginary part (or both) are plus or minus infinity, false 
otherwise. */
template<class T>
inline bool rsIsInfinite(std::complex<T> z)
{
  if(  z.real() == RS_INF(T) || z.real() == -RS_INF(T) 
    || z.imag() == RS_INF(T) || z.imag() == -RS_INF(T) )
    return true;
  else
    return false;
}

/** Returns the number of finite values in the complex array "a" of length "N". */
template<class T>
int rsGetNumFiniteValues(std::complex<T> *a, int N);

/** Given an array "z" of "N" complex numbers, this function copies those values from "z" to "zF"
that are finite. The return value returns the number of copied values, i.e. the effective (used)
length of the returned array "zF". The remaining values in "zF" are left as is. */
template<class T>
int rsCopyFiniteValues(const std::complex<T> *z, std::complex<T> *zF, int N);

/** Returns the product of finite values in the complex array "a" of length "N". */
template<class T>
std::complex<T> rsProductOfFiniteFactors(std::complex<T> *a, int N);

/** Given an array "z" of "N" complex numbers, this function copies those values from "z" to "zL"
that are in the left half-plane (including the imaginary axis). The return value returns the
number of copied values, i.e. the effective (used) length of the returned array "zL". The
remaining values in "zL" are left as is. */
template<class T>
int rsOnlyLeftHalfPlane(std::complex<T> *z, std::complex<T> *zL, int N);

/** Similar to onlyLeftHalfPlane - copies upper halfplane values (including real axis). */
template<class T>
int rsOnlyUpperHalfPlane(std::complex<T> *z, std::complex<T> *zU, int N);

/** Zeros the imaginary parts in the passed complex array "z" when the absolute value of the
imaginary part is below the given threshold. */
template<class T>
void rsZeroNegligibleImaginaryParts(std::complex<T> *z, int length, T threshold);

/** Applies complex conjugation to all values in the buffer. */
template<class T>
void rsConjugate(std::complex<T> *z, int length);

/** Compares two complex numbers for a less-than condition by first comparing real parts and if
they are equal, comparing imaginary parts. */
template<class T>
bool rsComplexLessByReIm(const std::complex<T>& left, const std::complex<T>& right);

/** Like rsComplexLessByReIm but comparing imaginary parts first */
template<class T>
bool rsComplexLessByImRe(const std::complex<T>& left, const std::complex<T>& right);

/** Sorts an array of complex numbers according to the less-than criterion defined by the
function complexLessByReIm. */
template<class T>
void rsSortComplexArrayByReIm(std::complex<T> *z, int length);

/** Returns true, when for all values in the passed z-array:
|z[i].im| <= |relativeTolerance*z[i].re|. That means that the absolute value of the imaginary part
must be a factor "relativeTolerance" smaller than the real part for a number to be considered
purely real. */
template<class T>
bool rsAreAllValuesReal(std::complex<T> *z, int length, T relativeTolerance = T(0));

/** Returns true, if z[i] and z[i+1] are complex conjugates for all i from 0..length-2. */
template<class T>
bool rsAreNeighborsConjugates(std::complex<T> *z, int length, T absoluteTolerance = T(0));


/* Computes the 2 roots (zeros) of the quadratic equation x^2 + p*x + q = 0. */
template<class T>
inline void solveQuadraticEquation(std::complex<T> p, std::complex<T> q,
  std::complex<T> &root1, std::complex<T> &root2)
{
  std::complex<T> tmp = sqrt(0.25*p*p - q);
  root1 = -0.5*p + tmp;
  root2 = -0.5*p - tmp;
}

#endif
