#ifndef RAPT_COMPLEXFLOAT64X2_H_INCLUDED
#define RAPT_COMPLEXFLOAT64X2_H_INCLUDED

/** In this file are some operattors and functions for complex variables where the real and 
imaginary parts are each a SIMD vector. For some reason, the standard library functions of 
std::complex dont work anymore when the template parameter to std::complex is a SIMD type, so we 
provide explicit specializations here. */

/** Divides two complex numbers. */
inline std::complex<rsFloat64x2> operator/(
  const std::complex<rsFloat64x2>& a, const std::complex<rsFloat64x2>& b) 
{ 
  double* reA = a.real().asArray();
  double* imA = a.imag().asArray();
  double* reB = b.real().asArray();
  double* imB = b.imag().asArray();

  double  s0  = 1.0 / (reB[0]*reB[0] + imB[0]*imB[0]);
  double  re0 = s0  * (reA[0]*reB[0] + imA[0]*imB[0]);
  double  im0 = s0  * (imA[0]*reB[0] - reA[0]*imB[0]);

  double  s1  = 1.0 / (reB[1]*reB[1] + imB[1]*imB[1]);
  double  re1 = s1  * (reA[1]*reB[1] + imA[1]*imB[1]);
  double  im1 = s1  * (imA[1]*reB[1] - reA[1]*imB[1]);

  return std::complex<rsFloat64x2>(rsFloat64x2(re0, re1), rsFloat64x2(im0, im1)); 
}
// \todo: provide optimized versions when left or right operand is real

/** Divides a complex number in place by another complex number. */
inline std::complex<rsFloat64x2>& std::complex<rsFloat64x2>::operator/=(
  const std::complex<rsFloat64x2>& a) 
{ 
  *this = *this / a;
  return *this;
}

/** Computes the complex exponential of z. */
inline std::complex<rsFloat64x2> rsExp(std::complex<rsFloat64x2> z)
{
  // e^z = e^(a + i*b) = e^a * e^(i*b) = e^a * (cos(b) + i*sin(b))
  double* re = z.real().asArray();  // real parts
  double* im = z.imag().asArray();  // imag parts
  double r0  = exp(re[0]);          // radius of 1st complex number
  double r1  = exp(re[1]);          // radius of 2nd complex number
  double re0 = r0 * cos(im[0]);     // real part of 1st complex number
  double im0 = r0 * sin(im[0]);     // imag part of 1st complex number
  double re1 = r1 * cos(im[1]);     // real part of 2nd complex number
  double im1 = r1 * sin(im[1]);     // imag part of 2nd complex number
  rsFloat64x2 vre(re0, re1);        // vector of resulting real parts
  rsFloat64x2 vim(im0, im1);        // vector of resulting imag parts
  return std::complex<rsFloat64x2>(vre, vim);
}

inline std::complex<rsFloat64x2> operator+(const std::complex<rsFloat64x2>& z) { return z; }

/*
inline std::complex<rsFloat64x2> operator==(
  const std::complex<rsFloat64x2>& a, const std::complex<rsFloat64x2>& b) 
{ 
  double* are = a.real().asArray();
  double* aim = a.imag().asArray();
  double* bre = b.real().asArray();
  double* bim = b.imag().asArray();
  return (are[0] == bre[0]) && (aim[0] == bim[0]) && (are[1] == bre[1]) && (aim[1] == bim[1]);
}
*/

#endif