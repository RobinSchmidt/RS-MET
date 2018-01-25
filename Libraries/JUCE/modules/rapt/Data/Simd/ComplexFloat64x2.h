#ifndef RAPT_COMPLEXFLOAT64X2_H_INCLUDED
#define RAPT_COMPLEXFLOAT64X2_H_INCLUDED

/** In this file are functions that implement some math functions for complex variables where the
real and imaginary parts are each a SIMD vector. For some reason, the standard library functions of
std::complex dont work anymore when the template parameter to std::complex is a SIMD type, so we 
provide explicit specializations here. */

/** Computes the complex exponential of z. */
inline std::complex<rsFloat64x2> exp(std::complex<rsFloat64x2> z)
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

#endif