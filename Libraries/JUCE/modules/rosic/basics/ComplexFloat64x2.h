#ifndef RAPT_COMPLEXFLOAT64X2_H_INCLUDED
#define RAPT_COMPLEXFLOAT64X2_H_INCLUDED

/** In this file are some operattors and functions for complex variables where the real and
imaginary parts are each a SIMD vector. For some reason, the standard library functions of
std::complex dont work anymore when the template parameter to std::complex is a SIMD type, so we
provide explicit specializations here. */


/*
// copied/edited from complex on mac (i commented out most of the declarations - it may become
// necessarry to uncomment and provide implementations later...
template<>
class std::complex<rsFloat64x2>
{
public:
  //typedef rsFloat64x2 value_type;

  constexpr complex(rsFloat64x2 re = 0.0f, rsFloat64x2 im = 0.0f);
  //explicit constexpr complex(const complex<rsFloat64x2>&);

  constexpr rsFloat64x2 real() const;
  constexpr rsFloat64x2 imag() const;
  //void real(rsFloat64x2);
  //void imag(rsFloat64x2);

  complex<rsFloat64x2>& operator= (rsFloat64x2 z)
  {
    return *this = z;
  }

  //complex<rsFloat64x2>& operator+=(rsFloat64x2);
  //complex<rsFloat64x2>& operator-=(rsFloat64x2);
  //complex<rsFloat64x2>& operator*=(rsFloat64x2);
  complex<rsFloat64x2>& operator/=(rsFloat64x2);

  //complex<rsFloat64x2>& operator=(const complex<rsFloat64x2>&);
  //template<class X> complex<rsFloat64x2>& operator= (const complex<X>&);
  //template<class X> complex<rsFloat64x2>& operator+=(const complex<X>&);
  //template<class X> complex<rsFloat64x2>& operator-=(const complex<X>&);
  //template<class X> complex<rsFloat64x2>& operator*=(const complex<X>&);
  //template<class X> complex<rsFloat64x2>& operator/=(const complex<X>&);
};

// from https://en.cppreference.com/w/cpp/numeric/complex:
// For any object z of type complex<T>, 
// reinterpret_cast<T(&)[2]>(z)[0] is the real part of z and 
// reinterpret_cast<T(&)[2]>(z)[1] is the imaginary part of z.
// may be relevant for implementing real() and imag() in an explicit specialization

*/





/** Returns the first (index 0) complex number in the complex of simd vectors. */
inline std::complex<double> get0(std::complex<rsFloat64x2> z)
{
  double* re = z.real().asArray();
  double* im = z.imag().asArray();
  return std::complex<double>(re[0], im[0]);
}

/** Returns the second (index 1) complex number in the complex of simd vectors. */
inline std::complex<double> get1(std::complex<rsFloat64x2> z)
{
  double* re = z.real().asArray();
  double* im = z.imag().asArray();
  return std::complex<double>(re[1], im[1]);
}

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

#ifdef _MSC_VER // doesn't compile with gcc
/** Divides a complex number in place by another complex number. */
inline std::complex<rsFloat64x2>& std::complex<rsFloat64x2>::operator/=(
  const std::complex<rsFloat64x2>& a)
{
  return *this = *this / a;
}
#endif
// this does compile with gcc, but there it compiles also, when leaving it out entirely:
//inline std::complex<rsFloat64x2>& operator/=(std::complex<rsFloat64x2>& a,
//                                             const std::complex<rsFloat64x2>& b)
//{
//  return a = a / b; // THIS NEEDS TESTING!
//}


inline std::complex<rsFloat64x2> operator+(const std::complex<rsFloat64x2>& z) { return z; }



/** Computes the complex exponential of z. */
inline std::complex<rsFloat64x2> rsExp(std::complex<rsFloat64x2> z)
{
  // e^z = e^(a + i*b) = e^a * e^(i*b) = e^a * (cos(b) + i*sin(b))
  double* re = z.real().asArray();  // real parts
  double* im = z.imag().asArray();  // imag parts
  double r0  = exp(re[0]);          // radius of 1st complex result
  double r1  = exp(re[1]);          // radius of 2nd complex result
  double re0 = r0 * cos(im[0]);     // real part of 1st complex result
  double im0 = r0 * sin(im[0]);     // imag part of 1st complex result
  double re1 = r1 * cos(im[1]);     // real part of 2nd complex result
  double im1 = r1 * sin(im[1]);     // imag part of 2nd complex result
  rsFloat64x2 vre(re0, re1);        // vector of resulting real parts
  rsFloat64x2 vim(im0, im1);        // vector of resulting imag parts
  return std::complex<rsFloat64x2>(vre, vim);
}



/** Preliminary, to satisfy compiler on mac  */
#ifndef _MSC_VER
inline bool isnan(std::complex<rsFloat64x2> z)
{
  return false;
}
inline bool isinf(std::complex<rsFloat64x2> z)
{
  return false;
}
inline std::complex<rsFloat64x2> copysign(std::complex<rsFloat64x2> z, std::complex<rsFloat64x2> w)
{
  return w;
}

//template<>
//inline std::complex<rsFloat64x2>& std::complex<rsFloat64x2>::operator=(
//  std::complex<rsFloat64x2> z) { return z; }

//template<>
//inline std::complex<rsFloat64x2>& std::complex<rsFloat64x2>::operator=(
//  const std::complex<rsFloat64x2>& z) { return z; }

#endif



//inline std::complex<rsFloat64x2>& std::complex<rsFloat64x2>::operator=(std::complex<rsFloat64x2>& a,                                             const std::complex<rsFloat64x2>& b)
//{
//  return a = b;
//}

//inline void std::complex<rsFloat64x2>::operator=(const std::complex<rsFloat64x2>& a)
//{
//  return *this = a;
//}
//inline std::complex<rsFloat64x2>& std::complex<rsFloat64x2>::operator=(const std::complex<rsFloat64x2>& a)
//{
//  return *this = a;
//}






#endif
