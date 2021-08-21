#ifndef RAPT_COMPLEX_H
#define RAPT_COMPLEX_H

/** Class for representing complex numbers. The primary reason for this class to exist is that
std::complex only admits float, double and long double for the underlying real type, but we 
sometimes need other types of complex numbers, too. Especially important are SIMD types, which is
the reason why we avoid branches in the implementation, because they do not tend to play nicely 
with SIMD. Avoiding branches may at times scarifice numerical accuracy (todo: elaborate). */

template<class T>
class rsComplex
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Standard constructor. Leaves real and imaginary part uninitialized. */
  rsComplex() {}

  /** Constructor. Initializes real part with given argument and imaginary part to zero via their 
  parametrized constructors. */
  rsComplex(T re_) : re(re_), im(T(0)) {}

  /** Constructor. Initializes real and imaginary part with given arguments via their parametrized 
  constructors. */
  rsComplex(T re_, T im_) : re(re_), im(im_) {}

  /** Constructor for implicit type-conversion of the complex number with T2 into a complex number 
  with type T for its components. */
  template<class T2>
  rsComplex(T2 re_, T2 im_) : re(T(re_)), im(T(im_)) {}


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  bool operator==(const rsComplex& z) const { return re == z.re && im == z.im; }
  bool operator!=(const rsComplex& z) const { return !(re == im); }

  rsComplex operator+() const { return rsComplex(+re, +im); }
  rsComplex operator-() const { return rsComplex(-re, -im); }

  rsComplex& operator+=(const rsComplex &z) { re += z.re; im += z.im; return *this; }
  rsComplex& operator-=(const rsComplex &z) { re -= z.re; im -= z.im; return *this; }


  rsComplex& operator*=(const rsComplex &z) 
  { 
    T tmp = re*z.re - im*z.im; 
    im    = re*z.im + im*z.re; 
    re = tmp; return *this; 
  }

  rsComplex& operator/=(const rsComplex &z)
  {
    T s = T(1) / (z.re*z.re + z.im*z.im);
    T r = s * (re*z.re + im*z.im);
    T i = s * (im*z.re - re*z.im);
    re = r; im = i; return *this;
  }
  // needs test

  rsComplex& operator+=(const T &r) { re += r; return *this; }
  rsComplex& operator-=(const T &r) { re -= r; return *this; }
  rsComplex& operator*=(const T &r) { re *= r; im *= r; return *this; }
  rsComplex& operator/=(const T &r) { T s = T(1)/r; re *= s; im *= s; return *this; }

  //-----------------------------------------------------------------------------------------------
  /** \name Access */

  T real() const { return re; }
  T imag() const { return im; }

  void real(const T& val) { re = val; }
  void imag(const T& val) { im = val; }



//protected:

  T re, im; // Don't change order or add further fields! Some algorithms rely on the data layout.

};

template<class T> 
inline rsComplex<T> rsInverse(const rsComplex<T> &z) 
{ 
  T s = T(1) / (z.re*z.re + z.im*z.im); 
  return rsComplex<T>(s*z.re, -s*z.im);
}
// maybe rename to rsInv, specify that rsInv should, by convention, always mean the left inverse for 
// types for which it matters

template<class T>
inline rsComplex<T> operator+(const rsComplex<T> &z, const rsComplex<T> &w)
{
  return rsComplex<T>(z.re+w.re, z.im+w.im);
}

template<class T>
inline rsComplex<T> operator+(const T &r, const rsComplex<T> &z)
{
  return rsComplex<T>(z.re+r, z.im);
}

template<class T>
inline rsComplex<T> operator+(const rsComplex<T> &z, const T &r)
{
  return rsComplex<T>(z.re+r, z.im);
}

template<class T>
inline rsComplex<T> operator-(const rsComplex<T> &z, const rsComplex<T> &w)
{
  return rsComplex<T>(z.re-w.re, z.im-w.im);
}

template<class T>
inline rsComplex<T> operator-(const rsComplex<T> &z, const T &r)
{
  return rsComplex<T>(z.re-r, z.im);
}

template<class T>
inline rsComplex<T> operator-(const T &r, const rsComplex<T> &z)
{
  return rsComplex<T>(r-z.re, -z.im);
}

template<class T>
inline rsComplex<T> operator*(const rsComplex<T> &z, const rsComplex<T> &w)
{
  return rsComplex<T>(z.re*w.re-z.im*w.im, z.re*w.im+z.im*w.re);
}

template<class T>
inline rsComplex<T> operator*(const rsComplex<T> &z, const T &r)
{
  return rsComplex<T>(z.re*r, z.im*r);
}

template<class T>
inline rsComplex<T> operator*(const T &r, const rsComplex<T> &z)
{
  return rsComplex<T>(z.re*r, z.im*r);
}

template<class T>
inline rsComplex<T> operator/(const rsComplex<T> &z, const rsComplex<T> &w)
{
  // Division causes problems with nested complex numbers: 1 nesting level works, but 2 don't. But
  // we don't care about dividing doubly nested complex numbers anyway (do they even form a 
  // division algebra? howeve, see RSLib for a fix that sacrifices performance):
  T s = T(1) / (w.re*w.re + w.im*w.im);
  return rsComplex<T>(s*(z.re*w.re + z.im*w.im), s*(z.im*w.re - z.re*w.im));
}

template<class T>
inline rsComplex<T> operator/(const rsComplex<T> &z, const T &r)
{
  T s = T(1) / r;
  return rsComplex<T>(s*z.re, s*z.im);
}

template<class T>
inline rsComplex<T> operator/(const T &r, const rsComplex<T> &z)
{
  T s = r / (z.re*z.re + z.im*z.im);
  return rsComplex<T>(s*z.re, -s*z.im);
}

//=================================================================================================
// Elementary math functions for complex numbers:

template<class T>
//rsComplex<T> rsAbs(rsComplex<T> z)
T rsAbs(rsComplex<T> z)
{
  return rsSqrt(z.re*z.re + z.im*z.im); 
  // try hypot, measure performance...but first make sure that hypot doesn't branch, if it does,
  // maybe try to implement a branchless rsHypot
}

template<class T>
//rsComplex<T> rsArg(rsComplex<T> z)
T rsArg(rsComplex<T> z)
{
  return rsAtan2(z.im, z.re);
}

template<class T>
rsComplex<T> rsExp(rsComplex<T> z)
{
  rsComplex<T> w;
  T r = rsExp(z.re);             // compute radius of w
  rsSinCos(z.im, &w.im, &w.re);  // set angle of w
  w.re *= r; w.im *= r;          // set radius of w
  return w;
}

template<class T>
rsComplex<T> rsLog(rsComplex<T> z)
{
  rsComplex<T> tmp;
  tmp.re = rsLog(rsAbs(z));
  tmp.im = rsArg(z);
  return tmp;
}

template<class T>
rsComplex<T> rsPow(rsComplex<T> basis, rsComplex<T> exponent)
{
  return rsExp(exponent * rsLog((rsComplex<T>)basis));
}


// ToDo: 
// -Do tests comparing it with T = float or double to std::complex with regard to numerical 
//  accuracy, handling of edge cases (inf, nan, etc.) and performance. Also instantiate it for int,
//  rsFraction, rsSimdVector and check that it works. Try it also with a type that uses heap memory
//  like rsMatrix.
// -Implement rsConj, rsSin, rsCos, etc. and also a couple of special functions (elliptic, etc.). 
//  But that requires them to be implemented in a branchless way.
// -Switch rsEngineersFilter back to using rsComplex to facilitate vectorization...basically, 
//  switch all uses of std::complex to rsComplex - while keeping an eye on the impact on numerical
//  accuracy.
// -Maybe allow implicit conversion to/from std::complex for convenience
// -Maybe let the elementary functions take their arguments by const reference


#endif
