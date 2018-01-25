#ifndef RAPT_FLOAT64X2_H_INCLUDED
#define RAPT_FLOAT64X2_H_INCLUDED

//=================================================================================================

/** A class for convenient handling of SIMD optimizations for a vector of two double-precision
floating point numbers. Requires the SSE2 instruction set. */

class rsFloat64x2
{
public:

  /** \name Construction */

  /** Standard constructor. Leaves both elements uninitialized. */
  inline rsFloat64x2() { /*v = _mm_setzero_pd();*/ }

  /** Constructor to copy an existing __m128d pair of values. */
  inline rsFloat64x2(const __m128d& rhs) : v(rhs) {}

  /** Constructor that initializes both elements to the given value. */
  inline rsFloat64x2(double a) : v(_mm_set1_pd(a)) {}

  /** Constructor that initializes the elements from two doubles. */
  inline rsFloat64x2(double a, double b) : v(_mm_setr_pd(a, b)) {}

  /** Constructor that initializes the elements from a 2-value array of doubles. */
  inline rsFloat64x2(double* values) 
  { 
    //*this = _mm_load_pd(values); // doesnt work
    v = _mm_setr_pd(values[0], values[1]);
  }


  /** \name Setup */

  /** Sets both elements to a. */
  inline void set(double a) { v = _mm_set1_pd(a); }
  // what's the difference between _mm_set1_pd and _mm_load1_pd? ...the latter takes a pointer?

  /** Sets the first element to a and the second element to b. */
  inline void set(double a, double b) { v = _mm_setr_pd(a, b); }


  /** \name Inquiry */

  // extract vector elements:
  inline double get0() { double d; _mm_storel_pd(&d, v); return d; }  // lower (index 0)
  inline double get1() { double d; _mm_storeh_pd(&d, v); return d; }  // upper (index 1)


  /** \name Operators */

  // arithmetic operators:
  inline rsFloat64x2& operator+=(const rsFloat64x2& b) { v = _mm_add_pd(v, b); return *this; }
  inline rsFloat64x2& operator-=(const rsFloat64x2& b) { v = _mm_sub_pd(v, b); return *this; }
  inline rsFloat64x2& operator*=(const rsFloat64x2& b) { v = _mm_mul_pd(v, b); return *this; }
  inline rsFloat64x2& operator/=(const rsFloat64x2& b) { v = _mm_div_pd(v, b); return *this; }

  //inline rsFloat64x2 operator+(const rsFloat64x2& b) { return _mm_add_pd(v, b.v); }
    // nope - we define it outside the class because that allows lhs or rhs to be of scalar type

  // unary minus:
  inline rsFloat64x2 operator-() 
  { 
    static const __m128d zero = _mm_setzero_pd();
    return _mm_sub_pd(zero, v);
  }

  /** Assignment from __m128d. */
  inline rsFloat64x2& operator=(const __m128d& rhs) { v = rhs; return *this; }

  /** Conversion to __m128d. */
  inline operator __m128d() const { return v; }


protected:

  __m128d v; // the value
  //__declspec(align(16)) __m128d v; // the value (define and ALIGN(N) macro for gcc/msc)

  // Note: Subclassing __m128d (instead of having a member of that type) works with the MS compiler 
  // but not with gcc. Apparently, MS allows primitive datatypes to be seen as classes while gcc 
  // doesn't.
};

// binary arithmetic operators:
inline rsFloat64x2 operator+(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_add_pd(a, b); }
inline rsFloat64x2 operator-(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_sub_pd(a, b); }
inline rsFloat64x2 operator*(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_mul_pd(a, b); }
inline rsFloat64x2 operator/(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_div_pd(a, b); }
// the binary operators with a scalar for the left or right hand side do not have to be defined due 
// to implicit conversions

// functions:
inline rsFloat64x2 rsMin(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_min_pd(a, b); }
inline rsFloat64x2 rsMax(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_max_pd(a, b); }
inline rsFloat64x2 rsSqrt(const rsFloat64x2& a) { return _mm_sqrt_pd(a); }

// implement clip (based on min/max), abs/sign (based on bit-masks), asArray

// todo: reduce_min, reduce_max, reduce_sum
// for more inspiration, see:
// https://msdn.microsoft.com/de-de/library/tyy88x2a(v=vs.90).aspx
// https://msdn.microsoft.com/de-de/library/9b07190d(v=vs.90).aspx
// http://johanmabille.github.io/blog/2014/10/10/writing-c-plus-plus-wrappers-for-simd-intrinsics-3/
// https://github.com/p12tic/libsimdpp 
// https://github.com/VcDevel/Vc
// https://github.com/NumScale/boost.simd

#endif