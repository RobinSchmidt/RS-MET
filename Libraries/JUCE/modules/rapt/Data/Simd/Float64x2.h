#ifndef RAPT_FLOATSIMDVECTORS_H_INCLUDED
#define RAPT_FLOATSIMDVECTORS_H_INCLUDED

//=================================================================================================

/** A class for convenient handling of SIMD optimizations for a vector of two double-precision
floating point numbers. */

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

// todo: reduce_min, reduce_max, reduce_sum
// for more inspiration, see:
// https://msdn.microsoft.com/de-de/library/tyy88x2a(v=vs.90).aspx
// https://msdn.microsoft.com/de-de/library/9b07190d(v=vs.90).aspx
// http://johanmabille.github.io/blog/2014/10/10/writing-c-plus-plus-wrappers-for-simd-intrinsics-3/
// https://github.com/p12tic/libsimdpp 
// https://github.com/VcDevel/Vc
// https://github.com/NumScale/boost.simd







//=================================================================================================

/** This is datatype to represent 4 32-bit floating point numbers at once.
\todo: 
-write some unit tests
-write some performance tests
-use the 128-bit SSE type internally to speed up the computations in the arithmetic operators
 -currently, this is only a prototype use 4 actual floats
-compare performance of SSE implementation to non SSE
*/

class rsFloat32x4
{

public:

  /** \name Construction */

  /** Constructor. Sets up the 4 elements to the given values. */
  rsFloat32x4(float v0 = 0.f, float v1 = 0.f, float v2 = 0.f, float v3 = 0.f)
  {
    setValues(v0, v1, v2, v3);
  }

  /** Constructor. Sets up the all 4 elements to the same give value. */
  rsFloat32x4(float value)
  {
    setValues(value, value, value, value);
  }

  /** Constructor for conversion from double to float and setting up all 4 elements to the same 
  value. */
  rsFloat32x4(double value)
  {
    float f = (float)value;
    setValues(f, f, f, f);
  }

  /** Conversion constructor for int. */
  rsFloat32x4(int value)
  {
    float f = (float)value;
    setValues(f, f, f, f);
  }

  /** \name Setup */

  /** Sets up all 4 values. */
  inline void setValues(float v0, float v1, float v2, float v3)
  {
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
  }

  //inline void setValues(float newValues[4]) { }


  /** \name Inquiry */

  /** Returns one of the 4 values. */
  inline float getValue(int index)
  {
    return v[index];
  }


  /** \name Operators */

  /** Unary minus */
  inline rsFloat32x4 operator-() const
  { 
    return rsFloat32x4(-v[0], -v[1], -v[2], -v[3]); 
  }

  /** Addition of 2 Float32x4 vectors. */
  inline rsFloat32x4 operator+(const rsFloat32x4& y)
  {
    return rsFloat32x4(v[0]+y.v[0], v[1]+y.v[1], v[2]+y.v[2], v[3]+y.v[3]);
  }

  /** Subtraction of 2 Float32x4 vectors. */
  inline rsFloat32x4 operator-(const rsFloat32x4& y)
  {
    return rsFloat32x4(v[0]-y.v[0], v[1]-y.v[1], v[2]-y.v[2], v[3]-y.v[3]);
  }

  /** Multiplication of 2 Float32x4 vectors. */
  inline rsFloat32x4 operator*(const rsFloat32x4& y)
  {
    return rsFloat32x4(v[0]*y.v[0], v[1]*y.v[1], v[2]*y.v[2], v[3]*y.v[3]);
  }

  /** Division of 2 Float32x4 vectors. */
  inline rsFloat32x4 operator/(const rsFloat32x4& y)
  {
    return rsFloat32x4(v[0]/y.v[0], v[1]/y.v[1], v[2]/y.v[2], v[3]/y.v[3]);
  }

  /** In place multiplication of this Float32x4 vector with another. */
  inline rsFloat32x4 operator*=(const rsFloat32x4& y)
  {
    v[0] *= y.v[0];
    v[1] *= y.v[1];
    v[2] *= y.v[2];
    v[3] *= y.v[3];
    return *this;
  }

  /** Less-or-equal comparison. Returns true, if all 4 values are less or equal. */
  inline bool operator<=(const rsFloat32x4& y)
  {
    return v[0] <= y.v[0] && v[1] <= y.v[1] && v[2] <= y.v[2] && v[3] <= y.v[3];
  }

  /** Greater-or-equal comparison. Returns true, if all 4 values are greater or equal. */
  inline bool operator>=(const rsFloat32x4& y)
  {
    return v[0] >= y.v[0] && v[1] >= y.v[1] && v[2] >= y.v[2] && v[3] >= y.v[3];
  }

  // the set of operators is still incomplete, we need
  // +, -, *, /;  +=, -=, *=, /=;  ==, !=, >=, <=, >, <; unary -
  // for the binary operators, we need also versions where either left or right operator can be a 
  // single float


protected:

  float v[4];  // vector of 4 float values - this will have to be replaced by the SSE type

};

#endif