#ifndef RAPT_FLOATSIMDVECTORS_H_INCLUDED
#define RAPT_FLOATSIMDVECTORS_H_INCLUDED

/** This is datatype to represent 4 32-bit floating point numbers at once.
\todo: 
-write some unit tests
-write some performance tests
-use the 128-bit SSE type internally to speed up the computations in the arithmetic operators
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

//=================================================================================================

/** A class for convenient handling of SIMD optimizations for a vector of two double-precision
floating point numbers. 

-not yet tested

// see https://msdn.microsoft.com/de-de/library/tyy88x2a(v=vs.90).aspx
// https://msdn.microsoft.com/de-de/library/9b07190d(v=vs.90).aspx

*/

class rsFloat64x2 : public __m128d
{

public:

  /** \name Construction/Destruction */

  /** Standard constructor. Initializes both elements to zero. */
  rsFloat64x2() { *this = rsFloat64x2(_mm_setzero_pd()); }

  /** Constructor to copy an existing pair of values. */
  rsFloat64x2(__m128d value) : __m128d(value) {}

  /** Constructor that initializes both elements to the given value. */
  rsFloat64x2(double value) { *this = _mm_set1_pd(value); }
  //rsFloat64x2(double value) { *this = _mm_load1_pd(&value); }

  /** Constructor that initializes the elements from a 2-value array of doubles. */
  rsFloat64x2(double* values) 
  { 
    //*this = _mm_load_pd(values); // doesnt work
    *this = _mm_setr_pd(values[0], values[1]);
  }

  /** Constructor that initializes the elements from two doubles. */
  rsFloat64x2(double a, double b) { *this = _mm_setr_pd(a, b); }


  /** \name Setup */

  /** Sets both elements to a. */
  inline void set(double a) { *this = _mm_set1_pd(a); }
    // what's the difference between _mm_set1_pd and _mm_load1_pd? ...the latter takes a pointer?

  /** Sets the first element to a and the second element to b. */
  inline void set(double a, double b) { *this = _mm_setr_pd(a, b); }


  /** \name Inquiry */

  // extract vector elements:
  inline double get0() { double d; _mm_storel_pd(&d, *this); return d; }  // lower (index 0)
  inline double get1() { double d; _mm_storeh_pd(&d, *this); return d; }  // upper (index 1)


  /** \name Operators */

  // arithmetic operators:
  inline rsFloat64x2& operator+=(const rsFloat64x2& b) { return *this = rsFloat64x2(_mm_add_pd(*this, b)); }
  inline rsFloat64x2& operator-=(const rsFloat64x2& b) { return *this = rsFloat64x2(_mm_sub_pd(*this, b)); }
  inline rsFloat64x2& operator*=(const rsFloat64x2& b) { return *this = rsFloat64x2(_mm_mul_pd(*this, b)); }
  inline rsFloat64x2& operator/=(const rsFloat64x2& b) { return *this = rsFloat64x2(_mm_div_pd(*this, b)); }

  // unary minus (can we do better than that?):
  inline rsFloat64x2& operator-()
  { 
    static const rsFloat64x2 zero; // standard constructor constructs a zero
    return *this = rsFloat64x2(_mm_sub_pd(zero, *this));
  }
};

inline rsFloat64x2 operator+(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_add_pd(a, b)); }
inline rsFloat64x2 operator-(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_sub_pd(a, b)); }
inline rsFloat64x2 operator*(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_mul_pd(a, b)); }
inline rsFloat64x2 operator/(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_div_pd(a, b)); }
  // declaring return values as rsFloat64x2& gives compiler warning (unit test passes anyway)

inline rsFloat64x2& rsMin(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_min_pd(a, b)); }
inline rsFloat64x2& rsMax(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_max_pd(a, b)); }

inline rsFloat64x2& rsSqrt(const rsFloat64x2& a) { return rsFloat64x2(_mm_sqrt_pd(a)); }
/*
for implementid SIMD vectors, see:
http://johanmabille.github.io/blog/2014/10/10/writing-c-plus-plus-wrappers-for-simd-intrinsics-3/
https://github.com/p12tic/libsimdpp
...maybe it's possible to make a template baseclass, templated on the vector and element types 
like __m128d and double in this case

template<class TElem, class TVec>
class rsSIMDVector : public TVec
{
  public:

  TVec convert(TElem x);

  inline TVec& add(const TVec& a, const TVec& b)

};

class rsFloat64x2 : public rsSIMDVector<double, __m128d>
{

};

 */
#endif