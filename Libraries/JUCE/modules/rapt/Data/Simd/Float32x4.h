#ifndef RAPT_FLOAT32X4_H_INCLUDED
#define RAPT_FLOAT32X4_H_INCLUDED

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

class rsFloat64x2 : public __m128d
{

public:

  rsFloat64x2(__m128d value) : __m128d(value) {}
  rsFloat64x2(double value) { *this = _mm_load1_pd(&value); }

  // extract vector elements:
  //double get0() {}
  //double get1() {}



};

inline rsFloat64x2 operator+(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_add_pd(a, b)); }
inline rsFloat64x2 operator-(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_sub_pd(a, b)); }
inline rsFloat64x2 operator*(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_mul_pd(a, b)); }
inline rsFloat64x2 operator/(const rsFloat64x2& a, const rsFloat64x2& b) { return rsFloat64x2(_mm_div_pd(a, b)); }

inline rsFloat64x2 rsMin(rsFloat64x2 a, rsFloat64x2 b) { return rsFloat64x2(_mm_min_pd(a, b)); }
inline rsFloat64x2 rsMax(rsFloat64x2 a, rsFloat64x2 b) { return rsFloat64x2(_mm_max_pd(a, b)); }



#endif