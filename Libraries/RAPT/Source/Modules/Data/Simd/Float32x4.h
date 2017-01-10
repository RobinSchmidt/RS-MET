#ifndef RAPT_FLOAT32X4_H_INCLUDED
#define RAPT_FLOAT32X4_H_INCLUDED

/** This is datatype to represent 4 32-bit floating point numbers at once.
\todo: 
-write some unit tests
-write some performance tests
-use the 128-bit SSE type internally to speed up the computations in the arithmetic operators
-compare performance of SSE implementation to non SSE
*/

class Float32x4
{

public:

  /** \name Construction */

  /** Constructor. Sets up the 4 elements to the given values. */
  Float32x4(float v0 = 0.f, float v1 = 0.f, float v2 = 0.f, float v3 = 0.f)
  {
    setValues(v0, v1, v2, v3);
  }

  /** Constructor. Sets up the all 4 elements to the same give value. */
  Float32x4(float value)
  {
    setValues(value, value, value, value);
  }

  /** Constructor for conversion from double to float and setting up all 4 elements to the same 
  value. */
  Float32x4(double value)
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
  inline Float32x4 operator-() const
  { 
    return Float32x4(-v[0], -v[1], -v[2], -v[3]); 
  }

  /** Addition of 2 Float32x4 vectors. */
  inline Float32x4 operator+(const Float32x4 &y)
  {
    return Float32x4(v[0]+y.v[0], v[1]+y.v[1], v[2]+y.v[2], v[3]+y.v[3]);
  }

  /** Subtraction of 2 Float32x4 vectors. */
  inline Float32x4 operator-(const Float32x4 &y)
  {
    return Float32x4(v[0]-y.v[0], v[1]-y.v[1], v[2]-y.v[2], v[3]-y.v[3]);
  }

  /** Multiplication of 2 Float32x4 vectors. */
  inline Float32x4 operator*(const Float32x4 &y)
  {
    return Float32x4(v[0]*y.v[0], v[1]*y.v[1], v[2]*y.v[2], v[3]*y.v[3]);
  }

  /** Division of 2 Float32x4 vectors. */
  inline Float32x4 operator/(const Float32x4 &y)
  {
    return Float32x4(v[0]/y.v[0], v[1]/y.v[1], v[2]/y.v[2], v[3]/y.v[3]);
  }


protected:

  float v[4];  // vector of 4 float values - this will have to be replaced by the SSE type

};

#endif