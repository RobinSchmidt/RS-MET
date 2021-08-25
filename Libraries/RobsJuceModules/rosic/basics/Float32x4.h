#ifndef RAPT_FLOAT32X4_H_INCLUDED
#define RAPT_FLOAT32X4_H_INCLUDED

//#include <emmintrin.h>
//#include <xmmintrin.h>

/** The code here is completely analogous to rsFloat64x2. */


//=================================================================================================
#if defined(RS_NO_SIMD_FLOAT32X4)

/** Fallback implementation. */

class rsFloat32x4
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction */

  /** Standard constructor. Leaves all 4 elements uninitialized. */
  inline rsFloat32x4() {}

  /** Constructor that initializes all elements to the given value. */
  inline rsFloat32x4(float a) { v[0] = v[1] = v[2] = v[3] = a; }

  /** Constructs a value from int (needed for implicit conversions). */
  inline rsFloat32x4(int a) { v[0] = v[1] = v[2] = v[3] = (float)a; }

  /** Constructs a value from double (needed for implicit conversions). */
  inline rsFloat32x4(double a) { v[0] = v[1] = v[2] = v[3] = (float)a; }

  /** Constructor that initializes the elements from four floats. */
  inline rsFloat32x4(float a, float b, float c, float d) { v[0] = a; v[1] = b; v[2] = c; v[3] = d; }

  /** Constructor that initializes the elements from a 4-value array of floats. */
  inline rsFloat32x4(float* p) { v[0] = p[0]; v[1] = p[1]; v[2] = p[2]; v[3] = p[3]; }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns our vector as array of 4 floats. */
  inline const float* asArray() const { return v; }

  /** Returns our vector as constant array of 4 floats. */
  inline float* asArray() { return v; }

  /** Returns the vector element with index i (valid indices are 0,1,2,3). */
  inline float get(size_t i) const { return asArray()[i]; }
  // redundant with [] ...but this is const

  /** Writes our vector into the 4-element float array p. (needs test, maybe implement a similar
  function for rsFloat64x2 - this has been added after copy/paste ) */
  inline void get(float* p) const { p[0] = v[0]; p[1] = v[1]; p[2] = v[2]; p[3] = v[3]; }

  /** Returns the sum of the values of both scalar elements in the vector. */
  inline float getSum() const { return v[0]+v[1]+v[2]+v[3]; }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets all 4 elements to the given number a. */
  inline void set(float a) 
  { 
    v[0] = v[1] = v[2] = v[3] = a;
  }

  /** Sets the 4 elements to the given numbers. */
  inline void set(float a, float b, float c, float d) { v[0] = a; v[1] = b; v[2] = c; v[3] = d; }

  /** Sets the 4 elements to the given integers, thereby converting them to float. */
  inline void set(int a, int b, int c, int d) { set((float) a, (float) b, (float) c, (float) d); }


  //-----------------------------------------------------------------------------------------------
  /** \name Constants */

  // the SSE implementation has code here...maybe copy it over and edit it or get rid of it in
  // the SSE version...it doesn't seem to be used


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Allows the four float values to be accessed (for reading and writing) as if this would be an
  array of four floats. Valid indices are 0,1,2,3. */
  inline float& operator[](const int i) { return asArray()[i]; }

  inline const float& operator[](const int i) const { return asArray()[i]; }

  // arithmetic operators:
  inline rsFloat32x4& operator+=(const rsFloat32x4& b) { v[0] += b.v[0]; v[1] += b.v[1]; v[2] += b.v[2]; v[3] += b.v[3]; return *this; }
  inline rsFloat32x4& operator-=(const rsFloat32x4& b) { v[0] -= b.v[0]; v[1] -= b.v[1]; v[2] -= b.v[2]; v[3] -= b.v[3]; return *this; }
  inline rsFloat32x4& operator*=(const rsFloat32x4& b) { v[0] *= b.v[0]; v[1] *= b.v[1]; v[2] *= b.v[2]; v[3] *= b.v[3]; return *this; }
  inline rsFloat32x4& operator/=(const rsFloat32x4& b) { v[0] /= b.v[0]; v[1] /= b.v[1]; v[2] /= b.v[2]; v[3] /= b.v[3]; return *this; }

  /** Comparison for equality. For two vectors to be considered equal, all scalar elements must be
  equal. */
  //inline bool operator==(const rsFloat32x4& b) const
  //{ return v[0] == b[0] && v[1] == b[1] && v[2] == b[2] && v[3] == b[3]; }


protected:

  float v[4];  // todo: maybe use an alignment specifier

};

// binary arithmetic operators:
inline rsFloat32x4 operator+(const rsFloat32x4& a, const rsFloat32x4& b) { return rsFloat32x4(a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]); }
inline rsFloat32x4 operator-(const rsFloat32x4& a, const rsFloat32x4& b) { return rsFloat32x4(a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3]); }
inline rsFloat32x4 operator*(const rsFloat32x4& a, const rsFloat32x4& b) { return rsFloat32x4(a[0]*b[0], a[1]*b[1], a[2]*b[2], a[3]*b[3]); }
inline rsFloat32x4 operator/(const rsFloat32x4& a, const rsFloat32x4& b) { return rsFloat32x4(a[0]/b[0], a[1]/b[1], a[2]/b[2], a[3]/b[3]); }


// limiting functions::
inline rsFloat32x4 rsMin(const rsFloat32x4& a, const rsFloat32x4& b) 
{ 
  return rsFloat32x4(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]), std::min(a[3], b[3]));
}
inline rsFloat32x4 rsMax(const rsFloat32x4& a, const rsFloat32x4& b) 
{ 
  return rsFloat32x4(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]), std::max(a[3], b[3]));
}

inline rsFloat32x4 rsAbs(const rsFloat32x4& a) 
{ 
  return rsFloat32x4(std::abs(a[0]), std::abs(a[1]), std::abs(a[2]), std::abs(a[3]));
}

inline rsFloat32x4 rsSqrt(const rsFloat32x4& a) 
{ 
  return rsFloat32x4(std::sqrt(a[0]), std::sqrt(a[1]), std::sqrt(a[2]), std::sqrt(a[3]));
}

inline float rsSign(float x) { return float(0.f < x) - float(x < 0.f); } // why needed?
inline rsFloat32x4 rsSign(rsFloat32x4 a)
{
  return rsFloat32x4(rsSign(a[0]), rsSign(a[1]), rsSign(a[2]), rsSign(a[3]));
}


#elif defined(RS_INSTRUCTION_SET_SSE)
//=================================================================================================

/** SSE implementation. */

class rsFloat32x4
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction */

  /** Standard constructor. Leaves both elements uninitialized. */
  inline rsFloat32x4() { /*v = _mm_setzero_pd();*/ }

  /** Constructor to copy an existing __m128 vector of values. */
  inline rsFloat32x4(const __m128 rhs) : v(rhs) {}

  /** Constructor that initializes both elements to the given value. */
  inline rsFloat32x4(float a) : v(_mm_set1_ps(a)) {}

  /** Constructs a value from int (needed for implicit conversions). */
  inline rsFloat32x4(int a) : v(_mm_set1_ps(float(a))) {}
  // gcc produces error (in 32-bit configuration);
  // inlining failed in call to always_inline '__m128 _mm_set1_ps(float)': target specific option mismatch|
  // unless the SSE instruction set is explicitly activated in the compiler settings

  /** Constructs a value from double (needed for implicit conversions). */
  inline rsFloat32x4(double a) : v(_mm_set1_ps(float(a))) {}

  /** Constructor that initializes the elements from four floats. */
  inline rsFloat32x4(float a, float b, float c, float d) : v(_mm_setr_ps(a, b, c, d)) {}

  /** Constructor that initializes the elements from a 4-value array of floats. */
  inline rsFloat32x4(float* p) { v = _mm_setr_ps(p[0], p[1], p[2], p[3]); }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns our vector as array of 4 floats. */
  inline float* asArray() const { return (float*) &v; }

  /** Returns the vector element with index i (valid indices are 0,1,2,3). */
  inline float get(size_t i) const { return asArray()[i]; }
  // redundant with [] ...but this is const

  /** Writes our vector into the 4-element float array p. (needs test, maybe implement a similar
  function for rsFloat64x2 - this has been added after copy/paste ) */
  inline void get(float* p) const { _mm_store_ps(p, v); }

  /** Returns the sum of the values of both scalar elements in the vector. */
  inline float getSum() const { float* a = asArray(); return a[0]+a[1]+a[2]+a[3]; }

  /** Returns the minimum of the values of both scalar elements in the vector. */
  //inline float getMin() const { float* a = asArray(); return (a[0] < a[1]) ? a[0] : a[1]; }

  /** Returns the maximum of the values of both scalar elements in the vector. */
  //inline float getMax() const { float* a = asArray(); return (a[0] > a[1]) ? a[0] : a[1]; }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets all 4 elements to the given number a. */
  inline void set(float a) { v = _mm_set1_ps(a); }

  /** Sets the 4 elements to the given numbers. */
  inline void set(float a, float b, float c, float d) { v = _mm_setr_ps(a, b, c, d); }

  /** Sets the 4 elements to the given integers, thereby converting them to float. */
  inline void set(int a, int b, int c, int d) { set((float) a, (float) b, (float) c, (float) d); }


  inline void set(const __m128& a) { v = a; }




  /** Sets the vector element with index i (valid indices are 0 and 1). */
  //inline void set(size_t i, float a)  { asArray()[i] = a; }
  // redundant with array access operator - maybe delete from rsFloat64x2, too

  //-----------------------------------------------------------------------------------------------
  /** \name Constants */

  /** Returns a vector that has a zero for all scalar elements. */
  inline static rsFloat32x4 zero() { static const __m128 z = _mm_setzero_ps(); return z; }

  /** Returns a vector that has a one for all scalar elements. */
  inline static rsFloat32x4 one()  { static const __m128 o = _mm_set1_ps(1.f); return o; }

  /** Returns a vector that has for all scalars a zero for the sign bit and the rest ones. This is
  useful for implementing the abs function. */
  inline static rsFloat32x4 signBitZero()
  {
    static const RAPT::rsInt32 i = 0x7fffffff;
    static const float   d = *((float*)(&i));
    static const __m128  r = _mm_set1_ps(d);
    return r;
  }

  /** Returns a vector that has for all scalars a one for the sign bit and the rest zeros. This is
  useful for implementing the sign function. */
  inline static rsFloat32x4 signBitOne()
  {
    static const RAPT::rsInt32 i = 0x80000000;
    static const float   d = *((float*)(&i));
    static const __m128  r = _mm_set1_ps(d);
    return r;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Allows the four float values to be accessed (for reading and writing) as if this would be an
  array of four floats. Valid indices are 0,1,2,3. */
  inline float& operator[](const int i) const { return asArray()[i]; }
  // maybe with this, we can get rid of set/get...hmm...but maybe not (or not just yet), they use
  // different operations - maybe they are faster? or safer (i.e. compatible with more
  // compilers? - test this first)

  // arithmetic operators:
  inline rsFloat32x4& operator+=(const rsFloat32x4& b) { v = _mm_add_ps(v, b); return *this; }
  inline rsFloat32x4& operator-=(const rsFloat32x4& b) { v = _mm_sub_ps(v, b); return *this; }
  inline rsFloat32x4& operator*=(const rsFloat32x4& b) { v = _mm_mul_ps(v, b); return *this; }
  inline rsFloat32x4& operator/=(const rsFloat32x4& b) { v = _mm_div_ps(v, b); return *this; }



  /** Comparison for equality. For two vectors to be considered equal, all scalar elements must be
  equal. */
  /*
  inline bool operator==(const rsFloat32x4& b) const
  {
    float* a = asArray();
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
    //return (this[0] == b[0]) && (this[1] == b[1]) && (this[2] == b[2]) && (this[3] == b[3]); // stack overflow - why?
    // todo: check, if there's an intrinsic vector comparison function...
    // check with rsFloat64x2, too - and why there is no stack overflow when using this[0]....
  }
  */


  // maybe these operators really should return an rsFloat32x4 - like the agner fog implementation

  /*
  // drag them out of the class - the implementation is the same in both cases
  inline bool operator<(const rsFloat32x4& b) const
  {
    float* a = asArray();
    return a[0] < b[0] && a[1] < b[1] && a[2] < b[2] && a[3] < b[3];
  }

  inline bool operator>(const rsFloat32x4& b) const
  {
    float* a = asArray();
    return a[0] > b[0] && a[1] > b[1] && a[2] > b[2] && a[3] > b[3];
  }

  inline bool operator>=(const rsFloat32x4& b) const
  {
    float* a = asArray();
    return a[0] >= b[0] && a[1] >= b[1] && a[2] >= b[2] && a[3] >= b[3];
  }

  inline bool operator<=(const rsFloat32x4& b) const
  {
    float* a = asArray();
    return a[0] <= b[0] && a[1] <= b[1] && a[2] <= b[2] && a[3] <= b[3];
  }
  */



  inline rsFloat32x4& operator=(const __m128& rhs) { v = rhs; return *this; }
  //inline rsFloat32x4& operator=(__m128 const & rhs) { v = rhs; return *this;  }


  inline operator __m128() const { return v; }


protected:

  __m128 v;
  //__declspec(align(16)) __m128 v; // the value (define and ALIGN(N) macro for gcc/msc)
  //__declspec(align(32)) __m128 v;
};

// binary arithmetic operators:
inline rsFloat32x4 operator+(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_add_ps(a, b); }
inline rsFloat32x4 operator-(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_sub_ps(a, b); }
inline rsFloat32x4 operator*(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_mul_ps(a, b); }
inline rsFloat32x4 operator/(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_div_ps(a, b); }

// limiting functions::
inline rsFloat32x4 rsMin(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_min_ps(a, b); }
inline rsFloat32x4 rsMax(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_max_ps(a, b); }

// bit-manipulations and related functions:
inline rsFloat32x4 rsBitAnd(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_and_ps(a, b); }
inline rsFloat32x4 rsBitOr( const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_or_ps( a, b); }
inline rsFloat32x4 rsBitXor(const rsFloat32x4& a, const rsFloat32x4& b) { return _mm_xor_ps(a, b); }
inline rsFloat32x4 rsAbs(const rsFloat32x4& a) { return rsBitAnd(a, rsFloat32x4::signBitZero()); }
inline rsFloat32x4 rsSign(const rsFloat32x4& a)
{
  rsFloat32x4 signOnly = rsBitAnd(a, rsFloat32x4::signBitOne());
  return rsBitOr(signOnly, rsFloat32x4::one());
}

inline rsFloat32x4 rsSqrt(const rsFloat32x4& a) { return _mm_sqrt_ps(a); }

#else
#error One of SIMD instruction sets or RS_NO_SIMD_FLOAT32X4 has to be defined.
#endif  // #if defined(RS_NO_SIMD_FLOAT32X4)

//=================================================================================================
// Functions and operators that do not depend on any particular instruction set:

inline rsFloat32x4 operator+(const rsFloat32x4& a) { return a; }                    // unary plus
inline rsFloat32x4 operator-(const rsFloat32x4& a) { return rsFloat32x4(0.f) - a; } // unary minus

// Comparisons:
inline bool operator< (const rsFloat32x4& a, const rsFloat32x4& b) { return a[0] <  b[0] && a[1] <  b[1] && a[2] <  b[2] && a[3] <  b[3]; }
inline bool operator<=(const rsFloat32x4& a, const rsFloat32x4& b) { return a[0] <= b[0] && a[1] <= b[1] && a[2] <= b[2] && a[3] <= b[3]; }
inline bool operator> (const rsFloat32x4& a, const rsFloat32x4& b) { return a[0] >  b[0] && a[1] >  b[1] && a[2] >  b[2] && a[3] >  b[3]; }
inline bool operator>=(const rsFloat32x4& a, const rsFloat32x4& b) { return a[0] >= b[0] && a[1] >= b[1] && a[2] >= b[2] && a[3] >= b[3]; }
inline bool operator==(const rsFloat32x4& a, const rsFloat32x4& b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]; }
inline bool operator!=(const rsFloat32x4& a, const rsFloat32x4& b) { return !(a == b); }

// Math:
inline rsFloat32x4 rsExp(const rsFloat32x4& x) { return rsFloat32x4(exp(x[0]), exp(x[1]), exp(x[2]), exp(x[3])); }
inline rsFloat32x4 rsLog(const rsFloat32x4& x) { return rsFloat32x4(log(x[0]), log(x[1]), log(x[2]), log(x[3])); }
inline rsFloat32x4 rsSin(const rsFloat32x4& x) { return rsFloat32x4(sin(x[0]), sin(x[1]), sin(x[2]), sin(x[3])); }
inline rsFloat32x4 rsCos(const rsFloat32x4& x) { return rsFloat32x4(cos(x[0]), cos(x[1]), cos(x[2]), cos(x[3])); }
inline rsFloat32x4 rsTan(const rsFloat32x4& x) { return rsFloat32x4(tan(x[0]), tan(x[1]), tan(x[2]), tan(x[3])); }

inline rsFloat32x4 rsClip(const rsFloat32x4& x, const rsFloat32x4& min, const rsFloat32x4& max) { return rsMax(rsMin(x, max), min); }


// maybe implement recriprocal and reciprocal sqrt (there are intrinsics for these)



//=================================================================================================
// Notes:

// see here for optimized vector math functions - free and open-source:
// http://gruntthepeon.free.fr/ssemath/
// code is available in the _third_party folder of rosic

/*

// reference:
// http://www.info.univ-angers.fr/pub/richer/ens/l3info/ao/intel_intrinsics.pdf

load:
__m128 _mm_load1_ps(float* p)    Loads a single SP FP value, copying it into all four words
__m128 _mm_load_ps(float* p)     Loads four SP FP values. The address must be 16-byte-aligned
__m128 _mm_loadu_ps(float* p)    Loads four SP FP values. The address need not be 16-byte-aligned.

set:
__m128 _mm_setzero_ps(void)                            Clears the four SP FP values.
__m128 _mm_set1_ps(float w )                           Sets the four SP FP values to w.
__m128 _mm_set_ps(float z, float y, float x, float w)  Sets the four SP FP values to the four inputs.
__m128 _mm_setr_ps(float z, float y, float x, float w) Sets the four SP FP values to the four inputs in reverse order.
..check which one actually reverses - the reference seems contradictory there

store:
void _mm_store_ps(float *p, __m128 a)    Stores four SP FP values. The address must be 16-byte-aligned.
void _mm_storeu_ps(float *p, __m128 a)   Stores four SP FP values. The address need not be 16-byte-aligned.

function that operate simultaneously on all elements:

math:
__m128 _mm_sqrt_ps(__m128 a)              square root
__m128 _mm_rcp_ps(__m128 a)               (approximate?) reciprocal
__m128 _mm_rsqrt_ps(__m128 a)             (approximate?) reciprocal square root
__m128 _mm_min_ps(__m128 a, __m128 b)     minimum
__m128 _mm_max_ps(__m128 a, __m128 b)     maximum

comparisons:
__m128 _mm_cmpeq_ps(__m128 a, __m128 b)    equality
__m128 _mm_cmplt_ps(__m128 a, __m128 b)    less than
__m128 _mm_cmple_ps(__m128 a, __m128 b)    less or equal
__m128 _mm_cmpgt_ps(__m128 a, __m128 b)    greater than
__m128 _mm_cmpge_ps(__m128 a, __m128 b)    greater or equal
__m128 _mm_cmpneq_ps(__m128 a, __m128 b)   inequality
__m128 _mm_cmpnlt_ps(__m128 a, __m128 b)   not less than
...reference page 38
the ps suffix operates on all elements, the same functions with the ss suffix operate only on the
1st and pass through all other values





*/


#ifdef RS_NEVER_DEFINED // code below obsolete - but maybe keep it as fall back when no SSE is available

/** This is datatype to represent 4 32-bit floating point numbers at once.
THIS IS NOT USABLE YET (it currently uses 4 actual floats)

\todo:
-write some unit tests
-write some performance tests
-use the 128-bit SSE type internally to speed up the computations in the arithmetic operators
-currently, this is only a prototype use 4 actual floats, see Float64x2 and model it after that
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

#endif
