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

  /** Constructs a value from int (needed for implicit conversions). */
  inline rsFloat64x2(int a) : v(_mm_set1_pd(double(a))) {}

  /** Constructor that initializes the elements from two doubles. */
  inline rsFloat64x2(double a, double b) : v(_mm_setr_pd(a, b)) {}

  /** Constructor that initializes the elements from a 2-value array of doubles. */
  inline rsFloat64x2(double* values)
  {
    //v = _mm_load_pd(values); // doesnt work
    v = _mm_setr_pd(values[0], values[1]);
  }

  /** \name Inquiry */

  /** Returns our vector as array of 2 doubles. */
  inline double* asArray() const { return (double*) &v; }

  /** Returns the vector element with index i (valid indices are 0 and 1). */
  inline double get(size_t i) const { return asArray()[i]; }

  /** Writes the two double into v0 and v1. */
  inline void get(double& v0, double& v1) const { _mm_storel_pd(&v0, v); _mm_storeh_pd(&v1, v); }

  /** Returns the 1st value (index 0). */
  //inline double get0() { double d; _mm_storel_pd(&d, v); return d; }
  //inline double get0() { return asArray()[0]; }
  inline double get0() const { return get(0); }
    // todo: performance test all versions
    // use get(0)

  /** Returns the 2nd value (index 1). */
  //inline double get1() { double d; _mm_storeh_pd(&d, v); return d; }
  //inline double get1() { return asArray()[1]; }
  inline double get1() const { return get(1); }

  /** Returns the sum of the values of both scalar elements in the vector. */
  inline double getSum() const { double* a = asArray(); return a[0]+a[1]; }

  /** Returns the minimum of the values of both scalar elements in the vector. */
  inline double getMin() const { double* a = asArray(); return (a[0] < a[1]) ? a[0] : a[1]; }

  /** Returns the maximum of the values of both scalar elements in the vector. */
  inline double getMax() const { double* a = asArray(); return (a[0] > a[1]) ? a[0] : a[1]; }


  /** \name Setup */

  /** Sets both elements to a. */
  inline void set(double a) { v = _mm_set1_pd(a); }
  // what's the difference between _mm_set1_pd and _mm_load1_pd? ...the latter takes a pointer?

  /** Sets the first element to a and the second element to b. */
  inline void set(double a, double b) { v = _mm_setr_pd(a, b); }

  /** Sets the vector element with index i (valid indices are 0 and 1). */
  inline void set(size_t i, double a)  { asArray()[i] = a; }


  inline void set0(double a) { set(size_t(0), a); }
  inline void set1(double a) { set(size_t(1), a); }

  //inline void set0(double a) { asArray()[0] = a; }
  //inline void set1(double a) { asArray()[1] = a; }



  /** \name Constants */

  /** Returns a vector that has a zero for both scalar elements. */
  inline static rsFloat64x2 zero() { static const __m128d z = _mm_setzero_pd(); return z; }

  /** Returns a vector that has a one for both scalar elements. */
  inline static rsFloat64x2 one()  { static const __m128d o = _mm_set1_pd(1.0); return o; }

  /** Returns a vector that has for both scalars a zero for the sign bit and the rest ones. This is
  useful for implementing the abs function. */
  inline static rsFloat64x2 signBitZero()
  {
    static const long long i = 0x7fffffffffffffff;
    static const double    d = *((double*)(&i));
    static const __m128d   r = _mm_set1_pd(d);
    return r;
  }

  /** Returns a vector that has for both scalars a one for the sign bit and the rest zeros. This is
  useful for implementing the sign function. */
  inline static rsFloat64x2 signBitOne()
  {
    static const long long i = 0x8000000000000000;
    static const double    d = *((double*)(&i));
    static const __m128d   r = _mm_set1_pd(d);
    return r;
  }


  /** \name Operators */

  /** Allows the two double values to be accessed (for reading and writing) as if this would be an
  array of two doubles. Valid indices are 0 and 1. */
  inline double& operator[](const int i) const { return asArray()[i]; }
    // maybe with this, we can get rid of set/get...hmm...but maybe not (or not just yet), they use
    // different operations - maybe they are faster? or safer (i.e. compatible with more
    // compilers? - test this first)

  // arithmetic operators:
  inline rsFloat64x2& operator+=(const rsFloat64x2& b) { v = _mm_add_pd(v, b); return *this; }
  inline rsFloat64x2& operator-=(const rsFloat64x2& b) { v = _mm_sub_pd(v, b); return *this; }
  inline rsFloat64x2& operator*=(const rsFloat64x2& b) { v = _mm_mul_pd(v, b); return *this; }
  inline rsFloat64x2& operator/=(const rsFloat64x2& b) { v = _mm_div_pd(v, b); return *this; }

  //inline rsFloat64x2 operator+(const rsFloat64x2& b) { return _mm_add_pd(v, b.v); }
    // nope - we define it outside the class because that allows lhs or rhs to be of scalar type

  /** Comparison for equality. For two vectors to be considered equal, both scalar elements must be
  equal. */
  inline bool operator==(const rsFloat64x2& b) const
  {
    //return (this[0] == b[0]) && (this[1] == b[1]);
    // gives warning: recursive on all control paths, will cause stack overflow

    double* A = asArray();
    double* B = b.asArray();
    return (A[0] == B[0]) && (A[1] == B[1]);
  }

  /** Comparison for equality. For two vectors to be considered unequal, it's sufficient when one
  of the scalar elements is unequal to the corresponding other. */
  inline bool operator!=(const rsFloat64x2& b) const
  { 
    return !(*this == b);
    //return (this[0] != b[0]) || (this[1] != b[1]); // produces compiler warning with msc
                                                     // about recursive call
  }

  /** Comparison for "less-than". For a vector to be considered less than another vector, both
  scalar elements must be less than the corresponding elements in the other vector. This
  definition is a bit arbitrary, but it makes sense for convergence tests in numerical algorithms
  (such as root finding) when the algorithm checks, whether a variable is within a given tolerance.
  These checks will evaluate to true only when both elements are within the tolerance, which is
  typically what is desired. */
  inline bool operator<(const rsFloat64x2& b) const
  {
    return (this[0] < b[0]) && (this[1] < b[1]);
    //double* A = asArray();
    //double* B = b.asArray();
    //return (A[0] < B[0]) && (A[1] < B[1]);
  }


  //// unary minus (commented out because the implementation outside the class turned out to be
  //// more efficient (suprisingly - maybe more tests are needed, why)
  //inline rsFloat64x2 operator-()
  //{
  //  return _mm_xor_pd(signBitOne(), v);

  //  //return _mm_sub_pd(zero(), v); // alternative implementation

  //  //static const __m128d zero = _mm_setzero_pd();
  //  //return _mm_sub_pd(zero, v);
  //}

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
inline rsFloat64x2 operator+(const rsFloat64x2& a) { return a; }                    // unary plus
inline rsFloat64x2 operator-(const rsFloat64x2& a) { return rsFloat64x2(0.0) - a; } // unary minus
// the binary operators with a scalar for the left or right hand side do not have to be defined due
// to implicit conversions

// limiting functions::
inline rsFloat64x2 rsMin(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_min_pd(a, b); }
inline rsFloat64x2 rsMax(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_max_pd(a, b); }
inline rsFloat64x2 rsClip(const rsFloat64x2& x, const rsFloat64x2& min, const rsFloat64x2& max)
{
  return rsMax(rsMin(x, max), min);
}
//inline rsFloat64x2 rsClip(rsFloat64x2 x, rsFloat64x2 min, rsFloat64x2 max)
//{
//  return rsMax(rsMin(x, max), min);
//}

// bit-manipulations and related functions:
inline rsFloat64x2 rsBitAnd(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_and_pd(a, b); }
inline rsFloat64x2 rsBitOr( const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_or_pd( a, b); }
inline rsFloat64x2 rsBitXor(const rsFloat64x2& a, const rsFloat64x2& b) { return _mm_xor_pd(a, b); }
inline rsFloat64x2 rsAbs(const rsFloat64x2& a) { return rsBitAnd(a, rsFloat64x2::signBitZero()); }
inline rsFloat64x2 rsSign(const rsFloat64x2& a)
{
  rsFloat64x2 signOnly = rsBitAnd(a, rsFloat64x2::signBitOne());
  return rsBitOr(signOnly, rsFloat64x2::one());
}

// reordering elements:
inline rsFloat64x2 rsSwap(const rsFloat64x2& x) { return rsFloat64x2(x[1], x[0]); }
// for rsFloat32x4, have a function rsReorder(rsFloat32x4 x, int to0, int to1, int to2, int to3)
// ...or maybe the to1, ...etc. can be template arguments, evaluated at compile time (might be faster)?


// math functions (except for sqrt, we need to fall back to the scalar versions):
inline rsFloat64x2 rsSqrt(const rsFloat64x2& a) { return _mm_sqrt_pd(a); }
inline rsFloat64x2 rsExp(const rsFloat64x2& x) { double* a = x.asArray(); return rsFloat64x2(exp(a[0]), exp(a[1])); }
inline rsFloat64x2 rsLog(const rsFloat64x2& x) { double* a = x.asArray(); return rsFloat64x2(log(a[0]), log(a[1])); }
inline rsFloat64x2 rsSin(const rsFloat64x2& x) { double* a = x.asArray(); return rsFloat64x2(sin(a[0]), sin(a[1])); }
inline rsFloat64x2 rsCos(const rsFloat64x2& x) { double* a = x.asArray(); return rsFloat64x2(cos(a[0]), cos(a[1])); }
inline rsFloat64x2 rsTan(const rsFloat64x2& x) { double* a = x.asArray(); return rsFloat64x2(tan(a[0]), tan(a[1])); }
// todo: asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, pow, fmod, floor, ceil
// get rid of the rs-prefix for all functions ...maybe we can get rid of that temporary a - we have the [] operator
// now - so x[0], x[1] should also work


// typecasting of pointers (to arrays) back and forth between double and rsFloat64x2:
inline double* rsCastPointer(rsFloat64x2* p) { return reinterpret_cast<double*> (p);      }
inline rsFloat64x2* rsCastPointer(double* p) { return reinterpret_cast<rsFloat64x2*> (p); }
// test this! it's quite hacky


// Explicit specializations of comparison functions that are used in RAPT. They return true, iff 
// *all* values in the 1st argument are less/greater/etc. than the corresponding values in the 
// 2nd argument because that's what's typically needed in numerical algorithms (for convergence
// tests etc.). They don't satisfy the trichotomy rule, though. This rule states that exactly 1
// of the 3 conditions must hold: x < y, x == y, x > y, but this rule is not useful for numerics.
// (it's useful for sorting, though...)
inline bool rsGreaterAbs(const rsFloat64x2& x, const rsFloat64x2& y)
{
  return fabs(x[0]) > fabs(y[0]) && fabs(x[1]) > fabs(y[1]);
}
// todo: rsLessAbs


//inline rsFloat64x2 exp(const rsFloat64x2& x) { double* a = x.asArray(); return rsFloat64x2(exp(a[0]), exp(a[1])); }
// without the rs-prefix, i get a shitload of compiler errors related to RealFunctions.h - perhaps because i
// invoke a non-templated version (for double) there - maybe it can be fixed by using std:: or :: there?


// for more inspiration, see:
// https://msdn.microsoft.com/de-de/library/tyy88x2a(v=vs.90).aspx
// https://msdn.microsoft.com/de-de/library/9b07190d(v=vs.90).aspx
// http://johanmabille.github.io/blog/2014/10/10/writing-c-plus-plus-wrappers-for-simd-intrinsics-3/
// https://github.com/p12tic/libsimdpp
// https://github.com/VcDevel/Vc
// https://github.com/NumScale/boost.simd

// for evaluating elementary math functions without resorting to the scalar versions, see:
// http://ito-lab.naist.jp/~n-sibata/pdfs/isc10simd.pdf

// reference:
// http://www.info.univ-angers.fr/pub/richer/ens/l3info/ao/intel_intrinsics.pdf

#endif
