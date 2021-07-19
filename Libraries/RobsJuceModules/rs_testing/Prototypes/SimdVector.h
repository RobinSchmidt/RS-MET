#pragma once

// Should eventually go into some configuration file:
#define RS_USE_SSE
//#define RS_USE_SSE2
//#define RS_USE_AVX
//#define RS_USE_AVX2

//-------------------------------------------------------------------------------------------------

/** A class for abstracting SIMD (single instruction, multiple data) processing. Client code can
work with vectors of, for example, 4 single precision floating point values as follows:

  rsSimdVector<float, 4> a, b, c;              // create vectors a, b, c
  a[0] =  2; a[1] =  3; a[2] =  5; a[3] =  7;  // assign elements of a
  b[0] = 11; b[1] = 13; b[2] = 17; b[3] = 19;  // assign elements of b
  c = a + b;                                   // compute element-wise sum

The scalar type T must be one of the primitive numerical datatypes in C++, like char, int, float, 
double, etc. and the vector size N must be a power of 2 and it must be >= 2.

Whether or not actual SIMD processing will take place on the hardware level will depend on the
compiler settings and preprocessor definitions. For example, if you compile for x86 and the macro
RS_USE_SSE is defined, the above code will actually compile to hardware SIMD instructions by 
virtue of having an explicit specialization of the template available (which would otherwise be
#ifdef'd out) that maps the code to operations on the __m128 SSE datatype. If that macro is not 
defined, the code will compile to 4 scalar operations for each of the vector elements. If you would 
have tried to use 8-vectors via rsSimdVector<float, 8>, you would need to have RS_USE_AVX defined 
to get operations mapped to __m256. If you don't have the AVX macro defined but you do have the 
SSE macro defined, the 8-vector code would not compile to 8 scalar instructions but to 2 SSE 
instructions on the lower and higher __m128 part of the __m256 vector. The class will always try 
to use as much hardware SIMD as is available as per the macro definitions and emulate the rest in
software, completely transparent to the client code. Client code will always look the same, 
regardless of the targeted instruction set.

Of course, this creates some caveats: client code may falsely assume that hardware SIMD processing 
takes place, when in fact, it's only emulated. Moreover, our compile-time dispatch here says 
nothing about the availability of the instructions at runtime. Client code will still need runtime 
dispatches to code with less SIMD or none at all to make sure that the code runs on machines with 
lesser instruction sets. */

template<class T, int N> // T: underlying scalar type, N: vector size (power of 2, >= 2)
class rsSimdVector
{

public:

  using V   = rsSimdVector<T, N>;
  using CV  = const V;
  using V2  = rsSimdVector<T, N/2>;
  using CV2 = const V2;

  V() {}
  V(T a) { lo() = a; hi() = a; }
  V(CV2& low, CV2& high) { lo() = low; hi() = high; }
  V(CV2&& low, CV2&& high) { lo() = low; hi() = high; } // do we need this?
  //V(CV2 low, CV2 high) { lo() = low; hi() = high; }
  // When we have the CV2& and CV2 versions both uncommented, we get a compile error:
  // "...constructor overload resolution was ambiguous" for rsSimdVector<double,16>, keeping only
  // one of them gives a compile error for rsSimdVector<float,4> when the new rsAbs implementation
  // is uncommented

  /** Write Access to the i-th element of the vector where i = 0,...,N-1. */
  inline T& operator[](const int i) { rsStaticAssert(i >= 0 && i < N); return v[i]; }

  /** Read Access to the i-th element of the vector where i = 0,...,N-1. */
  inline const T& operator[](const int i) const { rsStaticAssert(i >= 0 && i < N); return v[i]; }

  /** Returns the level of software emulation that takes place for this particular template 
  instantiation with the current compiler- and preprocessor settings. If the instantiation 
  provides full hardware SIMD, this function will return 0. If the hardware SIMD kicks in on the 
  level of half-vectors, it returns 1. If hardware SIMD kicks in on the level quarter-vectors, it 
  returns 2. And so on. For example, for rsSimdVector<float, 16> when you only have SSE 
  available, the 16 element SIMD vector will be emulated by 4 SSE __m128 variables, so the 
  emulation-level will be 2. */
  static int getEmulationLevel() { return V2::getEmulationLevel() + 1; }


  V operator+(CV b) const { return V(lo()+b.lo(), hi()+b.hi()); }
  V operator-(CV b) const { return V(lo()-b.lo(), hi()-b.hi()); }
  V operator*(CV b) const { return V(lo()*b.lo(), hi()*b.hi()); }
  V operator/(CV b) const { return V(lo()/b.lo(), hi()/b.hi()); }
  // ToDo: try, if it makes a difference (performance-wise), if we have CV& and/or CV&& versions

  //V operator+(CV& b) const { return V(lo()+b.lo(), hi()+b.hi()); }
  //V operator-(CV& b) const { return V(lo()-b.lo(), hi()-b.hi()); }
  //V operator*(CV& b) const { return V(lo()*b.lo(), hi()*b.hi()); }
  //V operator/(CV& b) const { return V(lo()/b.lo(), hi()/b.hi()); }

  //V operator+(CV&& b) const { return V(lo()+b.lo(), hi()+b.hi()); }
  //V operator-(CV&& b) const { return V(lo()-b.lo(), hi()-b.hi()); }
  //V operator*(CV&& b) const { return V(lo()*b.lo(), hi()*b.hi()); }
  //V operator/(CV&& b) const { return V(lo()/b.lo(), hi()/b.hi()); }


  /** Returns a reference to the lower half-vector. */
  inline V2& lo() { return *((V2*) &v[0]); }

  /** Returns a reference to the higher half-vector. */
  inline V2& hi() { return *((V2*) &v[N/2]); }

  /** Returns a const reference to the lower half-vector. */
  inline CV2& lo() const { return *((CV2*) &v[0]); }

  /** Returns a const reference to the higher half-vector. */
  inline CV2& hi() const { return *((CV2*) &v[N/2]); }
  // Do we need to consider endianness to decide which is the lower and which is the upper 
  // half-vector? Maybe it matters for the actual memory layout - but actually, the memory layout
  // should probably not need to bother us? I think, we don't really need to care which is which as
  // long as everything is consistent.


  T v[N]; // maybe we should specify an alignment?
};

// Some macros to shorten the boilerplate:
#define V rsSimdVector<T, N>
#define CV const V
#define TIV template<class T, int N> inline V

TIV operator+(const T& a, CV& b) { V c; c.lo()=a+b.lo(); c.hi()=a+b.hi(); return c; }
TIV operator-(const T& a, CV& b) { V c; c.lo()=a-b.lo(); c.hi()=a-b.hi(); return c; }
TIV operator*(const T& a, CV& b) { V c; c.lo()=a*b.lo(); c.hi()=a*b.hi(); return c; }
TIV operator/(const T& a, CV& b) { V c; c.lo()=a/b.lo(); c.hi()=a/b.hi(); return c; }
TIV operator+(const V& a) { return a; }        // unary plus
TIV operator-(const V& a) { return V(0) - a; } // unary minus - can we do better?
// ToDo: try passing arguments by value, check, if this incurs a performance hit

// Unary functions (int and float):
TIV rsAbs( V x) { return V(rsAbs( x.lo()), rsAbs( x.hi())); }
TIV rsSign(V x) { return V(rsSign(x.lo()), rsSign(x.hi())); }

// Unary functions (float only):
TIV rsFloor(V x) { return V(rsFloor(x.lo()), rsFloor(x.hi())); }
TIV rsCeil( V x) { return V(rsCeil( x.lo()), rsCeil( x.hi())); }
TIV rsRound(V x) { return V(rsRound(x.lo()), rsRound(x.hi())); }

TIV rsExp( V x) { return V(rsExp( x.lo()), rsExp( x.hi())); }
TIV rsLog( V x) { return V(rsLog( x.lo()), rsLog( x.hi())); }
TIV rsSqrt(V x) { return V(rsSqrt(x.lo()), rsSqrt(x.hi())); }

TIV rsSin(V x) { return V(rsSin(x.lo()), rsSin(x.hi())); }
TIV rsCos(V x) { return V(rsCos(x.lo()), rsCos(x.hi())); }
TIV rsTan(V x) { return V(rsTan(x.lo()), rsTan(x.hi())); }

TIV rsSinh(V x) { return V(rsSinh(x.lo()), rsSinh(x.hi())); }
TIV rsCosh(V x) { return V(rsCosh(x.lo()), rsCosh(x.hi())); }
TIV rsTanh(V x) { return V(rsTanh(x.lo()), rsTanh(x.hi())); }

// ToDo:
// asin, acos, atan, asinh, acosh, atanh, cbrt, expm1, exp2, log2, splitIntFrac, erf, tgamma
// maybe write a macros RS_VECTORIZE_1, so we just need to write RS_VECTORIZE_1(rsSin) etc. to 
// reduce boilerplate. the 1 is for "unary"

// Binary functions:
// min, max, pow, atan2, fmod, hypot, logN ...

// Ternary functions:
TIV rsClip(V x, V a, V b) { return V(rsClip(x.lo(),a.lo(),b.lo()), rsClip(x.hi(),a.hi(),b.hi())); }

// todo: lerp

#undef V
#undef CV
#undef TIV

//-------------------------------------------------------------------------------------------------

template<class T, int N>
class rsSimdVectorEmulated : public rsSimdVector<T, N>
{

};
// Goal: instantiating this class should always result in an emulated simd processing, regardless
// whether or not explicit specializations are available...i'm not yet sure, if subclassing can 
// achieve this - it's a trial that may result in an error...

//-------------------------------------------------------------------------------------------------
/** Base case of 1-element "vectors" which degenerate to scalars. If no explicit instantiation for
a particular combination of type and vector-size can be found, neither for its half-vectors, 
quarter-vectors, etc., this is the fallback case which the code will eventually compile to, 
implementing the SIMD as full software emulation operating on scalar values.  */

#define V  rsSimdVector<T, 1>
#define CV const V
#define TIV template<class T> inline V

template<class T>
class rsSimdVector<T, 1>
{

public:

  V() {}
  V(T a) { v[0] = a; }

  inline T& operator[](const int i) { rsStaticAssert(i == 0); return v[0]; }
  inline const T& operator[](const int i) const { rsStaticAssert(i == 0); return v[0]; }

  inline V operator+(CV& b) const { return V(v[0]+b.v[0]); }
  inline V operator-(CV& b) const { return V(v[0]-b.v[0]); }
  inline V operator*(CV& b) const { return V(v[0]*b.v[0]); }
  inline V operator/(CV& b) const { return V(v[0]/b.v[0]); }

  inline V operator+(T s) const { return V(v[0] + s); }
  inline V operator-(T s) const { return V(v[0] - s); }
  inline V operator*(T s) const { return V(v[0] * s); }
  inline V operator/(T s) const { return V(v[0] / s); }

  static int getEmulationLevel() { return 0; }

  T v[1];
};

TIV operator+(T s, CV b) { return V(s + b.v[0]); }
TIV operator-(T s, CV b) { return V(s - b.v[0]); }
TIV operator*(T s, CV b) { return V(s * b.v[0]); }
TIV operator/(T s, CV b) { return V(s / b.v[0]); }

// Unary functions (int and float):
TIV rsAbs( V x) { return V(rsAbs( x.v[0])); }
TIV rsSign(V x) { return V(rsSign(x.v[0])); }

// Unary functions (float only):
TIV rsFloor(V x) { return V(rsFloor(x.v[0])); }
TIV rsCeil( V x) { return V(rsCeil( x.v[0])); }
TIV rsRound(V x) { return V(rsRound(x.v[0])); }

TIV rsExp(  V x) { return V(rsExp(  x.v[0])); }
TIV rsLog(  V x) { return V(rsLog(  x.v[0])); }
TIV rsSqrt( V x) { return V(rsSqrt( x.v[0])); }

TIV rsSin(  V x) { return V(rsSin(  x.v[0])); }
TIV rsCos(  V x) { return V(rsCos(  x.v[0])); }
TIV rsTan(  V x) { return V(rsTan(  x.v[0])); }

TIV rsCosh(  V x) { return V(rsCosh(  x.v[0])); }
TIV rsSinh(  V x) { return V(rsSinh(  x.v[0])); }
TIV rsTanh(  V x) { return V(rsTanh(  x.v[0])); }

// Ternary functions:
TIV rsClip(V x, V a, V b) { return V(rsClip(x.v[0], a.v[0], b.v[0])); }

#undef V
#undef CV
#undef TIV

// Macros for function definitions:
//#define RS_VECTORIZE_FUNCTION_1(f) 


//=================================================================================================
//-------------------------------------------------------------------------------------------------
/** Explicit specialization for a vector of 4 floats. */

#ifdef RS_USE_SSE
#define V  rsSimdVector<float, 4>
#define CV const V

template<>
class rsSimdVector<float, 4>
{

public:

  using V2  = rsSimdVector<float, 2>;
  using CV2 = const V2;

  static int getEmulationLevel() { return 0; }

  // Constructors:
  V() {}
  V(CV& a) : v(a.v) {}
  //V(CV&& a) : v(a.v) {} // gives compile error - we need to also define the = operator then
  V(__m128 x) : v(x) {}
  V(float a) : v(_mm_set1_ps(a)) {}
  V(int a) : v(_mm_set1_ps(float(a))) {}
  V(double a) : v(_mm_set1_ps(float(a))) {}
  V(float a, float b, float c, float d) : v(_mm_setr_ps(a, b, c, d)) {}
  V(float* p) { v = _mm_setr_ps(p[0], p[1], p[2], p[3]); }
  V(CV2& low, CV2& high) { lo() = low; hi() = high; }

  // Setup:
  void set(float a, float b, float c, float d) { v = _mm_setr_ps(a, b, c, d); }
    // has no template yet...i think, the template needs to be variadic



  //// Arithmetic operators
  //V operator+(CV& w) const { return V(_mm_add_ps(v, w.v)); }
  //V operator-(CV& w) const { return V(_mm_sub_ps(v, w.v)); }
  //V operator*(CV& w) const { return V(_mm_mul_ps(v, w.v)); }
  //V operator/(CV& w) const { return V(_mm_div_ps(v, w.v)); }
  // moved outside the class

  // Access operators
  V& operator=(const __m128& rhs) { v = rhs; return *this; }
  float& operator[](const int i) const { return asArray()[i]; }

  // Ugly boilerplate, because we need lo()/hi() functions here, too:
  float* asArray() const { return (float*) &v; }                 // helper
  inline V2& lo() { return *((V2*) &(asArray()[0])); } 
  inline V2& hi() { return *((V2*) &(asArray()[2])); }           // 2 = 4/2 = N/2
  inline CV2& lo() const { return *((CV2*) &(asArray()[0])); }
  inline CV2& hi() const { return *((CV2*) &(asArray()[2])); }

  // The data:
  __m128 v;
};

// Arithmetic operators:
inline V operator+(float a, CV b) { return V(_mm_add_ps(_mm_set1_ps(a), b.v)); }
inline V operator-(float a, CV b) { return V(_mm_sub_ps(_mm_set1_ps(a), b.v)); }
inline V operator*(float a, CV b) { return V(_mm_mul_ps(_mm_set1_ps(a), b.v)); }
inline V operator/(float a, CV b) { return V(_mm_div_ps(_mm_set1_ps(a), b.v)); }
inline V operator+(CV a, CV b) { return V(_mm_add_ps(a.v, b.v)); }
inline V operator-(CV a, CV b) { return V(_mm_sub_ps(a.v, b.v)); }
inline V operator*(CV a, CV b) { return V(_mm_mul_ps(a.v, b.v)); }
inline V operator/(CV a, CV b) { return V(_mm_div_ps(a.v, b.v)); }
inline V operator+(const V a) { return a; }          // unary plus
inline V operator-(const V a) { return V(0.f) - a; } // unary minus - can we do better?
// We need to define them outside the class to enable automatic type conversion when the left 
// argument is a scalar
// Passing arguments by const reference has given (very weird!) access violations in the unit test.
// -> check out other simd libraries, how they pass arguments and check, if passing by value incurs
// a performance hit (done - doesn't seem to make a difference), check all other operators, maybe 
// switch to pass-by-value ther, too

#undef V
#undef CV
#endif

// ToDo: rsSimdVector<double, 2> (needs SSE2)

// copied form Float64x2 in rosic:
// for more inspiration, see:
// https://msdn.microsoft.com/de-de/library/tyy88x2a(v=vs.90).aspx
// https://msdn.microsoft.com/de-de/library/9b07190d(v=vs.90).aspx
// http://johanmabille.github.io/blog/2014/10/10/writing-c-plus-plus-wrappers-for-simd-intrinsics-3/
// https://github.com/p12tic/libsimdpp  SIMD library for x86, ARM and more
// https://github.com/VcDevel/Vc  ditto
// https://github.com/agenium-scale/nsimd

//
// https://github.com/simd-everywhere/simde

// for evaluating elementary math functions without resorting to the scalar versions, see:
// http://ito-lab.naist.jp/~n-sibata/pdfs/isc10simd.pdf

// reference:
// http://www.info.univ-angers.fr/pub/richer/ens/l3info/ao/intel_intrinsics.pdf

// ToDo: implement a fallback version to be used, if SSE2 is not available and also a special ARM
// version (needed for M1 processor, i guess):
// https://docs.microsoft.com/en-us/cpp/intrinsics/arm-intrinsics?view=msvc-160

// Maybe turn this into a SIMD library, split it into several files:
// -the base implementation with the fallback to scalar code
// -one file for each explicit specializations, named e.g. Float32x4_SSE, Float32x4_NEON, 
//  Float32x8_AVX, etc.
// -maybe create a template for new specializations that contains dummy instructions like:
//  _simd_add, _simd_set, etc. and a script that reads in the template file and a second file that
//  contains translations for the dummy instructions for a particular instruction set. _simd_add
//  would be translated to _mm_add_ps, when the SSE translation file is used, etc. - it's a sort of
//  code generator
// -name ideas: AMSAP: as much simd as possible, SIP: simd if possible, SIA: simd if available,
//  SIPER: simd if possible else recurse

// Also interesting:
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p1101r0.html Vector Length Agnostic SIMD