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
software, completely transparent to the client code. 

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

  /** Write Access to the i-th element of the vector where i = 0,...,N-1. */
  inline T& operator[](const int i) { rsStaticAssert(i >= 0 && i < N); return v[i]; }

  /** Read Access to the i-th element of the vector where i = 0,...,N-1. */
  inline const T& operator[](const int i) const { rsStaticAssert(i >= 0 && i < N); return v[i]; }


  static bool isSimdEmulated() { return true; }
  // explicit instantiations that actually do hardware simd, should return false

  //int getEmulationLevel() { }
  // idea: the function should return an information, at which level the actual hardware simd
  // kicks in. 0: full hardware implementation. 1: the half-vectors use hardware simd, 
  // 2: quarter-vectors use hardware simd, etc., i'm not yet sure, how to implement this


  V operator+(CV& w) const { V u; u.lo() = lo() + w.lo(); u.hi() = hi() + w.hi(); return u; }
  V operator-(CV& w) const { V u; u.lo() = lo() - w.lo(); u.hi() = hi() - w.hi(); return u; }
  V operator*(CV& w) const { V u; u.lo() = lo() * w.lo(); u.hi() = hi() * w.hi(); return u; }
  V operator/(CV& w) const { V u; u.lo() = lo() / w.lo(); u.hi() = hi() / w.hi(); return u; }


private:

  /** Returns a reference to the lower half-vector. */
  V2& lo() { return *((V2*) &v[0]);  }

  /** Returns a reference to the higher half-vector. */
  V2& hi() { return *((V2*) &v[N/2]); }

  /** Returns a const reference to the lower half-vector. */
  CV2& lo() const { return *((CV2*) &v[0]);   }

  /** Returns a const reference to the higher half-vector. */
  CV2& hi() const { return *((CV2*) &v[N/2]); }


  T v[N]; // maybe we should specify an alignment?
};

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

template<class T>
class rsSimdVector<T, 1>
{

public:

  using V  = rsSimdVector<T, 1>;
  using CV = const V;

  V operator+(CV& w) const { V u; u.v[0] = v[0] + w.v[0]; return u; }
  V operator-(CV& w) const { V u; u.v[0] = v[0] - w.v[0]; return u; }
  V operator*(CV& w) const { V u; u.v[0] = v[0] * w.v[0]; return u; }
  V operator/(CV& w) const { V u; u.v[0] = v[0] / w.v[0]; return u; }


private:

  T v[1];

};

//-------------------------------------------------------------------------------------------------
/** Explicit specialization for a vector of 4 floats. */

#ifdef RS_USE_SSE
template<>
class rsSimdVector<float, 4>
{

public:

  using V  = rsSimdVector<float, 4>;
  using CV = const V;

  static bool isSimdEmulated() { return false; }

  // Constructors:
  V() {};
  V(__m128 x) : v(x) {}
  V(float a) : v(_mm_set1_ps(a)) {}
  V(int a) : v(_mm_set1_ps(float(a))) {}
  V(double a) : v(_mm_set1_ps(float(a))) {}
  V(float a, float b, float c, float d) : v(_mm_setr_ps(a, b, c, d)) {}
  V(float* p) { v = _mm_setr_ps(p[0], p[1], p[2], p[3]); }


  //// Arithmetic operators
  //V operator+(CV& w) const { return V(_mm_add_ps(v, w.v)); }
  //V operator-(CV& w) const { return V(_mm_sub_ps(v, w.v)); }
  //V operator*(CV& w) const { return V(_mm_mul_ps(v, w.v)); }
  //V operator/(CV& w) const { return V(_mm_div_ps(v, w.v)); }
  // moved outside the class

  // Access operators
  V& operator=(const __m128& rhs) { v = rhs; return *this; }
  float& operator[](const int i) const { return asArray()[i]; }


  __m128 v;

private:

  float* asArray() const { return (float*) &v; }

};

#define V  rsSimdVector<float, 4>
#define CV const V

// Arithmetic operators
inline V operator+(CV& a, CV& b) { return V(_mm_add_ps(a.v, b.v)); }
inline V operator-(CV& a, CV& b) { return V(_mm_sub_ps(a.v, b.v)); }
inline V operator*(CV& a, CV& b) { return V(_mm_mul_ps(a.v, b.v)); }
inline V operator/(CV& a, CV& b) { return V(_mm_div_ps(a.v, b.v)); }
inline V operator+(const V& a) { return a; }          // unary plus
inline V operator-(const V& a) { return V(0.f) - a; } // unary minus - can we do better?
// We need to define them outside the class to enable automatic type conversion when the left 
// argument is a scalar

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
