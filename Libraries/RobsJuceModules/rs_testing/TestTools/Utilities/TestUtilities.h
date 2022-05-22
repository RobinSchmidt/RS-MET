#ifndef RS_TESTUTILITIES_H
#define RS_TESTUTILITIES_H

// todo: merge file with other utility files

//old:
//#include "../Common/Prototypes.h"

// new:
//#include "../RaptLibraryCode/RaptInstantiations.h"
#include "rosic/rosic.h"
#include "rs_testing/rs_testing.h"



bool runUnitTest(bool (*test)(), const std::string& name);



//bool detectMemoryLeaks();  // currently works only in MSVC

/** This function should be called on program startup when automatic detection of memory leaks
should be turned on. */
inline void checkForMemoryLeaksOnExit()
{
#if defined _MSC_VER
  int tmpFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG); // gets the current flag
  tmpFlag |= _CRTDBG_LEAK_CHECK_DF;                  // turns on leak checking
  //tmpFlag &= ~_CRTDBG_CHECK_CRT_DF;                  // turns off CRT block checking bit
  _CrtSetDbgFlag(tmpFlag);                           // set flag to the new value;
#endif
}

// helper functions to create some vectors useful for testing purposes (maybe move them to
// somewhere else):
std::vector<double> rsLinearRangeVector(     int N, double min, double max);
std::vector<double> rsExponentialRangeVector(int N, double min, double max);
std::vector<double> rsRandomVector(          int N, double min, double max, int seed = 0);
std::vector<double> rsRandomIntVector(       int N, int    min, int    max, int seed = 0);
std::vector<double> rsApplyFunction(const std::vector<double>& v, double p,
  double (*f) (double, double));

// conversions to std::string:
std::string toString(int n);

// replace with own prng:
inline double random(double min, double max)
{
  double tmp = (1.0/RAND_MAX) * rand();  // between 0...1
  return RAPT::rsLinToLin(tmp, 0.0, 1.0, min, max);
}

// returns x^2 = x*x, useful for testing application of a unary function using a function pointer
//double rsSquare(double x);

/** Fills the given vector a with all zeros. */
template<class T>
void rsZero(std::vector<T>& a) 
{ 
  RAPT::rsFill(a, T(0));
}
// maybe move to RAPT

//
template<class T>
T square(T x)
{
  return x*x;
}

template<class T>
void rsFillWithComplexRandomValues(std::complex<T>* x, size_t N, T min, T max,
  unsigned long seed = 0)
{
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(min, max);
  prng.setSeed(seed);
  for(size_t n = 0; n < N; n++)
    x[n] = std::complex<T>(prng.getSample(), prng.getSample());
}

template<class T> // convenience function for std::vector
void rsFillWithComplexRandomValues(std::vector<std::complex<T>>& x, T min, T max,
  unsigned long seed = 0)
{
  rsFillWithComplexRandomValues(&x[0], x.size(), min, max, seed);
}
template<class T> // yet more convenient function
std::vector<std::complex<T>> rsComplexRandomVector(int N, T min, T max, unsigned long seed = 0)
{
  std::vector<std::complex<T>> x(N);
  rsFillWithComplexRandomValues(x, min, max, seed);
  return x;
}

template<class T>
T rsMaxComplexError(std::complex<T>* target, std::complex<T>* actual, size_t N)
{
  T maxErr = T(0);
  for(size_t n = 0; n < N; n++)
    maxErr = RAPT::rsMax(maxErr, abs(target[n]-actual[n]));
  return maxErr;
}

template<class T>
bool rsAlmostEqual(std::vector<std::complex<T>>& x, std::vector<std::complex<T>>& y, T tolerance)
{
  RAPT::rsAssert(x.size() == y.size());
  T maxErr = rsMaxComplexError(&x[0], &y[0], x.size());
  return maxErr <= tolerance;
}

inline bool isIndexPermutation(int* b, int L)
{
  for(int i = 0; i < L; i++)
    if( !rsArrayTools::contains(b, L, i) )
      return false;
  return true;
}
// Returns true, iff b contains every number from 0 to L-1. Since b is of length L, this implies
// that every number is contained exactly once, so b is a permutation of the numbers 0...L-1.
// used here only for debug -  move elsewhere

/** Applies the inner function to the value x and then the outer function to the result of that
inner function and returns the final result. This is known as function composition in
mathematics. */
template<class T, class F1, class F2>
T applyComposedFunction(T x, F1 innerFunction, F2 outerFunction)
{
  return outerFunction(innerFunction(x));
}

/** Returns true, if the function f maps the given argument x to itself. */
template<class T, class F>
bool mapsToItself(T x, F f)
{
  return x == f(x);
}

/** Checks, if the 2nd function is the inverse function of the first for the given input argument
x. */
template<class T, class F1, class F2>
bool mapsBack(T x, F1 forwardFunction, F2 maybeInverseFunction)
{
  return x == applyComposedFunction(x, forwardFunction, maybeInverseFunction);
  //return x == maybeInverseFunction(forwardFunction(x));
}
// maybe rename to isFunctionLocallyInverse, isFunctionInverseAt

/** Checks, if the 2nd function is the inverse function of the first for a given range of input
arguments between minValue and maxValue with given. */
template<class T, class F1, class F2>
bool isInverseFunction(F1 forwardFunc, F2 maybeInverseFunc, T minValue, T maxValue, T increment)
{
  T value = minValue;
  while(value < maxValue) {
    if( !mapsBack(value, forwardFunc, maybeInverseFunc) )
      return false;
    value += increment;
  }
  return true;
}





// for testing the callback performance (this is actually in jura, but anyway):
#define JUCE_API
#include "jura_framework/control/jura_Callbacks.h"

inline bool detectMemoryLeaks()
{
#ifdef _MSC_VER
  return _CrtDumpMemoryLeaks() == 1;
#else
  return false;
#endif
}

// get rid of that:
inline void appendTestResultToReport(std::string &reportString, const std::string &nameOfTest,
  bool result)
{
  if( result == true )
    reportString += nameOfTest + ": OK \n";
  else
    reportString += nameOfTest + ": !!! FAILED !!!\n";
}

/** Comparison function that compares with a given error tolerance and also returns true when the
involved numbers are NaNs or infinities. */
bool areNumbersEqual(double x, double y, double relativeTolerance);

/** Convenience function to convert a string to a window-type.  options: rc,hn,hm,bm,bh */
RAPT::rsWindowFunction::WindowType stringToWindowType(const std::string& wt);

//=================================================================================================
// Convenience functions for vectors:

template<class T>
void rsNormalizeChunks(std::vector<T>& v, int chunkSize)
{
  using AT = rsArrayTools;
  int N = (int) v.size();
  rsAssert(N % chunkSize == 0);  // v must contain integer number of chunks
  int numChunks = N / chunkSize;
  for(int i = 0; i < numChunks; i++) {
    int j = i*chunkSize;
    T norm = sqrt(AT::sumOfSquares(&v[j], chunkSize));
    AT::scale(&v[j], chunkSize, T(1)/norm);  }
}

/** Extracts a chunk starting at "start" given "size" from vector v. */
template<class T>
std::vector<T> rsGetChunk(std::vector<T>& v, int start, int size)
{
  rsAssert(start >= 0 && start+size <= (int) v.size());
  std::vector<T> c(size);
  for(int i = 0; i < size; i++)
    c[i] = v[start+i];
  return c;
}

/** Checks if vector y is a permutation of x, up to some tolerance. */
template<class T>
bool rsIsPermutation(const std::vector<T>& x, const std::vector<T>& y, T tol)
{
  int N  = (int) x.size();
  if((int) y.size() != N)
    return false;
  std::vector<bool> done(N);   // flags to indicate that a value was used up
  for(int i = 0; i < N; i++)  {
    int j;
    for(j = 0; j < N; j++)
      if(!done[j] && rsAbs(x[i]-y[j]) <= tol)
        break;
    if(j == N)
      return false;
    done[j] = true;  }
  return true;
}

//=================================================================================================

/** Under construction.
A subclass of std::vector that can be used as drop-in replacement for it and performs some
additional logging of certain member function calls by intercepting them (via overriding), doing
the logging and then just forwarding the request to the baseclass implementation. It can be used to
verify that return value optimzation works as intended for classes that uses std::vector for the 
underlying data storage by also templatizing them on the vector type to use and pass this type here 
in the respective unit tests. Note that the overriding only gives compile-time polymorphism because
the overriden functions are not virtual in the baseclass - but this is good enough for the intended
purpose. */

template<class T>
class rsLoggingVector : public std::vector<T>
{

public:

  using Base = std::vector<T>;

  rsLoggingVector() : Base() 
  {
    // The standard constructor never allocates
  }

  rsLoggingVector(size_t s) : Base(s) 
  {
    // Constructing with a size may potentially allocate and will actually do so, iff size > 0
    numPotentialAllocs++;
    if(s > 0)
      numActualAllocs++;
  }
  // do we need to templatize on the size-type here like we do in resize()?

  rsLoggingVector(std::initializer_list<T> l) : Base(l) 
  {
    // Constructing with an initializer list may potentially allocate and will actually do so, 
    // iff the size of the list is > 0
    numPotentialAllocs++;
    if(l.size() > 0)
      numActualAllocs++;
  }

  template<class S> // S: size-type - we want to catch all calls with int, size_t, etc.
  void resize(S newSize) 
  { 
    numPotentialAllocs++;
    if(newSize > capacity())
      numActualAllocs++;
    Base::resize(newSize); 
  }

  template<class S> 
  void resize(S newSize, const T& val)
  {
    numPotentialAllocs++;
    if(newSize > capacity())
      numActualAllocs++;
    Base::resize(newSize, val);
  }
  // what does this val variable do? it doesn't seem to initialize all elements to val


  /*
  // std::vector::resize looks like this in MSVC:

  void resize(_CRT_GUARDOVERFLOW const size_type _Newsize) {
  // trim or append value-initialized elements, provide strong guarantee
  _Resize(_Newsize, _Value_init_tag{});
  }

  void resize(_CRT_GUARDOVERFLOW const size_type _Newsize, const _Ty& _Val) {
  // trim or append copies of _Val, provide strong guarantee
  _Resize(_Newsize, _Val);
  }
  */

  //void resize(size_t newSize) { numResizeCalls++; Base::resize(newSize); }
  //void resize(int    newSize) { numResizeCalls++; Base::resize(newSize); }

  // ToDo:
  // -log also calls to: reserve, shrink_to_fit, copy-assign, copy-construct, etc.
  //  ...anything that potentially (re)-allocates
  // -maybe also log move-assign/construct calls
  // -check that we override all relevant overloads (for size_t, int, const, etc.) to make
  //  sure that we really intercept all such calls.


  // These logging values are static for two reasons: (1) We want to be able to globally access 
  // them from the test code because the to-be-tested classes may not expose any interface to their 
  // underlying data storage vector. (2) We are actually interested in the number of resize calls on
  // all vectors combined - not just those on a specific one.

  //static size_t numResizeCalls;      // number of calls to resize()

  static size_t numPotentialAllocs;  
  // Number of potential (re)-allocations. This counts all calls to resize, reserve, shrink_to_fit,
  // constructors that may allocate, etc. regardless whether or not a re-allocation actually does
  // occur. For example, resizing below the capacity won't allocate but we'll count it as potential 
  // alloc anyway. That makes sense because allocation is the worst case behavior which is what we 
  // are usually interested in.

  static size_t numActualAllocs;
  // Number of actual (re)-allocations - std::vector reallocates if and only if the capacity grows
  // ..i think...hmm - but it may also re-allocate on shrink_to_fit ...maybe rename to 
  // numLikelyAllocs. We cannot know for sure because the re-allcoation behavior is not completely
  // specified in all cases: shrink_to_fit may be ignored. Maybe have a 3rd variable 
  // numCertainAllocs which counts all resize/reserve calls above capacity. Perhaps we can 
  // know the number of actual allocations, if we use a custom (logging) allocator that can itself 
  // be inquired for the number of allocations, it did. ....hmmm....

};

//template<class T> size_t rsLoggingVector<T>::numResizeCalls = 0;
template<class T> size_t rsLoggingVector<T>::numPotentialAllocs = 0;
template<class T> size_t rsLoggingVector<T>::numActualAllocs = 0;


//=================================================================================================
// Convenience functions for matrices:

template<class T>
std::vector<T> rsRowsToVector(const rsMatrix<T>& A)
{
  int M = A.getNumRows();
  int N = A.getNumColumns();
  std::vector<T> v(M*N);
  for(int m = 0; m < M; m++)
    for(int n = 0; n < N; n++)
      v[m*N + n] = A(m, n);
  return v;
}
template<class T>
std::vector<T> rsColumnsToVector(const rsMatrix<T>& A)
{
  int M = A.getNumRows();
  int N = A.getNumColumns();
  std::vector<T> v(M*N);
  for(int m = 0; m < M; m++)
    for(int n = 0; n < N; n++)
      v[n*M + m] = A(m, n);
  return v;
}
template<class T>
rsMatrix<T> rsToMatrixRowWise(const std::vector<T>& v, int numRows, int numCols)
{
  int M = numRows;
  int N = numCols;
  rsAssert((int)v.size() == M*N);
  rsMatrix<T> A(M, N);
  for(int m = 0; m < M; m++)
    for(int n = 0; n < N; n++) 
      A(m, n) = v[m*N + n];
  return A;
}
template<class T>
rsMatrix<T> rsToMatrixColumnWise(const std::vector<T>& v, int numRows, int numCols)
{
  int M = numRows;
  int N = numCols;
  rsAssert((int)v.size() == M*N);
  rsMatrix<T> A(M, N);
  for(int m = 0; m < M; m++)
    for(int n = 0; n < N; n++) 
      A(m, n) = v[n*M + m];
  return A;
}
// Maybe move into class rsMatrix(View). Maybe the "toVector" functions into matrixView to be
// called like A.toVectorRowWise() and the "toMatrix" functions as static member functions
// in rsMatrix to be called like rsMatrix<float>::toMatrixRowWise()

/** Given an array of N desired eigenvalues and an NxN matrix whose columns are the desired 
eigenvectors, this function returns the NxN matrix A that has this eigensystem. It basically 
computes vecs * diag(vals) * inv(vecs). */
template<class T>
rsMatrix<T> fromEigenSystem(const std::vector<T>& vals, const rsMatrix<T>& vecs)
{
  int N = (int) vals.size();
  rsAssert(vecs.hasShape(N, N));
  rsMatrix<T> tmp(vecs);
  for(int j = 0; j < N; j++)
    tmp.scaleColumn(j, vals[j]);
  return tmp * rsLinearAlgebraNew::inverse(vecs);
}
template<class T>
rsMatrix<T> fromEigenSystem(const std::vector<T>& vals, const std::vector<T>& vecs)
{
  int N = (int) vals.size();
  rsAssert((int)vecs.size() == N*N);
  return fromEigenSystem(vals, rsToMatrixColumnWise(vecs, N, N));
}
// maybe move to Prototypes and maybe into the rsLinearAlgebra

/** This function checks, if the eigenvalues "vals1" and eigenvectors "vecs1" are compatible with 
the ones in the vals2, vecs2 arrays. The equivalence relation between the eigensystems is that 
vals2 must be some permutation of vals1 and the eigenvectors that correspond to a given eigenvalue
must match up to a scalar factor. That means if some eigenvalue k is found at indices i,j in the 
vals1,vals2 arrays respectively, the j-th vector in vecs2 must match the i-th vector in vecs1 up to
scaling. ...and all equalities are taken up to some given tolerance level "tol".  */
template<class T>
bool checkEigensystem(
  const std::vector<T>& vals1, const std::vector<T>& vecs1, 
  const std::vector<T>& vals2, const std::vector<T>& vecs2, T tol)
{
  int N  = (int) vals1.size();
  rsAssert((int) vals2.size() == N);
  rsAssert((int) vecs1.size() == N*N);
  rsAssert((int) vecs2.size() == N*N);

  using ILA  = rsIterativeLinearAlgebra;
  std::vector<bool> done(N);              // flags to indicate that an eigenvalue was used up
  T val;
  bool match;
  int i, j;
  for(i = 0; i < N; i++)
  {
    for(j = 0; j < N; j++)
      if(!done[j] && rsAbs(vals1[i]-vals2[j]) <= tol)
        break;
    if(j == N)
      return false;

    // We have found a correspondence between vals1[i] and vals2[j], so we check, if the 
    // corresponding eigenvectors match up to a scalar factor:
    match = ILA::isScalarMultiple(&vecs1[i*N], &vecs2[j*N], N, tol, &val);
    if(!match)
      return false;
    done[j] = true;
  }
  return true;
}

//=================================================================================================
// Stuff for facilitating tests for the sampler engine

/** Helper function to add a single region for the given sample to the engine. The region is added
to the first group, which is added if not already there. */
void addSingleSampleRegion(rosic::Sampler::rsSamplerEngine* se,
  const std::vector<float>& sample, float keyCenter = 60.f, double sampleRate = 44100);


/** Helper function to set up the sampler engine with a single sinewave region. */
void setupForLoopedWave(rosic::Sampler::rsSamplerEngine* se, int N = 2048, int shape = 0);
// rename to setupForLoopedWave
// Shapes: 0: sine, 1: saw, 2: square, 3: triangle.


/** Sets up the sampler engine with a looped DC region. */
void setupForLoopedDC(rosic::Sampler::rsSamplerEngine* se, int N, float keyCenter, 
  double sampleRate);


/** Fills the outL, outR arrays with the output of the given sampler engine for the given note. The
optional noteOffAt parameter specifies the sample instant at which a note-off is triggered in the
engine, if any. By default, no note-off will be triggered at all. */
void getSamplerNote(rosic::Sampler::rsSamplerEngine* se, float key, float vel,
  std::vector<float>& outL, std::vector<float>& outR, int noteOffAt = -1);

/** Class for representing midi note events for use in some of the sampler tests.  */
struct rsTestNoteEvent
{
  int key    = 0;   // note number in 0..127
  int vel    = 0;   // velocity in 0..127
  int time   = 0;   // in samples
  int length = 0;   // in samples
};

/** Fills the outL, outR arrays with the output of the given sampler engine for the given sequence
of notes. */
void getSamplerNotes(rosic::Sampler::rsSamplerEngine* se, 
  const std::vector<rsTestNoteEvent>& notes,
  std::vector<float>& outL, std::vector<float>& outR);

/** Generates the samples that are used in the test patches. */
void generateTestSamples();


//=================================================================================================

/** Experimental - goal: resemble numpy/scipy/matplotlib functionality, so we may easily port such 
code to C++. */

template<class T>
class rsNumPy
{
public:

  std::vector<T> linspace(T min, T max, int N) { return RAPT::rsRangeLinear(T(min), T(max), N); }

  std::vector<T> sin(const std::vector<T>& x)
  {
    std::vector<T> y(x.size());
    RAPT::rsArrayTools::applyFunction(&x[0], &y[0], (int)x.size(), &::sin);
    return y;
  }

  T pi = PI;




};
// move this into the research repo
// https://docs.scipy.org/doc/numpy-1.11.0/numpy-ref-1.11.0.pdf





#endif
