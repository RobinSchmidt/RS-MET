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
// Maybe move to RAPT. There actually already is such a function: rsSetZero(). So maybe get rid of
// this one.

template<class T>
T square(T x)
{
  return x*x;
}

template<class T>
void rsFillWithRandomValues(T* x, size_t N, T min, T max, unsigned long seed = 0)
{
  RAPT::rsNoiseGenerator<T> prng;
  prng.setRange(min, max);
  prng.setSeed(seed);
  for(size_t n = 0; n < N; n++)
    x[n] = prng.getSample();
}


/** Fills the array z of complex values with random values. Can be used with std::complex and 
RAPT::rsComplex and any other complex number type for which a function rsSetComplex() is suitably
defined (see implementation for how such a function needs to look like). */
template<class TReal, class TComplex>
void rsFillWithComplexRandomValues(
  TComplex* z, size_t N, TReal min, TReal max, unsigned long seed = 0)
{
  RAPT::rsNoiseGenerator<TReal> prng;
  prng.setRange(min, max);
  prng.setSeed(seed);
  for(size_t n = 0; n < N; n++)
  {
    TReal re = prng.getSample();
    TReal im = prng.getSample();
    rsSetComplex(&z[n], re, im);
  }
}

// ToDo: Adapt the functions below so they can also work with rsComplex. Don't hardcode usage of 
// std::complex into them. Use a template parameter TComplex instead.

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
// that every number is contained exactly once, so b is a permutation of the numbers 0...L-1. We do
// not need to use a containsOnce() function (which would be even more expensive).
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

// Produces spectral magnitudes of given signal x.
template<class T>
std::vector<T> rsSpectralMagnitudes(const std::vector<T>& x)
{
  int N = (int)x.size();
  rsAssert(rsIsPowerOfTwo(N), "This function currently only works for powers of 2." );

  // Create and set up an FFT object:
  using FFT = rsFourierTransformerRadix2<T>;
  FFT fft;
  fft.setBlockSize(N);
  fft.setDirection(FFT::directions::FORWARD);
  fft.setNormalizationMode(FFT::normalizationModes::NEVER_NORMALIZE);

  std::vector<T> mags(N/2), phases(N/2);
  fft.getRealSignalMagnitudesAndPhases(&x[0], &mags[0], &phases[0]);

  return mags;
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

// Tests if a given impulse response h is allpass in nature - up to some tolerance because it's
// necessarily truncated to finite length.
template<class T>
bool isAllpass(const std::vector<T>& h, T tol)
{
  // Compute magnitude spectrum:
  std::vector<T> mags = rsSpectralMagnitudes(h);

  // Check if the maximum deviation from unit frequency response is within the tolerance:
  T maxErr = T(0);
  for(size_t k = 0; k < mags.size(); k++)
    maxErr = rsMax(maxErr, rsAbs(T(1) - mags[k]));
  bool ok = maxErr <= tol;

  //// This can be uncommented in debug sessions to investigate problems when the test fails:
  //if(!ok)
  //{
  //  rsError("Filter is not allpass!");
  //  rsPlotVectors(mags);  
  //}

  return ok;
}
// rename to rsIsAllpass()


/** Tests if the given signal x is a unit impulse. This can be useful when checking pairs of 
filters that are supposed to be inverses of one another. Applying thme both to a unit impulse 
should give back that unit impulse. There may be other use cases for such checks but that happens
to be the one for which I wrote it */
template<class T>
bool rsIsUnitImpulse(const std::vector<T>& x, T tol)
{
  rsAssert(x.size() > 0);
  T maxErr = rsAbs(x[0] - T(1));
  for(size_t n = 1; n < x.size(); n++)
    maxErr = rsMax(maxErr, rsAbs(x[n]));
  return maxErr <= tol;
}

template<class T>
bool rsIsShiftedUnitImpulse(const std::vector<T>& x, int shift, T tol)
{
  bool ok = true;
  for(int n = 0; n < (int) x.size(); n++)
  {
    if(n == shift)
      ok &= rsIsCloseTo(x[n], 1.0, tol);
    else
      ok &= rsIsCloseTo(x[n], 0.0, tol);
  }
  return ok;
}




template<class TSig, class TFlt>
TSig rsGetSample(TFlt& filter, TSig in)
{
  return filter.getSample(in);
}


template<class TSig, class TPar>
TSig rsGetSample(rsProtoFDN<TSig, TPar>& fdn, TSig in)
{
  std::vector<TSig> vIn(1), vOut(1);
  // ToDo: Use getNumIn/OutputChannels like so:
  //std::vector<TSig> vIn(fdn.getNumInputChannels()), vOut(fdn.getNumOuputChannels()); 

  vIn[0] = in;
  fdn.processFrame(vIn, vOut);
  return vOut[0];
}

template<class T>
T rsGetSample(rsStateSpaceFilter<T>& ssf, T in)
{
  //std::vector<T> vIn(1), vOut(1);
  std::vector<T> vIn(ssf.getNumInputs()), vOut(ssf.getNumOutputs());
  vIn[0] = in;
  ssf.processFrame(&vIn[0], &vOut[0]);
  return vOut[0];
}

/** Returns N samples of the impulse response of the passed filter as std::vector. It is necessary
for you to pass a scale factor of the type of the filter's output signal (for example: 1.0 for
double), such that the compiler can deduce the template parameter. We also use it to scale the
input impulse to the filter, so it actually gets some purpose besides satisfying the compiler. The
filter class must support the functions reset() and getSample(). We basically use duck typing here.
Any object that has these two member functions with the correct signature (and hopefully also the
correct semantics) can be passed as filter. */
template<class TSig, class TFlt>
inline std::vector<TSig> impulseResponse(TFlt& filter, int length, TSig scale)
{
  std::vector<TSig> y(length);
  filter.reset();

  //// Old:
  //y[0] = filter.getSample(scale);
  //for(int n = 1; n < length; n++)
  //  y[n] = filter.getSample(TSig(0));

  // New:
  y[0] = rsGetSample(filter, scale);
  for(int n = 1; n < length; n++)
    y[n] = rsGetSample(filter, TSig(0));

  return y;
}

/** Returns the response y[n] of the given filter to the input signal x[n]. */
template<class TSig, class TFlt>
inline std::vector<TSig> filterResponse(TFlt& filter, int length, std::vector<TSig> x)
{
  std::vector<TSig> y(length);
  filter.reset();
  for(int n = 0; n < length; n++)
    y[n] = filter.getSample(x[n]);
  return y;
}
// ToDo: Document the length parameter. It's for allowing y to be shorter than x. But why? It would 
// actually make more sense to allow y to be longer to allow for a ringout. This old documentation 
// text is wrong:
// The length determines the length of y which may be different from the length of x to allow the 
// filter to ring out, for example. WAIT - NO - this is wrong! The code below does not actually 
// implement a ringout phase! I think, the "length" parameter currently only serves to allow y to 
// be shorter than x. But why did I do it like that? Maybe change the implementation such that y
// may be shorter or longer than x. Or mayb just get rid of the length parameter and take x.size()
// for the length of y.



/** Computes the transfer function H(z) at the given z the hard way, i.e. as the z-transform of the 
impulse response. It's only an approximation though because we truncate the infinite sum at/ N-1. 
N should be large enough such that the impulse response has sufficiently decayed at the end. We 
compute  Ht = sum_{n=0}^{N-1} h[n] * z^{-n}  where h[n] is the impulse response of the filter. In 
the actual z-trafo, the upper limit of the sum would be infinity. */
template<class T, class TFlt>
inline rsComplex<T> rsEvaluateTransferFunctionNumerically(TFlt& filter, rsComplex<T> z, int N)
{
  filter.reset();
  rsComplex<T> z1 = T(1) / z;                         // z1 = z^-1
  rsComplex<T> zn = T(1);                             // zn = z^-n, initialized with n = 0
  //rsComplex<T> H  = filter.getSample(T(1)) * zn;      // H: acccumulator for H(z)
  rsComplex<T> H  = rsGetSample(filter, T(1)) * zn;     // H: acccumulator for H(z)
  for(int n = 1; n < N; n++)
  {
    zn *= z1;                                         // z^-n
    //H  += filter.getSample(T(0)) * zn;
    H  += rsGetSample(filter, T(0)) * zn;
  }
  return H;
}
// Maybe take and return a T rather than rsComplex<T>. Then T itself already will be the complex 
// type in instantiations. That makes it more flexible. Ah! Nope! That doesn't compile because
// filter.getSample() usually doesn't take a complex argument.

/** Helper function to test the result of filter.getTransferFunctionAt() against a naively computed
transfer function value. This is meant for unit testing the getTransferFunctionAt() member function
that I typically give to many of my filter classes. */
template<class T, class TFlt>
inline bool rsTestGetTransferFunctionAt(TFlt& filter, rsComplex<T> z, int N, T tol)
{
  //  Compute transfer function H(z) at the given z numerically:
  rsComplex<T> Ht = rsEvaluateTransferFunctionNumerically(filter, z, N);

  // Compute the transfer function using the filter's getTransferFunctionAt() method:
  rsComplex<T> H = filter.getTransferFunctionAt(z);

  // Compute the error and check if its absolute value is within the tolerance:
  rsComplex<T> err = H - Ht;
  T errAbs = rsAbs(err);
  return errAbs <= tol;
}

template<class T, class TFlt>
inline bool rsTestGetTransferFunction(TFlt& filter, rsComplex<T> z, int N, T tol)
{
  //  Compute transfer function H(z) at the given z numerically:
  rsComplex<T> Ht = rsEvaluateTransferFunctionNumerically(filter, z, N);

  // Let the filter produce its transfer function object via getTransferFunction() and evaluate the
  // returned function object at z:
  //rsSparseTransferFunction<T> H;   // Old - when the class had only one template parameter
  rsSparseTransferFunction<T, T> H;  // New - not sure if using T for TTol is always appropriate
  filter.getTransferFunction(&H);
  rsComplex<T> Hz = H(z);

  // Compute the error and check if its absolute value is within the tolerance:
  rsComplex<T> err = Hz - Ht;
  T errAbs = rsAbs(err);
  return errAbs <= tol;
}




template<class T>
inline std::vector<T> ampToDb(const std::vector<T>& x, T minDb)
{
  std::vector<T> y(x.size());
  for(size_t i = 0; i < x.size(); i++)
    y[i] = rsMax(rsAmpToDb(rsAbs(x[i])), minDb);
  return y;
}



//=================================================================================================
// Filtering

template<class T>
void rsBiquadResponse(const T* x, T* y, int N, T b0, T b1, T b2, T a1, T a2)
{
  RAPT::rsAssert(x != y, "Not suitable for in-place computation");

  if(N < 1) 
    return;
  y[0] = b0 * x[0];

  if(N < 2) 
    return;
  y[1] = b0 * x[1] + b1 * x[0] - a1 * y[0];

  for(int n = 2; n < N; n++)
    y[n] = b0 * x[n] + b1 * x[n-1] + b2 * x[n-2] - a1 * y[n-1] - a2 * y[n-2];
}
// Maybe move this into RAPT::rsArrayTools, maybe make a variant that allows in-place computation
// maybe distiguish it from this one by having only a single in/out array

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
  std::vector<bool> done(N);           // Flags to mark values that were already "used up"
  for(int i = 0; i < N; i++)  {
    int j;
    for(j = 0; j < N; j++)
      if(!done[j] && rsAbs(x[i]-y[j]) <= tol)
        break;
    if(j == N)
      return false;
    done[j] = true;  }
  return true;

  // ToDo: Document the algorithm. I think, we iterate through all entries of x and for each, we 
  // try to find a matching entry in y. But we consider only those elements of y that were not yet
  // "used up" by a previous iteration. That lets the algorithm do the right thing also in case on
  // arrays in which the same entry may occurr more than once - I think. Verify that! It was just 
  // my recollection. Check, if we have unit tests for this function, too. If not, write some.
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


// This is a kludge. Get rid! We need it to convert from rsMatrix into a format the GNUPlotter
// understands, i.e. matrices defined by a pointer-to-pointer structure. The best would be, if
// GNUPlotter would support to take the matrix in flat storage format.
template<class T>
T** createRowPointers(RAPT::rsMatrix<T>& M)
{
  int N = M.getNumRows();
  T** rp = new T*[N];
  for(int i = 0; i < N; i++)
    rp[i] = M.getRowPointer(i);
  return rp;
}
template<class T>
void deleteRowPointers(T** rowPointers, const RAPT::rsMatrix<T>& M)
{
  // New:
  delete[] rowPointers;

  // Old - wrong:
  //for(int i = 0; i < M.getNumRows(); i++)
  //  delete rowPointers[i];

  // This needs some testing!
}
// ToDo: maybe use rs-prefix





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
// Convenience functions for polynomials:

/** Randomizes the values of the coefficients of the polynomial p. Leaves the degree as is. */
template<class T>
void randomizeCoeffs(rsPolynomial<T>* p, T min, T max, int seed, bool roundToInt = false)
{
  rsNoiseGenerator<T> prng;
  prng.setRange(min, max);
  prng.setSeed(seed);
  for(int i = 0; i <= p->getDegree(); i++)  // <= is not a bug, numCoeffs is degree + 1
  {
    T c = prng.getSample(); 
    if(roundToInt)
      c = rsRound(c);
    p->setCoeff(i, c);
  }
}
// ToDo:
// -Check what happens when T is a complex type. It probably won't compile. Maybe make a version
//  that can handle complex polynomials as well. Maybe an explicit specialization for complex 
//  polynomials is needed?

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

/** Fills the outL, outR arrays with the output of the given sampler engine for the given sequence
of musicla events (note-ons, note-offs, control-changes, etc.). */
void getSamplerOutput(rosic::Sampler::rsSamplerEngine* se, 
  const std::vector<rosic::Sampler::rsMusicalEvent<float>>& events,
  float* outL, float* outR, int numFrames);
// needs tests

/** Generates the samples that are used in the test patches. */
void generateTestSamples();


//=================================================================================================
// Printing:




//=================================================================================================
// Conversions between different enums

int convertEnumMode_RBJ_to_SVF(int rbjMode);



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
