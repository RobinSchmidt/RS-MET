#ifndef RAPT_ARRAYTOOLS_H_INCLUDED
#define RAPT_ARRAYTOOLS_H_INCLUDED

/** A collection of functions that operate on 1-dimensional arrays. 


todo: 
-sort the functions in the header by purpose (fill, search, sort, arithmetic/combine (add,
 subtract,...), aggregation (sum, mean, product,...), filtering, etc.)
-add to the documentation, to which std::algorithm a function corresponds, if applicable
-maybe rename to rsArrayTools (done), class rsArray should be an actual dynamically allocated array
 with the same interface as std::vector - having my own implementation might be more efficient since
 std::vector initializes the memory - says this blog-post, at least
 https://lemire.me/blog/2012/06/20/do-not-waste-time-with-stl-vectors/
-we could make a baseclass rsArrayView which could also be used as baseclass for rsMatrix and 
 rsMultiArray
-make everything const that is possible (also by-value parameters, local variables, etc. - and use
 constexpr for compile-time constants)
 ->done up to copy
 ...maybe change the const by-value parameters to by-reference parameters
-inline, where it makes sense (trivial functions like copy/convert)
-use workspace pointers instead of heap allocation, where applicable, keep the functions with 
 heap-allocations as convenience functions (create workspace -> call worker -> delete workspace)
-maybe introduce another template parameters for the index-type such that it can be used with int 
 or size_t (the former is used throughout much of my own code but size_t is used by std::vector
 ...this incompatibility actually kinda sucks anyway). But be careful with functions that return an 
 index (like find...) where we use the convention that -1 indicates failure ("not found"). Maybe 
 define constants rsNone(T) that resolves to -1 or whatever else is appropriate...but maybe when we
 just keep using -1 and let it be converted to size_t, it will map to the maximum (that's what 
 2-complement does for unsigned integers, right?) which actually is a reasonable convention to 
 encode "not found" (although, std::find uses v.size() and not max(size_t)) - so it *may* just work 
 fine without any further ado (unit tests would be needed)  */

class rsArrayTools
{

public:

  /** Adds the elements of 'buffer1' and 'buffer2' - type must define operator '+'. The 'result'
  buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T>
  static inline void add(const T *buffer1, const T *buffer2, T *result, const int length);

  /** Adds the scalar 'valueToAdd' to the elements of 'buffer' - the type must define
  operator '+'. The 'result' buffer may be the same as 'buffer'. */
  template <class T>
  static inline void add(const T *buffer, const T valueToAdd, T *result, const int length);

  /** Adds a weighted, circularly shifted copy of the buffer to itself - the shift-offset may be
  non-integer in which case linear interpolation will be used. 
  \todo: maybe generalize such that a circularly shifted 2nd buffer can be added (which may or
  may not be the same buffer).  */
  template <class T>
  static void addCircularShiftedCopy(T *buffer, const int length, const double offset, 
    const T weight);
  // allocates heap memory - todo: use a workspace parameter
  // why is the offset a double and not T? is this a bug? if yes - fix, if not, document why

  /** Adds length-L array y into length-N array x starting at n (in x), taking care of not reading 
  beyond the limits of y and writing beyond the limits of x */
  template<class T>
  static void addInto(T *x, const int N, const T *y, int L, int n = 0);

  /** Adds the bufferToAdd, multiplied by some weight, into the inputAndResult. This is useful for
  implementing update rules of the form: x_new = x_old + weight * delta_x - but actually x_new and 
  x_old are the same array and it's more like x += w * delta_x. */
  template<class T>
  static void addWithWeight(T* inputAndResult, const int N, const T* bufferToAdd, const T weight);

  /** Applies the affine transformation y = a*x + b to all array elements. */
  template<class T>
  static void affineTrafo(const T* x, T* y, const int N, const T a, const T b);

  /** Allocates memory for a 2-dimensional array (i.e. a matrix) with equal dimensions in both
  directions. 
  \todo: remove - redundant with allocateMatrix  */
  template<class T>
  static void allocateSquareArray2D(T**& theArray, const int size);

  /** Applies the function f given by the function-pointer to all elements in inBuffer and stores
  the result in outBuffer (both buffers may be equal). */
  template <class T>
  static void applyFunction(const T* inBuffer, T* outBuffer, const int length, T (*f) (T));
  // why can't we make the function-pointer const? ...rosic doesn't compile when trying
  // maybe use a 2nd template parameter F for the function-type, such that we are not restricted to
  // C-style function pointers

  /** Checks, if the two buffers are elementwise approximately equal within the given tolerance. To 
  be considered within the tolerance, the absolute value of the difference must not be greater 
  than the given tolerance. In the case of perfect equality, it is considered to be within the 
  tolerance to make it work correctly also for zero tolerance. */
  template <class T>
  static inline bool almostEqual(const T *buffer1, const T *buffer2, 
    const int length, const T tolerance);

  /** Checks, if the two buffers are elementwise equal. */
  template <class T>
  static inline bool equal(const T *buffer1, const T *buffer2, const int length);

  /** Circularly shifts the content of the buffer by 'numPositions' to the right - for leftward
  shifts use negative values for numPositions. If the absolute value of 'numPositions' is greater
  than the length of the buffer, it will use numPositions modulo the length - so if the length is 6
  and numPositions is 8, it will whift by 2 positions. */
  template<class T>
  static void circularShift(T *buffer, const int length, const int numPositions);
  // allocates temporary heap memory - todo: make a version that doesn't and uses a workspace
  // that is passed into the function - but keep the one without workspace as convenience function

  /** Circularly shifts the content of the buffer by 'numPositions' to the right - for leftward
  shifts use negative values for numPositions. The function behaves analogous to
  circularShift(T*, int, int) but allows for non-integer shifts by using linear interpolation of
  the buffer. */
  template <class T>
  static void circularShiftInterpolated(T *buffer, const int length, const double numPositions);
  // // Allocates heap memory - todo: pass a workspace.

  /** Clears the buffer, i.e. fills it with all zeros. */
  template <class T>
  static void clear(T* buffer, const int length) { fillWithZeros(buffer, length); }

  /** Restricts the values in the buffer to the range between min and max for types that define the
  operators '<' and '>'. */
  template <class T>
  static void clip(T *buffer, const int length, const T min, const T max);
  // todo: maybe let it use only '<' -> less requirements on the type -> more generally applicable

  /** Returns -1 if a < b, 0 if a == b, +1 if a > b. The elements are compared succesively starting
  at index 0 and when an unequal element is encountered, the buffer with the greater element is
  considered to be greater. This is consistent with strcmp or numbers in a positional number
  system. */
  template <class T>
  static int compare(const T *a, const T *b, const int length);

  /** Similar to rsCompare(T *a, T *b, int length) but allows for the 2 buffers to have different
  lengths. If they match up to the length of the shorter buffer, the longer buffer is considered
  equal, iff all the remaining entries in the tail zero, otherwise the longer buffer is considered
  to be greater. */
  template <class T>
  static int compare(const T *a, int na, const T *b, const int nb);

  /** Searches the array for an element (via the '==' operator of type T) and returns true if the
  element was found. */
  template <class T>
  static inline bool contains(const T *buffer, const int length, const T elementToFind);

  /** Convolves an array x (seen as input signal) with another array h (seen as impulse response)
  and stores the result in the array y. The type T must define the operators: *, += and a 
  constructor that takes an int and initializes to zero when 0 is passed. The y array must have a 
  length equal to xLength+hLength-1. The pointer to the y array does not have to be distinct from 
  the x and/or h pointer - i.e. the function can be used for in-place convolution, for example for 
  overwriting a sequence with the convolution product of the sequence with some other sequence or 
  even with itself (such as when squaring polynomials). */
  template <class T>
  static void convolve(const T *x, const int xLength, const T *h, const int hLength, T *y);

  /** Convolves the array x with the two-element array h and stores the result in y. The y array 
  is allowed to alias to the x array. */
  template <class T>
  static inline void convolveWithTwoElems(const T* x, const int xLength, const T* h, T* y);

  /** Convolves the array x with the two elements [h0 h1] and stores the result in y. The y array 
  is allowed to alias to the x array. This special case is needed for multiplying in a linear 
  factor into an array of polynomial coefficients. */
  template <class T>
  static inline void convolveWithTwoElems(
    const T* x, const int xLength, const T h0, const T h1, T* y);


  /** Copies the data of one array into another one and converts the type if necessary. */
  template <class T1, class T2>
  static inline void convert(const T1 *source, T2 *destination, const int length);

  /** Convolves x with h and stores the result in x. The xLength parameter denotes the number of
  values in the x-array that will be considered as input signal. The actual array must be longer
  than that (namely xLength+hLength-1) to store the appended values. */
  template <class T>
  static void convolveInPlace(T *x, const int xLength, const T *h, const int hLength);
  // DEPRECATED - we can now do in-place covolution with the regular convolve function
  // ...but maybe keep it as convenience function

  /** Copies the data of one array into another one, converting the datatype, if necessarry. */
  template <class T>
  static inline void copy(const T *source, T *destination, const int length);
  // todo: maybe make a function "move" that allows for overlap between src and dst - can be 
  // implemented like: if(src < dst) copyBackward(src,dst,N) else copyForward(src,dst,N) where the
  // two copy functions traverse the source array from front to back or the other way around


  // old version:
  /** Copies the data of one array into another one. */
  //template <class T>
  //void copy(const T *source, T *destination, int length);


  /** Copies values from the source into the target buffer if they match (via the '==' operator of
  type T) one of the elements contained in the elementsToMatch buffer. The source and target
  buffers should have the same length, because potentially all elements in the source buffer must
  be copied. If not all elements match the criterion, the target buffer's extra elements will stay
  as they were. The return value informs, how many values have been copied. Target- and source
  buffer need not to be distinct, so the function can be used to strip off undesired elements from
  a buffer in place. Just remember that in any case in the extra elements in the targetBuffer are
  garbage. */
  template <class T>
  static int copyIfMatching(const T *sourceBuffer, T *targetBuffer, int sourceAndTargetLength,
                            const T *elementsToMatch, int matchLength);

  /** Similar to rsCopyIfMatching, but here the criterion is to NOT match one of the elements in the
  set of elementsToStrip. */
  template <class T>
  static int copyIfNotMatching(const T *sourceBuffer, T *targetBuffer, int sourceAndTargetLength,
                               const T *elementsToStrip, int stripLength);

  /** Copies the data of one array into another one where the lengths of the source- and target-
  arrays may be different - in this case, the target array will be filled by linearly interpolating
  the values in the source array. The type T must define a multiplication operator with a left
  operand of type double and an addition operator with both operands of type T. At the right border
  of the source buffer, a periodicity assumption is made. */
  template <class T>
  static void copyBufferWithLinearInterpolation(const T *source, int sourceLength, T *destination,
    int destinationLength);
  // rename to copyInterpolated

  /** Copies a section of length "copyLength" starting at "copyStart" from "source" to 
  "destination". If "copyStart" is less than 0 and/or the desired "copyLength" is such that the end 
  of the copied section is beyond the end of the source-array, the head and tail of "destination" 
  will be filled with zeros appropriately, i.e. we assume "source" to contain zero values for 
  indices < 0 and indices >= sourceLength. The types of the elements of the source and destination
  arrays can be different, in which case a type conversion will be performed. */
  template<class T1, class T2>
  static void copySection(const T1 *source, int sourceLength, T2 *destination, 
    int copyStart, int copyLength);
  // what if the destination array is too short?

  // old - without type conversion
  //template<class T>
  //void rsCopySection(T *source, int sourceLength, T *destination, int copyStart, int copyLength);

  /** Computes the cumulative sum y[n] = x[n] + y[n-1] of some signal.
  \todo comment the order, periodic parameters
  */
  //template <class T>
  //void rsCumulativeSum(T *buffer, int length, int order = 1, bool periodic = false);

  /** Computes the cumulative sum of x and stores it in y. Can also be used in place (i.e. y may 
  point to the same array as x). */
  template <class T>
  static void cumulativeSum(const T *x, T *y, int N);

  /** Computes a cumulative sum of arbirtry order of x and stores it in y. Can be used in place. */
  template <class T>
  static void cumulativeSum(const T *x, T *y, int N, int order);
  // why two functions? use default argument instead

  /** Frees memory allocated previously via rsAllocateSquareArray2D. */
  template<class T>
  static void deAllocateSquareArray2D(T**& theArray, int size);
  // move to rsMatrixTools

  /** Decimates the array x of length N by the given factor and writes the result into y, which
  must be of length N/factor (floor-division). This is a naive decimation without any 
  pre-filtering. */
  template <class T>
  static void decimate(const T* x, const int N, T* y, const int factor);

  /** Like decimate, but instead of just taking every "factor"-th sample of the input, it computes
  a mean over "factor" samples. Note that the mean is computed using samples in the forward 
  direction from the readout position in the input signal (todo: make a version that uses a mean
  centered at the input read position) */
  template <class T>
  static void decimateViaMean(const T* x, const int N, T* y, const int factor);

  /** Deconvolves the impulse response h out of the signal y resulting in the signal x which has a
  length of yLength-hLength+1. It's the inverse of convolve. */
  template <class T>
  static void deConvolve(const T *y, int yLength, const T *h, int hLength, T *x);

  /** De-interleaves a buffer of interleaved data. @see rsInterleave */
  template <class T>
  static void deInterleave(T *buffer, int numFrames, int numElementsPerFrame);
  // // Allocates heap memory - todo: pass a workspace.

  /** Computes the difference y[n] = x[n] - x[n-1] of some signal. The initial condition x[-1] is
  determined from the 'periodic' parameter - if true, the signal is assumed as periodic and the
  x[-1] sample will be assigned to he same value as the last in the buffer - otherwise, it will be
  assigned to zero. When multiplying the result with the interval between successive samples,
  this function can be used for numeric differentiation. If the order is greater than 1, the
  operation will be performed iteratively on the respective result of the previous pass.  */
  template <class T>
  static void difference(T *buffer, int length, int order = 1, bool periodic = false);
  // maybe rename to backwardDifference

  /** Divides the elements of 'buffer1' and 'buffer2' - type must define operator '/'. The
  'result' buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T1, class T2, class TR>
  static void divide(const T1 *buffer1, const T2 *buffer2, TR *result, int length);

  //template <class T>
  //static bool equals(const T *x, const T *y, int length, T tolerance);

  /** Returns the Euclidean norm (i.e. the length) of the N-dimensional vector x. This is the 
  square-root of the sum of the squares of the elements. */
  template<class T>
  static T euclideanNorm(const T* x, int N);

  /** Fills the array with values given by a function. For example, calling it like:
        fill(a, length, [](int i){ return 3*i+1; });
   gives each array element the value 3*i+1 where i is the array index. */
  template<class T, class F>
  static void fill(T* a, int N, F indexToValueFunction);

  /** Fills the passed array with an impulse, i.e. all samples except one are zero. You can 
  optionally set the value of the impulse and its sample location (by default, an impulse of 
  height one occurs at sample zero - that's the so called unit impulse in DSP jargon). */
  template <class T>
  static inline void fillWithImpulse(T *buffer, int length, T value = T(1), int position = 0);

  /** Fills the passed array with the values of the indices - the type T must have a
  constructor that takes an int and performs appropriate conversion. */
  template <class T>
  static void fillWithIndex(T *buffer, int length);

  /** Fills the buffer with NaN values. This may be useful as initialization for testing code that
  is supposed to initialize the buffer correctly - failing to do so will leave NaNs which are 
  easily discovered during debugging. */
  template <class T>
  static inline void fillWithNaN(T *buffer, int length);

  /** Fills the buffer with random values between min and max. */
  template <class T>
  static void fillWithRandomValues(T *buffer, int length, double min, double max, int seed);

  /** Fills the buffer with random integer values between min and max. */
  template <class T>
  static void fillWithRandomIntegers(T *buffer, int length, int min, int max, int seed);


  /** Fills the buffer with values ranging (exponentially scaled) from min to max (both end-values
  inclusive). The logs of the values are equidistant. */
  template <class T>
  static void fillWithRangeExponential(T *buffer, int length, T min, T max);

  /** Fills the buffer with equidistant values ranging from min to max (both end-values 
  inclusive). \todo: rename min/max into start/end or x0,x1  */
  template <class T>
  static void fillWithRangeLinear(T *buffer, int length, T min, T max);
  // corresponds to std::iota? and/or NumPy's linspace

  /** Generalizes fillWithRangeLinear (shape == 1) and fillWithRangeExponential (shape == 0) and 
  allows for shapes between and beyond 0 and 1. It uses a power rule as follows: If shape != 0, it 
  fills the array with equidistant values between min^shape and max^shape and then takes the 
  shape-th root of the array values. For shape == 0 (with some tolerance), this breaks down because
  we can't take the 0-th root, so we fill the array using fillWithRangeExponential, which is the 
  correct limit. The filled array has the property that its generalized mean with parameter 
  p == shape is equal to the same generalized mean of min and max. @see generalizedMean. */
  template <class T>
  static void fillWithRange(T *buffer, int length, T min, T max, T shape);
  // todo: make a unit test, avoid roundoff error at endpoints

  /** Fills the passed array with one value at all indices. */
  template <class T>
  static void fillWithValue(T *buffer, int length, T value);
  // rename to fill and/or make alias

  /** Fills the passed array with all zeros - the type must have a constructor that takes an int
  and initializes to the zero element when 0 is passed. */
  template <class T>
  static void fillWithZeros(T *buffer, int length);
  // rename to clear and/or make alias

  /** Filters the signal in input buffer x and stores the result in output buffer y. The filter
  realizes the difference equation:
  y[n] = b[0]*x[n] + b[1]*x[n-1] + b[2]*x[n-2] + ... + b[bOrder]*x[n-bOrder]
                   - a[1]*y[n-1] - a[2]*y[n-2] - ... - a[aOrder]*y[n-aOrder]
  so the arrays of coefficients must have a length of their order + 1. The output buffer may have
  a different length than the input buffer. If it is longer, it will be assumed that the missing
  samples in the input buffer are zero. The input and output buffers may also be identical (i.e.
  point to the same location), in which case the filtering will be done in place. */
  template <class T>
  static void filter(const T *x, int xLength, T *y, int yLength, 
    const T *b, int bOrder, const T *a, int aOrder);
  // Allocates heap memory - todo: pass a workspace.

  /** \todo check and comment this function  */
  template <class T>
  static void filterBiDirectional(const T *x, int xLength, T *y, int yLength, 
    const T *b, int bOrder, const T *a, int aOrder, int numRingOutSamples = 10000);
  // Allocates heap memory - todo: pass a workspace.

  /** Returns the index of the first value that matches the elementToFind, returns -1 when no 
  matching element is found. */
  template <class T>
  static inline int findIndexOf(const T *buffer, T elementToFind, int length);
  // rename to firstIndexOf make also a lastIndexOf

  /** Returns the index of the maximum absolute value in the buffer. */
  template <class T>
  static inline int findMaxAbs(const T *buffer, int length);

  /** Returns the index of a peak or value in the length-N array x right or equal to n0 or -1 if 
  none is found. */
  template<class T>
  static int findPeakOrValleyRight(const T *x, int N, int n0);

  /** Returns the index of a peak or value in the length-N array x left or equal to n0 or -1 if 
  none is found. */
  template<class T>
  static int findPeakOrValleyLeft(const T *x, int N, int n0);

  /** Returns the first index in the ascendingly sorted array "a", where the value is greater 
  than or equal to "splitValue". The idea is that you may split the array at that index i into
  two subarrays where in the left subarray a[0...i-1], all values are less than splitValue and in 
  the right subarray a[i...N-1], all values are greater-or-equal to splitValue. Note that if all 
  values in the array are actually less that the splitValue, N will be returned which is an index 
  after the last valid index - don't dereference it in this case! */
  template<class T>
  static int findSplitIndex(const T* a, int N, T splitValue);
  // T should be Ordered -> use concept Ordered/Sortable in c++20
  // maybe make a function findPositionOfValue that returns a floating point position, i.e. a 
  // fractional index, using linear interpolation (and extrapolation)

  /** Like splitIndex, but instead of just returning the first index i, where a[i] >= splitValue, 
  it checks, if a[i-1] is closer to the splitValue than a[i]. If it is, then it returns i-1
  instead of i. Also, if all entries are less than splitValue, it will return N-1 instead of N, so
  it will always return a valid index. */
  template<class T>
  static int findSplitIndexClosest(const T* a, const int N, const T splitValue);
  // used in rsEnvelopeExtractor<T>::fillSparseAreas

  /** Given a buffer of values, this function returns the first index where there's a nonzero value
  in the buffer. If there's no nonzero element at all, it returns -1.

  \todo use a different convention: return "length", when no index is found - this allows to return
  an unsigned int and extends the range where ths function can be used by a factor of (almost) 2

  */
  template <class T>
  static int firstIndexWithNonZeroValue(const T *buffer, int length);
  // todo: rename to something shorter - maybe 
  // firstNonZeroIndex            may not be that clear
  // indexOfFirstNonZero          seems a good compromise between length and clarity
  // firstIndexWithNonZeroValue   for length comparison reference

  /** Scales and offsets the passed buffer such that the minimum value hits 'min' and the
  maximum value hits 'max'. */
  template <class T>
  static void fitIntoRange(T *buffer, int length, T min, T max);

  /** Computes the geometric mean of the values in array x. This is the N-th root of the product of 
  all values. Due to the identity log(a*b) = log(a) + log(b) (which generalizes to products of 
  more than two factors in the obvious way), this can also be seen as taking the arithmetic average 
  of the logarithms of the values and then exponentiating that result (which is, in fact, the way 
  it is implemented because taking products of a great number of factors tends to create numerical 
  overflow or underflow problems). */
  template<class T>
  static T geometricMean(const T* x, int N);

  /** Computes the generalized mean with exponent p of the numbers in array x. The generalized mean
  is the p-th root of the arithmetic mean of all values raised to the p-th power. For p = 1, this 
  reduces to the usual arithmetic mean, for p = 2, it is the root-mean-square, for p = -1, it's the
  harmonic mean etc. The case p = 0 is a special case which, in the limit, becomes the geometric 
  mean. For p going to minus or plus infinity, it tends to the minimum or the maximum of all values
  respectively (...has not been tested yet for very large absolute values of p...use reasonable 
  values - like -10..+10 or something!). https://en.wikipedia.org/wiki/Generalized_mean */
  template<class T>
  static T generalizedMean(const T* x, int N, T p);

  /** Fills the array h with the impulse of the filter specified by the direct form coefficients
  given in b and a. \todo comment on the sign of the a-coeffs, whether or not an a0 is included,
  etc., maybe move to RSLib into the filter section */
  template <class T>
  static void impulseResponse(T *h, int hLength, const T *b, int bOrder, const T *a, int aOrder);

  /** Interleaves a buffer of non-interleaved data. */
  template <class T>
  static void interleave(T *buffer, int numFrames, int numElementsPerFrame);
  // Allocates heap memory - todo: pass a workspace.

  /** Returns a linearly interpolated value from the array at the given (non-integer) position. If
  the position is out of range, 0 is returned. */
  template<class T>
  static T interpolatedValueAt(const T *buffer, int length, double position);
    // Why is position a double? T may be bad if T is a vector or complex type - but maybe some 
    // generic TPos would be better?

  /** Returns a linearly interpolated value from the array at the given (non-integer) position. If
  the position is out of range, the output is clamped to the end values of the array. */
  template<class T>
  static T interpolateClamped(const T *buffer, int length, double position);

  /** Returns true, iff the passed buffer has only zero values. */
  template <class T>
  static inline bool isAllZeros(const T *buffer, int length);

  /** Returns true, iff the passed buffer has only values close to zero within the given 
  tolerance. */
  template <class T>
  static inline bool isAllZeros(const T *buffer, int length, T tolerance);

  /** Returns true, if the passed buffer has values equal to the given value, false otherwise. */
  template <class T>
  static inline bool isFilledWithValue(const T *buffer, int length, T value);

  /** Returns true if the value in array x at position n is greater than its left and right 
  neighbour (the caller must be sure, that n-1, n, n+1 are valid array indices). */
  template<class T>
  static bool isPeak(const T *x, int n);
  // inline this

  /** Similar to isPeak, but uses "greater-or-equal" comparison instead of "greater". */
  template<class T>
  static inline bool isPeakOrPlateau(const T *x, int n);

  /** Returns true if the value in array x at position n is smaller than its left and right 
  neighbour (the caller must be sure, that n-1, n, n+1 are valid array indices). */
  template<class T>
  static bool isValley(const T *x, int n);

  /** Returns true if the value in array x at position n is larger or smaller than its left and 
  right neighbour (the caller must be sure, that n-1, n, n+1 are valid array indices). */
  template<class T>
  static bool isPeakOrValley(const T *x, int n);

  /** Checks whether the buffer is sorted in ascending order, that is buffer[i] <= buffer[i+1] for
  all i, but it actually uses only the < operator in the implementation: it will return false, iff
  buffer[i+1] < buffer[i] for some value of i. */
  template <class T>
  static bool isSortedAscending(const T *buffer, int length);

  /** Checks whether the buffer is sorted in strictly ascending order, that is 
  buffer[i] < buffer[i+1] for all i. */
  template <class T>
  static bool isSortedStrictlyAscending(const T *buffer, int length);


  /** Shifts the content of the buffer numPlaces to the left, filling it up with zeros from the
  right. */
  template <class T>
  static void leftShift(T *buffer, int length, int numPlaces);

  /** Limits the passed value into the range between min and max (inclusive) and returns the
  result. */
  template <class T>
  static T limitToRange(T value, T min, T max);

  /** Fills the length-N array maxXY with the element-wise maximum of the arrays x and y. */
  template <class T>
  static void maxElementWise(const T *x, const T* y, const int N, T* maxXY);
  // todo: make a similar function for min

  /** Finds and returns the maximum absolute value of the buffer. */
  template <class T>
  static T maxAbs(const T *buffer, int length);

  template <class T>
  static T maxAbs(const std::complex<T> *buffer, int length);


  /** Finds and returns the index with the maximum absolute value of the buffer. */
  template <class T>
  static int maxAbsIndex(const T* const buffer, int length);

   /** Returns the maximum deviation (absolute value of the difference) between two buffers. */
  template <class T>
  static T maxDeviation(const T *buffer1, const T *buffer2, int length);

  /** Returns the index of the maximum deviation (see maxDeviation). */
  template <class T>
  static int maxDeviationIndex(const T *buffer1, const T *buffer2, int length);

  /** Returns the maximum value of the differences of adjacent array entries.  */
  template <class T>
  static T maxDifference(const T *buffer, int length);

  /** Returns the index of maximum value of the buffer (">"-operator must be defined). */
  template <class T>
  static int maxIndex(const T *buffer, int length);

  /** Returns the maximum value of the buffer (">"-operator must be defined). */
  template <class T>
  static T maxValue(const T *buffer, int length);

  /** Returns the index of minimum value of the buffer ("<"-operator must be defined). */
  template <class T>
  static int minIndex(const T *buffer, int length);

  /** Returns the minimum value of the buffer ("<"-operator must be defined). */
  template <class T>
  static T minValue(const T *buffer, int length);

  /** Computes the mean (i.e. the DC-component) from the passed buffer. The type must define
  operators: +=, / and a constructor which takes an int and initializes to zero when 0 is passed
  and a typecast from int. */
  template <class T>
  static T mean(const T *buffer, int length);

  /** Computes the mean of the differences of the array elements. */
  template <class T>
  static T meanDifference(const T* buffer, int length);

  /** Returns the mean of the squares of the values in the array. */
  template<class T>
  static T meanSquare(const T *x, int N);

  /** Returns the number of nonzero values in the vector array x. */
  template<class T>
  static int numNonZeros(const T *x, int N);

  /** Computes the mean of the absoulte values of the differences of arrays x and y */
  template<class T>
  inline static T meanOfAbsoluteDifferences(const T* x, const T* y, const int N)
  { return sumOfAbsoluteDifferences(x, y, N) / T(N); }

  /** Returns the median of the passed buffer.  */
  template <class T>
  static T median(const T *buffer, int length);
  // Allocates heap memory - todo: pass a workspace.

  /** Applies a 3-point moving average filter to the length-N array x and stores the result in y, 
  which may point to the same memory location, i.e. the filter may be used in place. The endpoints
  are either held fixed or handled using a 1-sided 2-point average, depending on the optional 
  endsFixed parameter (default: true, because fixed ends may be the more typical use-case - for 
  example, when a parameter trajectory should be smoothed). For smoothing, it may be useful to 
  apply the function iteratively multiple times (although, it may be more efficient and give 
  similar results to use a bidirectional IIR filter with a Gaussian impulse response in this 
  case). */
  template<class T>
  static void movingAverage3pt(const T* x, int N, T* y, bool endsFixed = true);
  // todo: maybe make a version that preserves the mean (calculate mean before and after and add
  // the difference
  // maybe have a version that leaves the endpoints alone - rationale: the array may represent a
  // trajectory that should be smoothed, but the start- and enpoints should not change - maybe
  // have a boolean parameter "fixEnds"
  // -maybe include a stride parameter - needed when this is used to filter nD arrays like images
  // there's a 5 point version under construction near the unit test of this, so if one is needed,
  // look there first

  template<class T>
  static void movingMedian3pt(const T* x, int N, T* y);


  /** Multiplies the elements of 'buffer1' and 'buffer2' - type must define operator '*'. The
  'result' buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T1, class T2, class TR>
  static void multiply(const T1 *buffer1, const T2 *buffer2, TR *result, int length);

  /** Writes the element-wise negation of the source buffer into the destination buffer, i.e flips 
  the signs. */
  template<class T>
  static void negate(const T *source, T *destination, int length);

  /** Like negate but flips only signs of elements with even indices. */
  template<class T>
  static void negateEven(const T *source, T *destination, int length);
  // needs test

  /** Like negate but flips only signs of elements with odd indices. */
  template<class T>
  static void negateOdd(const T *source, T *destination, int length);
  // needs test


  /** Normalizes the maximum absolute value of the passed array by multiplying the whole array 
  through by "maximum"/maxAbs(buffer) - where "maximum" is the passed argument and maxAbs(buffer)
  is the maximum absolute value in the buffer. It may optionally subtract the mean of the array
  before doing this normalization. The type must define: >, *=, / and a constructor that takes an 
  int and initializes to 0 when 0 is passed. Additionaly, it must be suitable for rsAbs - that 
  additionaly requires definition of unary '-' and '<'. */
  template <class T>
  static void normalize(T *buffer, int length, T maximum = T(1), bool subtractMean = false);

  /** Normalizes the mean of the array to the given value (default 1) - useful for window functions
  and lowpass filters that should have unit gain at DC. */
  template <class T>
  static void normalizeMean(T* x, int N, T newMean = 1);

  /** Rearranges/permutes and array of class T into bit-reversed order (in place). The 'length'
  MUST be the 'numBits' th power of two (this is not checked for). */
  //template <class T>
  //RS_INLINE void orderBitReversedInPlace(T* buffer, int length, int numStages);

  /** Reorders an array of class T and length N into bit-reversed order (in place). N must be a
  power of 2 and log2(N) must be passed in the argument log2N. log2(N) is the number of bits in the
  binary representation of the array index (verify this). */
  template <class T>
  static void orderBitReversed(T *buffer, int N, int log2N);

  /** Rearranges/permutes and array of type T into bit-reversed order. The 'length' MUST be the
  'numBits' th power of two (this is not checked for). */
  template <class T>
  static void orderBitReversedOutOfPlace(const T *inBuffer, T *outBuffer, int length, int numBits);

  /** Returns the product of the elements in the buffer for types which define the
  multiplication operator (the *= version thereof) and a constructor which can take an int
  paramater as argument and initializes to the multiplicative neutral element of that class when 1
  is passed . */
  template <class T>
  static T product(const T* const buffer, int length);

  /** Shifts the values in the 4-element array "a" one position back, discarding the last one a[3] 
  and inserting the new x at the front a[0]. */
  template<class T>
  static inline void pushFrontPopBack4(T x, T* a);

  /** Removes mean (i.e. the DC-component) from the passed buffer. The type must define operators:
  +=, -=, / and a constructor which takes an int and initializes to zero when 0 is passed and a
  typecast from int. */
  template <class T>
  static void removeMean(T *buffer, int length);

  /** Reverses the order of the elements the passed array. */
  template <class T>
  static void reverse(T *buffer, int length);

  /** Fills array y with the reversed content of array x. x and y should not be overlapping except
  if they point to the same array, in which case we fall back to reverse(T*, int) - the in-place 
  version with just one array parameter. */
  template <class T>
  static void reverse(const T* x, T* y, int length);

  /** Shifts the content of the buffer numPlaces to the right, filling it up with zeros from the
  left. */
  template <class T>
  static void rightShift(T *buffer, int length, int numPlaces);

  /** Returns the square-root of the mean of the squares (RMS value) of the values in the array. */
  template<class T>
  static inline T rootMeanSquare(const T *x, int N) { return rsSqrt(meanSquare(x, N)); }

  /** Scales the buffer by a constant factor. */
  template <class T1, class T2>
  static inline void scale(T1 *buffer, int length, T2 scaleFactor);

  /** Scales the "src" buffer by a constant factor and writes the result into the "dst" buffer. */
  template <class T1, class T2>
  static inline void scale(const T1 *src, T1 *dst, int length, T2 scaleFactor);

  /** Given the sequence y of length yLength, this function returns a sequence x which, when
  convolved with itself, gives y. yLength is assumed to be odd, the index of first nonzero value
  is assumed to be even and the first nonzero value in y is assumed to be positive (for real
  sequences), because this is what will happen, when convolving a (real) sequence with itself. To
  disambiguate the square-root, the function will return a sequence with its 1st nonzero value
  being positive. If the original sequence x (before it was convolved with itself to give y)
  started with a negative value, the result of taking the square-root of the squared sequence y
  will have all signs reversed with respect to the original sequence x (again, in the case of real
  sequences - \todo explain what happens in the complex case).
  The length of x will be (yLength+1)/2.

  \todo: this apparently works only if the y-sequence was actually constructed by squaring some
  sequence. For a general sequence y, terms of y and x^2 will match only up to the n-th term
  where n = (yLength+1)/2, i.e. the length of x. ...find out if this function is useless
  therefore and/or if there is a useful generalization - like using a sequence x which has the
  same length as y, convolve it with itslef and truncate the result to the length of y */
  template <class T>
  static void sequenceSqrt(const T *y, int yLength, T *x);

  /** Shifts the content of the buffer numPlaces to the right, filling it up with zeros from the
  left. If numPlaces is negative, the contents will be shifted to the left, filling up with zeros 
  from the right. */
  template <class T>
  static void shift(T *buffer, int length, int numPlaces);

  /** Shifts the values in the array x up or down such that the new minimum off all values will 
  become zero and writes the result to y. The return value is the minimum of the x-values - if you 
  add that value to the resulting y-values, you should get your old x-array back (up to roundoff 
  error). */ 
  template <class T>
  static T shiftToMakeMinimumZero(const T* x, int N, T* y);
  // could also be called subtractMinimum - but the current name better reveals the intention

  /** Subtracts the elements of 'buffer2' from 'buffer1' - type must define operator '-'. The
  'result' buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T>
  static inline void subtract(const T *buffer1, const T *buffer2, T *result, int length);

  /** Returns the sum of the elements in the buffer for types which define the
  addition operator (the += version thereof) and a constructor which can take an int
  paramater as argument and initializes to the additive neutral element of that class when 0
  is passed . */
  template <class T>
  static T sum(const T* buffer, int length);

  /** Given two arrays x and y of lengths N, this function computes the sum of the products 
  x[i]*y[i] where i runs from 0 to N-1. */
  template<class T>
  static T sumOfProducts(const T *x, const T *y, int N);

  /** Given an array x of length N, this function commputes the sum of the squares of the values
  in x. */
  template<class T>
  static T sumOfSquares(const T *x, int N);

  /** Computes the sum of the absoulte values of the differences of arrays x and y */
  template<class T>
  static T sumOfAbsoluteDifferences(const T *x, const T *y, const int N);

  /** Swaps the contents of of buffer1 and buffer2 using an auxiliary buffer bufferTmp. All buffers
  are assumed to have a size of sizeInBytes. */
  static void swapDataBuffers(void *buffer1, void *buffer2, void *bufferTmp, int sizeInBytes);
  // rename to swapContent, make a version that operates in place. maybe this function should not
  // be in rsArrayTools - the void-pointers and byte-size do not really follow the general 
  // pattern here

  /** Applies an affine tranform y = a*x + b to the array x, where a and b are chosen such that the
  transformed values span a range from targetMin to targetMax. Note that this works only, if the
  array is non-constant - if the current maximum equals the minimum, the formula will produce 
  NaNs. */
  template<class T>
  static void transformRange(const T* x, T* y, int N, T targetMin, T targetMax);

  /** Transposes a 2-dimensional square-array, such that the rows become columns and vice versa.
  The input-array is passed via "in" and the output array will be stored in "out". The function
  CANNOT be used in place. */
  template<class T>
  static void transposeSquareArray(T **in, T **out, int size);

  /** In-place version of rsTransposeSquareArray(T **in, T **out, int size). */
  template<class T>
  static void transposeSquareArray(T **theArray, int size);
  // move to MatrixFunctions, rename to transposeMatrix

  /** Given an array that contains values that have been subject to some kind wrap-around into some
  period p, this function (heuristically) unwraps them again. This is useful for plotting a 
  phase-response - if you have an array of angles in the interval -pi..pi (as you will get when 
  using arg(std::complex z) on a transfer function value), subtract off 2*pi, and unwrap with 
  period 2*pi. */
  template<class T>
  static void unwrap(T* a, int N, T p);
  // ToDo: check, if we have unit-test in place - it doesn't seem to always work...

  /** Returns the sum-over-i w[i]*x[i]. */
  template <class T>
  static T weightedSum(const T *w, const T *x, rsUint32 length);  // redundant with sumOfProducts

  /** Forms a weighted sum of the two buffers. */
  template <class T>
  static inline void weightedSum(const T *buffer1, const T *buffer2, T *result, 
    int length, T weight1, T weight2);

};

//=================================================================================================
// inlined implementations

template <class T>
inline void rsArrayTools::add(const T *buffer1, const T *buffer2, T *result, const int length)
{
  for(int i = 0; i < length; i++)
    result[i] = buffer1[i] + buffer2[i];
}

template <class T>
inline void rsArrayTools::add(const T *buffer, const T valueToAdd, T *result, const int length)
{
  for(int i = 0; i < length; i++)
    result[i] = buffer[i] + valueToAdd;
}

template<class T>
void rsArrayTools::addWithWeight(T* xy, const int N, const T* d, const T w)
{
  for(int n = 0; n < N; n++)
    xy[n] += w * d[n];
}

/*
template <class T>
inline bool rsArrayTools::almostEqual(const T *buffer1, const T *buffer2, 
  const int length, const T tolerance)
{
  for(int i = 0; i < length; i++) {
    if(rsAbs(buffer1[i]-buffer2[i]) > tolerance)
      return false; }
  return true;
}
// make a version that does not rely on taking the difference and then using abs - this will work 
// only for signed data types, unsigned types may wrap around to huge values
*/

/*
template <class T>
inline bool rsArrayTools::almostEqual(const T* x, const T* y, const int N, const T tol)
{
  for(int i = 0; i < N; i++) {
    T d;
    if(x[i] > y[i])
      d = x[i] - y[i];
    else
      d = y[i] - x[i];
    if(d > tol)
      return false; }
  return true;
}
// this is supposed to work without using an abs-function, so it may used with unsigned types
// that's still not suitable for rsPixelRGB - we need a function rsAlmostEquals for single 
// data-values and provide implementations for different types
*/

template <class T>
inline bool rsAlmostEqual(const T& x, const T& y, const T& tol)
{
  return rsAbs(x-y) <= tol;
}
// suitable for float and double - move to somewhere else

template <class T>
inline bool rsArrayTools::almostEqual(const T* x, const T* y, const int N, const T tol)
{
  for(int i = 0; i < N; i++)
    if(!rsAlmostEqual(x[i], y[i], tol))
      return false;
  return true;
}

template <class T>
inline void rsArrayTools::copy(const T *source, T *destination, const int length)
{
  for(int i = 0; i < length; i++)
    destination[i] = source[i];
}
// todo: switch at compile time between using memcopy for trivially copyable types and using the 
// assignment operator in a loop (such as now) otherwise - see:
// https://en.cppreference.com/w/cpp/types/is_trivially_copyable
// ...but maybe for short arrays, the overhead of calling memcpy may outweigh its efficiency, so 
// for very small arrays, the loop is actually faster? -> do benchmarks
// see here, around 50min:  https://www.youtube.com/watch?v=ZeU6OPaGxwM

template <class T>
inline bool rsArrayTools::contains(const T *buffer, const int length, const T elementToFind)
{
  for(int i = 0; i < length; i++) {
    if(buffer[i] == elementToFind)
      return true; }
  return false;
  //return (rsFindFirstOccurrenceOf(buffer, length, elementToFind) != -1);
}

template <class T1, class T2>
inline void rsArrayTools::convert(const T1 *source, T2 *destination, const int length)
{
  for(int i = 0; i < length; i++)
    destination[i] = (T2)source[i];
}

template <class T>
inline void rsArrayTools::convolveWithTwoElems(const T* x, const int xLength, const T* h, T* y)
{
  rsAssert(xLength > 0, "Input length must be > 0");
  y[xLength] = x[xLength-1]*h[1];
  for(int n = xLength-1; n > 0; n--)
    y[n] = x[n]*h[0] + x[n-1]*h[1];
  y[0] = x[0]*h[0];
}

template <class T>
inline void rsArrayTools::convolveWithTwoElems(
  const T* x, const int xLength, const T h0, const T h1, T* y) 
{
  rsAssert(xLength > 0, "Input length must be > 0");
  y[xLength] = x[xLength-1]*h1;
  for(int n = xLength-1; n > 0; n--)
    y[n] = x[n]*h0 + x[n-1]*h1;
  y[0] = x[0]*h0;
}

template <class T>
inline void rsArrayTools::decimate(const T* x, const int N, T* y, const int D)
{
  for(int i = 0; i < N/D; i++)
    y[i] = x[i*D];
}

template <class T>
inline void rsArrayTools::decimateViaMean(const T* x, const int N, T* y, const int D)
{
  const T s = T(1) / T(D);       // scaler
  for(int i = 0; i < N/D; i++)
    y[i] = s * sum(&x[i*D], D);  // s * sum(...) == mean(...)
}

template <class T>
inline bool rsArrayTools::equal(const T *buffer1, const T *buffer2, const int length)
{
  for(int i = 0; i < length; i++)
  {
    if(buffer1[i] != buffer2[i])
      return false;
  }
  return true;
}

template<class T>
T rsArrayTools::euclideanNorm(const T* x, int N)
{
  return sqrt(rsArrayTools::sumOfSquares(x, N));
}

template<class T, class F>
inline void rsArrayTools::fill(T* a, int N, F indexToValueFunction)
{
  for(int i = 0; i < N; i++)
    a[i] = indexToValueFunction(i);
}

template <class T>
inline void rsArrayTools::fillWithImpulse(T *buffer, int length, T value, int position)
{
  rsAssert(position >= 0 && position < length, "Impulse location out of range");
  fillWithZeros(buffer, length);
  buffer[position] = value;
}

template <class T>
inline void rsArrayTools::fillWithNaN(T* x, int N)
{
  for(int i = 0; i < N; i++) 
    x[i] = RS_SIGNALING_NAN(T); // or maybe it should be the quiet nan?
}

template <class T>
inline int rsArrayTools::findIndexOf(const T *buffer, T elementToFind, int length)
{
  for(int i = 0; i < length; i++)
  {
    if( buffer[i] == elementToFind )
      return i;
  }
  return -1;
}

template <class T>
inline int rsArrayTools::findMaxAbs(const T *buffer, int length)
{
  int maxIndex = 0;
  T maxValue   = T(0);
  for(int i=0; i<length; i++)
  {
    if( absT(buffer[i]) > maxValue )
    {
      maxValue = absT(buffer[i]);
      maxIndex = i;
    }
  }
  return maxIndex;
}

template<class T>
int rsArrayTools::findSplitIndex(const T* A, int N, T key)
{
  rsAssert(isSortedAscending(A, N), "array must be sorted");
  int imin = 0;
  int imax = N-1;

  if(A[imax] < key) // new - needs tests
    return N;

  while( imin < imax ) {
    int imid = imin/2 + imax/2;
    //rsAssert(imid < imax); // only for debug
    if( A[imid] < key )
      imin = imid + 1;
    else
      imax = imid;
  }
  return imin;
}
// compare to this: https://en.wikipedia.org/wiki/Binary_search_algorithm
// what about RSLib? look, if we have something like hat there already

template <class T>
inline bool rsArrayTools::isAllZeros(const T *buffer, int length)
{
  return isFilledWithValue(buffer, length, T(0));
}

template <class T>
inline bool rsArrayTools::isAllZeros(const T* a, int N, T tol)
{
  for(int i = 0; i < N; i++)
    if( rsGreaterAbs(a[i], tol) )
      return false;
  return true;
}

template <class T>
inline bool rsArrayTools::isFilledWithValue(const T *buffer, int length, T value)
{
  for(int n = 0; n < length; n++)
  {
    if(buffer[n] != value)
      return false;
  }
  return true;
}

template<class T>
bool rsArrayTools::isPeakOrPlateau(const T *x, int n)
{
  if(x[n] >= x[n-1] && x[n] >= x[n+1])
    return true;
  return false;
}

template <class T>
bool rsArrayTools::isSortedAscending(const T *buffer, int length)
{
  for(int i = 0; i < length-1; i++) 
  {
    //if(!(buffer[i] <= buffer[i+1])) return false; // old, requires, <=
    if(buffer[i+1] < buffer[i]) return false; // new, requires only < to be implemented
  }
  return true;
}

template<class T>
inline void rsArrayTools::pushFrontPopBack4(T x, T* a)
{
  a[3] = a[2];  // old a[3] is discarded
  a[2] = a[1];  // values in between..
  a[1] = a[0];  // ...are shifted
  a[0] = x;     // input x becomes the new a[0]
}
// todo: make versions for 1,2,3,N (using a loop or memmove) -> make performance tests, 
// which version is fastest for what range of lengths maybe rename to updateFifoBuffer4
// maybe make a version that returns the (old) a[3] element - the caller may be interested in it

template <class T1, class T2>
inline void rsArrayTools::scale(T1 *buffer, int length, T2 scaleFactor)
{
  for(int n = 0; n < length; n++)
    buffer[n] *= (T1)scaleFactor;
}

template <class T1, class T2>
inline void rsArrayTools::scale(const T1 *src, T1 *dst, int length, T2 scaleFactor)
{
  for(int n = 0; n < length; n++)
    dst[n] = scaleFactor * src[n];
}

template <class T>
inline void rsArrayTools::subtract(const T *buffer1, const T *buffer2, T *result, int length)
{
  for(int i = 0; i < length; i++)
    result[i] = buffer1[i] - buffer2[i];
}

template <class T>
inline void rsArrayTools::weightedSum(const T *buffer1, const T *buffer2, T *result, int length, 
  T weight1, T weight2)
{
  for(int n = 0; n < length; n++)
    result[n] = weight1 * buffer1[n] + weight2 * buffer2[n];
}


#endif