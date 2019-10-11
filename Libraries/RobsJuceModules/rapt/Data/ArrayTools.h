#ifndef RAPT_ARRAYTOOLS_H_INCLUDED
#define RAPT_ARRAYTOOLS_H_INCLUDED

/** A collection of functions that operate on 1-dimensional arrays. 

todo: 
-make everything const that is possible (also by-value parameters, local variables, etc. - and use
 constexpr for compile-time constants)
 ->done up to copyBuffer
 ...maybe change the const by-value parameters to by-reference parameters
-inline, where it makes sense (trivial functions like copy/convert)
-maybe turn into an actual class (with members) implementing a dynamically sized array

*/

class rsArray
{

public:

  /** Adds the elements of 'buffer1' and 'buffer2' - type must define operator '+'. The 'result'
  buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T>
  static void add(const T *buffer1, const T *buffer2, T *result, const int length);

  /** Adds the scalar 'valueToAdd' to the elements of 'buffer' - the type must define
  operator '+'. The 'result' buffer may be the same as 'buffer'. */
  template <class T>
  static void add(const T *buffer, const T valueToAdd, T *result, const int length);

  /** Adds a weighted, circularly shifted copy of the buffer to itself - the shift-offest may be
  non-integer in which case linear interpolation will be used. 
  \todo: maybe generalize such that a circularly shifted 2nd buffer can be added (which may or
  may not be the same buffer).  */
  template <class T>
  static void addCircularShiftedCopy(T *buffer, const int length, const double offset, 
    const T weight);
  // allocates heap memory - todo: use a workspace parameter

  /** Adds length-L array y into length-N array x starting at n (in x), taking care of not reading 
  beyond the limits of y and writing beyond the limits of x */
  template<class T>
  static void addInto(T *x, const int N, const T *y, int L, int n = 0);

  /** Applies the affine transformation y = a*x + b to all array elements. */
  template<class T>
  static void affineTrafo(const T* x, T* y, const int N, const T a, const T b);

  /** Allocates memory for a 2-dimensional array (i.e. a matrix) with equal dimensions in both
  directions. 
  \todo: remove - redundant with rsAllocateMatrix  */
  template<class T>
  static void allocateSquareArray2D(T**& theArray, const int size);

  /** Applies the function f given by the function-pointer to all elements in inBuffer and stores
  the result in outBuffer (both buffers may be equal). */
  template <class T>
  static void applyFunction(const T* inBuffer, T* outBuffer, const int length, T (*f) (T));
  // why can't we make the function-pointer const? ...rosic doesn't compile when trying

  /** Checks, if the two buffers are elementwise approximately equal within the given tolerance. */
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
  // that is passed inot the function - but keep the one without workspace as convenience function

  /** Circularly shifts the content of the buffer by 'numPositions' to the right - for leftward
  shifts use negative values for numPositions. The function behaves analogous to
  circularShift(T*, int, int) but allows for non-integer shifts by using linear interpolation of
  the buffer. */
  template <class T>
  static void circularShiftInterpolated(T *buffer, const int length, const double numPositions);

  /** Restricts the values in the buffer to the range between min and max for types that define the
  operators '<' and '>'. */
  template <class T>
  static void clip(T *buffer, const int length, const T min, const T max);

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
  static bool contains(const T *buffer, const int length, const T elementToFind);

  /** Convolves an array x (seen as input signal) with another array h (seen as impulse response)
  and stores the result in the array y. The type must define the operators: *, += and a constructor
  that takes an int and initializes to zero when 0 is passed. The y array must have a length equal 
  to xLength+hLength-1. The pointer to the y array does not have to be distinct from the x and/or h
  pointer - i.e. the function can be used for in-place convolution, for example for overwriting a
  sequence with the convolution product of the sequence with some other sequence or even with
  itself (such as when squaring polynomials). */
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
  static inline void convertBuffer(const T1 *source, T2 *destination, const int length);
  // rename to convert

  /** Convolves x with h and stored the result in x. The xLength parameter denotes the number of
  values in the x-array that will be considered as input signal. The actual array must be longer
  than that (namely xLength+hLength-1) to store the appended values. */
  template <class T>
  static void convolveInPlace(T *x, const int xLength, const T *h, const int hLength);
  // DEPRECATED - we can now do in-place covolution with the regular convolve function
  // ...but maybe keep it as convenience function

  /** Copies the data of one array into another one, converting the datatype, if necessarry. */
  template <class T1, class T2>
  static void copyBuffer(const T1 *source, T2 *destination, const int length);
  // rename to copy

  // old version:
  /** Copies the data of one array into another one. */
  //template <class T>
  //void copyBuffer(const T *source, T *destination, int length);


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

  /** Copies a section of length "copyLength" starting at "copyStart" from "source" to 
  "destination". If "copyStart" is less than 0 and/or the desired "copyLength" is such that the end 
  of the copied section is beyond the end of the source-array, the head and tail of "destination" 
  will be filled with zeros appropriately, i.e. we assume "source" to contain zero values for 
  indices < 0 and indices >= sourceLength. */
  template<class T1, class T2>
  static void copySection(const T1 *source, int sourceLength, T2 *destination, 
    int copyStart, int copyLength);

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

  /** Deconvolves the impulse response h out of the signal y resulting in the signal x which has a
  length of yLength-hLength+1. It's the inverse of convolve. */
  template <class T>
  static void deConvolve(const T *y, int yLength, const T *h, int hLength, T *x);

  /** De-interleaves a buffer of interleaved data. @see rsInterleave */
  template <class T>
  static void deInterleave(T *buffer, int numFrames, int numElementsPerFrame);

  /** Computes the difference y[n] = x[n] - x[n-1] of some signal. The initial condition x[-1] is
  determined from the 'periodic' parameter - if true, the signal is assumed as periodic and the
  x[-1] sample will be assigned to he same value as the last in the buffer - otherwise, it will be
  assigned to zero. When multiplying the result with the interval between successive samples,
  this function can be used for numeric differentiation. If the order is greater than 1, the
  operation will be performed iteratively on the respective result of the previous pass.  */
  template <class T>
  static void difference(T *buffer, int length, int order = 1, bool periodic = false);

  /** Divides the elements of 'buffer1' and 'buffer2' - type must define operator '/'. The
  'result' buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T1, class T2, class TR>
  static void divide(const T1 *buffer1, const T2 *buffer2, TR *result, int length);

  /** Fills the passed array with a unit impulse. */
  template <class T>
  static inline void fillWithImpulse(T *buffer, int length);

  /** Fills the passed array with the values of the indices - the type T must have a
  constructor that takes an int and performs appropriate conversion. */
  template <class T>
  static void fillWithIndex(T *buffer, int length);

  /** Fills the buffer with NaN values. This may be useful as initialization for testing code that
  is supposed to initialize the buffer correctly - failing to do so will leave NaNs which are 
  easily discovered druing debugging. */
  template <class T>
  static inline void fillWithNaN(T *buffer, int length);

  /** Fills the buffer with random values between min and max. */
  template <class T>
  static void fillWithRandomValues(T *buffer, int length, double min, double max, int seed);

  /** Fills the buffer with values ranging (exponentially scaled) from min to max (both end-values
  inclusive). */
  template <class T>
  static void fillWithRangeExponential(T *buffer, int length, T min, T max);

  /** Fills the buffer with values ranging from min to max (both end-values inclusive). 
  \todo: rename min/max into start/end  */
  template <class T>
  static void fillWithRangeLinear(T *buffer, int length, T min, T max);

  /** Fills the passed array with one value at all indices. */
  template <class T>
  static void fillWithValue(T *buffer, int length, T value);

  /** Fills the passed array with all zeros - the type must have a constructor that takes an int
  and initializes to the zero element when 0 is passed. */
  template <class T>
  static void fillWithZeros(T *buffer, int length);

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

  /** \todo check and comment this function - maybe move it to RSLib */
  template <class T>
  static void filterBiDirectional(const T *x, int xLength, T *y, int yLength, 
    const T *b, int bOrder, const T *a, int aOrder, int numRingOutSamples = 10000);

  /** Returns the index of the first value that matches the elementToFind, return -1 when no 
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

  /** Given a buffer of values, this function returns the first index where there's a nonzero value
  in the buffer. If there's no nonzero element at all, it returns -1.

  \todo use a different convention: return "length", when no index is found - this allows to return
  an unsigned int and extends the range where ths function can be used by a factor of (almost) 2

  */
  template <class T>
  static int firstIndexWithNonZeroValue(const T *buffer, int length);

  /** Scales and offsets the passed buffer such that the minimum value hits 'min' and the
  maximum value hits 'max'. */
  template <class T>
  static void fitIntoRange(T *buffer, int length, T min, T max);

  /** Fills the array h with the impulse of the filter specified by the direct form coefficients
  given in b and a. \todo comment on the sign of the a-coeffs, whether or not an a0 is included,
  etc., maybe move to RSLib into the filter section */
  template <class T>
  static void impulseResponse(T *h, int hLength, const T *b, int bOrder, const T *a, int aOrder);

  /** Interleaves a buffer of non-interleaved data. */
  template <class T>
  static void interleave(T *buffer, int numFrames, int numElementsPerFrame);

  /** Returns a linearly interpolated value from the array at the given (non-integer) position. If
  the position is out of range, 0 is returned. */
  template<class T>
  static T interpolatedValueAt(const T *buffer, int length, double position);

  /** Returns a linearly interpolated value from the array at the given (non-integer) position. If
  the position is out of range, the output is clamped to the end values of the array. */
  template<class T>
  static T interpolateClamped(const T *buffer, int length, double position);

  /** Returns true, if the passed buffer has only zero values, false otherwise. */
  template <class T>
  static inline bool isAllZeros(const T *buffer, int length);

  /** Returns true, if the passed buffer has values equal to the given value, false otherwise. */
  template <class T>
  static inline bool isFilledWithValue(const T *buffer, int length, T value);

  /** Returns true if the value in array x at position n is larger than its left and right 
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

  /** Shifts the content of the buffer numPlaces to the left, filling it up with zeros from the
  right. */
  template <class T>
  static void leftShift(T *buffer, int length, int numPlaces);

  /** Limits the passed value into the range between min and max (inclusive) and returns the
  result. */
  template <class T>
  static T limitToRange(T value, T min, T max);

  /** Finds and returns the maximum absolute value of the buffer. */
  template <class T>
  static T maxAbs(const T *buffer, int length);

  /** Finds and returns the index with the maximum absolute value of the buffer. */
  template <class T>
  static int maxAbsIndex(const T* const buffer, int length);

   /** Returns the maximum deviation (absolute value of the difference) between two buffers. */
  template <class T>
  static T maxDeviation(const T *buffer1, const T *buffer2, int length);

  /** Returns the index of the maximum deviation (see maxDeviation). */
  template <class T>
  static int maxDeviationIndex(const T *buffer1, const T *buffer2, int length);

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

  /** Returns the median of the passed buffer. */
  template <class T>
  static T median(const T *buffer, int length);

  /** Multiplies the elements of 'buffer1' and 'buffer2' - type must define operator '*'. The
  'result' buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T1, class T2, class TR>
  static void multiply(const T1 *buffer1, const T2 *buffer2, TR *result, int length);

  /** Writes the element-wise negation of the source buffer into the destination buffer. */
  template<class T>
  static void negate(const T *source, T *destination, int length);

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

  /** Shifts the values in the 4-element array a one position back, discarding the last one a[3] 
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

  /** Fills array y with the reversed content of array x. x and y must be distinct and 
  non-overlapping. */
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
  static void scale(T1 *buffer, int length, T2 scaleFactor);

  /** Scales the "src" buffer by a constant factor and writes the result into the "dst" buffer. */
  template <class T1, class T2>
  static void scale(const T1 *src, T1 *dst, int length, T2 scaleFactor);

  /** Given the sequence y of length yLength, this function returns a sequence x which, when
  convolved with itself, gives y. yLength is assumed to be odd, the index of first nonzero value
  is assumed to be even and the first nonzero value in y is assumed to be positive (for real
  sequences), because this is what will happen, when covolving a (real) sequence with itself. To
  disambiguate the square-root, the function will return a sequence with its 1st nonzero value
  being positive. If the original sequence x (before it was convolved with itself to give y)
  started with a negative value, the result of taking the square-root of the squared sequence y
  will have all signs reversed with respect to the original sequence x (again, in the case of real
  sequences - \todo explain what happens in the complex case).
  The length of x will be (yLength+1)/2.

  \todo: this apparently works only if the y-sequence was actually constructed by squaring some
  sequence, for a general sequence y, terms of y and x^2 will match only up to the n-th term
  where n = (yLength+1)/2, i.e. the length of x. ...find out if this function is useless
  therefore and/or if there is a usefule generalization - like using a sequence x which has the
  same length as y, convolve it with itslef and truncate the result to the length of y */
  template <class T>
  static void sequenceSqrt(const T *y, int yLength, T *x);

  /** Shifts the content of the buffer numPlaces to the right, filling it up with zeros from the
  left. If numPlaces is negative, the contents will be shifted to the left, filling up with zeros 
  from the right. */
  template <class T>
  static void shift(T *buffer, int length, int numPlaces);

  /** Subtracts the elements of 'buffer2' from 'buffer1' - type must define operator '-'. The
  'result' buffer may be the same as 'buffer1' or 'buffer2'. */
  template <class T>
  static void subtract(const T *buffer1, const T *buffer2, T *result, int length);

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

  /** Swaps the contents of of buffer1 and buffer2 using an auxiliary buffer bufferTmp. All buffers
  are assumed to have a size of sizeInBytes. */
  static void swapDataBuffers(void *buffer1, void *buffer2, void *bufferTmp, int sizeInBytes);

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
  // move to MatrixFunctions

  /** Given an array that contains values that have been subject to some kind wrap-around into some
  period p, this function (heuristically) unwraps them again. This is useful for plotting a 
  phase-response - if you have an array of angles in the interval -pi..pi (as you will get when 
  using arg(std::complex z) on a transfer function value), subtract off 2*pi, and unwrap with 
  period 2*pi. */
  template<class T>
  static void unwrap(T* a, int N, T p);

  /** Returns the sum-over-i w[i]*x[i]. */
  template <class T>
  static T weightedSum(const T *w, const T *x, rsUint32 length);  // redundant with sumOfProducts

  /** Forms a weighted sum of the two buffers. */
  template <class T>
  static  void weightedSum(const T *buffer1, const T *buffer2, T *result, int length, T weight1, T weight2);

};

//=================================================================================================
// inlined implementations

template <class T>
inline bool rsArray::almostEqual(const T *buffer1, const T *buffer2, 
  const int length, const T tolerance)
{
  for(int i = 0; i < length; i++)
  {
    if(rsAbs(buffer1[i]-buffer2[i]) > tolerance)
      return false;
  }
  return true;
}

template <class T>
inline bool rsArray::equal(const T *buffer1, const T *buffer2, const int length)
{
  for(int i = 0; i < length; i++)
  {
    if(buffer1[i] != buffer2[i])
      return false;
  }
  return true;
}

template <class T1, class T2>
inline void rsArray::convertBuffer(const T1 *source, T2 *destination, const int length)
{
  for(int i = 0; i < length; i++)
    destination[i] = (T2)source[i];
}

template <class T>
inline void rsArray::convolveWithTwoElems(const T* x, const int xLength, const T* h, T* y)
{
  y[xLength] = x[xLength-1]*h[1];
  for(int n = xLength-1; n > 0; n--)
    y[n] = x[n]*h[0] + x[n-1]*h[1];
  y[0] = x[0]*h[0];
}

template <class T>
inline void rsArray::convolveWithTwoElems(
  const T* x, const int xLength, const T h0, const T h1, T* y) 
{
  y[xLength] = x[xLength-1]*h1;
  for(int n = xLength-1; n > 0; n--)
    y[n] = x[n]*h0 + x[n-1]*h1;
  y[0] = x[0]*h0;
}

template <class T>
inline void rsArray::fillWithImpulse(T *buffer, int length)
{
  fillWithZeros(buffer, length);
  buffer[0] = T(1);
}

template <class T>
inline void rsArray::fillWithNaN(T* x, int N)
{
  for(int i = 0; i < N; i++) 
    x[i] = RS_SIGNALING_NAN(T); // or maybe it should be the quiet nan?
}

template <class T>
inline int rsArray::findIndexOf(const T *buffer, T elementToFind, int length)
{
  for(int i = 0; i < length; i++)
  {
    if( buffer[i] == elementToFind )
      return i;
  }
  return -1;
}

template <class T>
inline int rsArray::findMaxAbs(const T *buffer, int length)
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

template <class T>
inline bool rsArray::isAllZeros(const T *buffer, int length)
{
  return isFilledWithValue(buffer, length, T(0));
}

template <class T>
inline bool rsArray::isFilledWithValue(const T *buffer, int length, T value)
{
  for(int n = 0; n < length; n++)
  {
    if(buffer[n] != value)
      return false;
  }
  return true;
}

template<class T>
bool rsArray::isPeakOrPlateau(const T *x, int n)
{
  if(x[n] >= x[n-1] && x[n] >= x[n+1])
    return true;
  return false;
}

template<class T>
inline void rsArray::pushFrontPopBack4(T x, T* a)
{
  a[3] = a[2];  // old a[3] is discarded
  a[2] = a[1];  // values in between..
  a[1] = a[0];  // ...are shifted
  a[0] = x;     // input x becomes the new a[0]
}
// todo: make versions for 1,2,3,N (using a loop or memmove) -> make performance tests, 
// which version is fastest for what range of lengths maybe rename to updateFifoBuffer4


#endif