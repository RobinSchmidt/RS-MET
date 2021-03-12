//-------------------------------------------------------------------------------------------------
// preliminary - to fix compiler error on linux/gcc (move to somewhere else)

inline int rsFloorInt(double x)
{
  return int(x);
}

inline unsigned long rsBitReverse(unsigned long number, unsigned long numBits)
{
  unsigned long result = 0;
  for(unsigned long n=0; n<numBits; n++)
  {
    // leftshift the previous result by one and accept the new LSB of the current number on the
    // right:
    result   = (result << 1) + (number & 1);

    // rightshift the number to make the second bit from the right to the new LSB:
    number >>= 1;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------

//template <class T>
//void rsArrayTools::add(const T *buffer1, const T *buffer2, T *result, const int length)
//{
//  for(int i = 0; i < length; i++)
//    result[i] = buffer1[i] + buffer2[i];
//}
//
//template <class T>
//void rsArrayTools::add(const T *buffer, const T valueToAdd, T *result, const int length)
//{
//  for(int i = 0; i < length; i++)
//    result[i] = buffer[i] + valueToAdd;
//}

template <class T>
void rsArrayTools::addCircularShiftedCopy(
  T *buffer, const int length, const double offset, const T weight)
{
  T *tmp = new T[length];
  copy(buffer, tmp, length);
  circularShiftInterpolated(tmp, length, offset);
  scale(tmp, length, weight);
  add(buffer, tmp, buffer, length);
  delete[] tmp;
}

template<class T>
void rsArrayTools::addInto(T *x, const int N, const T *y, int L, int n)
{
  int r = 0;                // read start
  if(n < 0)
  {
    L += n;
    r -= n;
    n  = 0;
  }
  const int d = n + L - N;        // number of overhanging values
  if(d > 0)
    L -= d;
  for(int i = 0; i < L; i++)
    x[n+i] += y[r+i];
}
// the index manipulation code can be factored out

template<class T>
void rsArrayTools::affineTrafo(const T* x, T* y, const int N, const T a, const T b)
{
  for(int i = 0; i < N; i++)
    y[i] = a * x[i] + b;
}
// maybe inline this

template<class T>
void rsArrayTools::allocateSquareArray2D(T**& theArray, const int size)
{
  theArray = new T*[size];
  for(int i = 0; i < size; i++)
    theArray[i] = new T[size];
}

template <class T>
void rsArrayTools::applyFunction(const T *inBuffer, T *outBuffer, const int length, T (*f) (T))
{
  for(int i = 0; i < length; i++)
    outBuffer[i] = f(inBuffer[i]);
}

template <class T>
void rsArrayTools::circularShift(T *buffer, const int length, const int numPositions)
{
  int na = abs(numPositions);
  while(na > length)
    na -=length;
  T *tmp = new T[na];
  if(numPositions < 0)
  {
    memcpy(tmp, buffer, na*sizeof(T));
    memmove(buffer, &buffer[na], (length-na)*sizeof(T));
    memcpy(&buffer[length-na], tmp, na*sizeof(T));
  }
  else if(numPositions > 0)
  {
    memcpy(tmp, &buffer[length-na], na*sizeof(T));
    memmove(&buffer[na], buffer, (length-na)*sizeof(T));
    memcpy(buffer, tmp, na*sizeof(T));
  }
  delete[] tmp;
}
// this article describes an algo ("rotate") based on reversals that doesn't need auxiuliary
// memory: https://en.wikipedia.org/wiki/Block_sort - there's an implementation in Prototypes.cpp

template <class T>
void rsArrayTools::circularShiftInterpolated(T *buffer, const int length, const double numPositions)
{
  const double read = rsWrapAround(numPositions, (double)length);
  int w = 0;                       // write position
  int r = rsFloorInt(read);        // integer part of read position
  const double f  = read-r;        // fractional part of read position
  const double f2 = 1.0-f;
  T *tmp = new T[length];
  copy(buffer, tmp, length);
  while(r < length-1)
  {
    buffer[w] = f2*tmp[r] + f*tmp[r+1];
    r++;
    w++;
  }
  buffer[w] = f2*tmp[r] + f*tmp[0];
  r=0;
  w++;
  while(w < length)
  {
    buffer[w] = f2*tmp[r] + f*tmp[r+1];
    r++;
    w++;
  }
  delete[] tmp;
}

template <class T>
void rsArrayTools::clip(T *buffer, const int length, const T min, const T max)
{
  for(int i = 0; i < length; i++) {
    if(buffer[i] < min)
      buffer[i] = min;
    else if(buffer[i] > max)
      buffer[i] = max;
  }
}

template <class T>
int rsArrayTools::compare(const T *a, const T *b, const int length)
{
  for(int i = 0; i < length; i++)
  {
    if(a[i] < b[i])
      return -1;
    if(a[i] > b[i])
      return +1;
  }
  return 0;
}

template <class T>
int rsArrayTools::compare(const T *a, int na, const T *b, const int nb)
{
  const int nMin = rsMin(na, nb);
  for(int i = 0; i < nMin; i++)
  {
    if(a[i] < b[i])
      return -1;
    if(a[i] > b[i])
      return +1;
  }
  if(na > nb)
  {
    if(!rsIsAllZeros(&a[nMin], na-nMin))
      return +1;
  }
  else if(nb > na)
  {
    if(!rsIsAllZeros(&b[nMin], nb-nMin))
      return -1;
  }
  return 0;
}

/*
template <class T>
bool rsArrayTools::contains(const T *buffer, const int length, const T elementToFind)
{
  for(int i = 0; i < length; i++)
  {
    if(buffer[i] == elementToFind)
      return true;
  }
  return false;
  //return (rsFindFirstOccurrenceOf(buffer, length, elementToFind) != -1);
}
*/

template <class T>
void rsArrayTools::convolve(const T *x, const int xLength, const T *h, const int hLength, T *y)
{
  for(int n = xLength+hLength-2; n >= 0; n--) {
    T s = T(0);
    for(int k = rsMax(0, n-xLength+1); k <= rsMin(hLength-1, n); k++)
      s += h[k] * x[n-k];
    y[n] = s;
  }
}
// maybe optimize by getting rid of calling rsMin/rsMax in the head of the inner loop by splitting
// the outer loop into 3 partial loops....maybe

//template <class T1, class T2>
//void rsArrayTools::copy(const T1 *source, T2 *destination, const int length)
//{
//  for(int i = 0; i < length; i++)
//    destination[i] = (T2)source[i];
//}

// old version without type conversion:
//template <class T>
//void copy(const T *source, T *destination, int length)
//{
//  for(int i = 0; i < length; i++)
//    destination[i] = source[i];
//  // \todo: use memcpy - or maybe not - memcpy won't make deep copies
//}

template <class T>
void rsArrayTools::convolveInPlace(T *x, const int xLength, const T *h, const int hLength)
{
  convolve(x, xLength, h, hLength, x);
}

template <class T>
void rsArrayTools::copyBufferWithLinearInterpolation(const T *source, int sourceLength, T *destination,
  int destinationLength)
{
  double increment    = (double)sourceLength / (double)destinationLength;
  double readPosition = 0.0;
  for(int n = 0; n < destinationLength; n++)
  {
    int    iL = rsFloorInt(readPosition);   // integer part (left index)
    double f  = readPosition-iL;          // fractional part
    int    iR = iL+1;                     // right index

    // wraparound indices, if necesarry:
    while(iL >= sourceLength)
      iL -= sourceLength;
    while(iR >= sourceLength)
      iR -= sourceLength;

    // do the linear interpolation:
    destination[n]  = (1.0-f)*source[iL] + f*source[iR];

    readPosition   += increment;
  }
}

template <class T>
int rsArrayTools::copyIfMatching(const T *sourceBuffer, T *targetBuffer, int sourceAndTargetLength,
  const T *elementsToMatch, int matchLength)
{
  int numWritten = 0;
  for(int i = 0; i < sourceAndTargetLength; i++)
  {
    if(rsArrayTools::contains(elementsToMatch, matchLength, sourceBuffer[i]))
    {
      targetBuffer[numWritten] = sourceBuffer[i];
      numWritten++;
    }
  }
  return numWritten;
}

template <class T>
int rsArrayTools::copyIfNotMatching(const T *sourceBuffer, T *targetBuffer, int sourceAndTargetLength,
  const T *elementsToStrip, int stripLength)
{
  int numWritten = 0;
  for(int i = 0; i < sourceAndTargetLength; i++)
  {
    if(!rsArrayTools::contains(elementsToStrip, stripLength, sourceBuffer[i]))
    {
      targetBuffer[numWritten] = sourceBuffer[i];
      numWritten++;
    }
  }
  return numWritten;
}

template<class T1, class T2>
void rsArrayTools::copySection(const T1 *source, int sourceLength, T2 *destination, int copyStart,
  int copyLength)
{
  int cl, pl1, pl2;  // actual copy-, pre-padding-, post-padding-lengths
  if(copyStart >= 0)
  {
    // copying:
    cl = rsMin(copyLength, sourceLength-copyStart);
    convert(&source[copyStart], destination, cl);

    // post-padding:
    pl2 = copyLength-cl;
    fillWithZeros(&destination[cl], pl2);
  }
  else
  {
    // pre-padding:
    pl1 = rsMin(-copyStart, copyLength);
    fillWithZeros(destination, pl1);

    // copying:
    cl = rsMin(copyLength-pl1, sourceLength);
    convert(source, &destination[pl1], cl);

    // post-padding:
    pl2 = copyLength-cl-pl1;
    fillWithZeros(&destination[pl1+cl], pl2);
  }
}

//template<class T>
//void rsArrayTools::rsCopySection(T *source, int sourceLength, T *destination, int copyStart, int copyLength)
//{
//  int cl, pl1, pl2;  // actual copy-, pre-padding-, post-padding-lengths
//  if( copyStart >= 0 )
//  {
//    // copying:
//    cl = rsMin(copyLength, sourceLength-copyStart);
//    copy(&source[copyStart], destination, cl);

//    // post-padding:
//    pl2 = copyLength-cl;
//    fillWithZeros(&destination[cl], pl2);
//  }
//  else
//  {
//    // pre-padding:
//    pl1 = rsMin(-copyStart, copyLength);
//    fillWithZeros(destination, pl1);

//    // copying:
//    cl = rsMin(copyLength-pl1, sourceLength);
//    copy(source, &destination[pl1], cl);

//    // post-padding:
//    pl2 = copyLength-cl-pl1;
//    fillWithZeros(&destination[pl1+cl], pl2);
//  }
//}

template <class T>
void rsArrayTools::cumulativeSum(const T *x, T *y, int N)
{
  y[0] = x[0];
  for(int n = 1; n < N; n++)
    y[n] = x[n] + y[n-1];
}

template <class T>
void rsArrayTools::cumulativeSum(const T *x, T *y, int N, int order)
{
  copy(x, y, N);
  for(int i = 1; i <= order; i++)
    cumulativeSum(y, y, N);
}

template<class T>
void rsArrayTools::deAllocateSquareArray2D(T**& theArray, int size)
{
  for(int i = 0; i < size; i++)
    delete[] theArray[i];
  delete[] theArray;
}

template <class T>
int rsArrayTools::firstIndexWithNonZeroValue(const T *buffer, int N)
{
  for(int i = 0; i < N; i++)
  {
    if(buffer[i] != T(0))
      return i;
  }
  return -1;
}

template <class T>
void rsArrayTools::fillWithZeros(T *buffer, int length)
{
  for(int i = 0; i < length; i++)
    buffer[i] = T(0);
}

template <class T>
void rsArrayTools::deConvolve(const T *y, int yLength, const T *h, int hLength, T *x)
{
  int m = firstIndexWithNonZeroValue(h, hLength);
  if(m == -1)
  {
    // h is all zeros - return an all-zero x-signal:
    fillWithZeros(x, yLength-hLength+1);
    return;
  }
  T scaler = T(1) / h[m];
  x[0]     = scaler * y[m];
  for(int n = 1; n < yLength-hLength+1; n++)
  {
    x[n] = y[n+m];
    for(int k = m+1; k <= rsMin(hLength-1, n+m); k++)
      x[n] -= h[k] * x[n-k+m];
    x[n] *= scaler;
  }
  // Maybe this can be generalized such that the result x may be of arbitrary length? Here, we
  // assume xLength = yLength-hLength+1 but actually (i think), the sequence x is only finite if
  // y is indeed a convolution product of some signal x with h (or, if the polynomial
  // represented by the y array is divisible by the polynomial in the h array). If y and h are
  // arbitrary, the sequence x might be infinite. That means, we could want to calculate more
  // values. On the other hand, sometimes we may be interested only in the first few values of
  // x because the later ones are irrelevant for our problem, in which case we want to compute
  // less values. So, let's introduce an optional parameter xLength that defaults to zero (in
  // which case it's taken to be calculated as yLength-hLength+1). Then the loop becomes:
  // if( xLength == 0 )
  //   xLength = yLength-hLength+1;
  // for(int n = 1; n < xLength; n++)
}

template <class T>
void rsArrayTools::deInterleave(T *buffer, int numFrames, int numElementsPerFrame)
{
  T *tmp = new T[numFrames*numElementsPerFrame];
  int i, j;
  for(i = 0; i < numFrames*numElementsPerFrame; i++)
    tmp[i] = buffer[i];  // \todo use copy
  for(j = 0; j < numElementsPerFrame; j++)
  {
    int k = numFrames*j;
    for(i = 0; i < numFrames; i+=1)
      buffer[k+i] = tmp[numElementsPerFrame*i+j];
  }
  delete[] tmp;
}

template <class T>
void rsArrayTools::difference(T *buffer, int length, int order, bool periodic)
{
  T x, x1; // for temporary storage of the x, x[n-1] samples
  for(int o = 1; o <= order; o++) {
    if(periodic) x1 = buffer[length-1];
    else         x1 = T(0);
    for(int n = 0; n < length; n++) {
      x         = buffer[n];
      buffer[n] = x - x1;    // y[n] = x[n] - x[n-1]
      x1        = x;  }}
}

template <class T1, class T2, class TR>
void rsArrayTools::divide(const T1 *buffer1, const T2 *buffer2, TR *result, int length)
{
  for(int i = 0; i < length; i++)
    result[i] = buffer1[i] / buffer2[i];
}

template <class T>
void rsArrayTools::fillWithIndex(T *buffer, int length)
{
  for(int i = 0; i < length; i++)
    buffer[i] = T(i);
}

template <class T>
void rsArrayTools::fillWithRandomValues(T *buffer, int length, double min, double max, int seed)
{
  rsRandomUniform(min, max, seed);
  for(int i = 0; i < length; i++)
    buffer[i] = (T)rsRandomUniform(min, max, -1);
}

template <class T>
void rsArrayTools::fillWithRandomIntegers(T *buffer, int length, int min, int max, int seed)
{
  fillWithRandomValues(buffer, length, T(min)-T(0.5), T(max)+T(0.5), seed);
  for(int i = 0; i < length; i++)
  {
    buffer[i] = rsRound(buffer[i]);
    if(buffer[i] < T(min)) buffer[i] = T(min);
    if(buffer[i] > T(max)) buffer[i] = T(max);
  }
  // The +-0.5 is for making the min/max values equally likely to the "inner" values - but then, 
  // the rounding could produce numbers min-1, max+1 (in very unlikely edge cases, but possible 
  // anyway), so we need to catch these.
}

template <class T>
void rsArrayTools::fillWithValue(T *buffer, int length, T value)
{
  for(int i = 0; i < length; i++)
    buffer[i] = value;
}

template <class T>
void rsArrayTools::fillWithRangeExponential(T *x, int N, T min, T max)
{
  rsAssert(N >= 2);
  if(min == max)
    fillWithValue(x, N, min);
  else {
    for(int i = 0; i < N; i++)
      x[i] = (T)rsLinToExp((T)i, (T)0, (T)(N-1), min, max); }
}
// ...can(?) be optimized to avoid calling rsLinToExp for each i by precomputing some variables

template <class T>
void rsArrayTools::fillWithRangeLinear(T *x, int N, T min, T max)
{
  rsAssert(N >= 2);
  if(min == max)
    fillWithValue(x, N, min);
  else {
    T factor = (max-min) / (T)(N-1);
    for(int i = 0; i < N; i++)
      x[i] = (T)(factor * T(i) + min); }
}

template<class T>
void rsArrayTools::fillWithRange(T* x, int N, T min, T max, T p)
{
  static const T tol = (T)pow(2, 10) * RS_EPS(T); // found empirically
  if(rsAbs(p) <= tol) {
    fillWithRangeExponential(x, N, min, max);
    return; }
  fillWithRangeLinear(x, N, pow(min, p), pow(max, p));
  for(int i = 0; i < N; i++)
    x[i] = pow(x[i], T(1)/p);
}
// maybe we should set x[0] = min and x[N-1] = max to avoid roundoff error at the endpoints - the
// loop should run only for(i = 1; i < N-1; i++)...also in the inner functions...but this needs
// tests - especially when called with N=0 and N=1


template <class T>
void rsArrayTools::filter(const T *x, int xLength, T *y, int yLength, const T *b, int bOrder, const T *a, int aOrder)
{
  // allocate and intitialize memory for the filters internal state:
  int i, n;
  T *xOld = new T[bOrder+1];
  T *yOld = new T[aOrder+1];
  for(i=0; i<=bOrder; i++) xOld[i] = T(0);
  for(i=0; i<=aOrder; i++) yOld[i] = T(0);

  // compute the part of the signal where both buffers have values:
  int length = rsMin(xLength, yLength);
  T tmp;
  for(n = 0; n < length; n++)
  {
    // compute y[n]:
    tmp = b[0] * x[n];
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];

    // update state buffers:
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = x[n];
    yOld[1] = tmp;

    // store y[n]:
    y[n] = tmp;
  }

  // possibly compute the tail-part of y when the y-buffer is longer than x-buffer:
  for(int n = length; n < yLength; n++)
  {
    tmp = T(0);
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = T(0);
    yOld[1] = tmp;
    y[n]    = tmp;
  }

  // todo: we should divide the final value by a[0]

  // clean up memory:
  delete[] xOld;
  delete[] yOld;
}

template <class T>
void rsArrayTools::filterBiDirectional(const T *x, int xLength, T *y, int yLength, const T *b, int bOrder, const T *a,
  int aOrder, int numRingOutSamples)
{
  // allocate and intitialize memory for the filters internal state:
  int i, n;
  T tmp;
  T *xOld    = new T[bOrder+1];
  T *yOld    = new T[aOrder+1];
  T *ringOut = new T[numRingOutSamples];
  for(i=0; i<=bOrder; i++) xOld[i] = T(0);
  for(i=0; i<=aOrder; i++) yOld[i] = T(0);

  /*
  // backward pass through a portion of the x-buffer to warm up the filter:
  for(n=rsMin(xLength-1, 10000); n>=0; n--)
  {
    // compute y[n]:
    tmp = b[0] * x[n];
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];

    // update state buffers:
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = x[n];
    yOld[1] = tmp;
  }
  */

  // compute the part of the signal where both buffers have values:
  int length = rsMin(xLength, yLength);
  for(n = 0; n < length; n++)
  {
    tmp = b[0] * x[n];
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = x[n];
    yOld[1] = tmp;
    y[n]    = tmp;
  }

  // possibly compute the tail-part of y when the y-buffer is longer than x-buffer:
  for(int n = length; n < yLength; n++)
  {
    tmp = T(0);
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = T(0);
    yOld[1] = tmp;
    y[n]    = tmp;
  }

  // compute the ringout tail:
  for(n = 0; n < numRingOutSamples; n++)
  {
    tmp = T(0);
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1]    = y[n];
    yOld[1]    = tmp;
    ringOut[n] = tmp;
  }

  // backward pass through the ringout tail:
  for(n = numRingOutSamples-1; n>=0; n--)
  {
    tmp = b[0] * ringOut[n];
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = y[n];
    yOld[1] = tmp;
  }

  // backward pass through the y-buffer:
  for(n = yLength-1; n >= 0; n--)
  {
    tmp = b[0] * y[n];
    for(i=1; i<=bOrder; i++) tmp += b[i] * xOld[i];
    for(i=1; i<=aOrder; i++) tmp -= a[i] * yOld[i];
    for(i=bOrder; i>=2; i--) xOld[i] = xOld[i-1];
    for(i=aOrder; i>=2; i--) yOld[i] = yOld[i-1];
    xOld[1] = y[n];
    yOld[1] = tmp;
    y[n]    = tmp;
  }

  // \todo: get rid of the code duplication - the loop could be factored out into a function
  //

  // clean up memory:
  delete[] xOld;
  delete[] yOld;
  delete[] ringOut;
}

template <class T>
void rsArrayTools::fitIntoRange(T *buffer, int length, T min, T max)
{
  T range        = max-min;
  T currentMin   = minValue(buffer, length);
  T currentMax   = maxValue(buffer, length);
  T currentrsRange = currentMax - currentMin;
  if(currentrsRange == T(0))
    return; // avoid divide-by-zero
  T scaleFactor  = range/currentrsRange;
  T offset       = min - scaleFactor*currentMin;
  scale(buffer, length, scaleFactor);
  add(buffer, offset, buffer, length);
}

template<class T>
T rsArrayTools::geometricMean(const T* x, int N)
{
  T m = T(0);
  for(int i = 0; i < N; i++)
    m += log(x[i]);
  m /= T(N);
  m = exp(m);
  return m;
}
// I think, it's numerically better behaved to sum up the logarithms and take the exp at the and 
// rather than accumulating multiplicatively and extracting the N-th root - but that should be 
// verified experimentally.

template<class T>
T rsArrayTools::generalizedMean(const T* x, int N, T p)
{
  static const T tol = (T)pow(2, 10) * RS_EPS(T); // found empirically
  if(rsAbs(p) <= tol) 
    return geometricMean(x, N);
  T m = T(0);
  for(int i = 0; i < N; i++)
    m += pow(x[i], p);
  m /= T(N);
  m = pow(m, T(1)/p);
  return m;
}
// Maybe use m += exp(p * log(x[i]) - might be more effient than pow - but also less precise. Maybe
// take the min or max for abs(p) > thresh for some threshold...but that threshold must also be
// obtained empirically...and reasonable values may depend on the array values

template <class T>
void rsArrayTools::impulseResponse(T *h, int hLength, const T *b, int bOrder, const T *a, int aOrder)
{
  T x = T(1);
  filter(&x, 1, h, hLength, b, bOrder, a, aOrder);
}

template <class T>
void rsArrayTools::interleave(T *buffer, int numFrames, int numElementsPerFrame)
{
  T *tmp = new T[numFrames*numElementsPerFrame];
  int i, j;
  for(i=0; i<numFrames*numElementsPerFrame; i++)
    tmp[i] = buffer[i];  // \todo use copy
  for(j = 0; j < numElementsPerFrame; j++) {
    int k = numFrames*j;
    for(i = 0; i < numFrames; i+=1)
      buffer[numElementsPerFrame*i+j] = tmp[k+i]; }
  delete[] tmp;
}

template<class T>
T rsArrayTools::interpolatedValueAt(const T *x, int N, double n)
{
  if(n < 0.0 || n > N-1)
    return 0.0;
  else {
    int ni = (int)floor(n);            // integer part of n
    if(ni == N-1)
      return x[ni];
    double nf = n - ni;                 // fractional part of n
    return (1.0-nf)*x[ni] + nf*x[ni+1]; // linear interpolation
  }
}

template<class T>
T rsArrayTools::interpolateClamped(const T *y, int N, double n)
{
  if(n <= 0.0)
    return y[0];
  if(n >= N-1)
    return y[N-1];
  int ni = (int)n;
  double nf = n - ni;
  return (1-nf)*y[ni] + nf*y[ni+1];
}

template<class T>
bool rsArrayTools::isPeak(const T *x, int n)
{
  if(x[n] > x[n-1] && x[n] > x[n+1])
    return true;
  return false;
}

template<class T>
bool rsArrayTools::isValley(const T *x, int n)
{
  if(x[n] < x[n-1] && x[n] < x[n+1])
    return true;
  return false;
}

template<class T>
bool rsArrayTools::isPeakOrValley(const T *x, int n)
{
  return isPeak(x, n) || isValley(x, n);
}

template<class T>
int rsArrayTools::findPeakOrValleyRight(const T *x, int N, int n0)
{
  int nR = -1;
  for(int n = rsMax(1, n0); n < N-1; n++) {
    if(isPeakOrValley(x, n)) {
      nR = n;
      break; }}
  return nR;
}

template<class T>
int rsArrayTools::findPeakOrValleyLeft(const T *x, int N, int n0)
{
  int nL = -1;
  for(int n = rsMin(N-2, n0); n > 0; n--) {
    if(isPeakOrValley(x, n)) {
      nL = n;
      break; }}
  return nL;
}

/*
// moved to .h file
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
*/

template<class T>
int rsArrayTools::findSplitIndexClosest(const T* a, const int N, const T val)
{
  int i = findSplitIndex(a, N, val);

  if(i == N)     // new - needs tests
    return N-1;

  if(i > 0 && rsAbs(a[i]-val) > rsAbs(a[i-1]-val))
    i--;
  return i;
}

/* moved to .h
template <class T>
bool rsArrayTools::isSortedAscending(const T *buffer, int length)
{
  for(int i = 0; i < length-1; i++) {
    if(!(buffer[i] <= buffer[i+1]))
      return false; }
  return true;
}
*/

template <class T>
bool rsArrayTools::isSortedStrictlyAscending(const T *buffer, int length)
{
  for(int i = 0; i < length-1; i++) {
    if(!(buffer[i] < buffer[i+1]))
      return false; }
  return true;
}

template <class T>
void rsArrayTools::leftShift(T *buffer, int length, int numPlaces)
{
  //rsAssert(numPlaces >= 0 && numPlaces <= length, "we require 0 <= numPlaces <= length ");
  int i;
  for(i = 0; i < length-numPlaces; i++)
    buffer[i] = buffer[i+numPlaces];
  fillWithZeros(&buffer[i], numPlaces);
}

template <class T>
T rsArrayTools::limitToRange(T value, T min, T max)
{
  if(     value < min)  return min;
  else if(value > max)  return max;
  else                  return value;
}

template <class T>
void rsArrayTools::maxElementWise(const T* x, const T* y, const int N, T* maxXY)
{
  for(int n = 0; n < N; n++)
    maxXY[n] = rsMax(x[n], y[n]);
}

template <class T>
T rsArrayTools::maxAbs(const T *buffer, int length)
{
  T max = T(0);
  for(int i = 0; i < length; i++) {
    if(rsAbs(buffer[i]) > max)
      max = rsAbs(buffer[i]); }
  return max;
}
// maybe optimize such that rsAbs is called only once

template <class T>
T rsArrayTools::maxAbs(const std::complex<T>* buffer, int length)
{
  T maxSquared = T(0);
  for(int i = 0; i < length; ++i) {
    T absSquared = rsAbsSquared(buffer[i]);
    if(absSquared > maxSquared)
      maxSquared = absSquared; }
  return sqrt(maxSquared);
}


template <class T>
int rsArrayTools::maxAbsIndex(const T* const buffer, int length)
{
  T max = T(0);
  int maxIndex = 0;
  for(int i = 0; i < length; i++) {
    if(rsAbs(buffer[i]) > max) {
      max      = rsAbs(buffer[i]);
      maxIndex = i; }}
  return maxIndex;
}

template <class T>
T rsArrayTools::maxDeviation(const T *buffer1, const T *buffer2, int length)
{
  T max = T(0);
  for(int i = 0; i < length; i++)
  {
    T absDiff = rsAbs(buffer1[i]-buffer2[i]);
    if(absDiff > max)
      max = absDiff;
  }
  return max;
}

template <class T>
int rsArrayTools::maxDeviationIndex(const T *x, const T *y, int N)
{
  T maxErr = T(0);  // rename to maxDev
  int maxIdx = 0;
  for(int i = 0; i < N; i++) {
    T err = rsAbs(x[i]-y[i]);
    if(err > maxErr) {
      maxErr = err;
      maxIdx = i; }}
  return maxIdx;
}

template <class T>
T rsArrayTools::maxDifference(const T* x, int N)
{
  T x1   = T(0);
  T dMax = -RS_INF(T);
  for(int i = 0; i < N; i++) {
    T d = x[i] - x1;;
    if(d > dMax)
      dMax = d;
    x1 = x[i];
  }
  return dMax;
}

template <class T>
int rsArrayTools::maxIndex(const T *buffer, int length)
{
  T   value = buffer[0];
  int index = 0;
  for(int i = 0; i < length; i++) {
    if(buffer[i] > value) {
      value = buffer[i];
      index = i; }}
  return index;
}

template <class T>
T rsArrayTools::maxValue(const T *buffer, int length)
{
  return buffer[maxIndex(buffer, length)];
}

template <class T>
int rsArrayTools::minIndex(const T *buffer, int length)
{
  T   value = buffer[0];
  int index = 0;
  for(int i = 0; i < length; i++) {
    if(buffer[i] < value) {
      value = buffer[i];
      index = i; }}
  return index;
}

template <class T>
T rsArrayTools::minValue(const T *buffer, int length)
{
  return buffer[minIndex(buffer, length)];
}

template <class T>
T rsArrayTools::mean(const T *buffer, int length)
{
  return sum(buffer, length) / (T)length;
}

template <class T>
T rsArrayTools::meanDifference(const T *x, int N)
{
  T s = 0;              // sum (of differences)
  for(int i = 1; i < N; i++)
    s += x[i] - x[i-1];
  return s / T(N-1);    // for N values, there are N-1 differences
}

template<class T>
int rsArrayTools::numNonZeros(const T* x, int N)
{
  int n = 0;
  for(int i = 0; i < N; i++)
    if(x[i] != T(0))
      n++;
  return n;
}

template<class T>
T rsArrayTools::meanSquare(const T *x, int N)
{
  return sumOfSquares(x, N) / T(N);
}

template <class T>
T rsArrayTools::median(const T *buffer, int length)
{
  T* tmpBuffer = new T[length];
  copy(buffer, tmpBuffer, length);

  std::sort(tmpBuffer, &tmpBuffer[length]);
  T med;
  if(rsIsOdd(length))
    med = tmpBuffer[(length-1)/2];
  else
  {
    //med = (T)(0.5 * (tmpBuffer[length/2] + tmpBuffer[length/2-1]));
    med = T(0.5) * (tmpBuffer[length/2] + tmpBuffer[length/2-1]);
  }
    // todo: use T(0.5) * (tmpBuffer[length/2] + tmpBuffer[length/2-1]) - don't enforce conversion
    // to double in intermeidate result

  delete[] tmpBuffer;
  return med;
}
// maybe use std::n_thelement, see here:
// https://www.youtube.com/watch?v=sWgDk-o-6ZE 33:25

template <class T>
void rsArrayTools::movingAverage3pt(const T* x, int N, T* y, bool endsFixed)
{
  rsAssert(N >= 0);
  if(N == 0) return;
  if(N == 1) { y[0] = x[0]; return; }
  T t1 = x[0];
  T t2 = x[1];

  if(endsFixed)
    y[0] = t1;
  else
    y[0] = T(1/2.) * (t1 + t2);

  for(int n = 2; n < N; n++) {
    y[n-1] = T(1/3.) * (t1 + t2 + x[n]);
    t1 = t2;
    t2 = x[n]; }

  if(endsFixed)
    y[N-1] = t2;
  else
    y[N-1] = T(1/2.) * (t1 + t2);
}

template<class T>
void rsArrayTools::movingMedian3pt(const T* x, int N, T* y)
{
  T x1 = x[0];  // x[n-1]
  for(int n = 1; n < N-1; n++) {
    T xn = x[n];
    y[n] = rsMedian(x1, xn, x[n+1]);
    x1   = xn; }
  y[0]   = y[1];    // does this make sense?
  y[N-1] = y[N-2];
  // i think, it would use the 3pt forward median for y[0] and the 3pt backward median for y[N-1]
  // the two-point mediat would be the same as the two point average (even order medians use the
  // average of the middle two samples)...hmmm...
  // or maybe we should use linear extrapolation?
}


template <class T1, class T2, class TR>
void rsArrayTools::multiply(const T1 *buffer1, const T2 *buffer2, TR *result, int length)
{
  for(int i = 0; i < length; i++)
  {
    result[i] = TR(buffer1[i] * buffer2[i]);
    // does not compile when T1 == std::complex<float> and T2 == double. This occurs in
    // rsPolynomial<T>::evaluateWithDerivatives when we call 
    // rsArrayTools::multiply(&results[2], &rsFactorials[2], &results[2], numDerivatives-1);
    // becasue the rsFactorials array is an array of doubles...

    //result[i] = TR(buffer1[i]) * TR(buffer2[i]);
    // ...so we need to do the conversion to the result type for both operands separately, not just
    // for the result itself
  }
}

template<class T>
void rsArrayTools::negate(const T *x, T *y, int N)
{
  for(int i = 0; i < N; i++)
    y[i] = -x[i];
}

template<class T>
void rsArrayTools::negateEven(const T *x, T *y, int N)
{
  for(int i = 0; i < N; i += 2) y[i] = -x[i];
  for(int i = 1; i < N; i += 2) y[i] =  x[i];
}

template<class T>
void rsArrayTools::negateOdd(const T *x, T *y, int N)
{
  for(int i = 0; i < N; i += 2) y[i] =  x[i];
  for(int i = 1; i < N; i += 2) y[i] = -x[i];
}

template <class T>
void rsArrayTools::normalize(T *buffer, int length, T maximum, bool subtractMean)
{
  if(subtractMean == true)
  {
    //T mean = mean(buffer, length);
    add(buffer, -mean(buffer, length), buffer, length);
  }
  T max   = maxAbs(buffer, length);;
  T scale = maximum / max;
  for(int i = 0; i < length; i++)
    buffer[i] *= scale;
}

template<class T>
void rsArrayTools::normalizeMean(T* x, int N, T newMean)
{
  T m = mean(x, N);
  scale(x, N, newMean/m);
}

template <class T>
void rsArrayTools::orderBitReversed(T *buffer, int N, int log2N)
{
  int n, nr; // index and bit-reversed index
  for(n = 0; n < N; n++)
  {
    nr = (int)rsBitReverse(n, log2N);
    if(n < nr)
      rsSwap(buffer[n], buffer[nr]);
  }
}

template <class T>
void rsArrayTools::orderBitReversedOutOfPlace(const T *inBuffer, T *outBuffer, int length, int numBits)
{
  // gather up the values from the bit-reversed positions in the input array:
  for(int n = 0; n < length; n++)
    outBuffer[n] = inBuffer[rsBitReverse(n, numBits)];
}

template <class T>
T rsArrayTools::product(const T* const buffer, int length)
{
  T accu = T(1); // constructor call with 1 should initilize to multiplicative neutral element
  for(int n = 0; n < length; n++)
    accu *= buffer[n];
  return accu;
}

template <class T>
void rsArrayTools::removeMean(T *buffer, int length)
{
  T m = mean(buffer, length);
  for(int i = 0; i < length; i++)
    buffer[i] -= m;
}

/*
// old:
template <class T>
void rsArrayTools::reverse(T *buffer, int length)
{
  T tmp;
  int lengthMinus1 = length-1;
  for(int i = 0; i <= (length-2)/2; i++)
  {
    tmp                    = buffer[lengthMinus1-i];
    buffer[lengthMinus1-i] = buffer[i];
    buffer[i]              = tmp;
  }
}
*/

// new:
template <class T>
void rsArrayTools::reverse(T* x, int N)
{
  for(int i = 0; i < N/2; i++)
    rsSwap(x[i], x[N-i-1]);
}


template <class T>
void rsArrayTools::reverse(const T* x, T* y, int N)
{
  if(x == y) { reverse(y, N); return; }  // reverse in place
  for(int i = 0; i < N; i++)
    y[i] = x[N-1-i];
}

template <class T>
void rsArrayTools::rightShift(T *buffer, int length, int numPlaces)
{
  //rsAssert(numPlaces >= 0 && numPlaces <= length, "we require 0 <= numPlaces <= length ");
  for(int i = length-1; i >= numPlaces; i--)
    buffer[i] = buffer[i-numPlaces];
  fillWithZeros(buffer, numPlaces);
}

template <class T>
void rsArrayTools::sequenceSqrt(const T *y, int yLength, T *x)
{
  int m2 = firstIndexWithNonZeroValue(y, yLength);
  if(m2 == -1)
  {
    // y is all zeros - return an all-zero x-sequence:
    fillWithZeros(x, (yLength+1)/2);
    return;
  }
  int m = m2/2;          // 1st index of nonzero value in x
  fillWithZeros(x, m);
  x[m] = rsSqrt(y[m2]);    // maybe use a complex generalization later
  T scaler = T(1) / (2*x[m]);
  for(int n = m+1; n < (yLength+1)/2; n++)
  {
    x[n] = y[m+n];
    for(int k = 1; k <= n-m-1; k++)
      x[n] -= x[m+k] * x[n-k];
    x[n] *= scaler;
  }
}

template <class T>
void rsArrayTools::shift(T *buffer, int length, int numPlaces)
{
  if(numPlaces >= 0)
    rightShift(buffer, length, numPlaces);
  else
    leftShift(buffer, length, -numPlaces);
}

template <class T>
T rsArrayTools::shiftToMakeMinimumZero(const T* x, int N, T* y)
{
  T minVal = minValue(x, N);
  add(x, -minVal, y, N);
  return minVal;
}

template <class T>
T rsArrayTools::sum(const T* buffer, int length)
{
  T accu = T(0); // constructor call with 0 should initilizes to additive neutral element
  for(int n = 0; n < length; n++)
    accu += buffer[n];
  return accu;
}

template<class T>
T rsArrayTools::sumOfProducts(const T *x, const T *y, int N)
{
  T s = T(0);
  for(int n = 0; n < N; n++)
    s += x[n]*y[n];
  return s;
}

template<class T>
T rsArrayTools::sumOfSquares(const T *x, int N)
{
  return sumOfProducts(x, x, N);
}

template<class T>
T rsArrayTools::sumOfAbsoluteDifferences(const T* x, const T* y, const int N)
{
  T s(0);
  for(int n = 0; n < N; n++)
    s += rsAbs(x[n] - y[n]);
  return s;
}

inline void rsArrayTools::swapDataBuffers(void *buffer1, void *buffer2, void *bufferTmp, int sizeInBytes)
{
  memcpy(bufferTmp, buffer1, sizeInBytes);
  memcpy(buffer1, buffer2, sizeInBytes);
  memcpy(buffer2, bufferTmp, sizeInBytes);
}

template<class T>
void rsArrayTools::transformRange(const T* x, T* y, int N, T targetMin, T targetMax)
{
  T currentMin = minValue(x, N);
  T currentMax = maxValue(x, N);
  T a = (targetMin - targetMax) / (currentMin - currentMax);
  T b = (currentMax*targetMin - currentMin*targetMax) / (currentMax - currentMin);
  affineTrafo(x, y, N, a, b);
}

template<class T>
void rsArrayTools::transposeSquareArray(T **in, T **out, int size)
{
  for(int i = 0; i < size; i++)
  {
    for(int j = 0; j < size; j++)
      out[i][j] = in[j][i];
  }
}

template<class T>
void rsArrayTools::transposeSquareArray(T **A, int N)
{
  int k = 1;
  for(int i = k-1; i < N-1; i++) {
    for(int j = k; j < N; j++)
      rsSwap(A[i][j], A[j][i]);
    k++; }
}

template<class T>
void rsArrayTools::unwrap(T* a, int N, T p)
{
  int k = 0; // k-factor in newValue = oldValue + k * period
  for(int n = 1; n < N; n++) {
    k = rsUnwrapFactor(a[n], a[n-1], p, k);
    a[n] += k*p; }
}

template <class T>
T rsArrayTools::weightedSum(const T *w, const T *x, rsUint32 length) // use int
{
  T s(0);
  for(rsUint32 i = 0; i < length; i++)
    s += w[i] * x[i];
  return s;
}

/*
template <class T>
void rsArrayTools::weightedSum(const T *buffer1, const T *buffer2, T *result, int length, T weight1, T weight2)
{
  for(int n = 0; n < length; n++)
    result[n] = weight1 * buffer1[n] + weight2 * buffer2[n];
}
*/


/*

todo:
isFreeOf - like contains but inverse
containsOnce - finds element like contains and then uses isFreeOf on rest of array
isPermutationOf - uses containsOnce on the first array for each element of a second array

maybe make a class rsVectorTools that just contains convenience functions for the functions from
rsArrayTools, such that we don't need the ugly &b[0] syntax and maybe can get rid of the length
parameters (because vectors know their lengths)...but maybe it should also include


maybe for more ideas what could be useful, see:


https://www.youtube.com/watch?v=h4Jl1fk3MkQ

https://www.youtube.com/watch?v=2olsGf6JIkU
https://www.youtube.com/watch?v=pUEnO6SvAMo


*/
