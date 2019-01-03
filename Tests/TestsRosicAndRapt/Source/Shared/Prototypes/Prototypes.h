#ifndef RS_PROTOTYPES_H
#define RS_PROTOTYPES_H

//#include "rapt/rapt.h"
#include "rosic/rosic.h"

// new implementation of classic IIR filter design:
#include "FilterDesign/PoleZeroPrototype.h"
#include "FilterDesign/PoleZeroMapper.h"
#include "FilterDesign/PoleZeroDesignerAnalog.h"
#include "FilterDesign/PoleZeroDesignerDigital.h"
#include "FilterDesign/ComplementaryFilters.h"

#include "Probability.h"
#include "Projection3Dto2D.h"
#include "Polygon.h"
#include "Drawing.h"

#include "SinusoidalModeling.h"

/** This file contains prototypical implementations of algorithms. These prototypes are not meant 
to be used for production code but are useful for a more readable proof-of-concept (because of lack 
of optimizations), for tweaking an algorithm's internal parameters which might not be even exposed 
in the production-code versions, and to create reference output for the unit-tests for production 
code. */

/** Prototype for rsResampler::signalValueViaSincAt(). It provides as additional parameters for 
tweaking: 
-pointer to a window-function
-parameter for the window (if applicable)
-switch for normalizing the output by the sum of the tap weights 
*/
double signalValueViaSincAt(double *x, int N, double t, double sincLength, double stretch,
  //FunctionPointer3DoublesToDouble windowFunction = rsExactBlackmanWindow, 
  double (*windowFunction)(double,double,double) = &RAPT::rsWindowFunction::exactBlackman,
  double windowParameter = 0.0, bool normalizeDC = true);

/** Generates polynomial coefficients of the polynomial used in Halpern filters. It's the T^2(w) 
polynomial in Eq. 8.18 in Paarmann: Design and Analysis of Analog Filters. */
void halpernT2(double *c, int N);

/** Generates polynomial coefficients of the polynomial used in Papoulis filters. It's the L^2(w) 
polynomial in Eq. 8.14 in Paarmann: Design and Analysis of Analog Filters */
void papoulisL2(double *c, int N);

//=================================================================================================

template<class TSig, class TPar>
class rsStateVectorFilter
{
  typedef const TSig& CRSig;
  typedef const TPar& CRPar;

public:

  /** Sets up the filter coefficients to simulate a biquad filter with given coeffs. */
  void setupFromBiquad(CRPar b0, CRPar b1, CRPar b2, CRPar a1, CRPar a2);

  /** Sets up the two poles of this filter. You need to pass real and imaginary parts of both 
  poles separately. If there are two real poles, the imaginary parts p1im, p2im should both be zero 
  and if there's a complex pair, the imaginary parts should be negatives of each other, i.e p2im 
  should be -p1im. The poles determine the coefficients in the state update matrix. */
  void setPoles(CRPar p1re, CRPar p1im, CRPar p2re, CRPar p2im);

  /** Assuming the poles are already fixed, this function computes the mixing coefficients such 
  that the first 3 samples of the impulse response will equal what you pass to this function. This 
  is used to compute the mixing coefficients after the poles have been determined. */
  void setImpulseResponseStart(TPar h[3]);

  // maybe make a setZeros function, too


  /** Produces one output sample at a time. */
  inline TSig getSample(CRSig in)
  {
    updateState(in);
    return cx*x + cy*y + ci*in;
  }

  /** Resets the filter state. */
  void reset()
  {
    x = y = TSig(0);
  }

protected:

  /** Used internally in getSample to update the filter state. */
  inline void updateState(CRSig in)
  {
    TSig t = x;             // temporary
    x = xx*x + xy*y + in;   // update x
    y = yx*t + yy*y + in;   // update y
  }

  /** This is a function to fudge with the poles in cases where they are (almost) equal. Such a 
  case cannot be represented exactly by this filter structure (a singular matrix in the mixing 
  coefficient calculation would occur), so we use distinct poles close to the originally desired 
  poles. The effect is a slight misadjustment of the filter in these particular cases. */
  void makePolesDistinct();
   // maybe return a bool to inform, if the poles were modified, maybe also return a bool from
   // setPoles in order to be able to make client code aware of the fudging

  TPar xx = 0, xy = 0, yx = 0, yy = 0;  // matrix coeffs
  TPar cx = 0, cy = 0, ci = 1;          // mixing coeffs
  TSig x  = 0, y  = 0;                  // state vector

};

//=================================================================================================

/** A filter representing a mode that uses four parallel decaying sine filters with different decay
rates but otherwise the same parameters. The different decay rates mix to create interesting 
envelope shapes. It uses a single rsFloat32x4 variable for the signal, i.e. it computes four 
parallel filters at once using SSE2 vector instructions. */

class rsModalFilterFloatSSE2
{

public:

  /** Sets up the mode parameters. Omega is the radian frequency (2*pi*f/fs), the phase is in 
  radians and the decay time constants are in samples. See also rsDampedSineFilter. */
  void setParameters(
    double omega,   double energy,  double phase, 
    double attack1, double attack2, double attackBlend,
    double decay1,  double decay2,  double decayBlend);

  /** Produces the vector of the 4 outputs of the 4 individual decaying sine filters. The actual 
  scalar output sample would be the sum of these 4. */
  inline rsFloat32x4 getSampleVector(rsFloat32x4 in)
  {
    rsFloat32x4 y = b0*in + b1*x1 - a1*y1 - a2*y2; // todo: use all plusses (more efficient)
    x1 = in;  // maybe multiply by b0 at the output instead of input for better reponse to amplitude
    y2 = y1;  // modulation, try (transposed) direct form 2
    y1 = y;
    return y;
  }


  /** Produces a scalar output sample that adds up all the 4 decaying sines. Whne a single mode is
  synthesiszed, you can use this function. When many modes are added, it makes more sense to just
  call the vector function and accumulate the vectors and do just a single sum after the 
  accumulation. */
  //inline float getSample(float in) { return getSampleVector(rsFloat32x4(in)).getSum(); }



  // some test functions for performance measurements (not inlined, so we can look at the generated 
  // assembly code):
  rsFloat32x4 getSampleVectorTestDF1( rsFloat32x4 in);
  rsFloat32x4 getSampleVectorTestDF2( rsFloat32x4 in); 
  rsFloat32x4 getSampleVectorTestTDF2(rsFloat32x4 in);
  inline float getSample(float in) 
  { 
    static const rsFloat32x4 x = 1.f;
    rsFloat32x4 y = getSampleVectorTestDF1(x);
    //rsFloat32x4 y = getSampleVectorTestDF2(x);
    //rsFloat32x4 y = getSampleVectorTestTDF2(x);
    return 1.f; 
  }
  
  /** Resets the state variables to all zeros. */
  void reset() { x1 = y1 = y2 = 0; }

protected:

  rsFloat32x4 x1 = 0, y1 = 0, y2 = 0, b0 = 0, b1 = 0, a1 = 0, a2 = 0;
  //rsFloat32x4 b0 = 0, b1 = 0, a1 = 0, a2 = 0, x1 = 0, y1 = 0, y2 = 0;
  //rsFloat32x4 tmp;

  rsFloat32x4 x2 = 0, b2 = 0; // test - make a full biquad


};

class rsModalBank
{

protected:


  static const int maxNumModes = 1024;
  rsModalFilterFloatSSE2 modeFilters[maxNumModes];

  struct ModeParameters
  {
    float frequency, amplitude, phase, attack1, attack2, attackBlend, decay1, decay2, decayBlend;
  };
  ModeParameters modeParams[maxNumModes];


  double sampleRate;
  int numModes;


  //std::vector<rsModalFilterFloatSSE2> modeFilters;

};

//=================================================================================================

/** Another go at a general ordinary differential equation solver with a probably more convenient
interface than the old one (which required subclassing to define a concrete ODE (system)).

References
(1) Numerical Recipies
(2) Mathematik, Ahrens et al, 4.Aufl.

-maybe rename to rsExplicitInitialValueSolver (and provide also an implicit one)
-maybe factor out a solver that doesn't carry around x and where f only depends on y - or maybe
 subsume systems that depend explicitly on x by incorporating an identity function into the vector
 of functions f(y), i.e. f(y1, y2, y3, ...) = (y1, f2(y1,y2,y3..), f3(y1,y2,y3..), ...), the 
 derivative of y1 is always 1, so x += h translates to y1 += 1*h in the new notation. that would 
 simplify interface and implementation but requires more understanding from the user and does not 
 allow to have a different datatype for x
-maybe move the state variables to a subclass (rsMultiStepInitialValueSolver or something)
*/

template<class Tx, class Ty>
class rsInitialValueSolver 
{

public:


  /** \name Setup */

  void setStepSize(Tx newSize) { h = newSize; }




  /** \name Evaluation */

  Ty getSampleForwardEuler()
  {
    f0 = f(x, y);
    x += h;
    y += h*f0;
    return y;
  }

  Ty getSampleRungeKutta4()
  {
    Ty k1 = f(x,         y);
    Ty k2 = f(x + 0.5*h, y + 0.5*h*k1);
    Ty k3 = f(x + 0.5*h, y + 0.5*h*k2);
    Ty k4 = f(x +     h, y +     h*k3);
    x += h;
    y += (h/6.) * (k1 + 2*k2 + 2*k3 + k4); // optimize away division
    return y;
  }


  Ty getSampleAdamsBashforth2()
  {
    f0 = f(x, y);  // new evaluation
    x += h;
    y += 0.5*h*(3*f0 - f1);
    f1 = f0;
    return y;
  }
  // todo: make sure f0, f1, are initialized correctly (do this by doing an RK step)
  // see (2) page 481


  Ty getSampleAdamsBashforth4()
  {
    f0 = f(x, y);  // new evaluation
    x += h;
    y += (h/24.) * (55*f0 - 59*f1 + 37*f2 - 9*f3); // optimize away division
    f3 = f2;
    f2 = f1;
    f1 = f0;
    return y;
  }

  // rsCubicSpline<Tx, Ty> getSolution(Tx x0, Tx x1, double accuracy);
  // should produce an object of class rsCubicSpline that represents the solution



protected:

  Tx x = 0;
  Ty y = 0;
  Tx h = 1;  // step size

  std::function<Ty(Tx,Ty)> f; // this is the "f" in y'(x) = f(x,y)

  // state variables for multistep methods:
  Ty f0, f1, f2, f3, f4; // f(x[n],y[n]), f(x[n-1],y[n-1]), ...

};


//=================================================================================================
// stuff for sliding maximum filter:

template<class T>
class rsBuffer
{

public:

  /** The actual capacity will be the next power of two of given value. */
  rsBuffer(size_t capacity);

  /** Initializes all the buffer elements with given value (default is zero). */
  void initBufferValues(T value = T(0));

protected:

  /** Wraparound */
  inline size_t wrap(size_t i) const { return i & mask; } 

  std::vector<T> data;
  size_t mask;

};
// maybe this class should be named rsRingBuffer or rsCircularBuffer or rsRingBufferBase and the 
// class below should be named rsDelayLine

//-------------------------------------------------------------------------------------------------

template<class T>
class rsRingBuffer : public rsBuffer<T>
{

public:


  rsRingBuffer(size_t capacity) : rsBuffer(capacity) {}


  /** Sets up a new buffer length and adjusts the left index accodingly, leaving the right index
  where it is. The reason to do it that way and not the otherway around is that in a delayline,
  the right index represents the write-head and the left index represents the read-head and on a
  change of  the delay time, we want to move the read-head. */
  void setLength(size_t newLength) 
  {
    length = newLength;
    updateLeftIndex();
  }

  inline size_t getLength() const { return length; }



  inline T getSample(T in)
  {
    rightIndex = wrap(rightIndex+1); // incremenet before write/read to allow further operations
    leftIndex  = wrap(leftIndex+1);  // after getSample (indices must be still valid)
    data[rightIndex] = in;           // right index is write index
    T out = data[leftIndex];         // left index is read index
    return out;
  }
  // rename to pushRightPopLeft maybe make a subclass delayline that defines an alias function
  // name getSample...any maybe merge class with rsDoubleEndedQueue



  /** Returns the maximum value in the range between the two pointers. */
  T getMaximum()
  {
    size_t i = rightIndex;
    T ex = data[i];
    while(i != leftIndex) {
      i = wrap(i-1);
      if(data[i] > ex)  // or >= or maybe use a comparison function and rename function to
        ex = data[i];   // getOptimum/Extremum
    }
    return ex;
  }

  void reset();

protected:

  /** Updates the current left index according right index and length. */
  void updateLeftIndex()  { leftIndex = wrap(rightIndex - length); } 

  /** Updates the current right index according left index and length. */
  //void updateRightIndex() { rightIndex = wrap(leftIndex + length); } 
  // maybe uncomment if needed - it's not typically used

  size_t rightIndex = 0, leftIndex = 0; // rename back to writeIndex, readIndex
  size_t length = 0; // rename to capacity ..or use data.size() - nope it's the current length
  // maybe name them generally rightEnd, leftEnd or something ... rgt, lft. L,R - then we may
  // also make rsDoubleEndedQueue a subclass and inherit the data and wrap function
};

//-------------------------------------------------------------------------------------------------


template<class T>
bool rsGreater(const T& a, const T& b);

template<class T>
bool rsLess(const T& a, const T& b);  // merge with function in RAPT...SortAndSearch

/** Implements the double-ended queue data structure. This implementation provides a static, finite 
capacity that has to be passed at construction time. Due to implementation details, the actual 
capacity is always a power of two minus two, so when you pass an arbitrary number to the 
constructor that is not of the form 2^k - 2 the smallest possible k will be chosen such that 
2^k - 2 is greater than or equal to your requested capacity. 

The implementation is based on a circular buffer. Wikipedia says, the C++ std::deque class uses a 
multiple-array implementation which is probably not so good for realtime purposes:
https://en.wikipedia.org/wiki/Double-ended_queue */

template<class T>
class rsDoubleEndedQueue : public rsBuffer<T>
{

public:

  /** Constructor. Allocates enough memory for an internal buffer to hold the "capacity" number of 
  queued values. The memory allocated for the buffer will be the smallest power of two that is 
  greater or equal to the desired capacity plus 1 (two additional and unusable index value are 
  required to make the index arithemetic work right). */
  rsDoubleEndedQueue(size_t capacity) : rsBuffer<T>(capacity+2) {}


  /** \name Data Access */

  /** Appends a value to right/head side of the queue. */
  inline void pushFront(T value)
  {
    RAPT::rsAssert(!isFull(), "Trying to push onto full deque");
    data[head] = value;
    head = wrap(head+1);
  }

  /** Appends a value to left/tail side of the queue. */
  inline void pushBack(T value)
  {
    RAPT::rsAssert(!isFull(), "Trying to push onto full deque");
    data[tail] = value;
    tail = wrap(tail-1);
  }

  /** Removes a value from the right/head side of the queue and returns it. */
  inline T popFront()
  {
    RAPT::rsAssert(!isEmpty(), "Trying to pop from empty deque");
    head = wrap(head-1);
    return data[head];
  }

  /** Removes a value from the left/tail side of the queue and returns it. */
  inline T popBack()
  {
    RAPT::rsAssert(!isEmpty(), "Trying to pop from empty deque");
    tail = wrap(tail+1);
    return data[tail];
  }

  /** Returns the value at the head of the queue, i.e. the rightmost value. */
  inline T readHead() const 
  { 
    RAPT::rsAssert(!isEmpty(), "Trying to read from empty deque");
    return data[wrap(head-1)]; 
  }

  /** Returns the value at the tail of the queue, i.e. the leftmost value. */
  inline T readTail() const 
  { 
    RAPT::rsAssert(!isEmpty(), "Trying to read from empty deque");
    return data[wrap(tail+1)]; 
  }
  // or maybe readFirst/Last Front/Back or peekFront/peekBack

  /** Clears the queue and optionally resets the data in the underlying buffer to all zeros. The 
  clearing of the queue itself just resets some pointers/indices, turning the data into 
  inaccessible garbage - but it should have no effect for the functionality of the queue whether 
  this data buffer is uninitialized/garbage or all zeros. */
  void clear(bool cleanGarbage = false)
  {
    head = 1;
    tail = 0;
    if(cleanGarbage)
      initBufferValues(0);
  }


  /** \name Inquiry */

  /** Returns the number of values that are currently in the queue. */
  inline size_t getLength() const { return wrap(head - tail - 1); }

  /** Returns the maximum allowed length for the queue. Due to the way, the head- and tail pointers 
  operate, the usable length is 2 less than the size of the underlying data buffer. */
  //inline size_t getMaxLength() const { return data.size()-1; } //
  inline size_t getMaxLength() const { return data.size()-2; } //

  /** Returns true if the queue is empty. */
  inline bool isEmpty() const { return getLength() == 0; }

  /** Returns true, if the queue is full. */
  //inline bool isFull() const { return getLength() >= getMaxLength(); } 
  inline bool isFull() const { return getLength() > getMaxLength(); } 
  // shouldn't that be >= ? ...somehow it seems, it can actually take one more entry ...but only
  // under certain conditions? more tests needed....


protected:

  size_t head = 1, tail = 0; // chek out correct terminology in the literature
  // The head-pointer is always one position to the right of rightmost value and the tail pointer
  // is always one position to the left of the leftmost value (both with wraparound). This is 
  // required in order to compute the length of the queue from the had/tail indices. Placing them
  // directly ON the leftmost/rightmost values would not allow to distiguish between an empty
  // and one-element queue (in both cases, we would have head == tail). So the actual queued values 
  // are at indices i with: tail < i < head, not: tail <= i <= head and when head = tail + 1, the
  // queue is empty because no index i fits in between tail and head
};

/*
ToDo: 
-allow to dynamically resize the capacity at runtime - but such resize operations should be an
absolute exception in realtime code (only a last resort) 
-maybe provide push-functions that take a vector argument and pop-functions that pop a range of
values at once
-maybe provide some sort of random access functions:
T fromFront(size_t i) where i = 0 returns the front element, 1 the one after the front, etc.
or fromHead - analogously for fromBack/fromTail/fromRight
-maybe binary index search functions (at which index, counted from front (or back) does the value
in the que get greater than some value
size_t searchFromFront...this may become relevant for optimizing the moving maximum filter
*/

// maybe make convenience subclasses rsQueue and rsStack - both datastructures have different 
// subsets of functionality of the double ended queue actually, but the convenience functions could
// have more stacky or queuey names like push/pop/top, enqueue/dequeue/front/back

//-------------------------------------------------------------------------------------------------

/** Implements a realtime moving maximum filter, that is: a filter that looks at audio samples
x[n], x[n-1],..., x[n-k] and extracts the maximum value in this range of samples. The algorithm is 
based on a double ended queue and has an amortized complexity of O(1) per sample (a naive 
implementation would have complexity O(k) per sample). For an explanation of the algorithm, 
see: https://www.nayuki.io/page/sliding-window-minimum-maximum-algorithm  */

template<class T>
class rsMovingMaximumFilter
{

public:

  rsMovingMaximumFilter(size_t maxLength);

  /** Sets up the length of the filter, i.e. the number of samples within which a maximum is 
  searched. */
  void setLength(size_t newLength) { delayLine.setLength(newLength); }
  // todo: we may have to update the content of the maxDeque - when this is called during running 
  // the filter and new length is shorter than the old, some values that are currently in the deque
  // will never be removed - so we should remove them in the moment, when the length changes

  /** Sets the "greater-than" comparison function. Note that you can actually also pass a function
  that implements a less-than comparison in which case the whole filter turns into a moving-minimum
  filter. */
  void setGreaterThanFunction(bool (*greaterThan)(const T&, const T&))
  {
    greater = greaterThan;
  }
  //void setComparisonFunction(const std::function<bool(const T&, const T&)>& greaterThan) 
  //{ greater = greaterThan; }


  /** Returns up the length of the filter, i.e. the number of samples within which a maximum is 
  searched. */
  size_t getLength() const { return delayLine.getLength(); }



  /** \name Processing */

  /** Computes and returns an output sample.  */
  inline T getSample(T in)
  {
    // accept new incoming sample - this corresponds to 
    // Nayuki's Step 2 - "increment the array range’s right endpoint"
    while(!maxDeque.isEmpty() && greater(in, maxDeque.readTail()) )  
      maxDeque.popBack();
    maxDeque.pushBack(in);
    // maybe we could further reduce the worst case processing cost by using binary search to 
    // adjust the new tail-pointer instead of linearly popping elements one by one? search
    // between tail and head for the first element that is >= in (or not < in) - but maybe that 
    // would destroy the amortized O(1) cost? ...but actually, i don't think so
    // maybe using a "greater" function instead of a > operator is not such a good idea, 
    // performance-wise - we should do performance tests and maybe for production code, provide
    // a version of getSample that just uses >

    // update delayline (forget oldest sample and remember current sample for later), this 
    // corresponds to Nayuki's Step 3 - "increment the array range’s left endpoint":
    T oldest = delayLine.getSample(in);
    if(!maxDeque.isEmpty()) {            // happens when length is set to zero
      T maxVal = maxDeque.readHead();
      if(maxVal == oldest)
        maxDeque.popFront();          
      return maxVal;
    }
    else
      return in;
  }

  /** Computes an output sample using a naive algorithm that scans the whole ringbuffer for its 
  maximum value. The algorithm has a per sample complexity O(k) where k is the length of the 
  filter. The implementation is mainly for testing purposes and should probably not be used in 
  production code. */
  inline T getSampleNaive(T in)
  {
    delayLine.getSample(in);  // output of getSample not needed here
    T maxVal = delayLine.getMaximum(); // maybe pass greater function to the maximum finder
    return maxVal;
  }

  /** Resets the filter to its initial state. */
  void reset()
  {
    delayLine.reset();
    maxDeque.clear();
  }

protected:

  rsRingBuffer<T> delayLine;
  rsDoubleEndedQueue<T> maxDeque;

  bool (*greater)(const T&, const T&) = &rsGreater;
  //std::function<bool(const T&, const T&)> greater; // = &rsGreater;

};

// maybe call it rsMovingSelector, movingExtremum ...and even more general concept would be a 
// moving aggregator which could also be something like a moving median (or r-th quantile)
//
// https://www.nayuki.io/page/sliding-window-minimum-maximum-algorithm

/** A variant of the rsMovingMaximumFilter that also extracts the minimum. It's slightly more 
economic to extract min and max with a single filter object rather than using separate filters for
min and max (specifically, the delayline can be shared by both filters). */

template<class T>
class rsMovingMinMaxFilter : public rsMovingMaximumFilter<T>
{

public:

  rsMovingMinMaxFilter(size_t maxLength) : rsMovingMaximumFilter(maxLength), minDeque(maxLength) {}

  void setLessThanFunction(bool (*lessThan)(const T&, const T&)) { less = lessThan; }

  inline void getMinMax(T in, T* minVal, T* maxVal)
  {
    // update deque tails:
    while(!maxDeque.isEmpty() && greater(in, maxDeque.readTail()) )
      maxDeque.popBack();
    maxDeque.pushBack(in);
    while(!minDeque.isEmpty() && less(   in, minDeque.readTail()) )
      minDeque.popBack();
    minDeque.pushBack(in);

    // update delayline and init outputs:
    T oldest = delayLine.getSample(in);
    *maxVal = *minVal = in;

    // update deque heads:
    if(!maxDeque.isEmpty()) {
      *maxVal = maxDeque.readHead();
      if(*maxVal == oldest)
        maxDeque.popFront(); }
    if(!minDeque.isEmpty()) {
      *minVal = minDeque.readHead();
      if(*minVal == oldest)
        minDeque.popFront(); }
  }
  // hmm...maybe that shouldn't be inlined...we'll see

  void reset()
  {
    rsMovingMaximumFilter::reset();
    minDeque.clear();
  }

protected:

  rsDoubleEndedQueue<T> minDeque;
  bool (*less)(const T&, const T&) = &rsLess;

};

//-------------------------------------------------------------------------------------------------

template<class T>
class rsSlewRateLimiterLinear
{

public:


  /** \name Setup */

  /** Sets the maximum upward change from one sample to the next. */
  void setUpwardLimit(const T& newLimit)
  {
    upwardLimit = newLimit;
  }

  /** Sets the maximum downward change from one sample to the next. */
  void setDownwardLimit(const T& newLimit)
  {
    downwardLimit = newLimit;
  }

  /** Convenience function to set both, upward and downward limits, to the same value. */
  void setLimits(const T& newLimit)
  {
    setUpwardLimit(newLimit);
    setDownwardLimit(newLimit);
  }


  /** \name Processing */

  T getSample(T in)
  {
    y1 += rsClip(in-y1, -downwardLimit, upwardLimit);
    return y1;
  }

  void reset() { y1 = T(0); }

protected:

  T upwardLimit = RS_INF(T), downwardLimit = RS_INF(T);  
  T y1;  // previous output sample

};

//-------------------------------------------------------------------------------------------------

/** A class for applying a special kind of nonlinear smoothing algorithm. It is based on min/max
filtering over a certain number of samples and taking an adjustable linear combination of both and 
then applying a linear slew rate limiter to the result where the limit on the slew rate is 
adaptively updated according to the current values of min and max such that that it could gor from 
min to max in the same given number of samples (or some fraction or multiple thereof).

It is good for post-processing the raw output of an envelope follower to make it smoother. The raw
envelope follower output will typically show parasitic oscillations with the signal's frequency. 
Setting up a smopther with a length equal to the cycle-length of the (enveloped) input wave will
give optimal smoothing for the signal - in an ideal situation, the oscillations in the detected 
envelope will be completely flattened out. */

template<class T>
class rsMinMaxSmoother
{

public:

  rsMinMaxSmoother(size_t maxLength) : minMaxFilter(maxLength) {}


  /** \name Setup */

  /** Sets the smoothing length in samples. */
  void setLength(size_t newLength) 
  {
    minMaxFilter.setLength(newLength);
    updateSlewRateLimitingFactor();
  }

  /** Adjusts the scaling factor for the slew-rate. With a factor of one, the slew rate is set up
  such that a transition between the currently measured min and max can occur in between L samples
  where L is the length of the min-max filter. An amount of zero turn slew rate limitig off and 
  amounts higher than one make the whole thing sluggish and are probably not useful. */
  void setSlewRateLimiting(T newAmount)
  {
    slewLimitingAmount = newAmount;
    updateSlewRateLimitingFactor();
  }

  /** Sets the mix between min and max which is used as output signal (before it goes into the 
  adaptive, linear slew rate limiter) where 0 corresponds to min only, 1 corresponds to max only 
  and 0.5 to the arithmetic mean between min and max (todo: maybe a generalized mean could be 
  useful?). */
  void setMinMaxMix(T newMix)
  {
    mix = newMix;
  }


  /** \name Processing */

  /** Processes one sample at a time. */
  T getSample(T in)
  {
    T minVal, maxVal;
    minMaxFilter.getMinMax(in, &minVal, &maxVal);
    if(limitingFactor == RS_INF(T))
      slewLimiter.setLimits(limitingFactor); // avoids 0 * inf = NaN in case maxVal == minVal
    else
      slewLimiter.setLimits((maxVal-minVal) * limitingFactor);
    return slewLimiter.getSample( (1-mix)*minVal + mix*maxVal );
  }

  /** Resets the filter to its initial state. */
  void reset()
  {
    minMaxFilter.reset();
    slewLimiter.reset();
  }

protected:

  void updateSlewRateLimitingFactor()
  {
    limitingFactor = 1.0 / (slewLimitingAmount * minMaxFilter.getLength());
  }

  T limitingFactor = 0.0;
  T slewLimitingAmount = 1.0;
  T mix = 0.5; // mixing between minimum an maximum output, 0: min only, 1: max only
  rsMovingMinMaxFilter<T> minMaxFilter;
  ::rsSlewRateLimiterLinear<T> slewLimiter;

};


// the stuff below is just for playing around - maybe move code elsewhere:
//=================================================================================================

/** A class for representing a particular kind of string with which we can do some computations 
just like with numbers. The set of all such strings forms a group (see group theory). The group 
operation (which we call addition here) is to concatenate two strings and then delete all pairs of
equal characters, i.e. the string "aaab" would be reduced to "ab", one "aa" pair is deleted. The 
inverse element to each string is obtained by reversing it. Adding "cba" to "abc" like abc+cba 
results in abccba wich would subsequently be reduced to the empty string (the rule of deleting 
equal pairs is used as often as applicable). The additive neutral element is the empty string. */

class rsGroupString
{

public:

  rsGroupString() {}

  rsGroupString(const std::vector<unsigned int>& initialString) { s = initialString; }

  // define operator =, 
  // maybe: < (lexicographical order), * (i still have to invent a suitable multiplication rule)

  bool operator==(const rsGroupString& t) const { return t.s == s; }


  /** Adds two GroupString objects. This addition is the group operation and is (conceptually) 
  performed by concatenating two strings and then deleting all doublets (iteratively, as often as 
  necessary to eliminate all of them). */
  rsGroupString operator+(const rsGroupString &rhs) const;

  /** Unary minus. Returns the additive inverse */
  rsGroupString operator-() const { return inverse(); }

  /** Binary subtraction by adding the additive inverse. */
  rsGroupString operator-(const rsGroupString &rhs) { return *this + (-rhs); }






  /** Returns the (additive) inverse which is just the string in reversed order. */
  rsGroupString inverse() const;
   // maybe later (when we have multiplication), rename to additiveInverse and implement a 
   // multiplicativeInverse, too
   // then the class should be renamed to fieldString

   // maybe let integers 0 and 1 be used and implement 1/s = s.multiplicativeInverse, etc.


  std::vector<unsigned int> get() const { return s; }

  // bool hasDoublets();
  // void removeDoublets();


protected:



  std::vector<unsigned int> s;  // we represent the characters as unsigned integers

  //int modulus = 26;    
  // the modulus, we use the 26 lowercase letters, but that is tweakable...but we don't need that 
  // yet but maybe later when we do operations on individual characters

};

/** Subclass of rsGroupString that lets use more conveniently work with strings over the alphabet
a,b,c,..,x,y,z. The class provides the conversions from/to std::string comparison operators etc. 
But all these convenenience functions have nothing to do with the actual algebraic structure, which
is why they have been factored out to keep the baseclass pure. */

class rsGroupString2 : public rsGroupString
{

public:

  /** Convenience function meant to be used for strings over the aplhabet a,b,c,...,x,y,z. We 
  represent 'a' as 0 and then count up to 'z' = 25. */
  rsGroupString2(const char* initialString);

  //rsGroupString2(const std::string& initialString);

  /** Constructor for conversion from baseclass object (?) */
  rsGroupString2(const rsGroupString& gs);





  /** Converts to a std::string. */
  std::string toString() const;

  /** Comparison operator (for some reason, we need to override it or the compiler balks). */
  bool operator==(const rsGroupString2& t) const  { return t.s == s; }

  //bool operator==(const std::string& str) const { return str == toString(); }  // obsolete?

  //rsGroupString operator+(const rsGroupString &rhs) const;

  //rsGroupString operator+(const std::string& &rhs) const
  //{

  //}



  /** Checks, if the passed unsigned integer corresponds to one of the allowed characters, i.e. is
  from the alphabet. */
  static bool isLowerCaseLetter(unsigned int c) { return c >= 97 && c <= 122; } // 'a' = 97, 'z' = 122

  /** Checks, if all characters are inthe valid range, i.e. inside our restricted alphabet. */
  bool checkCharacters() const;

};


#endif
