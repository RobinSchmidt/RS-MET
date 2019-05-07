#ifndef RS_PROTOTYPES_H
#define RS_PROTOTYPES_H

// todo: move all prototyes into rs_testing juce module - and maybe much of the other code, too

//#include "rapt/rapt.h"
#include "rosic/rosic.h"
using namespace RAPT;

// new implementation of classic IIR filter design:
#include "FilterDesign/PoleZeroPrototype.h"
#include "FilterDesign/PoleZeroMapper.h"
#include "FilterDesign/PoleZeroDesignerAnalog.h"
#include "FilterDesign/PoleZeroDesignerDigital.h"
#include "FilterDesign/ComplementaryFilters.h"
#include "FilterDesign/NonUniformFilter.h"


#include "ParticleBouncer.h"
#include "Probability.h"
#include "Projection3Dto2D.h"
#include "Polygon.h"
#include "Drawing.h"



/** This file contains prototypical implementations of algorithms. These prototypes are not meant 
to be used for production code but are useful for a more readable proof-of-concept (because of lack 
of optimizations), for tweaking an algorithm's internal parameters which might not be even exposed 
in the production-code versions, and to create reference output for the unit-tests for production 
code. */


/** Solves a pentadiagonal linear system of equations with given diagonals and right-hand side 
using a simple algorithm without pivot-search. lowerDiag1 is the one directly below the main 
diagonal, lowerDiag2 the one below lowerDiag1 - and similarly for upperDiag1/upperDiag2. In the 
process of the computations, the right hand side vector is destroyed. the same is true for mainDiag
and the two inner sub/superdiagonals lowerDiag1, upperDiag1. Note also that you can't use the same 
array for lowerDiag1 and upperDiag1, even if your matrix is symmetric.

..What about lowerDiag2/upperDiag2? are these preserved and may these point to the same vector? 
It's probably safest to assume that everything may get messed up and all arrays should be 
distinct. */
std::vector<double> solvePentaDiagonalSystem(
  std::vector<double>& lowerDiag2, std::vector<double>& lowerDiag1,
  std::vector<double>& mainDiag, 
  std::vector<double>& upperDiag1, std::vector<double>& upperDiag2,
  std::vector<double>& righHandSide);

/** Multiplies a pentadiagonal matrix with a vector...  */
std::vector<double> pentaDiagMatVecMul(
  std::vector<double>& lowerDiag2, std::vector<double>& lowerDiag1,
  std::vector<double>& mainDiag, 
  std::vector<double>& upperDiag1, std::vector<double>& upperDiag2,
  std::vector<double>& input);

/** Minimizes the sum of squared differences between adjacent array elements under the constraint 
that the sums of adjacent array elements must be equal to given values. The input array s is the 
length N-1 array of the desired sums, the output array v is the length N value array, such that
v[i] + v[i+1] = s[i] and sum_i (v[i+1] - v[i])^2 = min. You may optionally pass an array of 
weights for the squared differences in the cost function - if you do, the w array must have the 
same length as s, if you don't, unit weights will be used for each squared difference. With 
weights, we will minimize sum_i w[i] * (v[i+1] - v[i])^2 subject to the (same) constraints that
v[i] + v[i+1] = s[i] for all i = 0,..,N-2 */
std::vector<double> rsMinSqrDifFixSum(
  const std::vector<double>& s, 
  const std::vector<double>& w = std::vector<double>() );

std::vector<double> rsMinSqrCrvFixSum(
  const std::vector<double>& s, 
  const std::vector<double>& w = std::vector<double>() );



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

/** Calculates a chebyshev window of size N, store coeffs in out as in Antoniou
  -out should be array of size N 
  -atten is the required sidelobe attenuation (e.g. if you want -60dB atten, use '60') 
Dolph-Chebychev window generation code from here:
http://practicalcryptography.com/miscellaneous/machine-learning/implementing-dolph-chebyshev-window/
not recommended for production use because the complexity is O(N^2) - instead use an iFFT 
approach
References:
[1] Lyons, R., "Understanding Digital Signal Processing", Prentice Hall, 2004.
[2] Antoniou, A., "Digital Filters", McGraw-Hill, 2000.  */
void cheby_win(double *out, int N, double atten);



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
adaptively updated according to the current values of min and max such that that it could go from 
min to max in the same given number of samples (or some fraction or multiple thereof).

It is good for post-processing the raw output of an envelope follower to make it smoother. The raw
envelope follower output will typically show parasitic oscillations with the signal's frequency. 
Setting up a smoother with a length equal to the cycle-length of the (enveloped) input wave will
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
  where L is the length of the min-max filter. An amount of zero turns slew rate limitig off and 
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












//-------------------------------------------------------------------------------------------------

/** This is currently just a vague idea.... */

template<class T>
class rsMovingMedianFilter
{

public:





  T getSample(T x)
  {
    insert(Node(x));        // O(log(N))
    remove(oldestNode);     // O(log(N))
    return nodeHeap[0].val;
  }

protected:


  struct Node
  {
    Node(T value) : val(value) {}
    T val;
    Node *newer; 
    //Node *older;  // we'll see, if we need this
  };

  /** ...Will also update our newestNode pointer */
  void insert(Node node)
  {
    //newestNode->newer = ...  // update "newer" pointer of current newest node
    // ...
  }

  /** ...Will also update our oldestNode pointer. */
  void remove(Node node)
  {
    Node* tmp = oldestNode->newer;
    oldestNode = tmp;
    // more to do
  }



  std::vector<Node> nodeHeap;

  Node *oldestNode;
  Node *newestNode;   // we'll see, if we need this

};

/*
Idea:
-a naive median filter would at each sample have to sort an array of delayed values and return the
 middle value of that sorted array
-sorting is an O(N*log(N)) process, so that would be the complexity per sample, which is bad
-a better implementation would insert the new incoming sample into an already sorted array - but 
 that involves shifting a lot of data around which is still O(N) - still bad
-instead, we keep the delayed samples organized as a max-heap, i.e. a binary tree that has the 
 property that for each node, the right child has the value >= and the left child a value <= the 
 value of the node in question
-inserting a new value and removing an old value from such a heap is O(log(N)) (i think, verify)
-the root of the tree/heap is always our desired output sample representing the current median
-when a new sample comes in, it has to be inserted into the heap and the oldest sample has to be 
 removed
-to do this, the filter keeps a pointer to the oldest node...


*/
















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
