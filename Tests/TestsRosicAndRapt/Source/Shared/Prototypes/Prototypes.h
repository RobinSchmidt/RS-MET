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
class rsRingBuffer
{

public:


  rsRingBuffer(size_t capacity);

  void setLength(size_t newLength);

  //inline int getLength() const { return wrap(writeIndex - readIndex); }

  inline T getSample(T in)
  {
    data[writeIndex] = in;
    T out = data[readIndex];
    writeIndex = wrap(writeIndex+1);
    readIndex  = wrap(readIndex+1);
    return out;
  }

  /** Returns the maximum value in the range between the two pointers. */
#undef min  // sucks to have that here! fuck winmindef.h!
  T getMaximum()
  {
    size_t i = writeIndex-1;
    T ex = RS_MIN(T);  // opposite extremum as initial value
    while(i != readIndex-1) {
      if(data[i] > ex)  // or >= or maybe use a comparison function and rename function to
        ex = data[i];   // getOptimum/Extremum
      i--;
    }
    return ex;
  }


  void reset();

protected:

  inline size_t wrap(size_t i) const { return i & mask; } 

  std::vector<T> data;
  size_t mask;
  size_t writeIndex = 0, readIndex = 0; 
  // maybe name them generally rightEnd, leftEnd or something ... rgt, lft. L,R - then we may
  // also make rsDoubleEndedQueue a subclass and inherit the data and wrap function

};



template<class T>
class rsDoubleEndedQueue
{

public:

  /** ...actual capacity will be a power of two */
  rsDoubleEndedQueue(size_t minCapacity);

  //void setMaximumLength(int newLength);


  inline void pushFront(T value)
  {
    head = wrap(head+1);
    data[head] = value;
  }

  inline void pushBack(T  value)
  {
    tail = wrap(tail-1);
    data[tail] = value;
  }

  inline T popFront()
  {
    T value = data[head];
    head = wrap(head-1);
    return value;
  }

  inline T popBack()
  {
    T value = data[tail];
    tail = wrap(tail+1);
    return value;
  }

  // this sort of generalizes a delayline - a regular integer delayline/ringbuffer would fill the 
  // queue initially with L zeros and then for each incoming sample do a pushFront/popBack so the 
  // ringbuffer has a fixed length L and can only write at the front and read at the back (and this
  // read/write is combined into a single operation). this dequeue here can read and write 
  // arbitarily from/to both ends...maybe factor out a ringbuffer...or...is it really exactly the 
  // same?

  inline int getLength() const { return head - tail; }

  inline T readHead() const { return data[head]; }

  inline T readTail() const { return data[tail]; }

  // or maybe readFirst/Last Front/Back



  //void writeData(T x);

  //T readData();

protected:

  inline int wrap(int i) const 
  { 
    int N = (int) data.size();
    while(i >= N) i -= N;
    while(i <  0) i += N;
    return i;
  } 
  // does this work for unsigned int as well? can we use a mask? i think so

  std::vector<T> data;  // maybe use pointer to T later,,,but bad for debugging...hmm
  //int mask;
  int head = 0, tail = 0;


  //size_t head = 0, tail = 0;
  //size_t mask;
  // or do we need int to correctly wrap around values below 0?
};


template<class T>
class rsMovingMaximumFilter
{

public:

  rsMovingMaximumFilter(size_t maxLength);


  inline T getSample(T in)
  {
    T oldest = rngBuf.getSample(in);

    T maxVal = rngBuf.getMaximum();

    return maxVal;
  
    //return oldest; // preliminary
  }

protected:

  rsRingBuffer<T> rngBuf;
  rsDoubleEndedQueue<T> dqueue;

};

// generalize to allow for movingMinimum or movingWhatever filter - it should take a comparison
// function as parameter, i.e. a function that takes two T-values as input and returns a bool 
// where the default is 
// bool greater(a, b) { return a > b; }  // or maybe greaterOrEqual?
// maybe call it rsMovingSelector ...and even more general concept would be a moving aggregator 
// which could also be something like a moving median (or r-th quantile)
//

// https://www.nayuki.io/page/sliding-window-minimum-maximum-algorithm



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
