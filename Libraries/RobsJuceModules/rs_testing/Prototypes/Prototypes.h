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

#include "ModalAnalyzer.h"
#include "ParticleBouncer.h"
#include "Probability.h"
#include "Projection3Dto2D.h"
#include "Polygon.h"
#include "Drawing.h"
#include "QuantumSystems.h"
#include "Relativity.h"
#include "SineParameterEstimator.h"

/** This file contains prototypical implementations of algorithms. These prototypes are not meant 
to be used for production code but are useful for a more readable proof-of-concept (because of lack 
of optimizations), for tweaking an algorithm's internal parameters which might not be even exposed 
in the production-code versions, and to create reference output for the unit-tests for production 
code. */


// todo: move to somewhere else - could be useful in other contexts:
static constexpr int allBits = -1;                                      // all bits are 1
static constexpr int allBitsButFirst = std::numeric_limits<int>::max(); // only 1st bit is 1
static constexpr int firstBitOnly = allBits ^ allBitsButFirst;          // only 1st bit is 0

// for unsiged int types, the bit twiddling is different:
//static size_t allBits = std::numeric_limits<size_t>::max();
//static size_t firstBitOnly = allBits - (allBits >> 1);
//static size_t allBitsButFirst= allBits ^ firstBitOnly;



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

/** Retrieves damped sine filter design parameters from its coefficients. See
@see rsDampedSineFilter for meaning of parameters. The phase p is returned in the interval
0...2pi. */
void rsDampedSineFilterAnalysis(double b0, double b1, double a1, double a2, double* w, double* A,
  double* d, double* p);

void rsDampedSineFilterAnalysis2(double b0, double b1, double a1, double a2, double* w, double* A,
  double* d, double* p); // other algorithm for the same thing

void rsDampedSineFilterResidueAndPole(double b0, double b1, double a1, double a2, 
  std::complex<double>* residue, std::complex<double>* pole);


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

/** Simulates the dynamics of a rotating rigid body around its three pricipal axes of intertia. If
they are all different, when it initially rotates around the axis of the middle/intermediate moment
of inertia with some small perturbation of having a rotational component around any of the other 
two axes, the rotation axis periodically flips over. This is known as the "tennis racket effect" 
because it also occurs when throwing up a tennis racket in a particular way. It is due to the 
rotation around the intermediate axis being an unstable equilibrium of the dynamic equations that
describe the rotation. Rotation around any of the other two principal axes (those with maximum and 
minimum moment of inertia) are stable equlibria.


// see:
https://en.wikipedia.org/wiki/Tennis_racket_theorem
https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)
https://en.wikipedia.org/wiki/Moment_of_inertia
https://en.wikipedia.org/wiki/Poinsot%27s_ellipsoid
https://en.wikipedia.org/wiki/Polhode
https://www.youtube.com/watch?v=1VPfZ_XzisU
https://arxiv.org/pdf/1606.08237.pdf

*/

template<class T>
class rsTennisRacket
{


public:

  /** \name Setup */

  /** Sets the ratio of the moments of inertia. This determines the frequency of the flips....
  todo: have a function setFrequency */
  void setInertiaRatio(T newRatio)
  {
    // maybe wrap into if(newRatio != ratio):
    ratio = newRatio;
    I1 = ratio;
    // I2 = 1;  // always
    I3 = T(1) / ratio;
  }

  /** Sets the current state consisting of the angular velocities along the 3 principal axes. */
  void setState(T w1, T w2, T w3)
  {
    this->w1 = w1;
    this->w2 = w2;
    this->w3 = w3;
  }

  /** Sets the step size for the numerical integration scheme */
  void setStepSize(T newSize)
  {
    h = newSize;
  }


  /** \name Inquiry */

  T getW1() const { return w1; }
  T getW2() const { return w2; }
  T getW3() const { return w3; }



  /** \name Processing */

  void updateState(T M1, T M2, T M3)
  {
    // compute angular acceleration vector:
    T a1 = (M1 - (I3 - I2)*w2*w3) / I1;
    T a2 = (M2 - (I1 - I3)*w3*w1) / I2;
    T a3 = (M3 - (I2 - I1)*w1*w2) / I3;
    // formula from https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)

    // todo: add damping terms...


    // update angular velocities:
    w1 += h*a1;
    w2 += h*a2;
    w3 += h*a3;

    // optionally renormalize rotational energy:
    if(normalizeEnergy) {
      T E = (I1*w1*w1 + I2*w2*w2 + I3*w3*w3); // is actually twice the energy
      T s = sqrt(T(1)/E);
      w1 *= s; w2 *= s; w3 *= s;
    }
    // maybe factor this out - client code may call it after calling setState to ensure an
    // energy normalized state - maybe include a target energy
  }


  T getSample(T in)
  {
    updateState(0, in, 0); // todo: use injection vector
    return w2;             // todo: use output vector
  }

protected:

  // user parameters:
  T ratio = 2;


  // algo parameters:
  T I1 = 2, I2 = 1, I3 = T(0.5); // moments of inertia along principal axes
  T h = T(0.01);                 // step size
  bool normalizeEnergy = true;

  // state:
  T w1 = 0, w2 = 1, w3 = 0;      // angular velocities along principal axes

};




//=================================================================================================

/** A second order (2 poles, 2 zeros) filter implementation, whose internal state is represented as
a 2D vector, i.e. a point in the xy-plane. The state is updated by multiplying the current state
vector by a matrix (and adding the input value to both components). The output is computed as a 
linear combination of the state-vector's coordinates and the input. The state update matrix will 
have one of these two general forms:

  |p1 0 |     or:     r * |cos(w)  -sin(w)| 
  |0  p2|                 |sin(w)   cos(w)|

where in the first case, p1 and p2 are the two real poles and the x and y states decay 
exponentially and independently from each other when the input is switched off. In the second case,
the numbers r*cos(w), r*sin(w) are the real and imaginary parts of a pair of complex conjugate 
poles and we will see a spiraling/decaying rotation of the states when there's no input (anymore).
The filter structure can realize almost any biquad transfer function - the singular problematic 
case is when there are two equal real poles - in this case, they will be made artificially distinct 
by fudging them a bit. The effect of that fudging on the transfer function will be miniscule. The
advanatge of that filter structure is that it (presumably) responds well to modulation. */

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


//=================================================================================================

/** Baseclass for rsBinaryHeap and rsBinarySearchTree... */

template<class T>
class rsBinaryTree
{

public:

  rsBinaryTree(T* newData = nullptr, int newSize = 0, int newCapacity = 0)
  {
    setData(newData, newSize, newCapacity);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the data that should be treated as tree. The object does not take ownership of the data.
  Instead, it acts pretty much like a "view" (as in rsMatrixView) - it just operates on an existing
  data array whose lifetime is managed elsewhere. */
  void setData(T* newData, int newSize, int newCapacity)
  {
    data = newData;
    size = newSize;
    capacity = newCapacity;
  }

  void setData(std::vector<T>& newData)
  {
    setData(&newData[0], (int)newData.size(), (int)newData.size());
  }

  /** Sets the comparison function to be used. If it implements "less-than", you'll get a max-heap
  and if it implements "greater-than", you'll get a min-heap. By default, it's assigned to a 
  less-than function based on the < operator. */
  void setCompareFunction(const std::function<bool(const T&, const T&)>& newFunc) 
  { less = newFunc; }
  // question: what happens, if we use a less-or-equal or greater-or-equal function for this?

  /** Sets the function used to swap two elements. By default, it just uses rsSwap but you may want
  to do something extra whenever a swap takes place in certain circumstances, for example, to keep 
  track of when items get moved around (see rsMovingQuantileFilter for an example). */
  void setSwapFunction(const std::function<void(T&, T&)>& newFunc)
  { swap = newFunc; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the number of nodes in the tree. */
  int getSize() const { return size; }
  // maybe make a version that takes a node-index as parameter

  /** Returns the capacity, i.e. the maxumum number of nodes that can be stored in the underlying 
  data array. */
  int getCapacity() const { return capacity; }

  /** Read/write access to elements. Warning: overwriting elements may destroy the defining 
  property (heap, binary-search, etc.) of the tree, so it should be used only if you know what you 
  are doing. */
  T& operator[](int i)
  {
    rsAssert(i >= 0 && i < size, "Index out of range");
    return data[i]; 
  }

  T& operator[](size_t i)
  {
    rsAssert((int) i < size, "Index out of range");
    return data[i]; 
  }


protected:

  /** Index of parent of node i. */
  static inline int parent(int i) { return (i-1) / 2; }

  /** Index of left child of node i. */
  static inline int left(int i)   { return 2*i + 1; }

  /** Index of right child of node i. */
  static inline int right(int i)  { return 2*i + 2; }

  /** Returns true, iff index i is the index of a left child node. */
  static inline bool isLeft(int i) { return i == left(parent(i)); }
  // needs test - i think, we can do it simpler: the odd indices are left children and the even 
  // indices are right children...verify that

  /** Returns true, iff index i is the index of a right child node. */
  static inline bool isRight(int i) { return i == right(parent(i)); }
  // needs test

  // The actual data array:
  T* data = nullptr;
  int size = 0;
  int capacity = 0;

  // Comparison and swapping functions:
  std::function<bool(const T&, const T&)> 
    less = [](const T& a, const T& b)->bool { return a < b;  };
  std::function<void(T&, T&)> 
    swap = [](T& a, T& b) { rsSwap(a, b); };
  // maybe these should be references?

  // Some related classes need direct acces to our members:
  template<class U> friend class rsDoubleHeap;

};


//=================================================================================================

/** Class for representing an array A of data as binary heap with functions for establishing and 
maintaining the heap-property. To understand what that property means, we must first interpret the
flat array A as a binary tree in the following way:

  left(i)   = 2*i + 1          left child index of index i
  reight(i) = 2*i + 2          right child index of index i
  parent(i) = (i-1) / 2        parent index of index i

Given that, the heap property says that for every node with index i (except the root), it holds 
that: 

  A[parent(i)] >= A[i]

Instead of >=, we could have also used <=. In the former case, we are dealing with a max-heap in 
the latter with a min-heap.

References:
  (1) Introduction to Algorithms, 2nd Ed. (Cormen, Leiserson, Rivest, Stein)  */

template<class T>
class rsBinaryHeap : public rsBinaryTree<T>
{

public:

  using rsBinaryTree::rsBinaryTree;  // inherit constructors


  void buildHeap()
  {
    for(int i = size/2-1; i >= 0; i--)  // or should we use (size-1)/2 ?
      floatDown(i);
  }
  // Runs in O(N). From the code, it would appear as having O(N*log(N)) complexity because we call
  // an O(log(N)) function inside the loop. However, for heaps and binary search trees, the runtime
  // of floatDown depends on the argument i in such a way to give an overall O(N) behavior (see 
  // reference (1)). The memory complexity is O(1).


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Return true, iff the underyling data array currently satisfies the heap property. The 
  function is meant mostly for testing and debugging purposes and is currently implemented 
  recursively (i.e. inefficiently). */
  bool isHeap(int i = 0) const
  {
    if(i >= size)
      return true;
    bool result = true;
    int l = left(i);
    int r = right(i);
    if(l < size) result &= !less(data[i], data[l]) && isHeap(l);
    if(r < size) result &= !less(data[i], data[r]) && isHeap(r);
    return result;
  }
  // maybe convert recursion to iteration
  // maybe factor out to baseclass and rename to isTree, maybe we need to call a virtual function
  // satisfiesTreeProperty(int i) that take a (parent) node index and then iterate through half of
  // size of the tree...but that places a burden on subclasses having to implement the virtual
  // function - not good!

  //-----------------------------------------------------------------------------------------------
  /** \name Data Manipulation */

  /** Replaces the element at index i with the new given element x and rebalances the heap to 
  maintain the heap-property which amounts to floating the new element x up or down. The return 
  value is the array index, where the new element actually ended up. */
  int replace(int i, const T& x)
  {
    data[i] = x;
    return floatIntoPlace(i);
  }

  /** Inserts the element x into the heap at the correct position to maintain the heap-property and 
  returns the array-index where it was inserted. */
  int insert(const T& x)
  {
    rsAssert(size < capacity, "Capacity exceeded");
    if( size == capacity ) return -1; 
    // Trying to insert an item when the heap is already full is a bug on client code side. We return
    // -1 early here, to avoid an access violation, if client code has this bug in a release version.

    data[size] = x;
    size++;
    return floatUp(size-1);
  }

  /** Removes the element at given index i from the heap and re-orders the remaining elements to
  maintain the heap-property. */
  void remove(int i)
  {
    swap(data[i], data[size-1]);
    //rsSwap(data[i], data[size-1]);  // test
    size--;
    floatIntoPlace(i);
  }
  // what if the last element is removed? should we treat that as special case? because otherwise, 
  // in this case, floatIntoPlace would be called with an i that is invalid after size-- ... check
  // this case in a unit test!
  // could the return value from floatIntoPlace be interesting for the caller? It is the position
  // where the fomerly last heap element ended up after all the re-ordering business

  void removeFirst()
  {
    remove(0);
  }
  // needs test
  // maybe we should have a special implementation to remove the front element? i think, in this case,
  // we may simply call floatDown instead of floatIntoPlace
  // or maybe we need a function getAndRemoveFirst()...or popFirst(), extractMax, extractFirst

  T extractFirst()
  {
    T first = data[0];
    removeFirst();
    return first;
  }
  // needs test


  /** Lets the node at index i float up or down into its appropriate place and returns the 
  resulting new array index. */
  int floatIntoPlace(int i) { return floatUp(floatDown(i)); }
  // -maybe rename to rebalanceNode
  // -why is this public? can we move it into the protected section? hmm, 
  //  rsMovingQuantileFilterCore calls it directly - can this be avoided?


protected:

  /** Functions to establish or maintain the heap-property of the underlying data array. */

  /** Lets the node at index i float up the tree if its parent violates the heap property. Returns 
  the new array index. */
  int floatUp(int i)
  {
    while(i > 0) {
      int p = parent(i);
      if(less(data[p], data[i]))  {
        swap(data[i], data[p]);
        i = p; }
      else
        return i; }
    return i;
  }

  /** Assuming that the subtrees rooted at left(i) and right(i) satisfy the heap property, this 
  function makes sure that node/index i also satifies the heap property. If it doesn't already 
  satisfy it, the function lets the value at i float down the subtree rooted at i. It has a time
  complexity of O(log(N)) and memory complexity of O(1). The return value is the new array index of
  the value that was previously located at index i. */
  int floatDownRec(int i)
  {
    int l = left(i);
    int r = right(i);
    int b = i;         // b for "big"
    if(l < size && less(data[i], data[l])) b = l;
    if(r < size && less(data[b], data[r])) b = r; 
    if(b != i) { swap(data[i], data[b]); return floatDown(b); }
    return i;
  }
  // a.k.a. maxHeapify
  // That's the recursive implementation from (1) page 130. When the iterative version is ready,
  // move it to the rsBinaryHeapTest subclass - we don't need it anymore in production code, then. 
  // But it may be interesting to figure out, if the recursion actually incurs an overhead since 
  // it's tail recursion and smart compilers might be able to translate it to iteration themselves.
  // We also may want to keep it as reference for unit tests (to test, if the iterative version 
  // really does the same thing).
  // maybe move to subclass rsBinaryHeapTest in the unit test section - we don't need it in 
  // production code, when the iterative version works (which seems to be the case)


  int floatDownIt(int i)
  {
    while(i < size-1)   // check if we should use size
    {
      int l = left(i);
      int r = right(i);  // == l+1  ->  optimize (but maybe parallel is actually better than serial?)
      int b = i; 
      if(l < size && less(data[i], data[l])) b = l;
      if(r < size && less(data[b], data[r])) b = r;
      if(b != i) { 
        swap(data[i], data[b]);
        i = b;  }
      else
        return i;
        // really? hmm..without, we get a hang. i think it's safe to return here, because when 
        // data[i] >= data[l] and data[i] >= data[r], i.e. if the heap condition holds at node i, 
        // it will automatically hold for all its children (because that's what we assume as 
        // precondition), so have nothing more to do here. The heap condition can only be violated 
        // at a child node of i, if we actually *did* a swap - if we didn't have to do a swap, we 
        // are done.
    }
    return i;
  }
  // iterative (i.e. non-recursive) implementation - needs tests
  // code adapted from here:
  // https://sites.math.rutgers.edu/~ajl213/CLRS/Ch6.pdf
  // ..but it seems buggy and i had to make some changes - it seems to work now but needs more 
  // thorough testing

  int floatDown(int i)
  {
    //return floatDownRec(i);  // calls recursive version - change that!
    return floatDownIt(i);
  }

  // todo: implement functions: int insert(T*), void remove(int), extractMax,
  // increaseKey, heapMax

  template<class U> friend class rsDoubleHeap;
  template<class U> friend class rsDoubleHeap2;
  template<class U> friend class rsDoubleHeap3;
  // get rid of these - we have them to get access to the less and swap functions - maybe store 
  // references to them there - that's more elegant anyway. ...or maybe use references in 
  // rsBinaryTree and the actual function object in rsDoubleHeap/2/3

};
// -maybe make also a class rsBinarySearchTree where the left child is <= and the right child is >=
//  the parent
// -pehaps we could implement a general function: needsSwap(int parentIndex, int childIndex) - in 
//  the case of a search tree, it would first have to figure out, if the childIndex is a left or 
//  right child and order the arguments of the comparison according to that
// -in the old RSLib codebase, i did some sort of linked-tree - maybe that could be dragged in as
//  rsLinkedTree or rsDynamicTree or something like that - all the different trees could be in
//  a file Trees.h/cpp
// -i think, currently, the order of the children of a node is undefined - both children are <=
//  the parent, but maybe we could impose additional useful structure, if we also have 
//  left <= right - would that be easy to implement and/or useful?

//=================================================================================================

/** Works similar to rsBinaryHeap just that the property satified by the nodes is different. 
Here, it holds that:

   A[left(i)]  <= A[i]    and
   A[right(i)] >= A[i]

As a reminder, for a (max)heap, the property was A[left(i)] <= A[i] and A[right(i)] <= A[i], i.e.
both children must have values less-or-equal than the parent. Here, the left node must be 
less-or-equal and the right node must be greater-or-equal. */

template<class T>
class rsBinarySearchTree : public rsBinaryTree<T>
{

public:

  using rsBinaryTree::rsBinaryTree;  // inherit constructors


  bool isSearchTree(int i = 0) const
  {
    if(i >= size)
      return true;
    bool result = true;
    int l = left(i);
    int r = right(i);
    if(l < size) result &= !less(data[i], data[l]) && isSearchTree(l);
    if(r < size) result &= !less(data[r], data[i]) && isSearchTree(r);
    return result;
  }
  // actually, this test is not enough - it says yes to the "pseudotree" 
  // 50,20,80,10,60,30,90 - it satisfies the property at each node with respect to direct 
  // parents/children but not for the full subtrees - however - maybe for maintaining the
  // property, it's enough to check that? can replace create pseudotrees? ..oh - yes - that seems
  // to be the case! damn! maybe the whole idea of using the same structure as for the binary heap
  // does not work out as i hoped... :-(

  void buildSearchTree()
  {
    for(int i = size-1; i >= 0; i--)
    {
      //floatUp(i);
      //floatIntoPlace(i);
      //floatDown(i);
      fixNode(i);
    }
    // is it vital to do this in reverse order?
    // i think, we can't just let the nodes i float into place - we must call some sort of
    // fixNode function that either swaps p,l or p,r or l,r, or nothing
  }


  void fixNode(int i)
  {
    //floatIntoPlace(i);  // does not work

    int l = left(i);
    if(l >= size) return;             // node is leaf
    int r = right(i);
    if(r >= size) {                   // node has only left child - make sure, the data at the
      if(less(data[i], data[l]))      // parent i is not less than at left child l - if it is: swap
        swap(data[i], data[l]);
      return; }

    // ok - we have a node that has both children - we must figure out which of the 3 nodes i,l,r
    // is the middle element m and if m is not already i, then do a swap. we must also make sure 
    // that data[l] <= data[r]. we have 3 possible swapping scenarios: swap(i,l), swap(i,r), 
    // swap(l,r) ...plus, of course, the no-swap scenario

    //int m = i;                           // index of mid element, preliminary set to i
    //if(less(data[i], data[l])) m = l;
    //if(less(data[r], data[m])) m = r;

    // fix order of children:
    if(less(data[r], data[l]))
      swap(data[l], data[r]);

    // fix order of i,l:
    if(less(data[i], data[l]))
      swap(data[l], data[i]);

    // fix order of i,r:
    if(less(data[r], data[i]))
      swap(data[r], data[i]);

    // try to do this only with 2 comparisons and one swap

  }



  int replace(int i, const T& x)
  {
    // data[i] = x; return floatIntoPlace(i); // this is what the heap needs

    // we insert it (preliminarily) either at i or at the sibling of i:
    int p = parent(i);
    if(isLeft(i))
      if(less(data[p], x)) 
        i = right(p);   // right sibling
    else
      if(less(x, data[p])) 
        i = left(p);    // left sibling
    data[i] = x; 
    return floatIntoPlace(i);
  }


protected:

  int floatIntoPlace(int i) { return floatUp(floatDown(i)); }

  int floatUp(int i)
  {
    while(i > 0) {
      int p = parent(i);
      if(isLeft(i)) {
        if(less(data[p], data[i]))
          swap(data[i], data[p]);
        else
          return i; }
      else {
        if(less(data[i], data[p]))      // i and p are swapped compared to left nodes
          swap(data[i], data[p]);
        else
          return i;   }
      i = p; }
    return i;
  }

  int floatDown1(int i)
  {
    int l = left(i);
    if(l >= size) return i;
    if(less(data[i], data[l])) {
      swap(data[i], data[l]);
      return floatDown(l); }
    int r = right(i);
    if(r >= size) return i;
    if(less(data[r], data[i])) {
      swap(data[i], data[r]);
      return floatDown(r); }
    return i;
  }
  // recursive implementation - not working


  int floatDown2(int i)
  {
    while(i < size-1) 
    {
      int l = left(i);
      if(l >= size) return i;              // node i is a leaf
      int m = i;                           // index of mid element, preliminary set to i
      if(less(data[i], data[l])) m = l;    // left child is bigger than i, so mid index goes to l
      int r = right(i);

      /*
      if(r >= size)
      {
        //rsError("not yet implemented");

        if(m != i)
        {
          swap(data[m], data[i]);
          i = m;
        }
        continue;
      }
      else
      {
        if(less(data[r], data[m]))
          m = r;
      }
      */

      if(r < size && less(data[r], data[m])) m = r;

      if(m != i) {
        swap(data[m], data[i]);
        i = m;  }
      else
        return i;

      /*
      if(r < size && less(data[b], data[r])) b = r;

      if(b != r) {              // difference to heap: there it was: if(b != i)
        swap(data[i], data[b]);
        i = b;  }
      else
        return i;
        */


      //// no - this is wrong: we need to consider 3 cases: 
      ////   swap l with i, swap r with i, swap r with l
      //// see  the heap implementation
      //if(l < size && less(data[i], data[l])) {
      //  swap(data[i], data[l]);
      //  i = l; }
      //else if(r < size && less(data[r], data[i])) {
      //  swap(data[i], data[r]);
      //  i = r; }
      //else
      //  return i; 

    
    }
    return i;

    // hmm..i think, we have to figure out the index of the biggest and smallest element b,s among
    // i,l,r. if s==l and b==r, the node i should stay where it is
  }

  int floatDown(int i)
  {
    return floatDown1(i);
    //return floatDown2(i);
  }



};

//=================================================================================================

/** Data structure that combines a max-heap with a min-heap for the purpose of splitting a bunch
of values into two groups: the small values and the large values. The two heaps satisfy the 
property that all values in the small heap are less-or-equal to all values in the large heap. The 
small values are held in a max-heap such that the largest of them can be extracted in O(1) time and 
the large values are held in a min-heap such that the smallest of them can also be extracted in 
O(1) time. This facilitates to exctract the median or other quantiles. tbc....  */

template<class T>
class rsDoubleHeap
{

public:

  rsDoubleHeap() 
  {
    large.less = [](const T& a, const T& b)->bool 
    { 
      return b < a;
    };
    // large.less is actually a greater-than function which we obtain from the less-than operator 
    // by swapping the arguments - this turns large into a min-heap rather than a max-heap. todo: 
    // use the more generic name comp for "compare" in rsBinaryTree for the comparison function
  }

  void setData(T* newSmallData, int newSmallSize, int newSmallCapacity,
               T* newLargeData, int newLargeSize, int newLargeCapacity)
  {
    small.setData(newSmallData, newSmallSize, newSmallCapacity);
    large.setData(newLargeData, newLargeSize, newLargeCapacity);
  }

  void setData(std::vector<T>& newSmall, std::vector<T>& newLarge)
  {
    setData(&newSmall[0], (int) newSmall.size(), (int) newSmall.size(),
            &newLarge[0], (int) newLarge.size(), (int) newLarge.size());
  }

  
  void setSwapFunction(const std::function<void(T&, T&)>& newFunc)
  {
    small.setSwapFunction(newFunc);
    large.setSwapFunction(newFunc);
  }


  /** Replaces the element at given key  with the given new value and return the the key where the
  new element actual ended up after doing all the floating up/down and potential swapping business.

  out of date:
  If nS is the size of the small heap, we use the convention that indices i < nS are direct indices 
  into the small heap, and when i >= nS, we use i-nS as heap index into the large heap (as seen 
  from the class rsBinaryHeap). The same holds for the return value. */
  int replace(int key, const T& newValue)  // rename index to key
  {
    rsAssert(isIndexValid(key), "Key out of range");

    int i  = key;  // maybe get rid of index and use only i
    int nS = small.getSize();

    // The actual replacement, taking place in one of our two heaps:
    if(key < small.getSize())
      i = small.replace(i,    newValue);
    else
      i = large.replace(i-nS, newValue) + nS;

    // The potential swap of the front elements (i.e. max and min) of both heaps:
    if( small.less(large[0], small[0]) )  // means: if(large[0] < small[0])
    {
      small.swap(small[0], large[0]);
      int is = small.floatDown(0);
      int il = large.floatDown(0);
      if(key < nS)     // we replaced in small heap, then swapped, so newValue is now in large heap
        i = il + nS;
      else
        i = is;
    }

    return i;
  }
  // we use the convention that when i < nS, the datum gets inserted into the small heap and when 
  // i >= nS, the datum gets inserted into the large heap at i-nS, tbc....

  /** Transfers the smallest of the large values to the heap of small values and returns the index
  where it ended up. This will shrink the large heap and grow the large heap by one element. */
  int moveFirstLargeToSmall()
  {
    T e = large.extractFirst();
    return small.insert(e);
  }
  // hmm - return value should be always 0 anyway

  int moveFirstSmallToLarge()
  {
    T e = small.extractFirst();
    return large.insert(e) + small.getSize();
  }
  // return value should always be small.getSize(), size measured after the operation


  //T& operator[](int i) { return atIndex(i); }
  // get rid or make private - client code should use either atIndex or atKey

  /** Element access via an integer index. If nS is the number of values in the small heap, indices
  i < nS refer directly to samples in the small heap with that same small-heap-index i, whereas 
  indices >= nS are interpreted as large-heap-index i-nS into the large heap. */
  T& atIndex(int i)
  {
    rsAssert(isIndexValid(i), "Index out of range");
    int nS = small.getSize();
    if(i < nS)
      return small[i];
    else
      return large[i-nS];
  }

  /** Element access via an integer key. In this class, the key and index are the same thing, but 
  not in subclass rsDoubleHeap2. We implement the atKey function here too, to have the same 
  interface as the subclass. */
  T& atKey(int k) { return atIndex(k); }

  /** Conversion from index to key. */
  int indexToKey(int i) { return i; }

  /** Conversion from key to index. */
  int keyToIndex(int k) { return k; }



  bool isIndexValid(int i)
  {
    return i >= 0 && i < small.getSize() + large.getSize();
  }

  int getNumSmallValues() const { return (int) small.getSize(); }
  int getNumLargeValues() const { return (int) large.getSize(); }
  int getNumValues()      const { return (int) (small.getSize() + large.getSize()); }

  T getLargestSmallValue()   { return small[0]; }
  T getSmallestLargeValue()  { return large[0]; }
  T getLastSmallValue()      { return small[small.getSize()-1]; }
  T getLastLargeValue()      { return large[large.getSize()-1]; }
  // make them const
  // rename to getFirstSmallValue getLast.. the last is not necessarily the largest or smallest
  // the children are both larger/smaller but their order is not defined...hmm - can we define 
  // it to give a heap even more structure?


  /*
  T getMean()
  {


  }
  */

  /*
  T getMedian()
  {
    int L = small.getSize() + large.getSize();
    if(rsIsOdd(L))
    {
      return getSmallestLargeValue();
    }
    else
      return T(0.5) * (getLargestSmallValue() + getSmallestLargeValue());
  }
  // wait - no - this is only the median, if small.size() == large.size() or the difference is only
  // one..maybe
  */

  /** Returns true, iff this object satisfies the double-heap property. Meant mostly for testing 
  and debugging (costly!). */
  bool isDoubleHeap()
  {
    return small.isHeap() && large.isHeap()   
      && !(small.less(large[0], small[0]));  // last condition means small[0] <= large[0]
  }
  // maybe we should store references to the less and swap functions here too and also provide a
  // setLessFunction function - then we may get rid of the friend declarations in rsBinaryHeap


protected:

  rsBinaryHeap<T> small, large;  // the two heaps for the small and large numbers

  template<class U> friend class rsQuantileFilterCore; // temporary - try to get rid

};

/** Subclass that uses a different convention how the indices are interpreted. Instead of 
interpreting indices i >= nS (nS = number of small values) as i-nS into the large buffer, we use
unsigned indices and use the first bit as indicator, if an element from the large heap is meant,
in which case we use the remaining bits as actual index. This may be less natural and less 
convenient, but we need it when we want to deal with a dynamically changing nS - we need some
indicator that is independent from nS, because otherwise, the stored indices in the quantile 
filter become meaningless after a change of the heap-sizes. */
template<class T>
class rsDoubleHeap2 : public rsDoubleHeap<T> 
{

// maybe get rid of the subclass an implement everything in the baseclass - maybe have functions
// replaceByIndex, replaceByKey

public:

  int replace(int key, const T& newValue)
  {
    rsAssert(isKeyValid(key), "Key out of range");

    int k = key;
 
    // The actual replacement:
    if(isKeyInLargeHeap(k)) {
      k  = large.replace(toLargeHeapIndex(k), newValue);
      k |= firstBitOnly; }  
    // should we do  (k - largeIndexOffset) | firstBitOnly?
    else
      k = small.replace(k, newValue);

    // The potential swap:
    if(small.less(large[0], small[0]))
    {
      small.swap(small[0], large[0]);
      int ks = small.floatDown(0);
      int kl = large.floatDown(0);
      if(isKeyInLargeHeap(key)) // new value was in large and is now in small heap
        k = ks;
      else
        k = kl | firstBitOnly;
        // should we do  (kl - largeIndexOffset) | firstBitOnly?
    }

    return k;
  }

  int insert(const T& newValue)
  {
    if(small.less(newValue, large[0]))
      return small.insert(newValue);
    else
      return large.insert(newValue) | firstBitOnly;
  }

  void remove(int key)
  {
    rsAssert(isKeyValid(key), "Key out of range");
    if(isKeyInLargeHeap(key))
      large.remove(toLargeHeapIndex(key));
    else
      small.remove(key);
  }

  T& atKey(int k) 
  { 
    rsAssert(isKeyValid(k), "Key out of range");
    if(isKeyInLargeHeap(k))
      return large[toLargeHeapIndex(k)];
    else
      return small[k];
  }

  /** Returns a preliminary key which is a key from the end of either the small or large heap. 
  This is needed in rsQuantileFilterCore for appending nodes */
  int getPreliminaryKey(const T& newValue)
  {
    if(small.less(newValue, large[0]))
      return small.getSize();
    else
      return large.getSize() | firstBitOnly;
  }

  bool isKeyValid(int k) const 
  {
    if(isKeyInLargeHeap(k)) 
      return toLargeHeapIndex(k) < large.getSize();
    else
      return k < small.getSize();
  }

  static inline bool isKeyInLargeHeap(int k) 
  { 
    return k & firstBitOnly; 
  }

  inline int toLargeHeapIndex(int k) const
  {
    return rawLargeHeapIndex(k) + largeIndexOffset;
  }

  int indexToKey(int i) 
  { 
    int nS = small.getSize();
    if(i < nS)
      return i;
    else
      return (i-nS) | firstBitOnly;  
    // do we need " + largeIndexOffset" here too? if so, is should be done before the "or". or 
    // maybe we need to subtract it? but no...
  }
  // takes the double-heap index from 0....nS+nL-1
  // needs test

  int keyToIndex(int k) 
  { 
    if(isKeyInLargeHeap(k))
    {
      //int i = rawLargeHeapIndex(k);  // no offset is needed here?
      int i = toLargeHeapIndex(k);  // no offset is needed here?
      return i + small.getSize();
    }
    else
      return k;
  }
  // returns the double-heap index from 0....nS+nL-1
  // needs test



  // maybe have functions to convert back and forth between key and indices

  // have static const members for the masks that separate the first bit and remove the first
  // bit

  // maybe use int instead of size_t because otherwise, when convertig indices to int, the 
  // indicator bit gets cut off ...but maybe we should take care to avoid such conversions




private:


  /** Returns the raw large-heap index, which is the key k with the first bit shaved off. */
  static inline int rawLargeHeapIndex(int k) // rename to largeHeapKeyToIndex
  {
    return k & allBitsButFirst;
  }
  // maybe make a function toLargeHeapIndex(int k) 

  int largeIndexOffset = 0; // get rid

};

/** Yet another variant. This uses again the convention to use i-nS as large-heap-index if 
i >= nS, but uses the offset as in rsDoubleHeap2 as well. ....i'm trying out some stuff to figure 
out what could work to make modulation of L,p possible in rsMovingQuantileFilterCore... 

OK - rsDoubleHeap2 works now - this renders this code obsolete - delete it*/

template<class T>
class rsDoubleHeap3 : public rsDoubleHeap<T>
{

  // we use private inheritance to make atIndex, etc. unavailable...maybe override them instead
  // ...

public:

  int replace(int key, const T& newValue)
  {
    rsAssert(isKeyValid(key), "Invalid key");
    int i  = key; 
    int nS = small.getSize();

    // replacement:
    if(isKeyInLargeHeap(i))
      i = large.replace(toLargeHeapIndex(i), newValue);
    else
      i = small.replace(i, newValue);

    // potential swap:
    if(small.less(large[0], small[0])) {
      small.swap(small[0], large[0]);
      int is = small.floatDown(0);
      int il = large.floatDown(0);
      if(isKeyInLargeHeap(key))
        i = is; // replacement took place in large, then we swapped, now it's in small
      else
        i = il + nS; }

    return i;   // return the new key/index
  }
  
  T& atIndex(int i)
  {
    rsAssert(isIndexValid(i), "Invalid index");
    if(isKeyInLargeHeap(i))  // works, because keys and indices are the same thing
      return large[i-small.getSize()];
    else
      return small[i];
  }

  T& atKey(int k) 
  { 
    rsAssert(isKeyValid(k), "Invalid key");
    if(isKeyInLargeHeap(k))
      return large[toLargeHeapIndex(k)];
    else
      return small[k];
  }

  bool isKeyValid(int k) const
  {
    if(isKeyInLargeHeap(k)) {
      k = toLargeHeapIndex(k);
      return k < large.getSize(); }
    else
      return k < small.getSize();
  }

  //int indexToKey(int i) { return i; } 
  //int keyToIndex(int k) { return k; }
  // not needed - the inherited versions do the same

  void incrementLargeIndexOffset() { largeIndexOffset++; }
  void decrementLargeIndexOffset() { largeIndexOffset--; }

protected:

  inline bool isKeyInLargeHeap(int k) const { return k >= small.getSize(); }
  inline int  toLargeHeapIndex(int k) const { return (k + largeIndexOffset) - small.getSize(); }
  int largeIndexOffset = 0;
  // the offset turned out to be not needed after all - delete it after the implementation of 
  // modulation of L,p is complete (maybe it will turn out that we will need it again, but i don't
  // think so

};
// maybe, if it works and is not too much extra code, absorb in baseclass..but maybe not - it's 
// kinda weird stuff and perhaps not always useful




//=================================================================================================


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



/** Class to exctract moving quantiles (such as the median) from a signal in realtime. If the 
quantile runs over N samples, the filter takes O(log(N)) operations per sample. It achieves this
by using an rsDoubleHeap together with a circular buffer of indices into that double-heap. The 
process that takes place in getSample is to replace the oldest sample in the double-heap with the 
new incoming sample. The potential re-ordering of the heaps due to such an replacement is kept 
track of by the circular buffer, such that it always points to the oldest sample in the 
double-heap.  */

template<class T>
class rsQuantileFilterCore
{

public:


  rsQuantileFilterCore(int maxLength = 2)
  { 
    rsAssert(maxLength >= 2, "A maxLength of at least 2 is needed");
    auto swapNodes = [&](Node& a, Node& b) { this->swapNodes(a, b); };
    dblHp.setSwapFunction(swapNodes);
    setMaxLength(maxLength);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the maximum length of the filter. This may re-allocate memory. */
  void setMaxLength(int newMaxLength)
  {
    // newMaxLength = rsNextPowerOfTwo(newMaxLength); // use later for optimizing wraparounds
    // mask = newMaxLength-1;

    small.resize(newMaxLength);
    large.resize(newMaxLength);
    keyBuf.resize(newMaxLength);
    setLengthAndReadPosition(L, p, true);
  }
  // maybe try to optimize the memory usage - we actually just need one nodes array of length
  // newMaxLength of which one part is used for the min-heap and the other for the max-heap - but 
  // maybe that would make modulating the length more difficult - if we use only one buffer, we may 
  // have to move around more data, when the length changes (i guess) - we'll see - hmm maybe not

  /** Sets the length of the filter, i.e. the number of past samples that we look at. */
  void setLength(int newLength, bool hard = false) 
  { setLengthAndReadPosition(newLength, p, hard); }

  /** Sets the read position in the sorted array of stored past values. This array does not exist 
  literally but only conceptually (in the naive implementation, this actually exists literally). In 
  practice it's the largest of the small values, i.e. the front element of the min-heap of large 
  values in our double-heap. So what this function actually does is to update the sizes of the 
  max-heap (of small values) and the min-heap (of large values) in the double-heap. But 
  conceptually, think about it simply as the readout index in an array of stored past values. */
  void setReadPosition(int newPosition, bool hard = false) 
  { setLengthAndReadPosition(L, newPosition, hard); }

  /** Sets length and read-position at once. */
  void setLengthAndReadPosition(int newLength, int newPosition, bool hard = false);

  /** When a higher level class wants to implement non-integer readout positions, we need to do a 
  linear interpolation between the values at sorted-array positions to the left and to the right of 
  actual requested non-integer position. This sets the weight w for the value to the right, the 
  weight for the value to the left is the 1-w. By default, w = 1, so we give full weight to the 
  sample to the right which is the smallest of the large values. */
  void setRightWeight(T newWeight) { w = newWeight; }

  /** Sets a pointer to a signal buffer that is driven by client code to facilitate artifact-free
  modulation of the length. If you leave this unassigned, you may see artifacts when the length
  of the filter is switched from a lower to a higher value, due to the fact that we don't know the 
  older samples here and just have to make up some values to fill up the heaps. When the true old 
  samples are available via such a buffer, these artifacts can be avoided. As said, client code 
  is supposed to feed/drive this buffer because this will also facilitate sharing of such buffers 
  between a bunch of parallel quantile filters or to create a highpass version. That means before
  calling setLength (and then getSample) on this object, client code should have called getSample
  on the buffer object. */
  void setModulationBuffer(rsDelayBuffer<T>* newBuffer) { sigBuf = newBuffer; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the length of the filter in samples. */
  int getLength() const { return L; }

  /** Returns the readout position in the (conceptual, not actually existent) array of sorted past
  input samples. */
  int getReadPosition() const { return p; }

  /** Returns the maximum length (in samples) that the filter currently supports, i.e. the capacity
  of the allocated memory for the buffers and heaps. */
  int getMaxLength() const { return (int) small.size(); } // large.size() == buf.size();


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes an output sample from a given input sample x. */
  T getSample(T x)
  {
    storeInput(x);       // updates double-heap and key-buffer
    return readOutput(); // sample readout
  }

  /** Resets the filter into its initial state. */
  void reset()
  {
    for(int n = 0; n < L; n++) {
      dblHp.atIndex(n).value  = T(0);
      int k = dblHp.indexToKey(n);
      dblHp.atIndex(n).bufIdx = n;
      keyBuf[n] = k; }
    bufIdx = 0;
  }


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Internals */

  /** Accepts a new input sample and updates our internal buffers accordingly. This will have 
  (conceptually) the effect that the oldest input sample will be remvoved from the double-heap and 
  the new input will be inserted (what actually happens is a replacement, but that's an 
  implementation detail). */
  void storeInput(T x)
  {
    //rsAssert(isStateConsistent(), "inconsistent state");
    int k = keyBuf[bufIdx];        // heap-key of oldest sample
    int w = wrap(bufIdx + L);      // (write) index of new node in keyBuf
    keyBuf[w] = k;                 // store preliminary key (== old node's key) in kexBuf at w
    dblHp.replace(k, Node(x, w));  // replace the old node, reshuffles heaps and keyBuf
    bufIdx = wrap(bufIdx + 1);     // update position in circular buffer of keys
    //rsAssert(isStateConsistent(), "inconsistent state");
  }

  /** Produces the current ouput sample by reading out the largest of the small values and the 
  smallest of the large values and forming a linear combination. Has been factored out from 
  getSample because it's also needed for modulating the length when no delay-buffer of old samples 
  is available (i.e. when sigBuffer == nullptr because the user has not assigned it). */
  T readOutput()
  {
    T yS = dblHp.getLargestSmallValue().value;   // smaller value
    T yL = dblHp.getSmallestLargeValue().value;  // larger value
    T y  = (T(1)-w)*yS + w*yL;                   // linear interpolation
    return y;
  }

  /** This is called from setLengthAndReadPosition when its "hard" parameter is false. Instead of
  setting up a new length and position immediately and doing a hard reset, this function adapts the 
  sizes of the two heaps by removing or inserting data into the heaps and/or moving data between
  the small and large heap. If N is the number of nodes that need to be removed/moved/inserted, 
  the cost of this function is O(N*log(N)). */
  void modulateLengthAndReadPosition(int newLength, int newPosition);

  /** Moves the first (smallest) value of the large heap over to the small heap, such that it 
  becomes the new largest element of the small heap. This decreases the size of the large heap by 
  one and increases the size of the small heap by one. The return value informs whether this was 
  successful - it will fail when the large heap has only one element in it because both heaps must
  always contain at least one element. */
  bool moveFirstLargeToSmall();

  /** Analog to moveFirstLargeToSmall. */
  bool moveFirstSmallToLarge();

  /** Discards the oldest sample. This will shrink one of the two heaps, wherever the oldest sample
  happens to be. This will also shrink the circular buffer by one. */
  void discardOldestSample();

  /** Adds a (dummy) sample to one of the heaps as new oldest sample. This will also grow the 
  circular buffer by one. */
  void addOlderSample();

  /** Conversion from a coneceptual buffer index in the range 0..L-1 to the actual index in the 
  range 0..buf.size()-1 */
  //int convertBufIndex(int indexFromOldest)
  //{ return (indexFromOldest + bufIdx) % (int)buf.size(); }
  // we don't really need this anymore - this was meant for scaffold code

  /** Wraparound for the bufIdx (handles only forward wraparound, i.e. input should be 
  non-negative). */
  int wrap(int i) { rsAssert(i >= 0); return i % (int)keyBuf.size(); }
  // todo: optimize using a bitmask -> restrict buffer-sizes to powers of two
 

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  /** A node stores an incoming signal value together with its index in the circular buffer. The 
  "less-than" comparison is based on the signal value. */
  struct Node
  {
    T value = T(0);
    int bufIdx = 0;   // maybe use size_t for better compatibility with rsDelayBuffer

    Node(T v = T(0), int i = 0) { value = v; bufIdx = i; }

    bool operator<(const Node& b) const 
    { return this->value < b.value; }

    bool operator==(const Node& b) const 
    { return this->value == b.value && this->bufIdx == b.bufIdx; }
  };

  /** The swapping for nodes must do the actual swap of a and b as usual but also let the circular
  buffer keep track of what gets swapped. */ 
  void swapNodes(Node& a, Node& b)
  {
    rsSwap(a, b);
    rsSwap(keyBuf[a.bufIdx], keyBuf[b.bufIdx]);
  }

  std::vector<Node>   small, large;     // storage arrays of the nodes
  rsDoubleHeap2<Node> dblHp;            // maintains large/small as double-heap
  std::vector<int>    keyBuf;           // circular buffer of heap keys - rename to keyBuf
  rsDelayBuffer<T>*   sigBuf = nullptr; // (possibly shared) buffer of delayed input samples

  int bufIdx = 0;    // index into keyBuf, mayb rename to keyIdx
  int L      = 2;    // total length of filter
  int p      = 1;    // readout position, 1 <= p <= L-1
  T   w      = T(1); // weight for smallest large value in the linear interpolation


  //-----------------------------------------------------------------------------------------------
  /** \name Scaffolding (may be removed when class is finished) */

  // some scaffolding code for development:
  //int getNodeKey(Node node);  // O(N)
  //void fixInconsistentBufferKeys(Node startNode);
  void makeBufferConsistent();

  // self-tests for debugging:
  bool isStateConsistent(); 
  bool isNodeConsistent(const Node& n);
  bool isBufferSlotConsistent(int i);
  // maybe don't delete them, these are nice for documentation purposes

};

//=================================================================================================

/** This is a more user-friendly, audio-oriented, plugin-ready wrapper around the core algorithm
of the moving quantile filter. It uses a more intuitive, audio-friendly parametrization and also 
adds features such as highpass mode, .... */

template<class T>
class rsQuantileFilter
{

public:

  rsQuantileFilter()
  {
    allocateResources();
    core1.setModulationBuffer(&delayLine);
    core2.setModulationBuffer(&delayLine);
    dirty = true;
    //updateInternals(); // hangs - why?
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Enumeration of availabel filter modes. */
  enum class Mode
  {
    lowpass,
    highpass,
    bandpass,
    bandstop  
    // todo: lowshelf, highshelf, bandshelf
  };

  /** Sets the mode of the filter. @see Mode */
  void setMode(Mode newMode) { mode = newMode; }

  void setSampleRate(T newSampleRate)
  {
    sampleRate = newSampleRate;
    allocateResources();
    dirty = true;
  }

  void setFrequency(T newFrequency)
  {
    setLength1(T(1) / newFrequency); 
    // maybe experiment with proportionality factors - figure out, where the cutoff frequency is
    // for a median filter by feeding sinusoidal inputs and measure the output amplitude (the 
    // output may be a distorted sine? or maybe it will be zero, when the length is exactly one
    // cycle?...yeah - that makes sense and the same is true for any symmetric waveform - that's 
    // going to be some weird ass nonlinear filter!)
  }


  void setMaxLength(T newMaxLength)
  {
    maxLength = newMaxLength;
    allocateResources();
    dirty = true;
  }

  void setSampleRateAndMaxLength(T newSampleRate, T newMaxLength)
  {
    sampleRate = newSampleRate;
    maxLength  = newMaxLength;
    dirty = true;
  }

  void setLength1(T newLength) { length1 = newLength; dirty = true; }
  // get rid - user should use setFrequency

  void setQuantile1(T newQuantile1) { quantile1 = newQuantile1; dirty = true; }

  //void setQuantile2(T newQuantile2) { quantile2 = newQuantile2; dirty = true; }


  // getDelay()...but what should it return in bandpass mode? or maybe have two functions
  // getDelay1, getDelay2

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  T getSample(T x)
  {
    delayLine.getSample(x);
    // Must be called *before* updateInternals because in updateInternals, we may modulate the
    // length and the artifact-avoidance strategy for the modulation assumes that we have already
    // called getSample on the delay buffer. At this point, we are not (yet) interested in the 
    // output of the delayline - we will retrieve ouputs later via random access using the [] 
    // operator

    if(dirty) 
      updateInternals();

    T y1 = core1.getSample(x);
    if(mode == Mode::lowpass)
      y = y1;
    else if(mode == Mode::highpass)
      y = delayLine[delay1] - y;      // highpass = delayed_input - lowpass
    else if(mode == Mode::bandpass)
    {
      y = delayLine[delay1] - y;      // use core1 as highpass...
      y = core2.getSample(y);         // ...feeding core2 which is the lowpass
    }
    else if(mode == Mode::bandstop)
    {
      T y2 = core2.getSample(x);
      y2 = delayLine[delay2] - y2;
      y = y1 + y2;
    }

    //y2 = core2.getSample(x);


    // experimental: highpass mode:
    //

    //y = y1;  // preliminary

    return y;
  }
  // maybe make a function to produce lowpass and highpass at the same time - client code may find
  // that useful

  void reset()
  {
    core1.reset();
    core2.reset();
    delayLine.reset();
    y = T(0);
  }


protected:

  /** Computes filter algorithm parameters length L, readout point p (both in samples) and weight 
  w for the linear interpolation from user parameters length (in seconds), quantile 
  (normalized 0..1) and sampleRate. ...q is the non-integer read-position which is equal to the
  introduced delay (verify!) */
  void convertParameters(T length, T quantile, T sampleRate, int* L, int* p, T* w, T* q)
  {
    rsAssert(quantile >= T(0) && quantile <= T(1), "Quantile needs to be between 0 and 1");
    *L  = (int) round(length * sampleRate);  // length of filter in samples
    *L  = rsMax(*L, 2);                      // ...needs to be at least 2
    *q  = quantile * sampleRate * (*L - 1);  // readout position in sorted array
    *p  = (int) floor(*q);                   // integer part (floor)
    *w  = *q - *p;                           // fractional part
    *p += 1;                                 // algo wants the next one
    if(*p > *L - 1) {                        // quantile == 1 (maximum) needs special care
      *p = *L - 1; *w = T(1);  }
  }

  /** Updates the internal algorithm parameters and embedded objects according to the user 
  parameters. */
  void updateInternals()
  {
    // compute internal core parameters:
    int L, p; T w;
    convertParameters(length1, quantile1, sampleRate, &L, &p, &w, &delay1);

    // Set the core up:
    core1.setLengthAndReadPosition(L, p);
    core1.setRightWeight(w);

    // todo: set up a second core (for bandpass mode) 

    dirty = false;
  }
  // when we want to have two filters for implementing bandpass, etc. we may need a function to 
  // compute L,p,w from length,sampleRate,quantile

  void allocateResources()
  {
    int mL = (int) ceil(maxLength * sampleRate);
    core1.setMaxLength(mL);
    core2.setMaxLength(mL);
    delayLine.setCapacity(mL);
  }

  // user parameters:
  //T feedback   = 0.0;
  T sampleRate = 44100;         // in Hz
  T maxLength  = 0.1;           // in seconds
  T length1    = 0.01;          // in seconds  ...move to algo parameters or get rid
  T quantile1  = 0.5;           // in 0..1, 0.5 is median
  T bandwidth  = 1.0;           // in octaves
  Mode mode    = Mode::lowpass;

  // filter state, algo parameters, infrastructure:
  T y = 0.0;                       // previous output sample (for feedback)
  T delay1 = 0.0;                  // required delay to implement highpass by subtraction
  T delay2 = 0.0;
  std::atomic<bool> dirty = true;  // flag to indicate that algo params must be recalculated

  // embedded objects:
  rsQuantileFilterCore<T> core1, core2;
  rsDelayBuffer<T> delayLine;

};
// todo: make a highpass version by subtracting the result from a (delayed) original signal - maybe
// in a subclass rsMultiModeQuantileFilter - contains two lowpasses with different cutoff-freq, one 
// of which may be turned into a highpass...or maybe have two cores here and a delayline



/** This is a naive implementattion of (the core of) a moving quantile filter and meant only for 
producing test outputs to compare the production version of the rsQuantileFilterCore against. 
It's horribly inefficient - the cost per sample is O(N*log(N)) whereas the production version 
should run in O(log(N)). ..maybe move to unit tests.. */

template<class T>
class rsQuantileFilterNaive
{

public:

  rsQuantileFilterNaive(int numSmaller = 20, int numLarger = 20)
    : buf(numSmaller+numLarger)
  {
    setLengths(numSmaller, numLarger);
  }

  void setLengths(int numSmaller, int numLarger)
  {
    nS = numSmaller;
    nL = numLarger;
    updateBufferLengths();
  }

  int getLength() const { return nS + nL; }


  void prepareSortedDelayBuffer(T x)
  {
    T y = buf.getSample(x);
    buf.copyTo(&tmp[0], true);
    rsHeapSort(&tmp[0], getLength());
  }

  T getSample(T x)
  {
    prepareSortedDelayBuffer(x);
    return tmp[nS];
    // pehaps, we should use interpolation here - provide functions setLength, setQuantile
  }

  T getSampleMedian(T x)
  {
    prepareSortedDelayBuffer(x);
    int L = getLength();
    T p = 0.5 * (L-1);
    int i = (int) floor(p);
    T   f = p - i;
    //return f*tmp[i] + (1-f)*tmp[i+1]; // this is wrong, but matches rsQuantileFilter
    return (1-f)*tmp[i] + f*tmp[i+1]; // this is correct and matches rsArrayTools::median
    // verify, if we use f and 1-f correctly
  }


  void reset()
  {
    buf.reset();
    //for(int i = 0; i < getLength(); i++)
    //  buf[i] = 0;
  }


protected:

  void updateBufferLengths()
  {
    buf.setLength((size_t)getLength());
    tmp.resize(   (size_t)getLength());
    //reset();
  }

  rsDelayBuffer<T> buf;   // circular buffer
  std::vector<T> tmp;
  int nS = 0, nL = 0;    // maybe use size_t

};

//=================================================================================================

template<class T>
class rsNoiseGeneratorTriModal
{

public:

  rsNoiseGeneratorTriModal()
  {
    selectorLowpass.setMode(selectorLowpass.LOWPASS_IIT);
    selectorLowpass.setSampleRate(1);
    selectorLowpass.setCutoff(0.1);
    selector.setRange(-1, +1);
    ng1.setRange(-1.0, -0.3);
    ng2.setRange(-0.3, +0.3);
    ng3.setRange(+0.3, +1.0);
    setOrder(7);
  }

  void setOrder(int newOrder)
  {
    ng1.setOrder(newOrder);
    ng2.setOrder(newOrder);
    ng3.setOrder(newOrder);
  }

  inline T getSample()
  {
    T s = selector.getSample();
    s = selectorLowpass.getSample(s);
    if(s < thresh1)
      return ng1.getSample();
    if(s < thresh2)
      return ng2.getSample();
    return ng3.getSample();
  }

  rsOnePoleFilter<T, T> selectorLowpass;
  // the selector is lowpassed such that successive samples tend to be selected from the same 
  // distribution


protected:

  rsNoiseGenerator<T> selector;
  rsNoiseGenerator2<T> ng1, ng2, ng3;
  T thresh1 = -0.2, thresh2 = +0.2;

};

//=================================================================================================

template<class T>
void fitQuadratic(T x1, T y1, T x2, T y2, T x3, T y3, T* a0, T* a1, T* a2)
{
  T k1 = y1 / ((x1-x2)*(x1-x3));
  T k2 = y2 / ((x2-x1)*(x2-x3));
  T k3 = y3 / ((x3-x1)*(x3-x2));
  T b1 = -k1*(x2+x3);
  T b2 = -k2*(x1+x3);
  T b3 = -k3*(x1+x2);
  T c1 = k1*x2*x3;
  T c2 = k2*x1*x3;
  T c3 = k3*x1*x2;
  *a2  = k1 + k2 + k3;  // coeff for x^2
  *a1  = b1 + b2 + b3;  // coeff for x^1
  *a0  = c1 + c2 + c3;  // coeff for x^0

  // Formulas were derived from setting up 3 polynomials in product form, where each has zeros at 
  // all but one of the datapoints, say xi, and to have value yi at xi and then adding them up 
  // (idea due to Lagrange):
  //   p1(x) = k1*(x-x2)*(x-x3)       p1 has zeros at at x2,x3
  //   p2(x) = k2*(x-x1)*(x-x3)       p2 has zeros at at x1,x3
  //   p3(x) = k3*(x-x1)*(x-x2)       p3 has zeros at at x1,x2
  // Require:
  //   p1(x1) = y1, p2(x2) = y2, p3(x3) = y3
  // Solve these for the ki, i.e. k1,k2,k3. For example, k1 = y1 / ((x1-x2)*(x1-x3)). Plug, for 
  // example, k1 back into the p1 equation and multiply it out to obtain its coeffs - do the same 
  // for p2 and p3 and then obtain the final polynomial coeffs by adding the corresponding  coeffs 
  // of each of the partial polynomials.

  // operations: add: 9, sub: 6, mul: 12, div: 3, neg: 3, ass: 12, tmp: 9
}

//=================================================================================================
// the stuff below is just for playing around - maybe move code elsewhere, like the research-repo:

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

  rsGroupString(int length) { s.resize(length); }

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

  /** Read/write access to i-th character. */
  unsigned int& operator[](int i) { return s[i]; }

  unsigned int last() const 
  { 
    if(s.size() > 0) return s[s.size()-1]; 
    else             return 0; // is this reasonable? or should we use a special "error" signal
  }

  void append(unsigned int x) { s.push_back(x); }

  void removeLast()
  {
    if(s.size() > 0)
      s.pop_back();
  }

  /** Length of the string. */
  int length() const { return (int) s.size(); }



  /** Returns the (additive) inverse which is just the string in reversed order. */
  rsGroupString inverse() const;
   // maybe later (when we have multiplication), rename to additiveInverse and implement a 
   // multiplicativeInverse, too
   // then the class should be renamed to fieldString

   // maybe let integers 0 and 1 be used and implement 1/s = s.multiplicativeInverse, etc.


  void resize(int newSize) { s.resize((int)newSize); }

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
