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

/** Baseclass for rsBinaryHeap and rsBinarySearchTree...

*/


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
    //buildTree();   // maybe make this call optional
  }

  void setData(std::vector<T>& newData)
  {
    setData(&newData[0], (int)newData.size(), (int)newData.size());
  }

  /** Sets the comparison function to be used. If it implements "less-than", you'll get a max-heap
  and if it implements "greater-than", you'll get a min-heap. By default, it's assigned to a 
  less-than function based on the < operator. */
  void setCompareFunction(const std::function<bool(const T&, const T&)>& newFunc)
  {
    less = newFunc;
  }
  // question: what happens, if we use a less-or-equal or greater-or-equal function for this?

  /** Sets the function used to swap two elements. By default, it just uses rsSwap but you may want
  to do something extra whenever a swap takes place in certain circumstances, for example, to keep 
  track of when items get moved around (see rsMovingQuantileFilter for an example). */
  void setSwapFunction(const std::function<void(T&, T&)>& newFunc)
  {
    swap = newFunc;
  }


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
  // should perhaps also be moved to rsBinaryHeap again - the replacement in the rsBinarySearchTree
  // needs a different approach (i think)

  /*
  // todo: 
  int insert(const T& x)
  {
  rsAssert(size < capacity, "Capacity exceeded");
  // ...maybe in this case, we should return -1 immediately

  data[size] = x;
  size++;
  return floatUp(size-1);
  }
  // not yet tested
  */

  // int insert(const T& x); // returns insertion index
  // void remove(int);       // removes item at index i
  // removing of element i should work as follows:
  // -check, which of the two child-nodes should get promoted to the position i
  // -move left(i) or right(i) up and apply the same process to its children (i.e check, which of 
  //  its children should get promoted up)
  // -and so on all the way to the bottom
  // -decrement size


  /** Lets the node at index i float up or down into its appropriate place and returns the 
  resulting new array index. */
  int floatIntoPlace(int i) { return floatUp(floatDown(i)); }
  // -maybe rename to rebalanceNode
  // -why is this public? can we move it into the protected section? hmm, rsMovingPercentileFilter
  //  calls it directly - can this be avoided?

  // get rid - the user should call buildHeap or buildSearchTree explicitly
  /*
  void buildTree()
  {
    for(int i = size/2-1; i >= 0; i--)  // or should we use (size-1)/2 ?
      floatDown(i);
  }
  */
  // Runs in O(N). From the code, it would appear as having O(N*log(N)) complexity because we call
  // an O(log(N)) function inside the loop. However, for heaps and binary search trees, the runtime
  // of floatDown depends on the argument i in such a way to give an overall O(N) behavior (see 
  // reference (1)). The memory complexity is O(1).

protected:

  /** Must be implemented by subclasses to let node with array index i float up, if necessarry. It 
  should return the new array index of the node. */
  virtual int floatUp(int i) = 0;

  /** Must be implemented by subclasses to let node with array index i float down, if 
  necessarry. It should return the new array index of the node. */
  virtual int floatDown(int i) = 0;

  /** Index of parent of node i. */
  static inline int parent(int i) { return (i-1) / 2; }

  /** Index of left child of node i. */
  static inline int left(int i)   { return 2*i + 1; }

  /** Index of right child of node i. */
  static inline int right(int i)  { return 2*i + 2; }

  /** Returns true, iff index i is the index of a left child node. */
  static inline bool isLeft(int i) { return i == left(parent(i)); }
  // needs test

  /** Returns true, iff index i is the index of a right child node. */
  static inline bool isRight(int i) { return i == right(parent(i)); }
  // needs test

  // The actual data array:
  T* data = nullptr;
  int size = 0;
  int capacity = 0;

  // Comparison and swapping functions:
  //std::function<bool(const T&, const T&)> less;  // rename to comp
  //std::function<void(T&, T&)> swap;


  // maybe these should be references?
  std::function<bool(const T&, const T&)> 
    less = [](const T& a, const T& b)->bool 
  { 
    return a < b; 
  };

  std::function<void(T&, T&)> 
    swap = [](T& a, T& b) { rsSwap(a, b); };


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


protected:

  /** Functions to establish or maintain the heap-property of the underlying data array. */

  /** Lets the node at index i float up the tree if its parent violates the heap property. Returns 
  the new array index. */
  int floatUp(int i) override
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

  int floatDown(int i) override
  {
    //return floatDownRec(i);  // calls recursive version - change that!
    return floatDownIt(i);
  }

  // todo: implement functions: int insert(T*), void remove(int), extractMax,
  // increaseKey, heapMax

  template<class U> friend class rsDoubleHeap;

};
// -maybe make also a class rsBinarySearchTree where the left child is <= and the right child is >=
//  the parent
// -pehaps we could implement a general function: needsSwap(int parentIndex, int childIndex) - in 
//  the case of a search tree, it would first have to figure out, if the childIndex is a left or 
//  right child and order the arguments of the comparison according to that
// -in the old RSLib codebase, i did some sort of linked-tree - maybe that could be dragged in as
//  rsLinkedTree or rsDynamicTree or something like that - all the different trees could be in
//  a file Trees.h/cpp

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

  int floatUp(int i) override
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

  int floatDown(int i) override
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


  /** Replaces the element at index i with the given new values and returns the the index where the
  new element actual ended up after doing all the floating up/down and potential swapping business.
  If nS is the size of the small heap, we use the convention that indices i < nS are direct indices 
  into the small heap, and when i >= nS, we use i-nS as heap index into the large heap (as seen 
  from the class rsBinaryHeap). The same holds for the return value. */
  int replace(int index, const T& newValue)
  {
    rsAssert(isIndexValid(index), "Index out of range");

    int i  = index;  // maybe get rid of index and use only i
    int nS = small.getSize();

    // The actual replacement, taking place in one of our two heaps:
    if(index < small.getSize())
      i = small.replace(i,    newValue);
    else
      i = large.replace(i-nS, newValue) + nS;

    // The potential swap of the front elements (i.e. max and min) of both heaps:
    if( small.less(large[0], small[0]) )  // means: if(large[0] < small[0])
    {
      small.swap(small[0], large[0]);
      int is = small.floatDown(0);
      int il = large.floatDown(0);
      if(index < nS)   // we replaced in small heap, then swapped, so newValue is now in large heap
        i = il + nS;
      else
        i = is;
    }

    return i;
  }
  // we use the convention that when i < nS, the datum gets inserted into the small heap and when 
  // i >= nS, the datum gets inserted into the large heap at i-nS, tbc....


  T& operator[](int i)
  {
    rsAssert(isIndexValid(i), "Index out of range");
    int nS = small.getSize();
    if(i < nS)
      return small[i];
    else
      return large[i-nS];
  }

  bool isIndexValid(int i)
  {
    return i >= 0 && i < small.getSize() + large.getSize();
  }

  T getLargestSmallValue()
  {
    return small[0];
  }

  T getSmallestLargeValue()
  {
    return large[0];
  }


protected:

  rsBinaryHeap<T> small, large;  // the two heaps for the small and large numbers

};

//=================================================================================================



/** Class to exctract moving quantiles (such as the median) from a signal in realtime. If the 
percentile runs over N samples, the filters takes O(log(N)) operations per sample. It achieves this
by using two heaps (a max-heap of the smaller-than-percentile values and a min-heap of the 
bigger-than-percentile values) and a circular buffer in a clever way....

*/

// other implementation - hopefully simpler - less indirections
template<class T>
class rsMovingQuantileFilter
{

public:


  rsMovingQuantileFilter()
  { 
    auto swapNodes = [&](Node& a, Node& b) { this->swapNodes(a, b); };
    heaps.setSwapFunction(swapNodes);
  }



  void setMaxLength(int newMaxLength)
  {
    small.resize(newMaxLength);
    large.resize(newMaxLength);
    buf.setCapacity(newMaxLength);
    updateBuffers();
  }

  void setLength(int newLength)
  {
    L = newLength;
    updateBuffers();
  }

  void setQuantile(int newQuantile)
  {
    q = newQuantile;
    updateBuffers();
  }

  int getLength() const
  {
    return L;
  }


  bool isIndexPermutation(int* b, int L)
  {
    for(int i = 0; i < L; i++)
      if( !rsArrayTools::contains(b, L, i) )
        return false;
    return true;
  }
  // Returns true, iff b contains every number from 0 to L-1. Since b is of length L, this implies
  // that every number is contained exactly once, so b is a permutation of the numbers 0...L-1.
  // used here only for debug -  move elsewhere

  
  bool isDelayBufferValid()
  {
    std::vector<int> tmp(L);
    buf.copyTo(&tmp[0]);
    return true;  // preliminary

    //return isIndexPermutation(&tmp[0], L); // i think, i have implemented such a function in 
                                           // ScratchPad.cpp for the matrix stuff
  }
  
  // for debugging

  T getSample(T x)
  {
    // under construction...i have no idea what i'm doing...the goal is to replace the oldest 
    // sample with the new incoming sample, thereby keeping track of where it goes in the heap, 
    // storing that number in the circular buffer

    rsAssert(isDelayBufferValid()); // for debug

    int hi = buf.getOldest();       // hi: heap-index of oldest sample
    int bi = heaps[hi].bufIndex;    // bi: buffer-index of oldest sample

    int hj = heaps.replace(hi, Node(x, bi));
    int bj = heaps[hj].bufIndex;   // needed?
    // create a Node from the new value, replace the oldest node with it and retrieve the 
    // heap-index where the new node went to - this may re-shuffle the buf as well

    buf.advancePointers();

    //hi = buf.getSample(hj);
    // store the heap-index of the new value in the buffer - the return value is not really 
    // interesting anymore - we are already done with it
    // maybe we should not use getSample but updatePointers

    //buf[bi] = hi; // verify! nope!

    //bi = heaps[hi].bufIndex;  // experimental - i don't know, what i'm doing
    //buf[bi] = hi;    

    // hi should not change here from what was returned from getOldest - but it does - but maybe
    // it's because of the swaps and actually ok? we don't really need it anymore anyway - it's
    // just for verification/debugging

    rsAssert(isDelayBufferValid()); // for debug


    T y = heaps.getSmallestLargeValue().value;
    // or should it be getSmallestLargeValue? in case of non-integer q, perhaps a linear 
    // combination of both? also for medians of even length - there, we need the average of both

    // maybe check, if buf contains a permutation of the values 0...L-1 - this should always be the
    // case, if everything is correct

    return y;


    //return heaps.getSmallestLargeValue();


    //return heaps[q];
    // or maybe q-1 or q+1? what about fractional q? we need a linear combination of values at
    // q and q+1
  }


  void reset()
  {
    for(int n = 0; n < L; n++)
    {
      heaps[n].value    = T(0);
      heaps[n].bufIndex = n;      // verify!

      buf[n]            = n;      
      // verify! i think, it should not matter - any permutation of the values 0..L-1 should do
    }
  }


protected:

  void updateBuffers()
  {
    rsAssert(L <= (int) small.size());
    int C = (int) small.size();
    heaps.setData(&small[0], q, C, &large[0], L-q, C);  // verify!
    buf.setLength(L);
    reset();  // is this needed?
  }

  /** A node stores an incomiong signal value together with its index in the circular buffer. */
  struct Node
  {
    Node(T v = T(0), int i = 0)
    {
      value    = v;
      bufIndex = i;
    }

    int bufIndex = 0;
    T value = T(0);
    bool operator<(const Node& b) const { return this->value < b.value; }
  };

  //std::vector<Node>  nodes;   // stores the incoming values together with their index
  // maybe keep separate vectors for the small and large heap - makes it easier


  std::vector<Node> small, large;

  rsDoubleHeap<Node> heaps;   // for keeping the data in the nodes array "semi-sorted"
  rsRingBuffer<int>  buf;     // circular buffer of indices into the nodes array


  // The swapping function must do the actual swap of a and b as usual but also let the circular
  // buffer keep trak of what gets swapped - the next index returned by buf.getOldest or 
  // buf.getSample should always hold the (double-heap-) index to the oldest value
  void swapNodes(Node& a, Node& b)
  {
    rsSwap(a, b);


    int i = a.bufIndex;
    int j = b.bufIndex;

    // todo: swap the corresponding indices in circular buffer:
    rsSwap(buf[a.bufIndex], buf[b.bufIndex]);  // is this correct?

    rsAssert(isDelayBufferValid());  // debug

    int dummy = 0;
  }


  int L = 0;  // total length of filter
  int q = 0;  // quantile as value 0 <= q < L


  // we need a circular buffer of heap indices and twe swap function should take care of updating
  // this circular buffer, too


};


template<class T>
class rsMovingQuantileFilterOld
{

public:


  rsMovingQuantileFilterOld(int numSmaller = 20, int numLarger = 20) 
    : buf(numSmaller+numLarger) // preliminary - we need to be able to adapt the capacity of ringbuffers at runtime
  {
    setLengths(numSmaller, numLarger);

    // assign the comparison and swapping functions in both heaps:
    auto less      = [&](const int& a, const int& b)->bool { return nodeLess(   a, b); };
    auto greater   = [&](const int& a, const int& b)->bool { return nodeGreater(a, b); };
    auto swapNodes = [&](int& a, int& b) { this->swapNodes(a, b); };
    //small.setCompareFunction(greater);
    //large.setCompareFunction(less);

    small.setCompareFunction(less);    // small is a max-heap
    large.setCompareFunction(greater); // large is a min-heap
    small.setSwapFunction(swapNodes);
    large.setSwapFunction(swapNodes);
  }


  void setLengths(int numSmaller, int numLarger)
  {
    nS = numSmaller;
    nL = numLarger;
    updateBufferLengths();
  }


  int getLength() const
  {
    return nS + nL;
  }


  T getSample(T x)
  {
    // under construction

    // -replace oldest sample in the heaps with the new incoming sample
    // -rebalance the heap in which the oldest sample resided
    // -if the largest (front) sample in the small heap is > than the smallest (front) sample in 
    //  the large heap: exchange them and rebalance both heaps
    // -the swaps that occurr during rebalancing the heaps should also lead to corresponding swaps 
    //  in our delay buffer, due to the way we have implemented our swap function
    //  ...err...i think, no, the delay buffer should remain as is?
    // -return the front sample of the "large" heap in case of an odd size and the average of both
    //  front samples for even sizes (size == getLength())

    //int L 

    int ni = buf.getOldest();      // node index of oldest node in nodes array
    int hi = nodes[ni].heapIndex;  // heap index of oldest node
    //int bi = nodes[ni].bufIndex;   // buffer index of oldest node - do we need this?

    // replace value of the oldest node with the new incoming value:
    nodes[ni].value = x;
    if(hi < nS) hi = small.floatIntoPlace(hi);          // node to be replaced is in small heap
    else        hi = large.floatIntoPlace(hi-nS) + nS;  // node to be replaced is in large heap

    // swap first elements of small/large heaps, if necessarry:
    if(heapsNeedExchange())
    {
      //rsSwap(nodes[0].heapIndex, nodes[nS].heapIndex);
      //rsSwap(nodes[0].value,     nodes[nS].value);       // is this correct?
      //rsSwap(nodes[0], nodes[nS]);
      // nope! the nodes array stays fixed - we don't swap anything there

      // the swap:
      rsSwap(nodes[0].heapIndex, nodes[nS].heapIndex); // update of our additional data
      rsSwap(heaps[0], heaps[nS]);                     // update of the actual heap elements

      //rsSwap(nodes[0].value,     nodes[nS].value);       // just a test - is this correct?

      // the update that may be required due to the swap:
      if(hi < nS)
      {
        // replacement took place in small heap, so after exchange, the new datum is in the large
        // heap
        small.floatIntoPlace(0);
        hi = large.floatIntoPlace(0) + nS;
      }
      else
      {
        // replacement took place in large heap, so after exchange, the new datum is in the small
        // heap
        large.floatIntoPlace(0);
        hi = small.floatIntoPlace(0);
      }
    }


    // debug:
    bool isSmallHeap = small.isHeap();
    bool isLargeHeap = large.isHeap();
    rsAssert(isSmallHeap && isLargeHeap);


    // hi now is the new heap index of the newly received value - we write it into the nodes 
    // array:
    nodes[ni].heapIndex = hi;  // is this correct? or do we need to do it after the next?

    ni = buf.getSample( (ni+1) % getLength() );  // ..and this?
    // this is the actual update of the buf - the returned ni value here should be equal to the one
    // we have already used all the time ...can we move this function call up? i think so




    //int idx = large.getElement(0);  // index to read from the heap
    int idx = large[0]; 
    //idx = nS;  // test..hmm..nope?
    //idx = ni;
    T val = nodes[idx].value;
    return val;
  }


  void reset()
  {
    for(int i = 0; i < getLength(); i++)
    {
      heaps[i] = i;
      buf[i] = i;
      nodes[i].heapIndex = i;
      nodes[i].value = T(0);
    }
    // does this initialization make sense?
  }





protected:

  T getLargestSmallValue()
  {
    return nodes[heaps[0]].value;
  }

  T getSmallestLargeValue()
  {
    return nodes[heaps[nS]].value;
  }

  bool heapsNeedExchange()
  {
    // for debug:
    T largestSmall  = getLargestSmallValue();  
    T smallestLarge =  getSmallestLargeValue();

    return getLargestSmallValue() > getSmallestLargeValue();
  }


  // function to be passed to the heap objects, to be used there for comparison and swapping:
  bool nodeLess(const int& left, const int& right)
  {
    return nodes[left].value < nodes[right].value;
  }
  bool nodeGreater(const int& left, const int& right)
  {
    return nodes[left].value > nodes[right].value;
  }
  void swapNodes(int& i, int& j)
  {
    rsSwap(nodes[i].heapIndex, nodes[j].heapIndex); // update of our additional data
    rsSwap(i, j);                                   // update of the actual heap elements
  }
  // expects two indices into the nodes array



  void updateBufferLengths()
  {
    nodes.resize(getLength());
    buf.setLength((size_t)getLength());
    heaps.resize(getLength());
    small.setData(&heaps[0],  nS, nS);
    large.setData(&heaps[nS], nL, nL);
    reset();
  }

  /** A node stores the value itself and the index in the heaps, where we use the convenetion that
  when a node it in the heap of larger values, the actualy index as seen from the "large" heap is
  given by heapIndex-nS - and this occurrs whenever heapIndex >= nS, when heapIndex < nS, it's an
  index inot the "smaller" heap */
  struct Node
  {
    int heapIndex = 0;
    //int bufIndex = 0;  // i think, this may be redundant - we'll see
    T value = T(0);
    // maybe we need a bool to indicate, if we are in the large heap - using the convention of
    // adding/subtracting nS is sometimes inconvenient...hmm...maybe not
  };


  std::vector<Node> nodes;


  rsRingBuffer<int> buf;   // circular buffer of indices into the nodes array
  std::vector<int> heaps;  // memory for the heaps
  // heaps and buf just store indices into the nodes array and the nodes-array itself contains the
  // actual data together with its current index inside the heap (and delay-buffer). We need this 
  // indirection in order to map directly to heap- and buffer-indices in O(1) time. We use the 
  // index in the nodes-array as key and update the data heapIndex,bufIndex of the array-entries
  // whenever we do a swap do to heap-rebalancing
  // -the heaps must be able to reference the data (for comparing)
  // -from buf.getOldest(), we must be able to retrieve the the position of the oldest datum in 
  //  the heap(s), so we can know, which heap-element we wnat to replace with the new incoming 
  //  datum


  rsBinaryHeap<int> small, large;
  // is it possible to use one binary serach tree instead of two binary heaps?


  // maybe heaps and nodes are not needed? can we just use buf to buffer the elements and let
  // small and large be indices into buf? then buf should be a buffer of T values - that would 
  // simplify the implementation a lot - less indirection - we actually don't need the two-way
  // association: the order of the circular buffer never changes

  // i think, maybe i should get rid of the nodes array and keeps only heaps - and a heap node 
  // should contain the catual value together with the buffer index and whenever two heap-nodes
  // are swapped, the corresponding buffer entries are also swapped (the buffer contains indices 
  // into the heaps)


  int nS = 0, nL = 0; // number of elements smaller than and larger than the percentile
  // maybe use size_t

};

// with a binary search tree, we may just need the tree and a circular buffer that always remembers
// which nod is the oldest, then on each sample:
// -replace the datum in node that holds the oldest value
// -rebalance the tree (when swapping nodes, the node index stored in the circular delay buffer 
//  must be updated along)
// -output is always the valzue stored in a specific, fixed node - if the tree is symmetric and we
//  we want the median, then use the root




template<class T>
class rsMovingQuantileFilterNaive
{

public:

  rsMovingQuantileFilterNaive(int numSmaller = 20, int numLarger = 20)
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

  int getLength() const
  {
    return nS + nL;
  }

  T getSample(T x)
  {
    T y = buf.getSample(x);
    buf.copyTo(&tmp[0]);
    rsHeapSort(&tmp[0], getLength());
    return tmp[nS];
  }

  void reset()
  {
    for(int i = 0; i < getLength(); i++)
      buf[i] = 0;
  }


protected:

  void updateBufferLengths()
  {
    buf.setLength((size_t)getLength());
    tmp.resize(   (size_t)getLength());
    reset();
  }



  rsRingBuffer<T> buf;   // circular buffer

  std::vector<T> tmp;

  int nS = 0, nL = 0;    // maybe use size_t

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
