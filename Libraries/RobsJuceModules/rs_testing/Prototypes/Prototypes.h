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
class rsBinaryHeap
{

public:

  rsBinaryHeap(T* newData = nullptr, int newSize = 0)
  {
    setData(newData, newSize);
    less = [](const T& a, const T& b)->bool { return a < b; };
    swap = [](      T& a,       T& b)       { rsSwap(a, b); };
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the data that should be treated as heap. The object does not take ownership of the data.
  Instead, it acts pretty much like a "view" (as in rsMatrixView) - it just operates on an existing 
  data array whose lifetime is managed elsewhere. */
  void setData(T* newData, int newSize)
  {
    data = newData;
    size = newSize;
    buildHeap(); // maybe make this call optional
  }
  // todo: also pass the capacity

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

  int getSize() const { return size; }

  //int getCapacity() const { return capacity; }



  /** Returns a const reference to the element at index i. */
  const T& getElement(int i) const 
  { 
    rsAssert(i >= 0 && i < size, "Index out of range");
    return data[i]; 
  }
  // maybe get rid - implement as [] operator returning a const-ref - simpler syntax in client code

  /** Read/write access to elements. Warning: overwriting elements may destroy the heap-property, 
  so it should be used only if you know what you are doing. */
  T& operator[](int i)
  {
    rsAssert(i >= 0 && i < size, "Index out of range");
    return data[i]; 
  }

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
  // needs test
  // maybe convert recursion to iteration




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

  /** Lets the node with array-index i float into its proper place and returns the resulting new 
  array index. */
  int floatIntoPlace(int i)
  {
    return floatUp(floatDown(i));
    // this calls the recursive version of floatDown - todo: use iterative version
  }
  // -maybe rename to rebalanceNode
  // -why is this public? can we move it into the protected section? hmm, rsMovingPercentileFilter
  //  calls it directly - can this be avoided?


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
  int floatDown(int i)
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


  int floatDown2(int i)
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


  // we should return the new index i, we also need floatUp(int i)


  void buildHeap()
  {
    for(int i = size/2-1; i >= 0; i--)  // or should we use (size-1)/2 ?
    {
      //floatDown(i);
      floatDown2(i);
    }
  }
  // runs in O(N). From the code, it would appear as having O(N*log(N)) complexity because we call
  // an O(log(N)) function inside the loop. However, the runtime of floatDown depends on the 
  // argument i in such a way to give an overall O(N) behavior (see reference (1)). The memory 
  // complexity is O(1).



  // todo: implement functions: int insert(T*), void remove(int), extractMax,
  // increaseKey, heapMax


  // maybe make them static:

  /** Index of parent of node i. */
  inline int parent(int i) const { return (i-1) / 2; }

  /** Index of left child of node i. */
  inline int left(int i)   const { return 2*i + 1; }

  /** Index of right child of node i. */
  inline int right(int i)  const { return 2*i + 2; }



  /** Returns true, iff index i is the index of a left child node. */
  inline bool isLeft(int i) const { return i == left(parent(i)); }

  /** Returns true, iff index i is the index of a right child node. */
  inline bool isRight(int i) const { return i == right(parent(i)); }

  // needs test





  T* data = nullptr;
  int size = 0;
  //int capacity = 0;

  //bool (*less)(const T& a, const T& b) = &RAPT::defaultLess;
  // comparison function used - this is currently a plain function pointer - maybe use 
  // std::function or a template parameter F later
  // maybe rename to compare or comp - it's not necessarily a less-than - can also be a 
  // greater-than comparison

  // Comparison and swapping functions:
  std::function<bool(const T&, const T&)> less;  // rename to comp
  std::function<void(T&, T&)> swap;
};
// -maybe make also a class rsBinarySearchTree where the left child is <= and the right child is >=
//  the parent
// -maybe a common baseclass can be factored out (implementing parent, left, right, storing swap 
//  maybe also storing the comparison function)
// -pehaps we could implement a general function: needsSwap(int parentIndex, int childIndex) - in 
//  the case of a search tree, it would first have to figure out, if the childIndex is a left or 
//  right child and order the arguments of the comparison according to that
// -in the old RSLib codebase, i did some sort of linked-tree - maybe that could be dragged in as
//  rsLinkedTree or rsDynamicTree or something like that - all the different trees could be in
//  a file Trees.h/cpp

//=================================================================================================

/** Class to exctract moving quantiles (such as the median) from a signal in realtime. If the 
percentile runs over N samples, the filters takes O(log(N)) operations per sample. It achieves this
by using two heaps (a max-heap of the smaller-than-percentile values and a min-heap of the 
bigger-than-percentile values) and a circular buffer in a clever way....

*/

template<class T>
class rsMovingQuantileFilter
{

public:


  rsMovingQuantileFilter(int numSmaller = 20, int numLarger = 20) 
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
    //  -> this will rebalance the heap in which the oldest sample resided
    // -if the largest (front) sample in the small heap is > than the smalles (front) sample in the
    //  large heap: exchange them and rebalance both heaps
    // -the swaps that occurr during rebalancing the heaps should also lead to corresponding swaps 
    //  in our delay buffer, due to the way we have implemented our swap function
    // -return the front sample of the "large" heap in case of an odd size and the average of both
    //  front samples for even sizes (size == getLength())

    //int L 

    int ni = buf.getOldest();      // node index of oldest node in nodes array
    int hi = nodes[ni].heapIndex;  // heap index of oldest node
    //int bi = nodes[ni].bufIndex;   // buffer index of oldest node - do we need this?


    nodes[ni].value = x;

    if(hi < nS)      // node to be replaced is in small heap
    {
      hi = small.floatIntoPlace(hi);   // it should float down to the bottom but doesn't

      int dummy = 0;
    }
    else
    {
      // node to be replaced is in large heap
      hi -= nS; // we need to offset the heap-index
      hi  = large.floatIntoPlace(hi);
      hi += nS;

      int dummy = 0;
    }

    // swap first elements of small/large heaps, if necessarry

    if(heapsNeedExchange())
    {
      rsSwap(nodes[0].heapIndex, nodes[nS].heapIndex);

      if(hi < nS)
      {
        // replacement took place in small heap, so after exchange, the new datum is in the large
        // heap

        small.floatIntoPlace(0);
        hi  = large.floatIntoPlace(0);
        hi += nS;
      }
      else
      {
        // replacement took place in large heap, so after exchange, the new datum is in the small
        // heap
        large.floatIntoPlace(0);
        hi = small.floatIntoPlace(0);
      }


      int dummy = 0;
    }

    // debug:
    bool isSmallHeap = small.isHeap();
    bool isLargeHeap = large.isHeap();
    rsAssert(isSmallHeap && isLargeHeap);



    nodes[ni].heapIndex = hi;  // is this correct?

    ni = buf.getSample( (ni+1) % getLength() );  // ..and this?




    int idx = large.getElement(0);  // index to read from the heap

    idx = nS;  // test
    T val = nodes[idx].value;


    return val;

    //return large.getElement(0);
    // preliminary - return the first (i.e. smallest) value of the large values
  }


  void reset()
  {
    for(int i = 0; i < getLength(); i++)
    {
      heaps[i] = i;
      buf[i] = i;
      nodes[i].heapIndex = i;
      nodes[i].bufIndex = i;
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

  bool nodeLess(const int& left, const int& right)
  {
    return nodes[left].value < nodes[right].value;
  }
  // is htis right? do we expect two indices into the nodes array?...yeah, i think so

  bool nodeGreater(const int& left, const int& right)
  {
    return nodes[left].value > nodes[right].value;
  }

  void swapNodes(int& i, int& j)
  {
    // ...something to do...
    //
    /*
    rsSwap(nodes[i].heapIndex, nodes[j].heapIndex);
    rsSwap(nodes[i].bufIndex,  nodes[j].bufIndex);
    rsSwap(nodes[i].value,     nodes[j].value);
    // is that correct? ..and/or complete...hmm...nooo
    */

    //rsSwap(nodes[i], nodes[j]);

    //rsSwap(buf[nodes[i].bufIndex], buf[nodes[j].bufIndex]);

    // i think, the nodes array stays fixed and we just swap the pointers to the nodes inside the 
    // heap and within the nodes array, the heap-indices

    rsSwap(nodes[i].heapIndex, nodes[j].heapIndex);
    rsSwap(i, j);

    // -could tha be all that is needed?
    // -i think, using the 

    return; 
  }
  // expects two indices into the nodes array



  void updateBufferLengths()
  {
    nodes.resize(getLength());
    buf.setLength((size_t)getLength());
    heaps.resize(getLength());
    small.setData(&heaps[0],  nS);
    large.setData(&heaps[nS], nL);
    reset();
  }

  struct Node
  {
    int heapIndex = 0;
    int bufIndex = 0;  // i think, this may be redundant - we'll see
    T value = T(0);
    // maybe we need a bool to indicate, if we are in the large heap - using the convention of
    // adding/subtracting nS is sometimes inconvenient
  };


  std::vector<Node> nodes;


  rsRingBuffer<int> buf;   // circular buffer
  std::vector<int> heaps;  // memory for the heaps
  // heaps and buf just store indices into the nodes array and the nodes-array itself contains the
  // actual data together with its current index inside the heap and delay-buffer. We need this 
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
