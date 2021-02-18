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

/*
moved to rapt:
static constexpr int allBits = -1;                                      // all bits are 1
static constexpr int allBitsButFirst = std::numeric_limits<int>::max(); // only 1st bit is 1
static constexpr int firstBitOnly = allBits ^ allBitsButFirst;          // only 1st bit is 0

// for unsiged int types, the bit twiddling is different:
//static size_t allBits = std::numeric_limits<size_t>::max();
//static size_t firstBitOnly = allBits - (allBits >> 1);
//static size_t allBitsButFirst= allBits ^ firstBitOnly;
*/

template<class T>
void weightedSum(const T* x1, int N1, T w1, const T* x2, int N2, T w2, T* y, int Ny)
{
  int n = 0;

  while(n < rsMin(N1, N2, Ny)) {
    y[n] = w1 * x1[n] + w2 * x2[n];
    n++; }

  if(n == Ny)
    return;

  // handle overhanging part in x1 or x2:
  if(N1 > N2) {
    while(n < rsMin(N1, Ny)) {
      y[n] = w1 * x1[n];
      n++; }}
  else if(N2 > N1) {
    while(n < rsMin(N2, Ny)) {
      y[n] = w2 * x2[n];
      n++; }}

  // handle trailing zeros in result, if necessarry:
  while(n < Ny) {
    y[n] = T(0);
    n++; }
}
// todo: document, write test, move into rsArrayTools

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


void rsCircularShift(int* a, int N, int k);
// circular shift without additional memory (using 3 reversals) - needs test


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

//=================================================================================================

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
// -maybe implement also a MinMaxHeap https://en.wikipedia.org/wiki/Min-max_heap

/** Works similar to rsBinaryHeap just that the property satisfied by the nodes is different.
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

  using rsBinaryTree<T>::rsBinaryTree;  // inherit constructors


  bool isSearchTree(int i = 0) const
  {
    if(i >= this->size)
      return true;
    bool result = true;
    int l = this->left(i);
    int r = this->right(i);
    if(l < this->size) result &= !this->less(this->data[i], this->data[l]) && isSearchTree(l);
    if(r < this->size) result &= !this->less(this->data[r], this->data[i]) && isSearchTree(r);
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
    for(int i = this->size-1; i >= 0; i--)
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

    int l = this->left(i);
    if(l >= this->size) return;       // node is leaf
    int r = this->right(i);
    if(r >= this->size) {             // node has only left child - make sure, the data at the
      if(this->less(this->data[i], this->data[l]))      // parent i is not less than at left child l - if it is: swap
        this->swap(this->data[i], this->data[l]);
      return; }

    // ok - we have a node that has both children - we must figure out which of the 3 nodes i,l,r
    // is the middle element m and if m is not already i, then do a swap. we must also make sure
    // that data[l] <= data[r]. we have 3 possible swapping scenarios: swap(i,l), swap(i,r),
    // swap(l,r) ...plus, of course, the no-swap scenario

    //int m = i;                           // index of mid element, preliminary set to i
    //if(less(data[i], data[l])) m = l;
    //if(less(data[r], data[m])) m = r;

    // fix order of children:
    if(this->less(this->data[r], this->data[l]))
      this->swap(this->data[l], this->data[r]);

    // fix order of i,l:
    if(this->less(this->data[i], this->data[l]))
      this->swap(this->data[l], this->data[i]);

    // fix order of i,r:
    if(this->less(this->data[r], this->data[i]))
      this->swap(this->data[r], this->data[i]);

    // try to do this only with 2 comparisons and one swap

  }



  int replace(int i, const T& x)
  {
    // data[i] = x; return floatIntoPlace(i); // this is what the heap needs

    // we insert it (preliminarily) either at i or at the sibling of i:
    int p = this->parent(i);
    if(this->isLeft(i))
      if(this->less(this->data[p], x))
        i = this->right(p);   // right sibling
    else
      if(this->less(x, this->data[p]))
        i = this->left(p);    // left sibling
    this->data[i] = x;
    return floatIntoPlace(i);
  }


protected:

  int floatIntoPlace(int i) { return floatUp(floatDown(i)); }

  int floatUp(int i)
  {
    while(i > 0) {
      int p = this->parent(i);
      if(this->isLeft(i)) {
        if(this->less(this->data[p], this->data[i]))
          this->swap(this->data[i], this->data[p]);
        else
          return i; }
      else {
        if(this->less(this->data[i], this->data[p]))      // i and p are swapped compared to left nodes
          this->swap(this->data[i], this->data[p]);
        else
          return i;   }
      i = p; }
    return i;
  }

  int floatDown1(int i)
  {
    int l = this->left(i);
    if(l >= this->size) return i;
    if(this->less(this->data[i], this->data[l])) {
      this->swap(this->data[i], this->data[l]);
      return floatDown(l); }
    int r = this->right(i);
    if(r >= this->size) return i;
    if(this->less(this->data[r], this->data[i])) {
      this->swap(this->data[i], this->data[r]);
      return floatDown(r); }
    return i;
  }
  // recursive implementation - not working


  int floatDown2(int i)
  {
    while(i < this->size-1)
    {
      int l = this->left(i);
      if(l >= this->size) return i;              // node i is a leaf
      int m = i;                           // index of mid element, preliminary set to i
      if(less(this->data[i], this->data[l])) m = l;    // left child is bigger than i, so mid index goes to l
      int r = this->right(i);

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

      if(r < this->size && less(this->data[r], this->data[m])) m = r;

      if(m != i) {
        swap(this->data[m], this->data[i]);
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

/** Extends rsQuantileFilter by a second core allowing to do more interesting stuff such as forming
linear combinations of lower an upper quantiles (such as min and max), etc. Using a second core
instead of just using two rsQuantileFilter objects is more economical because the delayline can be
shared between the cores. */

template<class T>
class rsDualQuantileFilter : public rsQuantileFilter<T>
{

public:

  rsDualQuantileFilter()
  {
    allocateResources();
    this->core.setDelayBuffer(&this->delayLine);
    this->core2.setDelayBuffer(&this->delayLine);
    this->dirty = true;
  }


  // setters for the parameters of the second core
  void setFrequency2(   T newFrequency) { setLength2(T(1) / newFrequency); }
  void setLength2(      T newLength)    { length2   = newLength;    this->dirty = true; }
  void setQuantile2(    T newQuantile1) { quantile2 = newQuantile1; this->dirty = true; }
  void setLowpassGain2( T newGain)      { loGain2   = newGain; }
  void setHighpassGain2(T newGain)      { hiGain2   = newGain; }

  void setCore2Equal()
  {
    length2   = this->length;
    quantile2 = this->quantile;
    loGain2   = this->loGain;
    hiGain2   = this->hiGain;
    delayScl2 = this->delayScl;
    this->dirty = true;
  }

  void setCore2Complementary()
  {
    setCore2Equal();
    quantile2 = T(1) - this->quantile;
    this->dirty = true;
  }


  // maybe have a setCore2Complementary which uses the same settings for core2 as core2 except the
  // quantile, for which we use: quantile2 = 1-quantile

  T getSample(T x)
  {
    this->delayLine.getSample(x);
    if(this->dirty)
      updateInternals();
    T yL1 = this->core.getSample(x);
    T yH1 = this->delayLine[this->delayScl*this->delay] - yL1;
    T yL2 = core2.getSample(x);
    T yH2 = this->delayLine[delayScl2*delay2] - yL2;
    return this->loGain * yL1 + this->hiGain * yH1 + loGain2 * yL2 + hiGain2 * yH2;
  }

  void reset() { this->core.reset(); core2.reset(); this->delayLine.reset(); }

  virtual void updateInternals() override
  {
    // rsError("has to be updated");
    // 
    /*
    // compute internal and set up core parameters:
    int L, p; T w;

    this->convertParameters(this->length, this->quantile, this->sampleRate, &L, &p, &w, &this->delay);
    this->core.setLengthAndReadPosition(L, p);
    this->core.setRightWeight(w);

    this->convertParameters(length2, quantile2, this->sampleRate, &L, &p, &w, &delay2);
    core2.setLengthAndReadPosition(L, p);
    core2.setRightWeight(w);
    */

    // update inherited baseclass members:
    T L = this->length * this->sampleRate;
    this->core.setLengthAndQuantile(L, this->quantile);
    this->delay = T(0.5)*(L-1);
    //this->delay = T(1.0) * (L-1) * this->quantile;  // experimental

    // update new subclass members:
    L = length2 * this->sampleRate;
    core2.setLengthAndQuantile(L, quantile2);
    delay2 = T(0.5)*(L-1);
    //delay2 = T(1.0) * (L-1) * quantile2;   // experimental

    this->dirty = false;
  }


protected:

  virtual void allocateResources() override
  {
    int mL = (int) ceil(this->maxLength * this->sampleRate);
    this->core.setMaxLength(mL);
    core2.setMaxLength(mL);
    this->delayLine.setCapacity(mL);
  }


  // the 2nd core and its set of parameters:
  rsQuantileFilterCore2<T> core2;
  T length2   = 0.0;
  T quantile2 = 0.5;
  T loGain2   = 0.0;
  T hiGain2   = 0.0;
  T delayScl2 = 1.0;
  T delay2    = 0.0;

};

//=================================================================================================

/** Extends rsQuantileFilter by producing a pseudo-resonance signal formed by... */

template<class T>
class rsQuantileFilterResonant : public rsQuantileFilter<T>
{

public:

  using Base = rsQuantileFilter<T>;

  rsQuantileFilterResonant()
  {
    //resoHighpass.setMode(rsOnePoleFilter<T,T>::HIGHPASS_MZT);
    resoHighpass.setMode(rsStateVariableFilter<T,T>::HIGHPASS);
    allocateResources();
    Base::dirty = true;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample rate at which this filter should operate. May re-allocate memory. */
  void setSampleRate(T newSampleRate)
  {
    sampleRate = newSampleRate;
    allocateResources();
    resoHighpass.setSampleRate(sampleRate);
    Base::dirty = true;
  }

  /** Sets the maximum length (in seconds) for this filter. May re-allocate memory. */
  void setMaxLength(T newMaxLength)
  {
    maxLength = newMaxLength;
    allocateResources();
    Base::dirty = true;
  }

  /** Sets sample rate and maximum length at the same time. May re-allocate memory. This may avoid
  some re-allocations compared to using setSampleRate followed by setMaxLength (or vice versa), 
  depending on the situation - so, if possible, it's recommended to set both at the same time. */
  void setSampleRateAndMaxLength(T newSampleRate, T newMaxLength)
  {
    sampleRate = newSampleRate;
    maxLength  = newMaxLength;
    allocateResources();
    Base::dirty = true;
  }

  void setResonanceMix(T newMix)
  {
    resoMix = newMix;
  }
  // maybe use a resoGain instead

  /*
  enum class ResoMode
  {


  };
  */


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Produces one output sample from a given input sample. */
  T getSample(T x)
  {
    delayLine.getSample(x);
    if(dirty) 
      updateInternals();
    T yL = core.getSample(x);                   // lowpass part
    T yD = delayLine[delayScl*delay];           // delayed input
    T yH = yD - yL;                             // highpass part
    T yF = loGain * yL + hiGain * yH;           // non-resonant filtered output
    T yR = getResonance(x, yL, yD);             // (pseudo) resonance
    return (T(1)-resoMix) * yF + resoMix * yR;  // crossfade between non-resonant and resonance
  }
  // maybe factor out a function to produce lowpass and highpass getSampleLoHi or something at the
  // same time - client code may find that useful - or maybe getOutputs to be consistent with
  // rsStateVariableFilter

  /** Resets the filter into its initial state. */
  void reset()
  {
    Base:reset();
    minCore.reset();
    maxCore.reset();
    bandpass.reset();
    resoHighpass.reset();
  }

  /** Updates the internal algorithm parameters and embedded objects according to the user
  parameters. This is called in getSample, if the state is dirty but sometimes it may be
  convenient to call it from client code, too. */
  virtual void updateInternals()
  {
    double L = length*sampleRate;
    core.setLengthAndQuantile(   L, quantile);

    minCore.setLengthAndQuantile(L, T(0));
    maxCore.setLengthAndQuantile(L, T(1));
    //minCore.setLengthAndQuantile(100, T(0));
    //maxCore.setLengthAndQuantile(100, T(1));
    //minCore.setLengthAndQuantile(rsMax(L, 10.0), T(0));
    //maxCore.setLengthAndQuantile(rsMax(L, 10.0), T(1));


    delay = T(0.5)*(L-1);

    T f = getFrequency();

    //T w = T(2*PI)*sampleRate/length; // == 2*PI*sampleRate*frequency
    T w = T(2*PI) * f  / sampleRate; // == 2*PI*frequency/sampleRate
    bandpass.setFrequencyAndAbsoluteBandwidth(w, T(0.00001));  // preliminary
    // todo: let the resonance frequency have an adjustable pitch-offset with respect to core 
    // filters...maybe the min/max cores could have an offset as well
    // and maybe we should use a power rule to generalize from constant absolute bandwidth (0) and
    // constant relative bandwidth (1) to anything between and beyond...and we generally need a 
    // reasonable user-parameter to control the bandwidth

    //T gain = bandpass.getMagnitudeAt(w);  // excluding the output gain
    //bandpass.setOutputGain(T(1)/gain);
    // for test/debug 
    // -gain tends to get huge for small bandwidths (seems to be in reciprocal 
    //  relationship)! 
    // -todo: plot impulse response

    //T fHp = rsMax(T(0), 0.9*f - T(1000));
    T fHp = rsMax(T(0), 0.5*f);
    //resoHighpass.setCutoff(fHp);
    resoHighpass.setFrequency(fHp);


    dirty = false;
  }



protected:

  /** Computes the pseudo-resonance signal. Inputs are the filter input x, the "lowpass" quantile
  filter output yL and the delayed input yD (rename to xD) which are used to generate the 
  resonance. */
  T getResonance(T x, T yL, T yD)
  {
    T min  = minCore.getSample(x);
    T max  = maxCore.getSample(x);

    // this computation should be factored out and we should provide various different ways to
    // compute wMin/wMax that can bew switched by the user using a setResonanceMode function:
    T wMin, wMax;
    wMin = 0;                              // preliminary.. maybe these should be members and 
    wMax = 0;                              // computed in updateinternals

    //if(yD >= yL) wMax = 1;
    //else         wMin = 1;

    // this is not a very good way to compute wMin/wMax - it should use some condition that is 
    // likely to cause a switch once per cycle, i.e. on the avarage switches once within a length
    // of the filter buffer. this condition here switches far too often - try several things
    // ...oh - but maybe that's just because the input is white noise - if it's somewhat lowpassed
    // noise, i may be better

    computeMinMaxWeights(&wMin, &wMax, x, yL, yD, min, max);
    T yR = wMin*min + wMax*max;

    // experimental - highpass and amplify the resonance to counteract resonance loss at high 
    // cutoff frequencies:
    yR  = resoHighpass.getSample(yR);
    T a = T(1);    // preliminary - should later depend on cutoff in some way

    // experimental:
    T p = 2*getFrequency()/sampleRate;  // 0..1 as f goes 0..fs
    //p = 1-p; p = p*p; a = 1 / (0.1 + p);
    //p = cos(0.5*PI*p); a = 1 / (0.1 + p);
    p = 1-sin(0.5*PI*p); a = 1 / (0.1 + p);
    // todo: figure out the best shape


    return a*yR;
  }

  void computeMinMaxWeights(T* wMin, T* wMax, T x, T yL, T yD, T yMin, T yMax)
  {
    T yB = bandpass.getSample(yD);  // or maybe feed x

    T thresh = T(0.0) * yL;
    //T thresh = T(0.2) * rsMax(rsAbs(yMin), rsAbs(yMax));
    // make user-adjustable, determines pulse-width, should perhaps be scaled by yL...and maybe we
    // should normalize the amplitude of the bandpass with respect to the bandwidth (we probably 
    // want a constant peak gain bandpass - maybe the RBJ bandpass is good for that)
    // ...it should take into account the magnitude of the bandpass output which may be a function
    // of the bandwidth - we shoould normalize everything such that thresholds between -1..+1 make
    // sense and have the same effect regardless of input volume and bandwidth setting

    T a = T(1);  // scaler to control the resonance amplitude
    T f = getFrequency();
    //a = 10 / length;
    //a = T(1) + T(4)*rsClip(f / T(20000), T(0), T(1)); 
    // yes, we need a factor, but we also need a highpass, otherwise we just amplify noise
    // maybe the highpass freq should be something like max(0, cutoff-2000) such that the highpass
    // sets in above 2kHz cutoff an is tuned 2kHz below the cutoff...of maybe use 
    // max(0, 0.5*cutoff-2000) or something - tweak the factor and offset to taste

    if(yB >= thresh) { *wMin = T(0); *wMax = a;    }
    else             { *wMin = a;    *wMax = T(0); }
    // maybe instead of hard-switching, we can make a sort of soft-switch?


    /*
    if(yD >= yL) { *wMin = T(0); *wMax = T(1); }
    else         { *wMin = T(1); *wMax = T(0); }
    */
  }

  virtual void allocateResources() override
  {
    Base::allocateResources();
    int mL = getMaxRequiredLengthInSamples();
    minCore.setMaxLength(mL);
    maxCore.setMaxLength(mL);
  }

  // filter cores to extract min and max:
  rsQuantileFilterCore2<T> minCore, maxCore;
  // can this be done more efficiently? i could perhaps use the rsMovingMinMaxFilter class but 
  // then the length would not be modulatable and i don't think it's easy to make it modulatable

  T resoMix = T(0);

  rsTwoPoleFilter<T, T> bandpass;
  // the two-pole is actually not really a bandpass but a resonator - maybe an RBJ bandpass is 
  // better?

  //rsOnePoleFilter<T, T> resoHighpass;

  rsStateVariableFilter<T, T> resoHighpass;

};



//=================================================================================================

/** This is a naive implementation of (the core of) a moving quantile filter and meant only for
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
    buf.copyTo(&tmp[0], true);          // O(N)
    rsHeapSort(&tmp[0], getLength());   // O(N*log(N))
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

/** This is a naive implementation of a "filter" that minimizes the distance between subsequent 
samples by re-ordering them. It always keeps a buffer of N samples from the past and when a new 
sample comes in, it replaces one of the buffered samples with the new one and the replaced sample 
from the buffer is used as output. Which one is replaced and used as output is one that is closest 
to the previous output sample. In the special case that the current input is closer to the previous
output than all buffered samples, it will be used as output and the buffer will left alone. */

template<class T>
class rsDistanceFilterNaive
{

public:

  void setLength(int newLength)
  {
    buf.resize(newLength);
    rsFill(buf, T(0));
  }
  // O(N)

  T getSample(T x)
  {
    T   dMin = rsAbs(y - x);
    int iMin = -1;
    for(int i = 0; i < (int) buf.size(); i++) {
      T d = rsAbs(y - buf[i]);
      if(d < dMin) {
        dMin = d;
        iMin = i; }}

    if(iMin != -1) {
      y = buf[iMin];
      buf[iMin] = x; }
    else
      y = x;

    return y;
  }
  // O(N)

  void reset()
  {
    rsFill(buf, T(0));
    y = T(0);
  }
  // O(N)

protected:

  std::vector<T> buf;     // buffer of stored samples
  T y = T(0);             // old output

};
// ..not yet tested
// todo:
// -implement this more efficiently using a heap: 
//  -we need a heap-search procedure which should run in log(N) time): int heap.findClosest(T x)
//  -this should be used to determine the sample to be replaced, replacement itself will also take
//   O(log(N)), so the overall complexity of getSample will be O(log(N))
//  ...oh - no, maybe a heap will not work and we need a binary search tree instead
// -use it to filter uniform white noise - it should preserve the uniform amplitude distribtion
//  while still imposing a correlation (that was actually the question that lead to the idea: how 
//  can we produce correlated noise with uniform amplitude distribution - because regular filters
//  tend to gaussianize it)
// -the naive version could generalize to other distance measures such as 2D Euclidean distance
//  for stereo signal - unfortunatley, i don't see, how this could be made efficient using a heap 
//  because there's no natural ordering for 2D vectors...but maybe one could be invented...but it 
//  should be related to distance between pairs of vectors

// ...how else could we impose correlations (and therefore, some sort of not-flat frequency 
// spectrum) without modifying the amplitude distribution? the amplitude distribution can be 
// modified by waveshaping and correlations can be introduced by filtering - so maybe a filtering
// process (like 2 or 3 point MA) can be used to spectrally shape the noise, which also turns a 
// uniform distribtuion into triangular or parabolic and then that can be followed by waveshaper 
// that counteracts this change (square? cube? sqrt?, cbrt?) ...or maybe it should be the 
// (inverse of?) the integral of the resulting distribution

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
// continued fraction stuff:

/** A class for generating the (integer) continued fraction expansion coefficients of a given
(floating point) number. You pass the number to the constructor and after that, you can pull out
successive cofficients via getNext(). */

template<class TInt, class TFloat>
class rsContinuedFractionGenerator
{

public:

  rsContinuedFractionGenerator(TFloat number) : state(number) {}

  TInt getNext()
  {
    TFloat flr = floor(state);
    state = TFloat(1) / (state - flr);
    return (TInt) flr;
  }
  // maybe rename to getNextCoeff ...or whatever these numbers are called


protected:

  TFloat state;

};
// https://www.youtube.com/watch?v=CaasbfdJdJg


template<class T>
rsFraction<T> rsContinuedFractionConvergent(T* a, int N)
{
  T p0 = 0, p1 = 1, p2 = 0;
  T q0 = 1, q1 = 0, q2 = 1;
  for(int i = 0; i < N; ++i) {
    p2 = a[i]*p1 + p0; p0 = p1; p1 = p2;
    q2 = a[i]*q1 + q0; q0 = q1; q1 = q2; }
  return rsFraction<T>(p2, q2);
}
// algorithm adapted from cfcv.c (by Hollos) - i don't really know, why it works
// this can actually be done directly using the generator, without the need for explicitly
// computing and storing the array a
// maybe move into class rsContinuedFractionGenerator...maybe as static method
// ..but i don't think that this continued fraction stuff should go into rsFraction - it's stuff
// on top of it
// -note the the convegents are not eqaul to the best approximants. there's some additional stuff
//  that needs to be done - maybe implement that in a function rsRationalApproximant:
//  https://en.wikipedia.org/wiki/Continued_fraction#Best_rational_approximations
//
// can we somehow figure out, how many of the CFE coeffs are correct without knowing the correct
// CFE? maybe by converting the convergents back to double and only add more coeffs as long as
// the back-converted number actually gets closer to the original number?
// -maybe rename to rsContinuedToRegularFraction - it just converts an array of (simple) continued
//  fraction coeffs to a normal, regular fraction - it doesn't really matter whether or not the
//  coefficient array is supposed to converge to some irrational number or if it's just another
//  representation of a fraction
// -maybe also implement a generalization that does not assume a *simple* continued fraction (i.e.
//  one where all numerators are 1 and there are no minusses)
// https://en.wikipedia.org/wiki/Continued_fraction
// https://en.wikipedia.org/wiki/Generalized_continued_fraction

// returns the (simple) continued fraction coeffs of the given rational number
template<class T>
std::vector<T> rsContinuedFraction(rsFraction<T> x)
{
  std::vector<T> c;
  T p = x.getNumerator();
  T q = x.getDenominator();
  T a = p/q, b;
  c.push_back(a);
  while(p > q*a) {
    b = p - q*a; p = q; q = b; a = p/q;
    c.push_back(a); }
  return c;
}
// algo adapted from from cfrat.c
// i think, this is some variation of the Euclidean algorithm

// Wikipedia says: "Even-numbered convergents are smaller than the original number, while
// odd-numbered ones are larger."
// https://en.wikipedia.org/wiki/Continued_fraction#Infinite_continued_fractions_and_convergents
// ...verify this. Knowing this could be useful

// other fun stuff that can be done with fractions:
// https://en.wikipedia.org/wiki/Egyptian_fraction


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

/** Sawtooth oscillator based on a PolyBlep. Convenience class that also demonstrates, how the 
bleps are supposed to be used in a quite minimal setting. */

template<class T>
class rsBlepSawOsc : protected rsBlepReadyOsc<T>
{

public:

  using Base = rsBlepReadyOsc<T>;

  void setPhaseIncrement(T newInc) { Base::setPhaseIncrement(newInc); }

  T getSample()
  {
    T stepDelay, stepAmp;
    T y = Base::getSampleSaw(&stepDelay, &stepAmp); // todo: maybe switch waveform later
    if(stepAmp != 0.0)
      blep.prepareForStep(stepDelay, stepAmp);
    return blep.getSample(y);
  }

  void reset(T start = T(0))
  {
    Base::resetPhase(start);
    blep.reset();
  }

protected:

  rsPolyBlep2<T, T> blep;

};

//=================================================================================================

/** A class for representing and performing computations with sparse matrices which are matrices in 
which most elements are zero. This implementation uses a std::vector of "elements" where each 
element stores its row- and column-indices and the actual value. */

template<class T>
class rsSparseMatrix
{


public:

  // todo: provide a (subset of) the same interface as rsMatrix(View)...maybe it's a good idea to 
  // keep the interface small to make it easier to provide different implementations with different
  // storage modes with the same interface that can be benchmarked against each other.

  rsSparseMatrix() {}

  rsSparseMatrix(int numRows, int numColumns)
  {
    rsAssert(numRows >= 1 && numColumns >= 1);
    this->numRows = numRows;
    this->numCols = numColumns;
  }

  //-----------------------------------------------------------------------------------------------

  /** Returns the number of nonzero elements in this matrix. */
  int getNumElements() const { return (int) elements.size(); }

  int getNumRows() const { return numRows; }

  int getNumColumns() const { return numCols; }



  bool isValidIndexPair(int i, int j) const 
  { return i >= 0 && i < numRows && j >= 0 && j < numCols; }


  /** Returns the diagonal part D of this matrix A, such that A = D + N (where N is the 
  non-diagonal part). Diagonal part means that only the diagonal elements are included, the others
  are set to zero, i.e. left out. */
  rsSparseMatrix<T> getDiagonalPart() const
  {
    rsSparseMatrix<T> D(numRows, numCols);
    for(size_t k = 0; k < elements.size(); k++)
      if(elements[k].i == elements[k].j)
        D.set(elements[k].i, elements[k].j, elements[k].value);
    return D;
  }
  // maybe this should return a vector instead of a sparse matrix, maybe make a function that 
  // splits A into diagonal and non-diagonal at once: splitOutDiagonalPart(Vec& D, Mat& N)

  /** Returns the diagonal part N of this matrix A, such that A = D + N (where D is the 
  diagonal part). */
  rsSparseMatrix<T> getNonDiagonalPart() const
  {
    rsSparseMatrix<T> N(numRows, numCols);
    for(size_t k = 0; k < elements.size(); k++)
      if(elements[k].i != elements[k].j)
        N.set(elements[k].i, elements[k].j, elements[k].value);
    return N;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Accessors. Element access via these is slow, so they should probably be only used, when 
  a matrix is built once and for all as a precomputation. When the matrix is used later e.g. in an
  iterative linear solver, you will probably want to use more efficient functions like product. */

  /** Read access. */
  T operator()(int i, int j) 
  { 
    Element e(i, j, T(0));
    if(elements.empty()) 
      return T(0);
    size_t k = (size_t) rsArrayTools::findSplitIndex(&elements[0], getNumElements(), e);
    if(k >= elements.size() || e < elements[k])
      return T(0);
    else
      return elements[k].value;
  }

  /** Sets the element at position (i,j) to the given new value. This may lead to insertion of a 
  new element (if there's no element yet at i,j) or removal of existing elements (if val is 
  zero). */
  void set(int i, int j, T val) 
  { 
    rsAssert(isValidIndexPair(i, j), "Index out of range");
    Element e(i, j, T(val));
    if(elements.empty() && val != T(0)) {
      elements.push_back(e);
      return;  }
    size_t k = (size_t) rsArrayTools::findSplitIndex(&elements[0], getNumElements(), e);
    if((k >= elements.size() || e < elements[k]) && val != T(0))
      rsInsert(elements, e, k);
    else {
      if(val != T(0))
        elements[k] = e;
      else
        rsRemove(elements, k); }
  }
  // todo: element removal needs tests


  //-----------------------------------------------------------------------------------------------

  /** Computes the matrix-vector product y = A*x where x must be numCols long and y must be numRows 
  long. The complexity is O(N+K) where N is the number of rows and K is the number of nonzero 
  elements. */
  void product(const T* x, T* y) const
  {
    rsAssert(x != y, "Can't be used in place");
    for(int j = 0; j < numRows; j++)
      y[j] = T(0);
    for(size_t k = 0; k < elements.size(); k++)
      y[elements[k].i] += elements[k].value * x[elements[k].j];
  }
  // -maybe include optional strideX, strideY parameters - or maybe implement a separate function with
  //  strides
  // -how can we implement the product with the transposed matrix? would it be
  //    y[elements[k].j] += elements[k].value * x[elements[k].i];


  /** Computes the matrix-vector product y = A*x and returns the maximum absolute value of the 
  change in y before and after the call. This is supposed to be used inside vector iterations of
  the form yNew = A * yOld, where y takes the role of yNew and x the role of yOld. This implies 
  that x and y must have the same dimensionality which in turn means the matrix A must be square. 
  The return value can be used in a convergence test. If x and y are the same, i.e. the function is
  applied in place, the returned y vector will actually not be the matrix-vector product A*x 
  because in the computations of the values with higher index, it will use already updated values 
  with lower index. That may seem undesirable from a mathematical cleanliness point of view, but in
  practice, that may actually speed up convergence (see Jacobi- vs Gauss-Seidel iteration - it's a 
  similar situation here). So, it can be used in place in such a context - just don't expect the 
  computed y to be the exact matrix-vector product then. (todo: test, if it really improves 
  convergence - i just assume it because of the similarity to Gauss-Seidel) */
  T iterateProduct(const T* x, T* y) const;
  // may not be needed


  // todo: maybe implement matrix-operators +,-,*. Such operations should result in another 
  // (possibly less) sparse matrix. A naive algorithm for addition can use element access and set. 
  // I think, this will have complexity O(N*M*(log(K1)+log(K2))) where N,M are number of rows and 
  // columns and K1, K2 are the numbers of nozero elements in operands 1 and 2. I think, it's 
  // possible to do in O(K1 + K2) or maybe O(K1*log(K1) + K2*log(K2))...or something. Maybe the 
  // naive algo should be part of the test-suite but not the class itself.


  // maybe move into class rsSparseLinearAlgebra to not overload this class

  /** Given the diagonal part D and non-diagonal part N of a matrix A, such that A = D + N, this 
  function solves the linear system A*x = (D+N)*x = b via Gauss-Seidel iteration and returns the 
  number of iterations taken. The resulting solution vector will be written into x. Whatever is 
  stored in the vector x before the call will be taken as initial guess. For the iteration to 
  converge, the matrix A must be strictly diagonally dominant. That means the diagonal element on 
  each row must have larger absolute value than the sum of the absolute values of the off-diagonal
  elements in the same row. 
  https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
  https://en.wikipedia.org/wiki/Diagonally_dominant_matrix   */
  static int solveGaussSeidel(const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& N, T* x, 
    const T* b, T tol);

  /** Convenience function that takes the matrix A and splits it internally into diagonal and 
  non-diagonal parts. This is slow, so it should be used only in testing. In production code, the 
  splitting can typically be done once and for all as pre-processing step. */
  static int solveGaussSeidel(const rsSparseMatrix<T>& A, T* x, const T* b, T tol)
  { return solveGaussSeidel(A.getDiagonalPart(), A.getNonDiagonalPart(), x, b, tol); }
  // maybe make it non-static, to be called like A.solveGaussSeidel


  /**
  https://en.wikipedia.org/wiki/Successive_over-relaxation
  needs a workspace of size N, becomes Gauss-Seidel for w = 1 and Jacobi for w = 0.  */
  static int solveSOR(const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& N, T* x, 
    const T* b, T tol, T* workspace, T w);

  static int solveSOR(const rsSparseMatrix<T>& A, T* x, const T* b, T tol, T* workspace, T w)
  { return solveSOR(A.getDiagonalPart(), A.getNonDiagonalPart(), x, b, tol, workspace, w); }


protected:

  /** Given a matrix position with row- and column-indices i,j, this function either returns the 
  flat index k at which the element is found in our elements array or, if it's not found, the index
  at which that element should be inserted. */
  /*
  int flatIndex(int i, int j)
  {

    return 0;  // preliminary
    // todo: do a binary search for element at position i,j, return its the flat index k if 
    // present, otherwise return the place where this element should be inserted
  }
  */

  struct Element
  {
    int i, j;   // row and column index
    T value;

    Element(int row, int col, T val) { i = row; j = col; value = val; }

    /** The less-than operator compares indices, i.e. a < b, iff a is supposed to be stored before 
    b in our elements array. The actual stored value plays no role in this comparison. This may 
    seem weird but it is just what is needed for the binary search (which uses the operator) that 
    is used in element access. ...todo: try to find a more elegant solution. */
    bool operator<(const Element& b) const
    {
      if(i   < b.i) return true;
      if(b.i <   i) return false;
      if(j   < b.j) return true;   // i == b.i, compare with respect to j
      return false;
    }
  };


  int numRows = 0, numCols = 0;  // number of rows and columns - maybe try to get rid
  
  std::vector<Element> elements;


};

// ToDo: 
// -Make another implementation (with the same interface) that stores rows. This saves one 
//  integer of storage space per entry because the row index is given implicitly. Maybe make a 
//  column-wise version, too - but that's less useful because with row-wise storage, it's more 
//  convenient and efficient to execute matrix-vector multiplications which is the most important 
//  operation in iterative linear solvers, which are the main application of sparse matrices.
// -Maybe templatize also on the integer type used for indices i,j. Users may want to use short
//  integers (16 bit) to save even more storage space, especially when T is float because then 
//  one entry is exactly 64 bits long). Maybe use TIdx, TVal for the two template parameters.
// -make a class rsSparseTensor

template<class T>
T rsSparseMatrix<T>::iterateProduct(const T* x, T* y) const
{
  rsAssert(numRows == numCols, "Can be used only for square matrices");
  size_t k = 0;
  T dMax = T(0);
  while(k < elements.size())
  {
    size_t i = elements[k].i;
    T yi = T(0);
    while(k < elements.size() && elements[k].i == i) {
      yi += elements[k].value * x[elements[k].j];
      k++; }
    T dyi = y[i] - yi;
    dMax = rsMax(dMax, rsAbs(dyi));
    y[i] = yi; 
  }
  return dMax;
}
// this function may not be that useful after all, maybe remove it (or put it into some sort of 
// code-attic)

template<class T>
int rsSparseMatrix<T>::solveGaussSeidel(
  const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& C, T* x, const T* b, T tol)
{
  size_t N = (size_t) D.numRows;
  rsAssert(D.numCols == N); // the matrices D,C must be square and have the same shape
  rsAssert(C.numRows == N);
  rsAssert(C.numCols == N);
  int numIts = 0;
  while(true) {

    // Perform one Gauss-Seidel step and record the maximum change in the value in any of the
    // elements in vector x:
    size_t k = 0;
    T dMax = T(0);
    for(size_t i = 0; i < N; i++) {
      T xi = b[i];  // xi will become the new, updated value for x[i]
      while(k < C.elements.size() && C.elements[k].i == i)  {
        xi -= C.elements[k].value * x[C.elements[k].j];
        k++;
      }
      xi /= D.elements[i].value;
      T dxi = x[i] - xi;
      dMax = rsMax(dMax, rsAbs(dxi));
      x[i] = xi; 
    }

    // Increment iteration counter and check convergence criterion:
    numIts++;
    if(dMax <= tol)
      break;
  }

  return numIts;
}
// todo: can the internal step be expressed using C.product(x, x) - this would be nice for 
// generalizing the algo to other implementations of sparse matrices...but it may be hard to keep
// track of the dMax

template<class T>
int rsSparseMatrix<T>::solveSOR(const rsSparseMatrix<T>& D, const rsSparseMatrix<T>& C, T* x,
  const T* b, T tol, T* wrk, T w)
{
  size_t N = (size_t) D.numRows;
  int numIts = 0;
  while(true)
  {
    // Perform one iteration of SOR and record the maximum change in the value in any of the
    // elements in vector x:
    for(size_t i = 0; i < N; i++)
      wrk[i] = x[i];

    size_t k = 0;
    T dMax = T(0);
    for(size_t i = 0; i < N; i++)
    {
      T xi = b[i];  // xi will become the new, updated value for x[i]
      while(k < C.elements.size() && C.elements[k].i == i) 
      {
        T xOld = wrk[C.elements[k].j];
        T xNew = x[  C.elements[k].j];
        xi -= C.elements[k].value * (w * xNew + (T(1)-w) * xOld);
        k++;
      }
      xi /= D.elements[i].value;
      T dxi = x[i] - xi;
      dMax = rsMax(dMax, rsAbs(dxi));
      x[i] = xi;
    }

    // Increment iteration counter and check convergence criterion:
    numIts++;
    if(dMax <= tol)
      break;
  }

  return numIts;
}

//=================================================================================================

/** A class that implements iterative algorithms for numerical linear algebra. It is written in
such a way that the same code can be used with dense or sparse matrices or any other kind of 
special matrix, as long as some appropriate low-level functions are added to this class for the
specific matrix-types. These low-level functions include things like retrieving the shape and for
computing matrix-vector or matrix-matrix products. These low-level special purpose implementations 
(one for each special matrix class) will then be used inside the actual computational algorithms 
which themselves will be the same for all the different matrix classes. If you want to use the 
algorithms for some new special matrix class, you will need to add just a small amount of 
boilerplate code for these operations - either here in the class or as global functions. */

class rsIterativeLinearAlgebra
{

public:

  template<class T, class TMat>
  static int largestEigenValueAndVector(const TMat& A, T* val, T* vec, T tol, T* workspace);
  // maybe rename to vonMisesIteration or eigenViaVonMises/eigenViaPowerIteration

  template<class T, class TMat>
  static int eigenspace(const TMat& A, T* vals, T* vecs, T tol, T* workspace);
  // each eigenvector is found in turn from the largest to the smallest via a variation of the von 
  // Mises iteration in which the projection of the iterates onto the already found eigenspace is
  // subtracted from the iterates
  // returns the maximum number of iterations, i.e. the number of iteration for the eigenvector 
  // that took the most iterations...hmm - or maybe it should return the sum of all iterations
  // -maybe it should take an additional parameter to specify, how many eigenpairs should be found
  //  and also a maximum number of iterations


  // Specializations of some low-level functions for sparse matrices (boilerplate):
  template<class T> static void product(const rsSparseMatrix<T>& A, const T* x, T* y) { A.product(x, y); }
  template<class T> static int numRows(const rsSparseMatrix<T>& A) { return A.getNumRows(); }
  template<class T> static int numColumns(const rsSparseMatrix<T>& A) { return A.getNumColumns(); }

  // todo: implement:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  // https://en.wikipedia.org/wiki/Biconjugate_gradient_method
  // https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

};

template<class T, class TMat>
int rsIterativeLinearAlgebra::largestEigenValueAndVector(
  const TMat& A, T* val, T* vec, T tol, T* wrk)
{
  rsAssert(numRows(A) == numColumns(A), "Can be used only for square matrices");
  using AT = rsArrayTools;
  int N = numRows(A);
  T L = AT::euclideanNorm(vec, N);
  AT::scale(vec, N, T(1) / L);
  int numIts = 0;
  while(true) {
    product(A, vec, wrk);
    L = AT::euclideanNorm(wrk, N);
    AT::scale(wrk, N, T(1) / L);
    T dMax = AT::maxDeviation(vec, wrk, N);
    if(dMax <= tol) {
      *val = L;                          // that's only the absolute value...
      int i = AT::maxAbsIndex(vec, N);
      if(vec[i] * wrk[i] < T(0))         // ...this figures out the sign
        *val = -(*val);
      AT::copy(wrk, vec, N);             // return the last iterate, not the 2nd to last
      break; }
    AT::copy(wrk, vec, N);
    numIts++; }
  return numIts;
}
// This implements the von Mises vector iteration.
// https://en.wikipedia.org/wiki/Power_iteration

// ToDo:
// -try using the maximum norm instead of the Euclidean - may be numerically more precise due to 
//  involving less arithmetic (none, actually)
// -maybe rename to MisesIteration
// -maybe include a maxNumIterations parameter
// -maybe factor out the stuff below product, such that the compiled version of this code can be 
//  re-used for dense matrices (to reduce binary size)..maybe into a normalizeAndTest function
//  that returns a boolean, the copy inside the if can be dragged outside the while(true) loop
// -try to come up with an in-place iteration like Gauss-Seidel - but it must take care to not
//  use the updated values as is but scale them by the current stimate of the eigenvalue, because 
//  the next iterate is expected to be scaled by the eigenvalue in the iteration - so the summation
//  loop need to be split inot two halfs: one with the scaling (using the updated values) and one 
//  without (using the old values)
// -try to find all eigenvalues and -vectors by subtracting the projection onto the already found
//  eigenspace from each iterate - this should be done right after product(vec, wrk); i think

template<class T, class TMat>
int rsIterativeLinearAlgebra::eigenspace(const TMat& A, T* vals, T* vecs, T tol, T* wrk)
{
  // Algorithm:
  // It works like the von Mises vector iteration for one eigenvector at the time, from large to 
  // small. To avoid converging to the same eigenvector again, instead of just forming the 
  // matrix-vector product of the matrix with the previous iterate, we subtract from that product 
  // the projection of it onto the space spanned by the already previously found eigenvectors.

  rsAssert(numRows(A) == numColumns(A), "Can be used only for square matrices");
  using AT = rsArrayTools;

  T (*norm)(const T*, int) = &AT::euclideanNorm;
  //T (*norm)(const T*, int) = &AT::maxAbs;
  // i think, any norm should work, but maxAbs apparently doesn't - why?

  int N = numRows(A);
  int numIts = 0;
  for(int n = 0; n < N; n++) {        // loop over the eigenvectors
    T* val = &vals[n];                // location of n-th eigenvalue
    T* vec = &vecs[n*N];              // start location of n-th eigenvector
    //T L = AT::euclideanNorm(vec, N);
    T L = norm(vec, N);
    AT::scale(vec, N, T(1) / L);
    while(true) {
      product(A, vec, wrk);
      for(int i = 0; i < n; i++) {
        T pi = T(0);
        for(int j = 0; j < N; j++) pi += wrk[j] * vecs[i*N + j];   // compute projection coeff
        for(int j = 0; j < N; j++) wrk[j] -= pi * vecs[i*N + j]; } // subtract projection
      //L = AT::euclideanNorm(wrk, N);
      L = norm(wrk, N);
      AT::scale(wrk, N, T(1) / L);
      T dMax = AT::maxDeviation(vec, wrk, N);
      if(dMax <= tol) {
        *val = L;
        int i = AT::maxAbsIndex(vec, N);
        if(vec[i] * wrk[i] < T(0)) 
          *val = -(*val);
        AT::copy(wrk, vec, N);
        break; }
      AT::copy(wrk, vec, N);
      numIts++; }}
  return numIts;
}
// -maybe get rid of the function that that computes only the largest - give this function another
//  parameter that determines, how many eigenvalues should be computed, and if it's 1, it just 
//  reduces to the function that computes the largest.
// -maybe have a boolean parameter to switch between finding eigenspace of A or A^-1 - in the 2nd
//  case, we may have to use solve(A, vec, wrk) or solve(A, wrk, vec) instead of 
//  product(A, vec, wrk) ...maybe we could also finde eigenvalues of A^T, A^-T
// -maybe try to improve convergence by using a matrix with suitably shifted eigenvalues, the 
//  shift should be such as to maximize the ratio between the largest and second-to-largest 
//  (remaining, absolute) eigenvalue. maybe the trace could help: it's the sum of all eigenvalues,
//  divide by N to get the average and subtract that average to center the eigenvalues around zero.
//  ...but for the 2nd eigenpair, we want to center the *remaining* eigenvalues around zero - maybe 
//  subtract all the already found eigenvalues from the trace before dividing by N -> experiments
//  needed. ..the determinant is the product of all eigenvalues btw - but i don't think that's 
//  helpful here (it's hard to compute anyway)

//=================================================================================================


template<class T> 
bool isHarmonic2(const rsBivariatePolynomial<T>& u, T tol = T(0))
{
  for(int m = 2; m <= u.getDegreeX(); m++) {
    for(int n = 2; n <= u.getDegreeY(); n++) {
      T c_xx = m * (m-1) * u.coeff(m  , n-2);   // coeff for x^(m-2) * y^(n-2) in u_xx
      T c_yy = n * (n-1) * u.coeff(m-2, n  );   // coeff for x^(m-2) * y^(n-2) in u_yy
      if(rsAbs(c_xx + c_yy) > tol)
        return false;   }}
  return true;
}
// -needs more tests
// -if it works, maybe replace the old implementation in rsBivariatePolynomial - this here is more
//  efficient
// -make a function makeHarmonicY that assigns:
//    u.coeff(m-2, n) = -m * (m-1) * u.coeff(m, n-2) / (n * (n-1))
//  and a makeHarmonicX that assigns:
//    u.coeff(m, n-2) = -n * (n-1) * u.coeff(m-2, n) / (m * (m-1))
//  and maybe some sort of symmetric version that assigns both - the goal is always that 
//  c_xx + c_yy = 0  ->  m * (m-1) * u[m, n-2] + n * (n-1) * u[m-2, n] = 0
//  m * (m-1) * u[m, n-2] = -n * (n-1) * u[m-2, n]
//  so:
//    u[m, n-2] / u[m-2, n] = -(n*(n-1)) / (m*(m-1))
//  i think, this is the general symmetry condition that the matrix for a harmonic polynomial must
//  satisfy

template<class T> 
rsPolynomial<std::complex<T>> getComplexFromHarmonicUV(
  const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v)
{
  rsAssert(rsBivariatePolynomial<T>::areHarmonicConjugates(u, v));
  // fails: seems like ux == -vy where we should have ux == vy - there's some sign confusion: i 
  // think P_y is not the harmonic conjugate of of P_x, -P_y is

  /*
  int M = rsMax(u.getDegreeX(), v.getDegreeX());
  rsPolynomial<std::complex<T>> p(M);
  for(int i = 0; i <= M; i++)
  p[i] = std::complex<T>(u.getCoeffPadded(i, 0), v.getCoeffPadded(i, 0));
  return p;
  // seems like if v.degX > u.degX, the last coeff is zero anyway - more tests needed
  */

  int M = rsMin(u.getDegreeX(), v.getDegreeX());
  rsPolynomial<std::complex<T>> p(M);
  for(int i = 0; i <= M; i++)
    p[i] = std::complex<T>(u.coeff(i, 0), v.coeff(i, 0));
  return p;
}
// needs more tests - especially: why should the loop go only up to M and not max(M,N), where 
// N = v.getDegreeX(), using a zero-padded accessor like getCoeffPadded(i, 0) for safe 
// out-of-range access returning 0? ...or maybe the loop should go only up to min(M,N)


template<class T> 
void firstFundamentalForm(const rsBivariatePolynomial<T>& x, const rsBivariatePolynomial<T>& y,
  const rsBivariatePolynomial<T>& z, rsBivariatePolynomial<T>& E, rsBivariatePolynomial<T>& F,
  rsBivariatePolynomial<T>& G)
{
  using BiPoly = rsBivariatePolynomial<T>;
  BiPoly x_u = x.derivativeX();     // dx/du
  BiPoly x_v = x.derivativeY();     // dx/dv
  BiPoly y_u = y.derivativeX();     // dy/du
  BiPoly y_v = y.derivativeY();     // dy/dv
  BiPoly z_u = z.derivativeX();     // dz/du
  BiPoly z_v = z.derivativeY();     // dz/dv
  E = x_u*x_u + y_u*y_u + z_u*z_u;
  F = x_u*x_v + y_u*y_v + z_u*z_v;
  G = x_v*x_v + y_v*y_v + z_v*z_v;
}
// see Weitz - Differentialgeometrie, p.154
// -maybe make a function that takes as input an rsVector3D of rsBivariatePolynomial instead of 
//  x,y,z separately
// -are E,F,G the entries of J^T * J where J is the Jacobian and E*G - 2*F is its determinant?

// todo: 
// -consider parametric surfaces given by a triple of bivariate polynomials:
//  x(u,v), y(u,v), z(u,v) (note that u,v are here the independent variables, i.e. the inputs to
//  the 3 component polynomials, in the context of complex analysis, they were used for the real 
//  and imaginary part of the function - don't get confused by this!)
// -write a function to compute the first fundamental form coeffs E,F,G of the surface
// -functions for principal curvatures, mean curvature and Gaussian curvature - these should
//  all be bivariate polynomials again (if i'm not mistaken)
// -circulation ond flux through a curve (path integrals over curl and divergence)
// -maybe integral functions that take univariate polynomials for the integration limits


// experimental - doesn't seem to work:
template<class T> 
void vectorPotential2(const rsTrivariatePolynomial<T>& f, const rsTrivariatePolynomial<T>& g,
  const rsTrivariatePolynomial<T>& h, const rsTrivariatePolynomial<T>& d, 
  rsTrivariatePolynomial<T>& F, rsTrivariatePolynomial<T>& G, rsTrivariatePolynomial<T>& H)
{
  using TP = rsTrivariatePolynomial<T>;
  TP fz   = f.integralZ();
  TP gz   = g.integralZ();
  TP fz_x = fz.derivativeX();
  TP gz_y = gz.derivativeY();
  TP a_x  = h + fz_x + gz_y;
  TP a    = a_x.integralX();

  TP fz_y = fz.derivativeY();
  TP gz_x = gz.derivativeX();
  TP a_y  = a.derivativeY();
  TP D    = d + fz_y - gz_x;
  TP b_x  = D - a_y;
  TP b    = b_x.integralX();

  F = b + gz;       // F(x,y,z) =  gz + b(x,y)
  G = a - fz;       // G(x,y,z) = -fz + a(x,y)
  H = TP(0,0,0);    // H(x,y,z) =  0
}


template<class T>
void potentialToDivergence(const rsBivariatePolynomial<T>& P, rsBivariatePolynomial<T>& D)
{
  int M = P.getDegreeX(); 
  int N = P.getDegreeY();
  D.initialize(M, N);
  for(int m = 0; m <= M; m++) {
    for(int n = 0; n <= N; n++) {
      if(m <= M-2 && n <= N-2)
        D.coeff(m, n) = (m+1)*(m+2)*P.coeff(m+2, n) + (n+1)*(n+2)*P.coeff(m, n+2);
      else if(m <= M-2 && n > N-2)
        D.coeff(m, n) = (m+1)*(m+2)*P.coeff(m+2, n);
      else if(m > M-2 && n <= N-2)
        D.coeff(m, n) = (n+1)*(n+2)*P.coeff(m, n+2);
      else
        D.coeff(m, n) = T(0); }}
}

template<class T>
void divergenceToPotential1(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int m = 0; m <= M; m++) 
    for(int n = 0; n <= N; n++) 
      P.coeff(m+2, n) = (D.coeff(m, n) - (n+1)*(n+2)*P.coeff(m, n+2)) / ((m+1)*(m+2));
}
// something is still wrong - now it produces too many nonzero entries - i think, it works only 
// when the bottom-right 2x2 section of the divergence is zero?
// seems like we could fix it by zeroing the last two rows ..or well - we need to do something
// that creates zeros in the last two rows in the computed divergence

template<class T>
void divergenceToPotential4(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int m = 0; m <= M; m++) 
    for(int n = 0; n <= N; n++) 
      P.coeff(m+2, n) = (D.coeff(m, n) - (n+1)*(n+2)*P.coeff(m, n+2)) / ((m+1)*(m+2));

  //for(int m = M; m <= M+2; m++) 
  //  for(int n = 0; n <= N+2; n++)
  //    P.coeff(m, n) = T(0);
}


template<class T>
void divergenceToPotential2(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int n = 0; n <= N; n++)
    for(int m = 0; m <= M; m++)
      P.coeff(m, n+2) = (D.coeff(m, n) - (m+1)*(m+2)*P.coeff(m+2, n)) / ((n+1)*(n+2));
}
// note that we had to exchange the loops compared to divergenceToPotential1 

// experimental:
template<class T>
void divergenceToPotential3(const rsBivariatePolynomial<T>& D, rsBivariatePolynomial<T>& P)
{
  int M = D.getDegreeX();
  int N = D.getDegreeY();
  P.initialize(M+2, N+2);
  for(int n = 0; n <= N; n++)
    for(int m = 0; m <= M; m++)
      if(m <= M-2 /* || n <= N-2*/)
        P.coeff(m, n+2) = (D.coeff(m, n) - (m+1)*(m+2)*P.coeff(m+2, n)) / ((n+1)*(n+2));
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
