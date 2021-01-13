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

/** Class for representing bivariate polynomials, i.e. polynomials in two variables x and y. They 
are represented by an MxN matrix A of coefficients which is supposed to be sandwiched between two 
vectors of powers of x and y as in X^T * A * Y, where X,Y denote the vectors constructed from 
powers of x and y and ^T denotes transposition, making X^T a row vector. For example, a polynomial
that has degree 2 in x and degree 3 in y looks like:

  p(x,y) = |x^0 x^1 x^2| * |a00 a01 a02 a03| * |y^0|
                           |a10 a11 a12 a13|   |y^1|
                           |a20 a21 a22 a23|   |y^2|
                                               |y^3|

         =   a00*x^0*y^0 + a01*x^0*y^1 + a02*x^0*y^2 + a03*x^0*y^3
           + a10*x^1*y^0 + a11*x^1*y^1 + a12*x^1*y^2 + a13*x^1*y^3
           + a20*x^2*y^0 + a21*x^2*y^1 + a22*x^2*y^2 + a23*x^2*y^3

so the shape of the matrix is 3x4 with M=3, N=4. In general, it's (degX+1)x(degY+1) where degX, 
degY are the degrees with respect to x and y. */

template<class T>
class rsBivariatePolynomial
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime 

  rsBivariatePolynomial() {}

  rsBivariatePolynomial(int degreeX, int degreeY) : coeffs(degreeX+1, degreeY+1) {}

  rsBivariatePolynomial(int degreeX, int degreeY, std::initializer_list<T> l) 
    : coeffs(degreeX+1, degreeY+1, l) {}


  //-----------------------------------------------------------------------------------------------
  // \name Evaluation (High Level)

  /** Evaluates the polynomial at the given (x,y). */
  T evaluate(T x, T y) const;

  /** Evaluates the polynomial for a given x. The result is a univariate polynomial in y. */
  rsPolynomial<T> evaluateX(T x) const;

  /** Evaluates the polynomial for a given y. The result is a univariate polynomial in x. */
  rsPolynomial<T> evaluateY(T y) const;

  /** Takes a univariate polynomial x = x(y) as input and replaces each occurrence of x in this
  bivariate polynomial p(x,y) by the polynomial expression x(y). This results in a univariate
  polynomial in y. */
  rsPolynomial<T> evaluateX(const rsPolynomial<T>& x) const;

  /** Takes a univariate polynomial y = y(x) as input and replaces each occurrence of y in this
  bivariate polynomial p(x,y) by the polynomial expression y(x). This results in a univariate
  polynomial in x. */
  rsPolynomial<T> evaluateY(const rsPolynomial<T>& y) const;

  //-----------------------------------------------------------------------------------------------
  // \name Evaluation (Low Level)

  /** Evaluates the polynomial for a given x. The result is a univariate polynomial in y whose 
  coefficients are stored in py. */
  void evaluateX(T x, T* py) const;
  // make function static, taking an rsMatrixView parameter, just like derivativeX

  /** Evaluates the polynomial for a given y. The result is a univariate polynomial in x whose 
  coefficients are stored in px. */
  void evaluateY(T y, T* px) const;

  /** Used internally by evaluateX(const rsPolynomial<T>& x). Let M,N be the degrees in x and y of
  this bivariate polynomial an K be the degree of the passed univariate polynomial. Then, the 
  degree of the result py will be given by K*M + N. The workspace must have a size of K*M+1.  */
  void evaluateX(const T* px, int xDeg, T* py, T* workspace) const;

  /** Used internally by evaluateY(const rsPolynomial<T>& y). Let M,N be the degrees in x and y of
  this bivariate polynomial an K be the degree of the passed univariate polynomial. Then, the 
  degree of the result px will be given by K*N + M. The workspace must have a size of K*N+1.  */
  void evaluateY(const T* py, int yDeg, T* px, T* workspace) const;


  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic

  /** Computes the coefficients of a bivariate polynomial r that is given as the product of two
  univariate polynomial p and q that are functions of x and y alone, respectively, such that 
  r(x,y) = p(x) * q(y). */
  static rsBivariatePolynomial<T> multiply(const rsPolynomial<T>& p, const rsPolynomial<T>& q);

  static void multiply(const T* p, int pDeg, const T* q, int qDeg, rsMatrixView<T>& r);

  static void weightedSum(const rsMatrixView<T>& p, T wp, const rsMatrixView<T>& q, T wq, 
    rsMatrixView<T>& r);


  /** Computes the coefficients of a bivariate polynomial r that is given as the composition of 
  a linear combination of x and y and a univariate polynomial p: r(x,y) = p(a*x + b*y). */
  static rsBivariatePolynomial<T> composeWithLinear(const rsPolynomial<T>& p, T a, T b);

  static rsBivariatePolynomial<T> composeWithLinearOld(const rsPolynomial<T>& p, T a, T b);

  /** Given a bivariate polynomial p(x,y) and two univariate polynomials x(t), y(t), this function 
  computes p(x(t),y(t)) which is a univariate polynomial in t.  */
  static rsPolynomial<T> compose(const rsBivariatePolynomial<T>& p,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y);



  /** Multiplies this bivariate polynomial with a univariate polynomial in y only and returns the 
  result which is again a bivariate polynomial. This amounts to convolving each row of our 
  coefficient matrix with the coefficient array of p. */
  rsBivariatePolynomial<T> multiplyY(const rsPolynomial<T>& p) const;


  void negate() { coeffs.negate(); }

  void scale(T s) { coeffs.scale(s); }

  // todo: make a similar multiplyX method - this needs to convolve the columns with p, so we will
  // need a convolution routine with strides


  /** Given a complex valued bivariate polynomial, this function splits it into two real-valued 
  bivariate polynomials...tbc... */
  static void splitRealImag(const rsBivariatePolynomial<std::complex<T>>& p,
    rsBivariatePolynomial<T>& pRe, rsBivariatePolynomial<T>& pIm);


  template<class T2>
  rsBivariatePolynomial<T2> convert(T2 dummy) const
  {
    rsBivariatePolynomial<T2> p(getDegreeX(), getDegreeY());
    for(int m = 0; m <= getDegreeX(); m++)
      for(int n = 0; n <= getDegreeY(); n++)
        p.coeff(m, n) = T2(coeffs(m, n));
    return p;
  }


  /** Given a complex polynomial (or more generally, a complex function), the associated Polya 
  vector field is the complex conjugate of the vector field that would result from just intepreting
  real and imaginary parts of the function's output as x- and y-coordinates of a 2D vector field. 
  That just means, the y-part is negated, i.e. fx(x,y) = Re(p(z)), fy(x,y) = -Im(p(z)) where 
  z = x + i*y. The reason for this seemingly unnatural negation is that the resulting vector field 
  will be conservative when the original function is analytic. This is a desirable property for 
  vector fields and it does not hold true without the negation (-> verify that).  */
  static void polyaVectorField(const rsPolynomial<std::complex<T>>& p,
    rsBivariatePolynomial<T>& px, rsBivariatePolynomial<T>& py);


  //-----------------------------------------------------------------------------------------------
  // \name Factory

  /** Creates the zero polynomial. */
  static rsBivariatePolynomial<T> zero() { return rsBivariatePolynomial<T>(0, 0, { 0 }); };


  //-----------------------------------------------------------------------------------------------
  // \name Calculus

  static void derivativeX(const rsMatrixView<T>& c, rsMatrixView<T>& d);
  rsBivariatePolynomial<T> derivativeX() const;

  static void derivativeY(const rsMatrixView<T>& c, rsMatrixView<T>& d);
  rsBivariatePolynomial<T> derivativeY() const;


  static void integralX(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c = T(0));
  rsBivariatePolynomial<T> integralX(T c = T(0)) const;
  // should the integration constant be a univariate polynomial in y instead of a fixed value?
  // ...yes, i think so ...or maybe get rid of the paramneter alltogether and just leave the matrix
  // entries corresponding to the terms not involving x as is

  static void integralY(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c = T(0));
  rsBivariatePolynomial<T> integralY(T c = T(0)) const;


  /** Computes the definite integral of the polynomial with respect to x with the integration 
  limits a,b. The result is a univariate polynomial in y whose coefficients are stored in py. */
  //static void integralX(T a, T b, T* py);

  /** Computes the definite integral of the polynomial with respect to y with the integration 
  limits a,b. The result is a univariate polynomial in x whose coefficients are stored in px. */
  //static void integralY(T a, T b, T* px);

  /** Computes the definite integral of the polynomial with respect to x with the integration 
  limits a,b. The result is a univariate polynomial in y. The types Ta, Tb can both be 
  independently either T (for a constant integration limit) or rsPolynomial<T> (for an integration
  limit that is a univariate polynomial in the variable that is not integrated over, here y). */
  template<class Ta, class Tb>
  rsPolynomial<T> integralX(Ta a, Tb b) const;
  // needs tests for when a and/or b are polynomials
  // if it's called with integer parameters a, b the compiler actually generates a version with
  // Ta,Tb = int even if T = double...hmm...that may be a bit bloatsome...maybe client code should
  // make sure, it does actually call it with double arguments in this case

  /** Computes the definite integral of the polynomial with respect to y with the integration 
  limits a,b. The result is a univariate polynomial in x. */
  template<class Ta, class Tb>
  rsPolynomial<T> integralY(Ta a, Tb b) const;
  // needs tests for when a and/or b are polynomials

  /** Given two bivariate polynomials px(x,y), py(x,y) that together constitute a 2D vector field 
  and assuming that this vector field is conservative (implying that a potential exists), this 
  function computes the potential. The potential of a 2D vector field given by fx(x,y), fy(x,y) is
  a single bivariate function P(x,y) whose partial derivatives with respect to x and y give the 
  original two functions fx, fy. Of course, one function gives in general less information than 
  two, so this works only for special kinds of vector fields, namely conservative vector fields. 
  The caller must ensure that px, py actually satisfy this condition - if they don't, the returned
  function is meaningless. 
  See: https://mathinsight.org/conservative_vector_field_find_potential  */
  static rsBivariatePolynomial<T> getPotential(
    const rsBivariatePolynomial<T>& px, const rsBivariatePolynomial<T>& py);
  // be consistent with regard to using get - either we should rename integralX to getIntegralX 
  // etc. or rename getPotenteial to potential - choose the variant that is consistent with 
  // rsPolynomial and rsNumericDifferentiator
  // -the convention is that the potential's *NEGATIVE* gradient should give the original 
  //  functions back - here we take just the gradient - change that...we may need to ripple the 
  //  negation through to the getPolyaPotential

  /** The Polya vector field of an analytic complex function (such as a polynomial) is 
  conservative, so a potential exists for such a Polya vector field. This function computes that 
  potential for a given (complex) Polynomial. The result is a bivariate real polynomial whose 
  partial derivatives with respect to x and y give the real part and the negative imaginary part 
  of the original complex polynomial. */
  static rsBivariatePolynomial<T> getPolyaPotential(const rsPolynomial<std::complex<T>>& p);


  rsBivariatePolynomial<T> getHarmonicConjugate() const;


  rsBivariatePolynomial<T> getLaplacian() const
  { return derivativeX().derivativeX() + derivativeY().derivativeY(); }


  static rsBivariatePolynomial<T> divergence2D(
    const rsBivariatePolynomial<T>& fx, const rsBivariatePolynomial<T>& fy)
  { return fx.derivativeX() + fy.derivativeY(); }
  // needs test, get rid of the 2D qualifier - it's clear that we need 2D divergence when we deal 
  // with bivariate polynomials

  static rsBivariatePolynomial<T> curl2D(
    const rsBivariatePolynomial<T>& fx, const rsBivariatePolynomial<T>& fy)
  { return fy.derivativeX() - fx.derivativeY(); }
  // needs test

  // https://en.wikipedia.org/wiki/Vector_calculus_identities#Divergence






  /** Computes the double integral of the polynomial over the given rectangle. This function 
  performs the integration over x first and then the integration over y. */
  T doubleIntegralXY(T x0, T x1, T y0, T y1) const;

  /** Computes the double integral of the polynomial over the given rectangle. This function 
  performs the integration over y first and then the integration over x. */
  T doubleIntegralYX(T x0, T x1, T y0, T y1) const;

  /** Given a scalar field p(x,y) and a parametric curve x(t), y(t) and start- and end-values a,b 
  for the parameter t, this function computes the path integral (a.k.a. line integral, curve 
  integral or contour integral) of the (polynomial) scalar field along the given (polynomial) path.
  Currently, this is implemented only for linear paths, i.e. the degrees of x(t),y(t) must be at 
  most 1. Path integrals over scalar fields can be thought of as computing the area of a curtain. 
  Imagine the path on the xy-plane being projected up to the surface landscape of p(x,y) and a 
  curtain hanging down from this projected line to the xy-plane. Of course, the "landscape" may
  also lie below the xy-plane, in which case the result would become negative. So it's more like
  a signed area difference. see: https://en.wikipedia.org/wiki/Line_integral */
  static T pathIntegral(const rsBivariatePolynomial<T>& p, 
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b);

  /** Given a vector field p(x,y), q(x,y) and a parametric curve x(t), y(t) and start- and 
  end-values a,b for the parameter t, this function computes the path integral (a.k.a. line 
  integral, curve integral or contour integral) of the (polynomial) vector field along the given
  (polynomial) path. Path integrals over vector fields can be thought of as computing the work that
  a force field does on a particle that moves through the field along the given path. It integrates
  over the scalar product of the vector field's direction vectors with the curve's "velocity" 
  vectors, i.e. over the component of the vector field that points along the curve's direction.
  see: https://en.wikipedia.org/wiki/Line_integral */
  static T pathIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b);
  // aka path integral of 2nd kind or vector path integral?
  // -maybe rename to lineIntegral or contourIntegral ...or maybe flowIntegral
  // http://ndp.jct.ac.il/tutorials/infitut2/node61.html

  /** ....
  It integrates over the 2D cross product (which is a scalar) of the vector field's direction 
  vectors with the curve's "velocity" vectors, i.e. over the component of the vector field that 
  points perpendicular to the curve's direction.
  */
  static T fluxIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b);

  /** Path integral over a vector field around the closed rectangular loop going along the 4 line 
  segments: (x0,y0) -> (x1,y0) -> (x1,y1) -> (x0,y1) -> (x0,y0). If x1 > x0 and y1 > y0, then we go
  right, up, left, down. By Green's theorem, this should be equal to the double integral over the 
  enclosed rectangle of the curl of the vector field. */
  static T loopIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    T x0, T x1, T y0, T y1);
  // aka circulation? maybe rename to circulationIntegral or loopFlow, flowAroundLoop

  /** Like loopIntegral but using the flux instead of the flow. The flux is the component of the 
  vector field is perpendicular to the curve. The integral measures, how much of a fluid flows
  out of the given rectangle. By Gauss' theorem, this should be equal to the double integral over 
  the enclosed rectangle of the divergence of the vector field. */
  static T outfluxIntegral(const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
    T x0, T x1, T y0, T y1);
  // i made up this name - figure out, if there is a more proper name for that, maybe loopFlux,
  // fluxThroughLoop

  // todo: implement flux integral around a rectangular loop

  // todo:
  // -implement scalar path integrals - they are possible only numerically, due to the square-root 
  //  in the integrand (for the speed) ...except when x(t), y(t) are both linear - in this case, the 
  //  speed is a constant number and we can just scale everything by it
  // -implement flux and circulation integrals, if possible - compare results to double integral
  //  (Gauss' and Green's theorem - the vector field must satisfy some conditions for them to hold)



  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void initialize(int degX, int degY)
  {
    coeffs.setShape(degX+1, degY+1);
    coeffs.setToZero();
  }
  // todo: write a setDegrees function that takes over the content of the old matrix and fills up
  // with zeros



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getDegreeX() const { return coeffs.getNumRows()-1; }

  int getDegreeY() const { return coeffs.getNumColumns()-1; }

  /** Implements a sort of weak equality comparison of this polynomial with a second polynomial q. 
  It has a tolerance for the comparisons of individual coefficients and in the second polynomial q
  may have a differently shaped coefficient matrix, if the coeffs that have no corresponding partner 
  in the respective other matrix are (close to) zero. So, the "overlapping" coeffs must be close to 
  their partner and the non-overlapping coeffs must be close to zero. */
  bool isCloseTo(const rsBivariatePolynomial<T>& q, T tol = T(0)) const
  {

    int M = rsMax(coeffs.getNumRows(),    q.coeffs.getNumRows());
    int N = rsMax(coeffs.getNumColumns(), q.coeffs.getNumColumns());
    for(int m = 0; m < M; m++) {
      for(int n = 0; n < N; n++) {
        T d = coeffs.getElementPadded(m, n) - q.coeffs.getElementPadded(m, n);
        if(rsAbs(d) > tol)
          return false;     }}
    return true;
  }
  // move out of class

  bool isHarmonic(T tol = T(0)) const { return getLaplacian().isCloseTo(zero(), tol); }

  static bool areHarmonicConjugates(
    const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v, T tol = T(0));


  //-----------------------------------------------------------------------------------------------
  // \name Operators

  bool operator==(const rsBivariatePolynomial<T>& rhs) const { return coeffs == rhs.coeffs; }


  rsBivariatePolynomial<T> operator+(const rsBivariatePolynomial<T>& q) const 
  { 
    int M = rsMax(coeffs.getNumRows(),    q.coeffs.getNumRows());
    int N = rsMax(coeffs.getNumColumns(), q.coeffs.getNumColumns());
    rsBivariatePolynomial<T> r(M-1, N-1);
    weightedSum(coeffs, T(1), q.coeffs, T(1), r.coeffs);
    return r;
  }

  rsBivariatePolynomial<T> operator-(const rsBivariatePolynomial<T>& q) const 
  { 
    int M = rsMax(coeffs.getNumRows(),    q.coeffs.getNumRows());
    int N = rsMax(coeffs.getNumColumns(), q.coeffs.getNumColumns());
    rsBivariatePolynomial<T> r(M-1, N-1);
    weightedSum(coeffs, T(1), q.coeffs, T(-1), r.coeffs);
    return r;
  }

  rsBivariatePolynomial<T> operator*(const rsBivariatePolynomial<T>& q) const
  {
    rsBivariatePolynomial<T> r(getDegreeX()+q.getDegreeX(), getDegreeY()+q.getDegreeY());
    rsMatrixView<T>::convolve(coeffs, q.coeffs, &r.coeffs);
    return r;
  }

  //rsBivariatePolynomial<T> operator*(T a, const rsBivariatePolynomial<T>& q) const;

  /** Unary minus. */
  rsBivariatePolynomial<T> operator-() const
  {
    rsBivariatePolynomial<T> r = *this;
    r.negate();
    return r;
  }

  /** Evaluation. */
  T operator()(T x, T y) const { return evaluate(x, y); }


  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Read and write access to the (i,j)th coefficient. */
  T& coeff(int i, int j) { return coeffs(i, j); }

  /** Read access to the (i,j)th coefficient. */
  const T& coeff(int i, int j) const { return coeffs(i, j); }

  T getCoeffPadded(int i, int j, T padding = T(0)) const 
  { return coeffs.getElementPadded(i, j, padding); }


protected:

  rsMatrix<T> coeffs;

};
// maybe implement a function getHarmonicConjugate which returns 2*x*y when p = x^2 - y^2, etc.
// ...are they obtained by a rotation? does every polynomial have such a conjugate?

template<class T>
rsBivariatePolynomial<T> operator*(const T& a, const rsBivariatePolynomial<T>& q)
{
  rsBivariatePolynomial<T> r = q; r.scale(a); return r;
}

template<class T>
T rsBivariatePolynomial<T>::evaluate(T x, T y) const
{
  T xm(1), yn(1), r(0);  // x^m, y^n, result
  for(int m = 0; m < coeffs.getNumRows(); m++) {
    yn = T(1);
    for(int n = 0; n < coeffs.getNumColumns(); n++) {
      r += coeffs(m, n) * xm * yn;
      yn *= y;  }
    xm *= x; }
  return r;
}
// todo: try find an algo that works like Horner's rule 

template<class T>
void rsBivariatePolynomial<T>::evaluateX(T x, T* py) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsArrayTools::fillWithZeros(py, N);
  T xm(1);
  for(int m = 0; m < M; m++) {
    for(int n = 0; n < N; n++)
      py[n] += coeffs(m, n) * xm;
    xm *= x; }
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateX(T x) const
{
  rsPolynomial<T> py(getDegreeY());
  evaluateX(x, py.getCoeffPointer());
  return py;
}

template<class T>
void rsBivariatePolynomial<T>::evaluateY(T y, T* px) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsArrayTools::fillWithZeros(px, M);
  T yn(1);
  for(int n = 0; n < N; n++) {
    for(int m = 0; m < M; m++)
      px[m] += coeffs(m, n) * yn;
    yn *= y; }
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateY(T y) const
{
  rsPolynomial<T> px(getDegreeX());
  evaluateY(y, px.getCoeffPointer());
  return px;
}

template<class T>
void rsBivariatePolynomial<T>::evaluateX(const T* x, int xDeg, T* r, T* xm) const
{
  using AT = rsArrayTools;
  int M = getDegreeX();
  int N = getDegreeY();
  int K = xDeg;
  AT::fillWithZeros(xm, K*M+1); 
  xm[0] = T(1);
  int Km = 1; 
  for(int n = 0; n <= N; n++)
    r[n] = coeffs(0, n);
  for(int m = 1; m <= M; m++) {
    AT::convolveInPlace(xm, Km, x, K+1);
    Km += K;
    for(int i = 0; i < Km; i++)
      for(int n = 0; n <= N; n++)
        r[i+n] += coeffs(m, n) * xm[i]; }
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateX(const rsPolynomial<T>& px) const
{
  int M = getDegreeX();
  int N = getDegreeY();
  int K = px.getDegree();
  int L = K*M + N;
  rsPolynomial<T> r(L);
  std::vector<T> xm(K*M+1);
  evaluateX(px.getCoeffPointerConst(), K, r.getCoeffPointer(), &xm[0]);
  return r;
}

template<class T>
void rsBivariatePolynomial<T>::evaluateY(const T* y, int yDeg, T* r, T* yn) const
{
  using AT = rsArrayTools;
  int M = getDegreeX();
  int N = getDegreeY();
  int K = yDeg;
  AT::fillWithZeros(yn, K*N+1); 
  yn[0] = T(1);                           // initially, y^n the constant 1, i.e. y^0
  int Kn = 1;                             // current effective length of yn
  for(int m = 0; m <= M; m++)
    r[m] = coeffs(m, 0);                  // copy coeffs from the left column
  for(int n = 1; n <= N; n++) {
    AT::convolveInPlace(yn, Kn, y, K+1);  // multiply y^n by y again to get the next power
    Kn += K;                              // length of yn has increased by K
    for(int i = 0; i < Kn; i++)
      for(int m = 0; m <= M; m++)
        r[i+m] += coeffs(m, n) * yn[i]; } // accumulate
}

template<class T>
rsPolynomial<T> rsBivariatePolynomial<T>::evaluateY(const rsPolynomial<T>& py) const
{
  int M = getDegreeX();
  int N = getDegreeY();
  int K = py.getDegree();
  int L = K*N + M;            // degree of result
  rsPolynomial<T> r(L);       // result
  std::vector<T> yn(K*N+1);   // workspace, holds coeff-array of powers of y(x)
  evaluateY(py.getCoeffPointerConst(), K, r.getCoeffPointer(), &yn[0]);
  return r;
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::multiply(
  const rsPolynomial<T>& p, const rsPolynomial<T>& q)
{
  rsBivariatePolynomial<T> r(p.getDegree(), q.getDegree());
  multiply(p.getCoeffPointerConst(), p.getDegree(), q.getCoeffPointerConst(), q.getDegree(),
    r.coeffs);
  return r;
}

template<class T>
void rsBivariatePolynomial<T>::multiply(const T* p, int pDeg, const T* q, int qDeg, 
  rsMatrixView<T>& r)
{
  rsAssert(r.hasShape(pDeg+1, qDeg+1));
  for(int m = 0; m <= pDeg; m++)
    for(int n = 0; n <= qDeg; n++)
      r(m, n) = p[m] * q[n];
}

template<class T>
void rsBivariatePolynomial<T>::weightedSum(
  const rsMatrixView<T>& p, T wp, const rsMatrixView<T>& q, T wq, rsMatrixView<T>& r)
{
  rsMatrixView<T>::weightedSum(p, wp, q, wq, r);
  /*
  int M = rsMax(p.getNumRows(),    q.getNumRows());
  int N = rsMax(p.getNumColumns(), q.getNumColumns());
  rsAssert(r.hasShape(M, N));
  for(int m = 0; m < M; m++)
    for(int n = 0; n < N; n++)
      r(m, n) = wp * p.getElementPadded(m, n) + wq * q.getElementPadded(m, n);
      */
}
// todo: move to rsMatrixView

template<class T>
void rsBivariatePolynomial<T>::derivativeX(const rsMatrixView<T>& c, rsMatrixView<T>& d)
{
  int M = c.getNumRows();
  int N = c.getNumColumns();
  rsAssert(d.hasShape(M-1, N));
  for(int m = 1; m < M; m++) {
    T s(m);
    for(int n = 0; n < N; n++)
      d(m-1, n) = s * c(m, n); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::derivativeX() const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  if(M == 1)
    return rsBivariatePolynomial<T>::zero();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M-1, N);
  derivativeX(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsBivariatePolynomial<T>::derivativeY(const rsMatrixView<T>& c, rsMatrixView<T>& d)
{
  int M = c.getNumRows();
  int N = c.getNumColumns();
  rsAssert(d.hasShape(M, N-1));
  for(int n = 1; n < N; n++) {
    T s(n);
    for(int m = 0; m < M; m++)
      d(m, n-1) = s * c(m, n); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::derivativeY() const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  if(N == 1)
    return rsBivariatePolynomial<T>::zero();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M, N-1);
  derivativeY(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsBivariatePolynomial<T>::integralX(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c)
{
  int M = p.getNumRows();
  int N = p.getNumColumns();
  rsAssert(pi.hasShape(M+1, N));
  for(int n = 0; n < N; n++)   // i think, this is wrong: just set pi(0,0) to c and the rest to 0
    pi(0, n) = c; 
  for(int m = 1; m <= M; m++) {
    T s = T(1) / T(m);
    for(int n = 0; n < N; n++)
      pi(m, n) = s * p(m-1, n); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::integralX(T c) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M+1, N);
  integralX(coeffs, q.coeffs, c);
  return q;
}

template<class T>
void rsBivariatePolynomial<T>::integralY(const rsMatrixView<T>& p, rsMatrixView<T>& pi, T c)
{
  int M = p.getNumRows();
  int N = p.getNumColumns();
  rsAssert(pi.hasShape(M, N+1));
  for(int m = 0; m < M; m++)   // i think, this is wrong: just set pi(0,0) to c and the rest to 0
    pi(m, 0) = c; 
  for(int n = 1; n <= N; n++) {
    T s = T(1) / T(n);
    for(int m = 0; m < M; m++)
      pi(m, n) = s * p(m, n-1); }
}

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::integralY(T c) const
{
  int M = coeffs.getNumRows();
  int N = coeffs.getNumColumns();
  rsBivariatePolynomial<T> q;
  q.coeffs.setShape(M, N+1);
  integralY(coeffs, q.coeffs, c);
  return q;
}

template<class T>
template<class Ta, class Tb>
rsPolynomial<T> rsBivariatePolynomial<T>::integralX(Ta a, Tb b) const
{
  rsBivariatePolynomial<T> P = integralX();
  rsPolynomial<T> Pb = P.evaluateX(b);
  rsPolynomial<T> Pa = P.evaluateX(a);
  return Pb - Pa;
}
// todo: make workspace-based version(s)

template<class T>
template<class Ta, class Tb>
rsPolynomial<T> rsBivariatePolynomial<T>::integralY(Ta a, Tb b) const
{
  rsBivariatePolynomial<T> P = integralY();
  rsPolynomial<T> Pb = P.evaluateY(b);
  rsPolynomial<T> Pa = P.evaluateY(a);
  return Pb - Pa;
}
// todo: make workspace-based version(s)

template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::getPotential(
  const rsBivariatePolynomial<T>& px, const rsBivariatePolynomial<T>& py)
{
  //rsAssert(px.derivativeY() == py.derivativeX(), "px, py are not a potential field");
  // we need a weaker notion of equality here: allow different formal shapes of the coeff matrices, 
  // allow tolerance for equality of coefficients, maybe have a function 
  // isPotentialField(px, py, tol) that calls px.isCloseTo(py, tol)

  rsBivariatePolynomial<T> Px, Px_y, gyp, gy;
  Px   = px.integralX();    // integrate px with respect to x
  Px_y = Px.derivativeY();  // differentiate the result with respect to y
  gyp  = py - Px_y;         // g'(y): derivative of integration "constant" g(y)..
  gy   = gyp.integralY();   // ..which is still a function of y
  return Px + gy;           // Px + gy is the desired potential function P(x,y)
}
// -maybe optimize: gyp has nonzero coeffs only for terms that are free of any x
// -maybe implement different algorithms, integrating py with respect to y first, etc. - they 
//  should all give the same result up to roundoff error


template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::getPolyaPotential(
  const rsPolynomial<std::complex<T>>& p)
{
  rsBivariatePolynomial<T> px, py;
  polyaVectorField(p, px, py);
  return getPotential(px, py);       // (px, py) is a potential field -> compute its potential
}

template<class T> 
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::getHarmonicConjugate() const
{
  rsAssert(isHarmonic());  // needs tolerance
  using BiPoly = rsBivariatePolynomial<T>;
  BiPoly u_dx_iy = derivativeX().integralY();
  BiPoly u_dy = derivativeY(); 
  for(int m = 0; m <= u_dy.getDegreeX(); m++)
    for(int n = 1; n <= u_dy.getDegreeY(); n++)
      u_dy.coeff(m, n) = T(0);
  return u_dx_iy - u_dy.integralX();
}
// -needs tests 

template<class T> 
T rsBivariatePolynomial<T>::doubleIntegralXY(T x0, T x1, T y0, T y1) const
{
  rsPolynomial<T>& ix = integralX(x0, x1);  // still a function of y
  return ix.definiteIntegral(y0, y1);       // just a number
}

template<class T> 
T rsBivariatePolynomial<T>::doubleIntegralYX(T x0, T x1, T y0, T y1) const
{
  rsPolynomial<T>& iy = integralY(y0, y1);  // still a function of x
  return iy.definiteIntegral(x0, x1);       // just a number
}
// make function names consistent maybe rename definiteIntegral to integral...but this could lead
// to ambiguities with the function that evaluates the indefinite integral at some x0 with 
// integration constant c

template<class T> 
T rsBivariatePolynomial<T>::pathIntegral(const rsBivariatePolynomial<T>& p, 
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b)
{
  rsAssert(x.getDegree() <= 1 && y.getDegree() <= 1, "Only implemented for linear paths");
  // ToDo: lift this restriction by switching to numerical integration, if any degree is larger.
  // When x(t) or y(t) is nonlinear, then xp or yp would still be functions of t and the sqrt would
  // appear in the integrand and we would have to evaluate x,y at t within the integral.

  using Poly   = rsPolynomial<T>;
  using BiPoly = rsBivariatePolynomial<T>;
  Poly pt = BiPoly::compose(p, x, y);     // p(t)
  Poly xp = x.derivative();               // x'(t) - optimize: it's just the coeff for x^1
  Poly yp = y.derivative();               // y'(t)   ...but take care: the degree may be zero
  T dx = xp(0);                           // dx/dt, constant bcs x is linear
  T dy = yp(0);                           // dy/dt, constant bcs y is linear
  T ds = sqrt(dx*dx + dy*dy);             // arc length (the "arc" is just a line element)
  return ds * pt.definiteIntegral(a, b);
}
// pt has 2 zero coeffs at the end in the test - check, why that happens - should we expect that or
// is that a bug?

template<class T> 
T rsBivariatePolynomial<T>::pathIntegral(
  const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b)
{
  using Poly   = rsPolynomial<T>;
  using BiPoly = rsBivariatePolynomial<T>;
  Poly ut = BiPoly::compose(u, x, y);  // u(t) - todo: support syntax: u.compose(x, y)
  Poly vt = BiPoly::compose(v, x, y);  // v(t)
  Poly xp = x.derivative();            // x'(t)
  Poly yp = y.derivative();            // y'(t)
  Poly P  = ut * xp + vt * yp;         // the scalar product in the integrand
  return P.definiteIntegral(a, b);
}
// todo: maybe have a function that leaves a,b, unspecified - this should return a single 
// univariate polynomial into which the limits can be inserted later, i.e. just return
// P.integral or p.indefiniteIntegral

template<class T> 
T rsBivariatePolynomial<T>::fluxIntegral(
  const rsBivariatePolynomial<T>& p, const rsBivariatePolynomial<T>& q,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, T a, T b)
{
  using Poly   = rsPolynomial<T>;
  using BiPoly = rsBivariatePolynomial<T>;
  Poly pt = BiPoly::compose(p, x, y);  // p(t)
  Poly qt = BiPoly::compose(q, x, y);  // q(t)
  Poly xp = x.derivative();            // x'(t)
  Poly yp = y.derivative();            // y'(t)
  Poly P  = pt * yp - qt * xp;
  return P.definiteIntegral(a, b);
}
// needs tests - this is on shaky grounds - especially with respect to computing and 
// (not) normalizing the normal vector (which is taken to be equal to (nx, ny) = (yp, -xp) here...
// but actually should be normalized...right?). ...in general, the function is very similar to
// pathIntegral with P = ut * xp + vt * yp  replaced by  P = pt * yp - qt * xp

template<class T> 
T rsBivariatePolynomial<T>::loopIntegral(const rsBivariatePolynomial<T>& p,
  const rsBivariatePolynomial<T>& q, T x0, T x1, T y0, T y1)
{
  using P = rsPolynomial<T>;
  P xt, yt;
  T r = T(0);
  xt = P({x0, x1-x0}); yt = P({y0}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // rightward
  xt = P({x1}); yt = P({y0, y1-y0}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // upward
  xt = P({x1, x0-x1}); yt = P({y1}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // leftward
  xt = P({x0}); yt = P({y1, y0-y1}); r += pathIntegral(p, q, xt, yt, T(0), T(1));  // downward
  return r;
}
// optimize: avoid creating so many temporary univariate polynomials xt, yt - reuse the existing 
// objects by re-assigning the coeffs

template<class T> 
T rsBivariatePolynomial<T>::outfluxIntegral(const rsBivariatePolynomial<T>& p, 
  const rsBivariatePolynomial<T>& q, T x0, T x1, T y0, T y1)
{
  using P = rsPolynomial<T>;
  P xt, yt;
  T r = T(0);
  xt = P({x0, x1-x0}); yt = P({y0}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // rightward
  xt = P({x1}); yt = P({y0, y1-y0}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // upward
  xt = P({x1, x0-x1}); yt = P({y1}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // leftward
  xt = P({x0}); yt = P({y1, y0-y1}); r += fluxIntegral(p, q, xt, yt, T(0), T(1));  // downward
  return r;
}
// maybe get rid of the duplication by using a pointer to fluxIntegral/pathIntegral

template<class T>
bool rsBivariatePolynomial<T>::areHarmonicConjugates(
  const rsBivariatePolynomial<T>& u, const rsBivariatePolynomial<T>& v, T tol)
{
  // We must have: u_x == v_y and u_y == -v_x:
  using BiPoly = rsBivariatePolynomial<T>;
  bool r = true;
  BiPoly ux = u.derivativeX();
  BiPoly uy = u.derivativeY();
  BiPoly vx = v.derivativeX(); vx.negate();
  BiPoly vy = v.derivativeY();
  r &= ux.isCloseTo(vy, tol);
  r &= uy.isCloseTo(vx, tol);
  return r;
}
// try to optimize: avoid creating the temporary BiPolys. Instead, compare coeffs in u, v directly.

// optimized version of composeWithLinear
template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::composeWithLinear(
  const rsPolynomial<T>& p, T a, T b)
{
  int N = p.getDegree();
  rsBivariatePolynomial<T> r(N, N);
  r.coeffs.setToZero();
  const T* c = p.getCoeffPointerConst();
  int N1 = N+1;
  std::vector<T> wrk(3*N1);             // workspace
  T* B  = &wrk[0*N1];                   // binomial coeffs
  T* an = &wrk[1*N1];                   // powers of a
  T* bn = &wrk[2*N1];                   // powers of b
  an[0] = T(1);                         // a^0
  bn[0] = T(1);                         // b^0
  for(int n = 1; n <= N; n++) {
    an[n] = a * an[n-1];                // a^n
    bn[n] = b * bn[n-1]; }              // b^n
  for(int n = 0; n <= N; n++) {
    rsNextPascalTriangleLine(B, B, n);  // B[k] is now B(n,k) the binomial coeff "n-choose-k"
    for(int k = 0; k <= n; k++)
      r.coeffs(k, n-k) += c[n] * B[k] * an[k] * bn[n-k]; } // += c[n] * B(n,k) * a^k * b^(n-k)
  return r;
}
// todo:
// -let the workspace be passed by the user
// -operate on a pre-allocated rsMatrixView
// -the polynomial p should be passed as raw coefficient array
// -maybe generalize to composeWithAffine(..., T a, T b, T c) that computes the composition with 
//  (a*x + b*y + c) - i think, we need trinomial coeffients and Pascal's pyramid for this:
//  https://en.wikipedia.org/wiki/Multinomial_theorem, 
//  https://en.wikipedia.org/wiki/Pascal%27s_pyramid

// old version of composeWithLinear, just for reference
template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::composeWithLinearOld(
  const rsPolynomial<T>& p, T a, T b)
{
  int N = p.getDegree();
  rsBivariatePolynomial<T> r(N, N);
  r.coeffs.setToZero();
  const T* c = p.getCoeffPointerConst();
  for(int n = 0; n <= N; n++)
    for(int k = 0; k <= n; k++)
      r.coeffs(k, n-k) += c[n] * (T) rsBinomialCoefficient(n, k) * pow(a, k) * pow(b, n-k);
  return r;
}
// This works, but is very inefficient: "Shlemiel the painter" strikes again in 
// rsBinomialCoefficient and pow is expensive! On the other hand, it does not allocate memory.
// ....maybe keep both versions...or well.. it actually does allocate for the returned object
// 

template<class T> 
rsPolynomial<T> rsBivariatePolynomial<T>::compose(const rsBivariatePolynomial<T>& p,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y)
{
  //   p(x,y) = \sum_m \sum_n a_{mn} x^m y^n
  // where: 
  //   x = x(t) = \sum_i b_i t^i
  //   y = y(t) = \sum_j c_j t^j
  // so the univariate p(t) is:
  //   p(t) = \sum_m \sum_n \left(  a_{mn} * (\sum_i b_i t^i)^m * (\sum_j c_j t^j)^n  \right)
  // so, the algo needs:
  //   -successive powers of the b,c arrays (iterated convolutions of the arrays with themselves)
  //   -the products of all possible combinations of those powers (more convolutions)
  //   -multiply these products by a coeff a_{mn} and accumulate the result into the output coeff
  //    array
  using Poly = rsPolynomial<T>;
  using AT   = rsArrayTools;
  int M = p.getDegreeX();
  int N = p.getDegreeY();
  int I = x.getDegree();
  int J = y.getDegree();
  int L = I+M + J+N;              // degree of result
  int strideX = M*(I+1)-1;        // number of columns in matrix xp (powers of x)
  int strideY = N*(J+1)-1;        // same for yp
  rsMatrix<T> xp(M+1, strideX);   // powers of x
  rsMatrix<T> yp(N+1, strideY);   // powers of y
  std::vector<T> tmp(L+1);        // holds coeffs for (x(t))^m * (y(t))^n as function of t
  Poly pt(L);                     // target univariate polynomial p(t)
  Poly::powers(x.getCoeffPointerConst(), I, xp.getDataPointer(), M, strideX);
  Poly::powers(y.getCoeffPointerConst(), J, yp.getDataPointer(), N, strideY);
  for(int m = 0; m <= M; m++) {
    for(int n = 0; n <= N; n++) {
      T*  xm  = xp.getRowPointer(m);            // coeffs for (x(t))^m
      T*  yn  = yp.getRowPointer(n);            // coeffs for (y(t))^n
      int lxm = m*I+1;                          // effective length of xm
      int lyn = n*J+1;                          // effective length of xm
      AT::convolve(xm, lxm, yn, lyn, &tmp[0]);
      for(int i = 0; i < lxm+lyn-1; i++)
        pt[i] += p.coeff(m, n) * tmp[i];  }}    // accumulation
  return pt;
}
// what about the composition where the inner polynomial is bivariate and the outer univariate?


template<class T>
rsBivariatePolynomial<T> rsBivariatePolynomial<T>::multiplyY(const rsPolynomial<T>& polyY) const
{
  int degX = getDegreeX();
  int degY = getDegreeY();
  int degP = polyY.getDegree();
  rsBivariatePolynomial<T> r(degX, degY+degP);
  const T* h = polyY.getCoeffPointerConst();     // "impulse response"
  for(int i = 0; i <= degX; i++) {
    const T* x = coeffs.getRowPointerConst(i);   // input
    T* y = r.coeffs.getRowPointer(i);            // output
    rsArrayTools::convolve(x, degY+1, h, degP+1, y); }
  return r;
}
// this assumes row-major storage of matrices

template<class T>
void rsBivariatePolynomial<T>::splitRealImag(const rsBivariatePolynomial<std::complex<T>>& p,
  rsBivariatePolynomial<T>& pRe, rsBivariatePolynomial<T>& pIm)
{
  int m = p.getDegreeX();
  int n = p.getDegreeY();
  pRe.initialize(m, n);
  pIm.initialize(m, n);
  for(int i = 0; i <= m; i++) {
    for(int j = 0; j <= n; j++) {
      pRe.coeff(i, j) = p.coeff(i, j).real();
      pIm.coeff(i, j) = p.coeff(i, j).imag(); }}
}

template<class T>
void rsBivariatePolynomial<T>::polyaVectorField(const rsPolynomial<std::complex<T>>& p,
  rsBivariatePolynomial<T>& px, rsBivariatePolynomial<T>& py)
{
  using Complex = std::complex<T>;
  using BiPolyC = rsBivariatePolynomial<Complex>;
  Complex one(1, 0), im(0, 1);
  BiPolyC bp = BiPolyC::composeWithLinear(p, one, im); // bp(x, y) = p(x + i*y) = p(z)
  splitRealImag(bp, px, py);                           // extract real and imaginary parts
  py.negate();                                         // apply complex conjugation
}

// ToDo:
// -compute Hessian and its determinant and maybe (squared) eigenvalues and -vectors
// -root finding: find points (x,y) for which p(x,y) = 0. i think, these are in general not 
//  isolated points but rather curves - for example, when p(x,y) = x^2 + y^2 - 1, the unit circle
//  gives the set of zeros and there we have y^2 = sqrt(1- x^2)...not sure how to deal with this
// -maybe implement constrained optimization via Lagarange multipliers - find extremum of p(x,y)
//  subject to c(x,y) = 0 where c is the constraint - the Lagrange function will be a trivariate 
//  polynomial
// -implement differentiation of implictit function (Kaprfinger, pg. 560)
// -implement differential operators in polar-coordinates (needs bivariate rational function for
//  the 1/r factor)


//=================================================================================================

template<class T>
class rsTrivariatePolynomial
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime 

  rsTrivariatePolynomial() {}

  rsTrivariatePolynomial(int degreeX, int degreeY, int degreeZ)
  {
    setDegrees(degreeX, degreeY, degreeZ);
  }

  rsTrivariatePolynomial(int degreeX, int degreeY, int degreeZ, std::initializer_list<T> l)
  {
    setDegrees(degreeX, degreeY, degreeZ);
    std::vector<T> vl(l);
    coeffs.setData(vl);  
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setDegrees(int degreeX, int degreeY, int degreeZ)
  {
    std::vector<int> shape({degreeX+1, degreeY+1, degreeZ+1});
    coeffs.setShape(shape);
    coeffs.setToZero();  // todo: take over old data
  }

  void fillRandomly(T min = T(0), T max = T(1), int seed = 0, bool roundToInt = false)
  {
    coeffs.fillRandomly(min, max, seed, roundToInt);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getDegreeX() const { return coeffs.getExtent(0)-1; }
  int getDegreeY() const { return coeffs.getExtent(1)-1; }
  int getDegreeZ() const { return coeffs.getExtent(2)-1; }

  //-----------------------------------------------------------------------------------------------
  // \name Evaluation

  T evaluate(T x, T y, T z) const;

  rsBivariatePolynomial<T> evaluateX(T x) const;

  // todo: partial evaluation for x,y,z (returning bivriate polynomials), xy,xz,yz (returning 
  // univariate polynomials)


  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic

  static void weightedSum(const rsTrivariatePolynomial<T>& p, T wp,
    const rsTrivariatePolynomial<T>& q, T wq, rsTrivariatePolynomial<T>& r)
  { 
    int L = rsMax(p.getDegreeX(), q.getDegreeX());
    int M = rsMax(p.getDegreeY(), q.getDegreeY());
    int N = rsMax(p.getDegreeZ(), q.getDegreeZ());
    r.setDegrees(L, M, N);
    rsMultiArray<T>::weightedSum(p.coeffs, wp, q.coeffs, wq, r.coeffs); 
  }

  rsTrivariatePolynomial<T> operator+(const rsTrivariatePolynomial<T>& p) const
  { rsTrivariatePolynomial<T> r; weightedSum(*this, T(1), p, T(1), r); return r; }

  rsTrivariatePolynomial<T> operator-(const rsTrivariatePolynomial<T>& p) const
  { rsTrivariatePolynomial<T> r; weightedSum(*this, T(1), p, T(-1), r); return r; }

  rsTrivariatePolynomial<T> operator*(const rsTrivariatePolynomial<T>& p) const;

  static rsPolynomial<T> compose(const rsTrivariatePolynomial<T>& p,
    const rsPolynomial<T>& x, const rsPolynomial<T>& y, const rsPolynomial<T>& z);

  static rsBivariatePolynomial<T> compose(const rsTrivariatePolynomial<T>& p,
    const rsBivariatePolynomial<T>& x, const rsBivariatePolynomial<T>& y,
    const rsBivariatePolynomial<T>& z);


  //-----------------------------------------------------------------------------------------------
  // \name Calculus

  // todo: use pointers for output parameters

  static void derivativeX(const rsMultiArray<T>& c, rsMultiArray<T>& d);
  rsTrivariatePolynomial<T> derivativeX() const;

  static void derivativeY(const rsMultiArray<T>& c, rsMultiArray<T>& d);
  rsTrivariatePolynomial<T> derivativeY() const;

  static void derivativeZ(const rsMultiArray<T>& c, rsMultiArray<T>& d);
  rsTrivariatePolynomial<T> derivativeZ() const;

  static void gradient(const rsTrivariatePolynomial<T>& f, rsTrivariatePolynomial<T>& f_x, 
    rsTrivariatePolynomial<T>& f_y, rsTrivariatePolynomial<T>& f_z)
  { f_x = f.derivativeX(); f_y = f.derivativeY(); f_z = f.derivativeZ(); }

  static rsTrivariatePolynomial<T> divergence(const rsTrivariatePolynomial<T>& fx, 
    const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz)
  { return fx.derivativeX() + fy.derivativeY() + fz.derivativeZ(); }

  static void curl(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, rsTrivariatePolynomial<T>& cx, 
    rsTrivariatePolynomial<T>& cy, rsTrivariatePolynomial<T>& cz);

  rsTrivariatePolynomial<T> laplacian()
  {
    return derivativeX().derivativeX() + derivativeY().derivativeY() + derivativeZ().derivativeZ();
  }
  // todo: implement computing the 2nd partial derivatives in one go -> optimization




  static void integralX(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c = T(0));
  rsTrivariatePolynomial<T> integralX(T c = T(0)) const;


  template<class Ta, class Tb>
  rsBivariatePolynomial<T> integralX(Ta a, Tb b) const;














  /** Computes the triple integral of the polynomial over the given cuboid. This function 
  performs the integration over x first, then over y, then over z. */
  T tripleIntegralXYZ(T x0, T x1, T y0, T y1, T z0, T z1) const;


  static T pathIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const rsPolynomial<T>& x, const rsPolynomial<T>& y, 
    const rsPolynomial<T>& z, T a, T b);


  static T pathIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const std::vector<rsVector3D<T>>& path);


  /** Computes the flux of a vector field given by 3 functions fx(x,y,z), fy(x,y,z), fz(x,y,z)
  through a parametric surface patch given by x(u,v), y(u,v), z(u,v) where u and v run from u0 to
  u1 and v0 to v1 respectively. If the vector field describes a fluid velocity, the flux integral 
  measures, how much of the fluid flows through the given surface patch per unit time. */
  static T fluxIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const rsBivariatePolynomial<T>& x, 
    const rsBivariatePolynomial<T>& y, const rsBivariatePolynomial<T>& z, T u0, T u1, T v0, T v1);

  /** Flux through a triangular patch between vertices P0, P1, P2. */
  static T fluxIntegral(const rsTrivariatePolynomial<T>& fx, const rsTrivariatePolynomial<T>& fy, 
    const rsTrivariatePolynomial<T>& fz, const rsVector3D<T>& P0, const rsVector3D<T>& P1, 
    const rsVector3D<T>& P2);

  /** Computes the flux of a vector field coming out of a cuboid bounded by the given coordinates. 
  By Gauss theorem, this should be equal to the triple integral of the divergence and i think, 
  computing it that way is more efficient and accurate. This function is mostly for the sake of 
  completeness and proof of concept. */
  static T outfluxIntegral(const rsTrivariatePolynomial<T>& fx, 
    const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz,
    T x0, T x1, T y0, T y1, T z0, T z1);




  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Read and write access to the (i,j)th coefficient. */
  T& coeff(int i, int j, int k) { return coeffs(i, j, k); }

  /** Read access to the (i,j)th coefficient. */
  //const T& coeff(int i, int j) const { return coeffs(i, j); }

  //T getCoeffPadded(int i, int j, T padding = T(0)) const 
  //{ return coeffs.getElementPadded(i, j, padding); }


protected:


  rsMultiArray<T> coeffs;

};

// evaluation:

template<class T>
T rsTrivariatePolynomial<T>::evaluate(T x, T y, T z) const
{
  T xl(1), ym(1), zn(1), r(0);  // x^l, y^m, z^n, result
  for(int l = 0; l < coeffs.getExtent(0); l++) {
    ym = T(1);
    for(int m = 0; m < coeffs.getExtent(1); m++) {
      zn = T(1);
      for(int n = 0; n < coeffs.getExtent(2); n++) {
        r += coeffs(l, m, n) * xl * ym * zn;
        zn *= z; }
      ym *= y; }
    xl *= x; }
  return r;
}
// have functions where x,y,z are uni- or bivariate polynomials - the result is then again a 
// polynomial (of the same kind)

template<class T>
rsBivariatePolynomial<T> rsTrivariatePolynomial<T>::evaluateX(T x) const
{
  int L = getDegreeX();
  int M = getDegreeY();
  int N = getDegreeZ();
  rsBivariatePolynomial<T> p_yz(M, N);
  T xl(1);   // x^l
  for(int l = 0; l <= L; l++) {
    for(int m = 0; m <= M; m++) {
      for(int n = 0; n <= N; n++) {
        p_yz.coeff(m, n) += coeffs(l, m, n) * xl; }}
    xl *= x; }
  return p_yz;
}

// arithmetic:

template<class T>
void rsConvolve3D(const rsMultiArray<T>& x, const rsMultiArray<T>& h, rsMultiArray<T>& y)
{
  int Lx = x.getExtent(0), Mx = x.getExtent(1), Nx = x.getExtent(2);
  int Lh = h.getExtent(0), Mh = h.getExtent(1), Nh = h.getExtent(2);
  int Ly = Lx + Lh - 1,    My = Mx + Mh - 1,    Ny = Nx + Nh - 1;
  y.setShape({Ly, My, Ny});
  for(int l = 0; l < Ly; l++) {
    for(int m = 0; m < My; m++) {
      for(int n = 0; n < Ny; n++) {
        T s = T(0);
        for(int i = rsMax(0, l-Lx+1); i <= rsMin(Lh-1, l); i++) {
          for(int j = rsMax(0, m-Mx+1); j <= rsMin(Mh-1, m); j++) {
            for(int k = rsMax(0, n-Nx+1); k <= rsMin(Nh-1, n); k++) {
              s += h(i, j, k) * x(l-i, m-j, n-k); }}}
        y(l, m, n) = s; }}}
  int dummy = 0;
}
// needs tests
// move to rsMultiArray ...how can this be implemented for general nD convolution?

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::operator*(
  const rsTrivariatePolynomial<T>& p) const
{
  rsTrivariatePolynomial<T> r;
  rsConvolve3D(coeffs, p.coeffs, r.coeffs);
  return r;
}

template<class T> 
rsPolynomial<T> rsTrivariatePolynomial<T>::compose(const rsTrivariatePolynomial<T>& p,
  const rsPolynomial<T>& x, const rsPolynomial<T>& y, const rsPolynomial<T>& z)
{
  using Poly = rsPolynomial<T>;
  Poly p_t;            // p(t) = p(x(t), y(t), z(t))
  //Poly one({1});       // constant polynomial one(t) = 1
  Poly one(0); one[0] = T(1);
  Poly xl, ym, zn;     // (x(t))^l, (y(t))^m, (z(t))^n
  xl = one;
  for(int l = 0; l < p.coeffs.getExtent(0); l++) {
    ym = one;
    for(int m = 0; m < p.coeffs.getExtent(1); m++) {
      zn = one;
      for(int n = 0; n < p.coeffs.getExtent(2); n++) {
        p_t = p_t + p.coeffs(l, m, n) * xl * ym * zn;
        zn = zn * z; }
      ym = ym * y; }
    xl = xl * x; }
  return p_t;
}
// try to avoid the duplication with the function below by templatizing - Poly and BiPoly both need
// a factory function to create the one-polynomial or maybe more generally a constant polynomial
// ...but maybe we later optimize in which case the code will look different in both cases

template<class T> 
rsBivariatePolynomial<T> rsTrivariatePolynomial<T>::compose(const rsTrivariatePolynomial<T>& p,
  const rsBivariatePolynomial<T>& x, const rsBivariatePolynomial<T>& y,
  const rsBivariatePolynomial<T>& z)
{
  //   p(x,y,z) = \sum_l \sum_m \sum_n a_{lmn} x^l y^m z^n
  // where: 
  //   x = x(u,v) = \sum_h \sum_i b_{hi} u^h v^i
  //   y = y(u,v) = \sum_j \sum_k c_{jk} u^j v^k
  //   z = z(u,v) = \sum_q \sum_r d_{qr} u^q v^r
  // so the bivariate p(u,v) is:
  //   p(u,v) = \sum_l \sum_m \sum_n  
  //            \left(   a_{lmn} 
  //                   * (\sum_h \sum_i b_{hi} u^h v^i)^l
  //                   * (\sum_j \sum_k c_{jk} u^j v^k)^m
  //                   * (\sum_q \sum_r d_{qr} u^q v^r)^n
  //            \right)
  // so, the algo needs:
  //   -successive powers of the b,c,d matrices (iterated 2D convolutions of the matrices with 
  //    themselves)
  //   -the products of all possible combinations of those powers (more 2D convolutions)
  //   -multiply these products by a coeff a_{lmn} and accumulate the result into the output coeff
  //    array
  //   -that's complicated - let's write an inefficient prototype first that just works like 
  //    evaluation but with bivariate polynomials instead of numbers


  using BiPoly  = rsBivariatePolynomial<T>;
  BiPoly p_uv;            // p(u,v) = p(x(u,v), y(u,v), z(u,v))
  BiPoly one(0, 0, {1});  // constant bivariate polynomial one(u,v) = 1
  BiPoly xl, ym, zn;      // (x(u,v))^l, (y(u,v))^m, (z(u,v))^n

  // This is quite inefficient (lots of temporary objects are created where we could potentially 
  // work in place) - but it's readable - may be optimized later:
  xl = one;
  for(int l = 0; l < p.coeffs.getExtent(0); l++) {
    ym = one;
    for(int m = 0; m < p.coeffs.getExtent(1); m++) {
      zn = one;
      for(int n = 0; n < p.coeffs.getExtent(2); n++) {
        p_uv = p_uv + p.coeffs(l, m, n) * xl * ym * zn;
        zn = zn * z; }
      ym = ym * y; }
    xl = xl * x; }
  return p_uv;
}

// calculus:

template<class T>
void rsTrivariatePolynomial<T>::derivativeX(const rsMultiArray<T>& c, rsMultiArray<T>& d)
{
  int L = c.getExtent(0);
  int M = c.getExtent(1);
  int N = c.getExtent(2);
  rsAssert(d.hasShape({ L-1, M, N }));
  for(int l = 1; l < L; l++) {
    T s(l);
    for(int m = 0; m < M; m++)
      for(int n = 0; n < N; n++)
        d(l-1, m, n) = s * c(l, m, n); }
}
// maybe don't require d to have the proper shape - instead, set up the shape of d

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::derivativeX() const
{
  rsTrivariatePolynomial<T> q(getDegreeX()-1, getDegreeY(), getDegreeZ());
  derivativeX(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsTrivariatePolynomial<T>::derivativeY(const rsMultiArray<T>& c, rsMultiArray<T>& d)
{
  int L = c.getExtent(0);
  int M = c.getExtent(1);
  int N = c.getExtent(2);
  rsAssert(d.hasShape({ L, M-1, N }));
  for(int m = 1; m < M; m++) {
    T s(m);
    for(int l = 0; l < L; l++)
      for(int n = 0; n < N; n++)
        d(l, m-1, n) = s * c(l, m, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::derivativeY() const
{
  rsTrivariatePolynomial<T> q(getDegreeX(), getDegreeY()-1, getDegreeZ());
  derivativeY(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsTrivariatePolynomial<T>::derivativeZ(const rsMultiArray<T>& c, rsMultiArray<T>& d)
{
  int L = c.getExtent(0);
  int M = c.getExtent(1);
  int N = c.getExtent(2);
  rsAssert(d.hasShape({ L, M, N-1 }));
  for(int n = 1; n < N; n++) {
    T s(n);
    for(int m = 0; m < M; m++)
      for(int l = 0; l < L; l++)
        d(l, m, n-1) = s * c(l, m, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::derivativeZ() const
{
  rsTrivariatePolynomial<T> q(getDegreeX(), getDegreeY(), getDegreeZ()-1);
  derivativeZ(coeffs, q.coeffs);
  return q;
}

template<class T>
void rsTrivariatePolynomial<T>::curl(const rsTrivariatePolynomial<T>& fx, 
  const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz, 
  rsTrivariatePolynomial<T>& cx, rsTrivariatePolynomial<T>& cy, rsTrivariatePolynomial<T>& cz)
{
  cx = fz.derivativeY() - fy.derivativeZ();
  cy = fx.derivativeZ() - fz.derivativeX();
  cz = fy.derivativeX() - fx.derivativeY();
}

template<class T>
void rsTrivariatePolynomial<T>::integralX(const rsMultiArray<T>& a, rsMultiArray<T>& ai, T c)
{
  int L = a.getExtent(0);
  int M = a.getExtent(1);
  int N = a.getExtent(2);
  rsAssert(ai.hasShape({ L+1, M, N }));

  ai(0, 0, 0) = c; 
  for(int m = 1; m < M; m++)
    for(int n = 1; n < N; n++)
      ai(0, m, n) = 0; 
  // i think, more generally, the integration "constant" c could be a bivariate polynomial in y,z 
  // and we would do: ai(0, m, n) = c(m, n) for (m,n) = (0,0)...(M-1,N-1)

  for(int l = 1; l <= L; l++) {
    T s = T(1) / T(l);
    for(int m = 0; m < M; m++)
      for(int n = 0; n < N; n++)
        ai(l, m, n) = s * a(l-1, m, n); }
}

template<class T>
rsTrivariatePolynomial<T> rsTrivariatePolynomial<T>::integralX(T c) const
{
  rsTrivariatePolynomial<T> q(getDegreeX()+1, getDegreeY(), getDegreeZ());
  integralX(coeffs, q.coeffs, c);
  return q;
}
// needs test

template<class T>
template<class Ta, class Tb>
rsBivariatePolynomial<T> rsTrivariatePolynomial<T>::integralX(Ta a, Tb b) const
{
  rsTrivariatePolynomial<T> P = integralX();
  rsBivariatePolynomial<T> Pb = P.evaluateX(b);
  rsBivariatePolynomial<T> Pa = P.evaluateX(a);
  return Pb - Pa;
}
// needs test


template<class T> 
T rsTrivariatePolynomial<T>::tripleIntegralXYZ(T x0, T x1, T y0, T y1, T z0, T z1) const
{
  rsBivariatePolynomial<T>& ix = integralX(x0, x1);  // still a function of y and z
  return ix.doubleIntegralXY(y0, y1, z0, z1);
}
// maybe have functions where only the innermost limits x0,x1 are constant, the middle limits 
// y0,y1 are univariate polynomials in x and the outermost limits z0, z1 are bivariate polynomials
// in x,y
template<class T> 
T rsTrivariatePolynomial<T>::pathIntegral(
  const rsTrivariatePolynomial<T>& fx, 
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz, 
  const rsPolynomial<T>& x, 
  const rsPolynomial<T>& y, 
  const rsPolynomial<T>& z, 
  T a, T b)
{
  using Poly    = rsPolynomial<T>;
  using TriPoly = rsTrivariatePolynomial<T>;
  Poly fxt = TriPoly::compose(fx, x, y, z);   // fx(t) = fx(x(t),y(t),z(t))
  Poly fyt = TriPoly::compose(fy, x, y, z);   // fy(t) = fy(x(t),y(t),z(t))
  Poly fzt = TriPoly::compose(fz, x, y, z);   // fz(t) = fz(x(t),y(t),z(t))
  Poly xp = x.derivative();                   // x'(t)
  Poly yp = y.derivative();                   // y'(t)
  Poly zp = z.derivative();                   // z'(t)
  Poly P  = fxt * xp + fyt * yp + fzt * zp;   // the scalar product in the integrand
  return P.definiteIntegral(a, b);
}

template<class T> 
T rsTrivariatePolynomial<T>::pathIntegral(
  const rsTrivariatePolynomial<T>& fx, 
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz,
  const std::vector<rsVector3D<T>>& path)
{
  T result = T(0);
  rsPolynomial<T> x(1), y(1), z(1);
  for(size_t i = 1; i < path.size(); i++) {
    x[0] = path[i-1].x; x[1] = path[i].x - x[0];
    y[0] = path[i-1].y; y[1] = path[i].y - y[0];
    z[0] = path[i-1].z; z[1] = path[i].z - z[0];
    result += pathIntegral(fx, fy, fz, x, y, z, T(0), T(1)); }
  return result;
}
// needs more tests

template<class T> 
T rsTrivariatePolynomial<T>::fluxIntegral(
  const rsTrivariatePolynomial<T>& fx,
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz,
  const rsBivariatePolynomial<T>& x,
  const rsBivariatePolynomial<T>& y,
  const rsBivariatePolynomial<T>& z,
  T u0, T u1, T v0, T v1)
{
  using BiPoly  = rsBivariatePolynomial<T>;
  using TriPoly = rsTrivariatePolynomial<T>;

  // vector field on the surface:
  BiPoly gx = TriPoly::compose(fx, x, y, z);  // gx(u,v) = fx(x(u,v), y(u,v), z(u,v))
  BiPoly gy = TriPoly::compose(fy, x, y, z);  // gy(u,v) = fy(x(u,v), y(u,v), z(u,v))
  BiPoly gz = TriPoly::compose(fz, x, y, z);  // gz(u,v) = fz(x(u,v), y(u,v), z(u,v))

  // partial derivatives of gx,gy,gz with respect to u,v:
  BiPoly ax = x.derivativeX();                // ax := d x(u,v) / du
  BiPoly ay = y.derivativeX();                // ay := d y(u,v) / du
  BiPoly az = z.derivativeX();                // az := d z(u,v) / du
  BiPoly bx = x.derivativeY();                // bx := d x(u,v) / dv
  BiPoly by = y.derivativeY();                // by := d y(u,v) / dv
  BiPoly bz = z.derivativeY();                // bz := d z(u,v) / dv 

  // components of the cross-product (B�rwolff, pg 355):
  BiPoly cx = ay*bz - az*by;
  BiPoly cy = az*bx - ax*bz;
  BiPoly cz = ax*by - ay*bx;

  // differential flux element and total flux through surface (B�rwollf, pg 600):
  BiPoly df = gx*cx + gy*cy + gz*cz;
  return df.doubleIntegralXY(u0, u1, v0, v1);
}

template<class T> 
T rsTrivariatePolynomial<T>::fluxIntegral(
  const rsTrivariatePolynomial<T>& fx,
  const rsTrivariatePolynomial<T>& fy,
  const rsTrivariatePolynomial<T>& fz,
  const rsVector3D<T>& P0, const rsVector3D<T>& P1, const rsVector3D<T>& P2)
{
  rsBivariatePolynomial<T> x(1,1), y(1,1), z(1,1);
  x.coeff(0, 0) = P0.x;
  y.coeff(0, 0) = P0.y;
  z.coeff(0, 0) = P0.z;

  x.coeff(1, 0) = P1.x - P0.x;
  y.coeff(1, 0) = P1.y - P0.y;
  z.coeff(1, 0) = P1.z - P0.z;

  x.coeff(1, 1) = P2.x - P1.x;
  y.coeff(1, 1) = P2.y - P1.y;
  z.coeff(1, 1) = P2.z - P1.z;

  return fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // Parametrization was derived by considering:
  //   p1(u) = P0 + u*(P1-P0), p2(u) = P0 + u*(P2-P0)
  // which leads to:
  //   p(u,v) = (1-v)*p1(u) + v*p2(u)  
  //          = P0 + (P1-P0)*u + (P2-P1)*u*v
}

template<class T> 
T rsTrivariatePolynomial<T>::outfluxIntegral(const rsTrivariatePolynomial<T>& fx,
  const rsTrivariatePolynomial<T>& fy, const rsTrivariatePolynomial<T>& fz,
  T x0, T x1, T y0, T y1, T z0, T z1)
{
  using BiPoly = rsBivariatePolynomial<T>;

  T dx = x1 - x0;
  T dy = y1 - y0;
  T dz = z1 - z0;
  BiPoly x, y, z;

  // flux through surface patch where z = z0:
  x = BiPoly(1, 1, { x0,   0, dx, 0 });
  y = BiPoly(1, 1, { y1, -dy,  0, 0 });
  z = BiPoly(0, 0, { z0             });
  T fz0 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where z = z1:
  //x = BiPoly(1, 1, { x0,  0, dx, 0  }); // has not changed
  y = BiPoly(1, 1, { y0, dy,  0, 0  });
  z = BiPoly(0, 0, { z1             });
  T fz1 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where y = y0:
  //x = BiPoly(1, 1, { x0, 0,  dx, 0 }); // has not changed
  y = BiPoly(0, 0, { y0            });
  z = BiPoly(1, 1, { z0,  dz, 0, 0 });
  T fy0 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where y = y1:
  //x = BiPoly(1, 1, { x0, 0,  dx, 0 }); // has not changed
  y = BiPoly(0, 0, { y1            });
  z = BiPoly(1, 1, { z1, -dz, 0, 0 });
  T fy1 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where x = x0:
  x = BiPoly(0, 0, { x0           });
  y = BiPoly(1, 1, { y0,  0, dy, 0 });
  z = BiPoly(1, 1, { z1, -dz, 0, 0 });
  T fx0 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  // flux through surface patch where x = x1:
  x = BiPoly(0, 0, { x1           });
  //y = BiPoly(1, 1, { y0, 0, dy, 0 }); // has not changed
  z = BiPoly(1, 1, { z0, dz, 0, 0 });
  T fx1 = fluxIntegral(fx, fy, fz, x, y, z, T(0), T(1), T(0), T(1));

  return fx0 + fx1 + fy0 + fy1 + fz0 + fz1;
}
// -optimize 
// -create the BiPoly objects only once and re-assign the coeffs -> less allocation
// -remove the fz0, etc variables - directly add intermediate results into an accumulator


// ToDo:
// -compute potential and vector-potential
// -transformed integrals (using the determinant of the Jacobian matrix)
// -constrained optimization via Lagrange multipliers (maybe)
// -Legendre transform (maybe)

//=================================================================================================

/** A class for representing and performing computations with functions that are defined as 
piecewise polynomials. */

template<class T>
class rsPiecewisePolynomial
{

public:



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of pieces. */
  int getNumPieces() const { return (int) pieces.size(); }

  /** Returns a constant reference to the piece at index i. */
  const rsPolynomial<T>& getPieceConstRef(int i) const 
  { rsAssert(i >= 0 && i < getNumPieces(), "Index out of range"); return pieces[i]; }

  /** Returns the index of the piece, where x belongs or -1, if x is out of range to the left or
  getNumPieces() if x is out of range to the right. */
  int getIndex(T x) const;

  /** Returns the left boundary of the domain, where the function is defined. */
  T getDomainMinimum() const
  {
    if(domains.empty())
      return T(0);
    return domains[0];
  }

  /** Returns the right boundary of the domain, where the function is defined. */
  T getDomainMaximum() const
  {
    if(domains.empty())
      return T(0);
    return rsLast(domains);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Evaluation

  /** Evaluates the function at the given x. If x is outside the range where we have defined 
  pieces, it returns zero. */
  T evaluate(T x) const;

  /** Evaluation as oprator. */
  T operator()(T x) const { return evaluate(x); }


  //-----------------------------------------------------------------------------------------------
  // \name Manipulations

  /** Adds another piece to the object, possibly with a scalar weight applied to the new piece. */
  void addPiece(const rsPolynomial<T>& p, T pL, T pU, T weight = T(1));

  /** Clears the polynomial, setting it back into a freshly constructed state. */
  void clear() { domains.clear(); pieces.clear(); }

  /** Scales the whole function in the y-direction by the given factor. */
  void scale(T factor);

  /** Stretches the whole function in the x-direction by the given factor. */
  void stretch(T factor);

  /** Integrates the function. The integration constant determines the function value at the left 
  boundary. */
  void integrate(T c = T(0));

  rsPiecewisePolynomial<T> integral(T c = T(0)) const
  {
    rsPiecewisePolynomial<T> p = *this;
    p.integrate(c);
    return p;
  }


  // todo: derivative

  /** Shifts the pieces up or down in the y-direction such that they match at the segment 
  boundaries. */
  void makeContinuous();


  //-----------------------------------------------------------------------------------------------
  // \name Combination

  rsPiecewisePolynomial<T>& operator+=(const rsPiecewisePolynomial<T>& q) 
  { 
    for(int i = 0; i < q.getNumPieces(); i++)
      addPiece(q.pieces[i], q.domains[i], q.domains[i+1]);
    return *this;
  }

  rsPiecewisePolynomial<T> operator+(const rsPiecewisePolynomial<T>& q) const 
  { rsPiecewisePolynomial<T> r = *this; r += q; return r; }

  rsPiecewisePolynomial<T>& operator-=(const rsPiecewisePolynomial<T>& q) 
  { 
    for(int i = 0; i < q.getNumPieces(); i++)
      addPiece(q.pieces[i], q.domains[i], q.domains[i+1], T(-1));
    return *this;
  }

  rsPiecewisePolynomial<T> operator-(const rsPiecewisePolynomial<T>& q) const 
  { rsPiecewisePolynomial<T> r = *this; r -= q; return r; }

  // todo: implement multiplication, maybe the addPiece function can be extended to also handle 
  // multiplication, such that we do not have to replicate the splitting logic


  /** Convolves two polynomial pieces p(x) and q(x) that are defined on the domains pL..pU and 
  qL..qU respectively and assumed to be zero outside these domains (L and U stand for lower and 
  upper boundaries of the domains). The result are 3 polynomial pieces rL,rM,rR that are adjacent 
  to each other (L,M,R stand for left, middle, right). These polynomials are output parameters and
  assigend by the function. These output pieces are defined on the domains rLL..rLU, rLU..rRL, 
  rRL..rRR respectively which are also output parameters. The left/right boundaries of the middle 
  segment are the same as the right/left boundaries of the adjacent left and right pieces, so there 
  are no output parameters for the middle section's domain boundaries because they would be 
  redundant. We get 3 pieces as output because for the left piece, the two input segments that are 
  convolved, are not yet fully overlapping as q gets shifted over p and for the right piece, they 
  are not fully overlapping anymore. Only in the middle segment, there's full overlap, i.e. q is 
  fully inside p or the other way around. */
  static void convolvePieces(
    const rsPolynomial<T>& p, T pL, T pU, const rsPolynomial<T>& q, T qL, T qU,
    rsPolynomial<T>& rL, T& rLL, T& rLU, rsPolynomial<T>& rM, rsPolynomial<T>& rR, T& rRL, T& rRU);

  /** Convolves this piecewise polynomial with another piecewise polynomial p and returns the 
  result which is again a piecewise polynomial. */
  rsPiecewisePolynomial<T> convolve(const rsPiecewisePolynomial<T>& p);


  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Creates a Irwin-Hall distribution of given order N between a and b. This is the uniform 
  distribution convolved with itself N times. It arises as amplitude distribution, when you add up 
  the outputs of N independent noise generators with uniform distributions between a..b. N=0 gives 
  the uniform distribution, N=1 gives a triangular distribution, N=2 a piecewise parabolic 
  distribution and so on. */
  static rsPiecewisePolynomial<T> irwinHall(int N, T a = T(0), T b = T(1));


protected:

  std::vector<T> domains;  // maybe rename to boundaries, ends, limits, borders
  std::vector<rsPolynomial<T>> pieces;
  // -the pieces are adjacent (no verlap, no gaps)
  // -piece[i] goes from domains[i] to domains[i+1]

  // todo: maybe have an evaluation mode parameter that determines what happens when the user wants 
  // evaluate outside the domain. possible values are: zero (as it is now), just use the left/right
  // polynomials for extrapolation, clamp output at whatever left/right polynomials give at the
  // boundaries - so modes could be named: zero, extrapolate, clamp ..maybe, the extrapolation 
  // could be restricted to use only the lower order coeffs
  // 

};

template<class T>
void rsPiecewisePolynomial<T>::addPiece(const rsPolynomial<T>& p, T pL, T pU, T weight)
{
  T tol(0);  // tolerance for (floating point) equality comparisons -> make member
  auto match = [&](T x, T y) -> bool { return rsAbs(x-y) <= tol; };

  rsAssert(pL < pU);
  if(pieces.empty()) {             // initialize with the first piece
    pieces.push_back(p);
    domains.push_back(pL);
    domains.push_back(pU);
    return;  }

  // I think, these cases can now be subsumed by the code below...but maybe it's neverless useful 
  // to handle them specially for efficiency - these cases here (especially the first) come up in 
  // the supposedly common case of building up a piecewise polynomial from scratch by appending one
  // segment after another:
  if(match(pL, rsLast(domains))) { // append piece at the right end
    pieces.push_back(p);
    domains.push_back(pU);
    return;  }
  if(match(pU, domains[0])) {      // prepend piece at the left end
    rsPrepend(pieces,  p);
    rsPrepend(domains, pL);
    return; }

  // Figure out start- and end indices for segment:
  int numPieces = getNumPieces();
  int iL = getIndex(pL);
  int iU;
  if(pU >= rsLast(domains))
    iU = numPieces;
  else
    iU = getIndex(pU);

  // (don't) handle gaps:
  if(iL == numPieces || iU == -1) {  // p starts after this or ends before this
    rsError("Gaps are not allowed"); // we can't handle them with the current implementation
    return; }

  // Function to split the piece at index i into two pieces at x0:
  auto split = [&](int i, T x0) 
  { 
    rsAssert(x0 > domains[i] && x0 < domains[i+1], "x0 outside domain of piece i");
    rsInsert(domains, x0, i+1);
    rsInsert(pieces, pieces[i], i);
  };
  // maybe move out and make it a member function, maybe also have a merge function.

  // Function to accumulate polynomial q into polynomial p:
  auto accumulate = [&](rsPolynomial<T>& p, const rsPolynomial<T>& q)
  { p.addWithWeight(q, weight); };
  // todo: use a std::function that is passed as parameter, so we can use the same function also 
  // for subtraction, multiplication, etc. - but for multiplication or division, we have to do 
  // something else in the case where we prepend or append the pieces here...we'll see...
  // maybe just implement addition, subtraction and scalar multiplication (as operators) - then
  // we have at least a vector-space. elementwise multiplication and division can come later. 
  // mulitplication may actually shrink the domain because the result is nonzero only where both
  // factors are nonzero, so we need some additional logic to update the domains - but maybe we 
  // don't need to mainpulate the domains - we can just let some sections be the zero polynomial
  // ...or we could cut off zero sections at start and end as post-processing
  // ...maybe this function should have an additional "mode" parameter, 0: add, 1: multiply, 
  // 2: divide...but no - divide doe not really make sense - from a mathematical perspective, we
  // would expect a piecewise rational function to come out and not whatever results from 
  // polynomial division (the results are equal only if this is divisible by p)


  // Handle case when left boundary of p is to the left of our left boundary:
  if(iL == -1) {
    rsPrepend(domains, pL);
    rsPrepend(pieces,  p);   // this works only for addition, not multiplication
    iL = 1;
    pL = domains[iL];        // so subsequent code can be used as if we had a match
    iU++;                    // because we have a piece more now
    numPieces++; }

  // Handle case when left boundary of p is in the middle of one of our existing pieces:
  if(!match(pL, domains[iL])) {
    split(iL, pL); iL++; iU++; numPieces++; }

  // Accumulate the new polynomial into the exisitng pieces that are now in full overlap:
  for(int i = iL; i < iU; i++)
    accumulate(pieces[i], p);

  // If end of new piece is aligned with end of an existing piece, we are done:
  if(match(pU, domains[iU]))
    return;

  // ...and if it's not aligned, we have two cases to consider:
  if(iU == numPieces) {          // End of new piece is beyond or right boundary 
    domains.push_back(pU);
    pieces.push_back(p);   }     // this works only for addition, not multiplication
  else  {                        // End of new piece is in the middle of some existing piece
    split(iU, pU);
    accumulate(pieces[iU], p); }
}

template<class T>
void rsPiecewisePolynomial<T>::scale(T factor)
{
  for(size_t i = 0; i < pieces.size(); i++)
    pieces[i].scale(factor);
}

template<class T>
void rsPiecewisePolynomial<T>::stretch(T factor)
{
  for(size_t i = 0; i < pieces.size(); i++)
    pieces[i].stretch(factor);
  for(size_t i = 0; i < domains.size(); i++)
    domains[i] *= factor;
}

template<class T>
void rsPiecewisePolynomial<T>::integrate(T c)
{
  for(size_t i = 0; i < pieces.size(); i++)
    pieces[i].integrate();
  T x  = domains[0];
  T yL = pieces[0](x);
  pieces[0].shiftY(c-yL);  // adjust start value
  makeContinuous();
}

template<class T>
void rsPiecewisePolynomial<T>::makeContinuous()
{
  for(size_t i = 1; i < pieces.size(); i++)
  {
    T x  = domains[i];
    T yL = pieces[i-1](x);
    T yR = pieces[i](x);
    pieces[i].shiftY(yL-yR);
  }
}

template<class T>
int rsPiecewisePolynomial<T>::getIndex(T x) const
{
  if(pieces.empty() || x < domains[0])
    return -1;
  if( x >= rsLast(domains))
    return getNumPieces();
  int i = rsArrayTools::findSplitIndex(&domains[0], (int) domains.size(), x);
  if(x < domains[i])
    i--;
  rsAssert(i >= -1 && i < getNumPieces());
  return i;
}

template<class T>
T rsPiecewisePolynomial<T>::evaluate(T x) const
{
  int i = getIndex(x);
  if(i < 0 || i >= getNumPieces()) 
    return T(0);
  else        
    return pieces[i](x);  // use evaluate function (needs to be written)
}

template<class T>
void rsPiecewisePolynomial<T>::convolvePieces(
  const rsPolynomial<T>& p, T pL, T pU, const rsPolynomial<T>& q, T qL, T qU,
  rsPolynomial<T>& rL, T& rLL, T& rLU, rsPolynomial<T>& rM, rsPolynomial<T>& rR, T& rRL, T& rRU)
{
  // Check, which one of p,q has higher degree and invoke the function recursively with swapped 
  // inputs if they are in disadvantageous order. We don't want q to be the higher degree 
  // polynomial because it's q that gets blown up to a bivariate polynomial. Mathematically, the 
  // results should be the same, but using the lower degree polynomial as q should be more 
  // efficient and maybe also more precise. The unit tests should also pass when we comment this 
  // optimization out:
  if(q.getDegree() > p.getDegree()) {
    convolvePieces(q, qL, qU, p, pL, pU, rL, rLL, rLU, rM, rR, rRL, rRU);
    return; }

  // Compute the domains for the 3 output segments:
  T wp = pU - pL;   // width of domain of p
  T wq = qU - qL;   // width of domain of q
  rLL  = pL + qL;
  rLU  = pL + qU;   // == rLL + wq
  rRL  = pU + qL;
  rRU  = pU + qU;   // == rRL + wq
  if(wq > wp)
    rsSwap(rLU, rRL);

  // Create the bivariate polynomial PQ(x,y) = p(y)*q(x-y) which is our integrand in the 
  // convolution integral:
  using BiPoly = rsBivariatePolynomial<T>;
  BiPoly Q  = BiPoly::composeWithLinear(q, T(1), T(-1)); //  Q(x,y) = q(x-y)
  BiPoly PQ = Q.multiplyY(p);                            // PQ(x,y) = p(y)*q(x-y)

  // Integrate out the dummy variable y (typically tau in literature about convolution, and our x 
  // is their t), leaving a polynomial only in x (aka t). We get 3 segments:
  using Poly = rsPolynomial<T>;
  Poly a({-qU, 1});                           // lower integration limit (in some cases)
  Poly b({-qL, 1});                           // upper integration limit (in some cases)
  rL                  = PQ.integralY(pL, b);  // left segment
  if(     wp > wq) rM = PQ.integralY(a,  b);  // middle segment when wp longer than wq
  else if(wq > wp) rM = PQ.integralY(pL, pU); // middle segment when wq longer than wp
  else             rM = Poly(0);              // p,q have same length -> no middle segment
  rR                  = PQ.integralY(a,  pU); // right segment

  // Adjust degrees by truncating the trailing zero coefficients. I think, they arise because the
  // matrix of coefficients of the bivariate polynomial Q is triangular. The integral could 
  // potentially produce higher order nonzero coeffs but doesn't due to the special structure of Q.
  int deg = 0;
  if(wp > wq) deg = p.getDegree();
  if(wp < wq) deg = q.getDegree();
  rM.setAllocatedDegree(deg);
  deg = p.getDegree() + q.getDegree() + 1;
  rL.setAllocatedDegree(deg);
  rR.setAllocatedDegree(deg);

  // old:
  // Cut off trailing zero coefficients in the produced segments:
  //rL.truncateTrailingZeros();
  //rM.truncateTrailingZeros();
  //rR.truncateTrailingZeros();

  // The middle section may still have close-to-zero trailing coeffs (i think, it happens only 
  // when wp > wq but this needs more tests)
  // i think, the degrees are degP+degQ+1 for the L/R sections and degP or deQ for the middle M
  // section and which one of the two it is, is determined by which one has the longer domain

  // Notes:
  // -The expressions for the integration limits were found by trial and error and need more tests,
  //  especially, when q has longer support than p. Perhaps, we should switch the roles of p and q
  //  in such a case, but maybe it's more advisable to select the roles of p and q by their 
  //  degrees. We should probably create the bivariate polynomial from whichever has lower degree.
  //  ...we'll see
  // ToDo:
  // -make a function that uses a workspace - if we later need to convolve many pieces in a 
  //  piecewise polynomial, we'll call this in a double-loop: each piece from one input is 
  //  convolved with each piece from the other (and then all the results get added up), so we want
  //  the operation to be efficient
}

template<class T>
rsPiecewisePolynomial<T> rsPiecewisePolynomial<T>::convolve(const rsPiecewisePolynomial<T>& q)
{
  rsPiecewisePolynomial<T> r;
  using Poly = rsPolynomial<T>;
  Poly rL, rM, rR;
  T rLL, rLU, rRL, rRU;
  for(int i = 0; i < getNumPieces(); i++) {
    for(int j = 0; j < q.getNumPieces(); j++) {
      const Poly& pi = getPieceConstRef(i);
      T pL = domains[i];
      T pU = domains[i+1];
      const Poly& qj = q.getPieceConstRef(j);
      T qL = q.domains[j];
      T qU = q.domains[j+1];
      convolvePieces(pi, pL, pU, qj, qL, qU, rL, rLL, rLU, rM, rR, rRL, rRU);
      r.addPiece(rL, rLL, rLU);
      if(rLU < rRL)
        r.addPiece(rM, rLU, rRL);
      r.addPiece(rR, rRL, rRU);   }}
  return r;
}

template<class T>
rsPiecewisePolynomial<T> rsPiecewisePolynomial<T>::irwinHall(int order, T a, T b)
{
  rsAssert(order >= 0);
  rsPiecewisePolynomial<T> p0; 
  p0.addPiece(rsPolynomial<T>({ T(1)/(b-a) }), a, b);  // our seed function
  rsPiecewisePolynomial<T> p = p0;
  for(int i = 1; i <= order; i++)
    p = p.convolve(p0);
  return p;
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
