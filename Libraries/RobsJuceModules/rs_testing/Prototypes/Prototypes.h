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

template<class T>
class rsQuantileFilterCore2 : public rsQuantileFilterCore<T>
{

public:

  /** Produces a sample that would have been produced, if the length L of the filter would be 
  longer by one sample, i.e. L+1. This is used to implement non-integer length filters by 
  crossfading between the outputs of two filters whose lengths differ by one. 
  
  wrong:
  It needs as input the
  same input sample that has been fed to getSample and it should be called after getSample. 

  correct again:
  From 
  the values returned by the regular getSample call and the call to this afterwards, a non-integer
  length filter sample can be computed by crossfading. */



  T readOutputLongerBy1()
  {
    rsAssert(sigBuf != nullptr, "To use this feature, the input buffer must be assigned.");
    T xL = (*sigBuf)[L];   // should be x[n-L], client code must assure this
    return readOutputWithOneMoreInput(xL);
  }


  T readOutputWithOneMoreInput(T xL)
  {
    T p1 = p * T(L) / T(L-1); // but wait - p is an integer - should we use p+w or p+(1-w)?
    T w1 = p1 - floor(p1);
    T yS, yL; // hmm...yL means yLarge but xL means x[n-L] - notational clash!
    Node nx(xL, 0); // we need to create a node

    if(dblHp.small.isLess(nx, dblHp.large[0]))  // means: if(x < large[0])
    {
      // x belongs in small heap
      yS = dblHp.get2ndLargestSmallValue().value;
      yS = rsMaxViaLess(yS, xL);  //
      yL = dblHp.getLargestSmallValue().value;
    }
    else
    {
      // x belongs in large heap
      yL = dblHp.get2ndSmallestLargeValue().value;
      yL = rsMin(yL, xL);
      yS = dblHp.getSmallestLargeValue().value;
    }
    T y = (T(1)-w1)*yS + w*yL;
    return y;
  }
  // needs test - compare against signal that has beed produced by a baseclass filter that actually
  // is one sample longer
  // this is wrong - it should not take as input the sample x[n] (that is stored already in the 
  // heaps after getSample). instead, it needs x[n-L] maybe make a function 
  // readOutputWithAdditionalInput

protected:

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
    core.setModulationBuffer(&delayLine);
    core2.setModulationBuffer(&delayLine);
    dirty = true;
  }


  // setters for the parameters of the second core
  void setFrequency2(   T newFrequency) { setLength2(T(1) / newFrequency); }
  void setLength2(      T newLength)    { length2   = newLength;    dirty = true; }
  void setQuantile2(    T newQuantile1) { quantile2 = newQuantile1; dirty = true; }
  void setLowpassGain2( T newGain)      { loGain2   = newGain; }
  void setHighpassGain2(T newGain)      { hiGain2   = newGain; }

  void setCore2Equal()
  {
    length2   = length;
    quantile2 = quantile;
    loGain2   = loGain;
    hiGain2   = hiGain;
    delayScl2 = delayScl;
    dirty = true;
  }

  void setCore2Complementary()
  {
    setCore2Equal();
    quantile2 = T(1) - quantile;
    dirty = true;
  }


  // maybe have a setCore2Complementary which uses the same settings for core2 as core2 except the
  // quantile, for which we use: quantile2 = 1-quantile

  T getSample(T x)
  {
    delayLine.getSample(x);
    if(dirty) 
      updateInternals();
    T yL1 = core.getSample(x);
    T yH1 = delayLine[delayScl*delay] - yL1;
    T yL2 = core2.getSample(x);
    T yH2 = delayLine[delayScl2*delay2] - yL2;
    return loGain * yL1 + hiGain * yH1 + loGain2 * yL2 + hiGain2 * yH2;
  }

  void reset() { core.reset(); core2.reset(); delayLine.reset(); }

  virtual void updateInternals() override
  {
    // compute internal and set up core parameters:
    int L, p; T w;

    convertParameters(length, quantile, sampleRate, &L, &p, &w, &delay);
    core.setLengthAndReadPosition(L, p);
    core.setRightWeight(w);

    convertParameters(length2, quantile2, sampleRate, &L, &p, &w, &delay2);
    core2.setLengthAndReadPosition(L, p);
    core2.setRightWeight(w);

    dirty = false;
  }


protected:

  virtual void allocateResources() override
  {
    int mL = (int) ceil(maxLength * sampleRate);
    core.setMaxLength(mL);
    core2.setMaxLength(mL);
    delayLine.setCapacity(mL);
  }


  // the 2nd core and its set of parameters:
  rsQuantileFilterCore<T> core2;
  T length2   = 0.0;
  T quantile2 = 0.5;
  T loGain2   = 0.0;
  T hiGain2   = 0.0;
  T delayScl2 = 1.0;
  T delay2    = 0.0; 

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

/** Class for representing fractions a.k.a. rational numbers, i.e. ratios of two integers. 
Numerator and denominator are kept as signed integers "num", "den". On construction and in 
arithmetic operations, fractions are always put into a canonical representation which is a reduced
form where the minus sign (if any) is put into the numerator. */

template<class T>  // T should be a signed int type
class rsFraction
{

public:

  rsFraction(T numerator = T(0), T denominator = T(1)) : num(numerator), den(denominator)
  { canonicalize(); }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void set(T numerator, T denominator) { num = numerator; den = denominator; canonicalize(); }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  T getNumerator()   const { return num; }
  T getDenominator() const { return den; }

  double toDouble() const { return double(num) / double(den); }

  //-----------------------------------------------------------------------------------------------
  // \name Arithmetic operators

  // unary minus:
  rsFraction operator-() const { return rsFraction(-num, den); }
  // optimization: this should avoid calling canonicalize()

  // +,-,*,/ where both arguments are fractions:
  rsFraction operator+(const rsFraction& b) const { return rsFraction(num*b.den + b.num*den, den * b.den); }
  rsFraction operator-(const rsFraction& b) const { return rsFraction(num*b.den - b.num*den, den * b.den); }
  rsFraction operator*(const rsFraction& b) const { return rsFraction(num * b.num, den * b.den); }
  rsFraction operator/(const rsFraction& b) const { return rsFraction(num * b.den, den * b.num); }

  // ..same for integer right arguments:
  rsFraction operator+(const T& b) const { return rsFraction(num + b*den, den); }
  rsFraction operator-(const T& b) const { return rsFraction(num - b*den, den); }
  rsFraction operator*(const T& b) const { return rsFraction(num * b, den); }
  rsFraction operator/(const T& b) const { return rsFraction(num, den * b); }

  // boilerplate for the +=, -=, *=, /= operators:
  rsFraction& operator+=(const rsFraction& b) { return *this = (*this) + b; }
  rsFraction& operator-=(const rsFraction& b) { return *this = (*this) - b; }
  rsFraction& operator*=(const rsFraction& b) { return *this = (*this) * b; }
  rsFraction& operator/=(const rsFraction& b) { return *this = (*this) / b; }
  rsFraction& operator+=(const T& b) { return *this = (*this) + b; }
  rsFraction& operator-=(const T& b) { return *this = (*this) - b; }
  rsFraction& operator*=(const T& b) { return *this = (*this) * b; }
  rsFraction& operator/=(const T& b) { return *this = (*this) / b; }


  //-----------------------------------------------------------------------------------------------
  // \name Comparison operators

  bool operator==(const rsFraction& b) const { return num == b.num && den == b.den; }
  bool operator!=(const rsFraction& b) const { return !(*this == b); }
  bool operator< (const rsFraction& b) const { return num * b.den <  b.num * den; }
  bool operator<=(const rsFraction& b) const { return num * b.den <= b.num * den; }
  bool operator> (const rsFraction& b) const { return num * b.den >  b.num * den; }
  bool operator>=(const rsFraction& b) const { return num * b.den >= b.num * den; }


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Reduces this number to lowest terms. */
  void reduce() { T gcd = rsGcd(num, den); num /= gcd; den /= gcd; }

  /** Reduces to lowest terms and ensures that denominator is nonnegative. */
  void canonicalize() { reduce(); if(den < 0) { num = -num; den = -den; }  }


  T num, den;  // numerator and denominator (they are always kept canonical)

};

// operators for integer left argument:
template<class T>
rsFraction<T> operator+(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getDenominator() + r.getNumerator(), r.getDenominator()); }

template<class T>
rsFraction<T> operator-(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getDenominator() - r.getNumerator(), r.getDenominator()); }

template<class T>
rsFraction<T> operator*(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getNumerator(), r.getDenominator()); }

template<class T>
rsFraction<T> operator/(const T& i, const rsFraction<T>& r)
{ return rsFraction<T>(i * r.getDenominator(), r.getNumerator()); }

// ToDo:
// -maybe use algorithms for the arithmetic operators that make overflow less likely (divide by gcd 
//  before computing products, use lcm in + and - instead of just computing products, etc.).
// -implement some functions like pow (with integer exponent)..maybe using the ^ operator - but 
//  care has to be taken to parenthesize expressions like (r^i) inside longer expressions due to 
//  C++ precendence rules
// -maybe detect if overflow will happen and trigger an assert
// -implement functions to convert to double or float
// -find best rational approximation of a double or float by using continued fractions 
//  maybe have a function rsFraction<T> findApproximant(TFloat x)
// -maybe there are more interesting things we can do with continued fractions......maybe a 
//  function getContinuedFractionExpansion that returns an array of integers?
//
// Notes:
// -maybe it's sometimes convenient to keep it in unreduced form - it may be easier to spot 
//  patterns in sequences of unreduced rational numbers that come from some computation
// -but this will be relevant only for research code, not production code
// -maybe we could introduce a compile-time switch (maybe a boolen template parameter) that 
//  controls if we canonicalize or not
// -enforced canonical representation is important for the == operator to work properly...maybe it 
//  should be implemented in a way that admits non-canonical representations? 
//  a/b == c/d  <->  a*d == b*c

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
