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


/** Implements a quantum system that describes the spin of a particle such as an electron. This is
the most simple and prototypical quantum system and can also be used as quantum bit (qubit). The 
user can set the system into various predefined states and has functions to manipulate the state. 
The state can also be measured in which case the system will - with a certain probability 
determined by the state - fall into one of two possible pure states corresponding to the measured
variable. An example of such a measured observable is the spin along the z-axis. The measured value
will be either +1 or -1 ("up" or "down") with probabilities determined by the current state. After 
the measurement, however, this state will have been changed into a pure state such that subsequent 
measurements of the same observable will always produce the same result (with probability one).

Unlike real quantum systems, we can look into the actual state which consists of two complex 
numbers. 

...

To specify any state as a ket vector |A>, we express it as a linear combination of two (somewhat 
arbitrarily) choosen basis ket vectors |u> and |d> for "up" and "down" spin:

  |A> = au*|u> + ad*|d>


maybe rename to rsQuantumSpin - what is used as the quantum "bit" is actually one of spin's 
components, such as the z-component (a pure "up" represents binary 1 and a pure "down" state
represents binary 0. Mixed states (with respect to the z-axis) represent a superposition of 0 and
1 (a pure left or right or in or out state is "mixed" with respect to the z axis)
done ..maybe make a simpler class rsQuantumBit that deals exclusively with the up/down direction
and doesn't consider the others at all - it should also not call the states "up" and "down" but
|1> and |0> respectively


References:
  (1) The Theoretical Minimum - Quantum Mechanics (Leonard Susskind, Art Friedman)

*/



// forward declarations:

//template<class T>
//class rsQuantumSpin;
//
//template<class T>
//inline std::complex<T> operator*(const rsQuantumSpin<T>& B, const rsQuantumSpin<T>& A);

template<class T>
class rsQuantumSpin
{

public:

  /** \name Construction */

  /** Default constructor. Creates a qubit in a pure "up" state. */
  rsQuantumSpin() { prepareUpState(); }

  /** Constructor to create a spin object with given up and down components. It does not verify, if
  these components specify a valid state. ...maybe do an assert... */
  rsQuantumSpin(const std::complex<T>& upComponent, const std::complex<T>& downComponent) 
  { 
    au = upComponent;
    ad = downComponent;
  }




  static rsQuantumSpin<T> up()    { rsQuantumSpin<T> s; s.prepareUpState();    return s; }
  static rsQuantumSpin<T> down()  { rsQuantumSpin<T> s; s.prepareDownState();  return s; }
  static rsQuantumSpin<T> right() { rsQuantumSpin<T> s; s.prepareRightState(); return s; }
  static rsQuantumSpin<T> left()  { rsQuantumSpin<T> s; s.prepareLeftState();  return s; }
  static rsQuantumSpin<T> in()    { rsQuantumSpin<T> s; s.prepareInState();    return s; }
  static rsQuantumSpin<T> out()   { rsQuantumSpin<T> s; s.prepareOutState();   return s; }



  /** \name Setup */

  void prepareUpState()    { au = 1; ad =  0;   }
  void prepareDownState()  { au = 0; ad =  1;   }
  void prepareRightState() { au = s; ad =  s;   }  // (1) Eq 2.5
  void prepareLeftState()  { au = s; ad = -s;   }  // (1) Eq 2.6
  void prepareInState()    { au = s; ad =  s*i; }  // (1) Eq 2.10
  void prepareOutState()   { au = s; ad = -s*i; }  // (1) Eq 2.10


  void setState(const std::complex<T>& newUpComponent, const std::complex<T>& newDownComponent)
  {
    au = newUpComponent;
    ad = newDownComponent;
  }
  // maybe it should automatically (optionally) renormalize the state?

  /*
  void randomizeState() 
  {
    au = std::polar(1, 2*PI*prng->getSample());
    ad = std::polar(1, 2*PI*prng->getSample());
    normalizeState();
  }
  */
  // maybe have an amount parameter between 0..1 - linearly interpolate between current state and
  // random new state - may be used to simulate decoherence

  /** Normalizes the state such that the total probability is unity - which it must be for a valid 
  state. */
  void normalize()
  {
    T r = sqrt(T(1) / getTotalProbability(*this));
    au *= r;
    ad *= r;
  }


  /** Sets a a pointer to a pseudo random number generator that is used in measurement operations.
  These operations destroy superposition, i.e. put the spin into on of two possible pure states. 
  Which of the two that is, is determined randomly according to the probabilities of the two 
  states. That random decision is what this PRNG is used for. Client code should make sure that the
  PRNG is correctly set up to produce numbers in the range 0...1. */
  void setRandomGenerator(rsNoiseGenerator<T>* newGenerator) { prng = newGenerator; }


  /** \name Inquiry */

  std::complex<T> getUpComponent()   const { return au; }
  std::complex<T> getDownComponent() const { return ad; }
  // todo: getLeft/Right/In/Out Component

  /** Returns the squared norm (or magnitude, length, radius) of a complex number. */
  static T getSquaredNorm(const std::complex<T>& z)
  {
    return z.real()*z.real() + z.imag()*z.imag(); // == conj(z) * z, (1) page 39
  }

  /** Returns the total probability for given ket A, i.e. the probability to be in any state at 
  all - which must, of course, always return unity for a valid state. The function can be used for 
  sanity checks and/or to (re)normalize random states. */
  static T getTotalProbability(const rsQuantumSpin& A)
  {
    return getSquaredNorm(A.getUpComponent()) + getSquaredNorm(A.getDownComponent()); // (1) Eq 2.4
  }


  /** Computes the up component of the given ket/state |A>. */
  static std::complex<T> getUpComponent(const rsQuantumSpin& A)
  {
    rsQuantumSpin u;
    u.prepareUpState();
    return u*A;         // (1), Eq 2.1
  }

  static std::complex<T> getDownComponent(const rsQuantumSpin& A)
  { rsQuantumBit d; u.prepareDownState(); return d*A; }

  // make similar functions for left,right,in,out components




  /** Returns the probability for the given state A to be measured in "up" configuration. */
  static T getUpProbability(const rsQuantumSpin& A)
  {
    rsQuantumSpin<T> u;
    u.prepareUpState();
    std::complex<T> r = (A*u) * (u*A); // (1), Eq 2.2
    return r.real();                   // imag should be zero
  }

  /** I'm not sure, if this is correct...this is supposed to return the probability for the given 
  state A to be measured in the target configuration t (i use this to generalize the 
  getUpProbability for getRightProbability and getInProbability)...but for other target states t,
  there may not even be measurement functiosn...so...dunno....todo: figure this stuff out. */
  static T getStateProbability(const rsQuantumSpin& A, const rsQuantumSpin& t)
  {
    std::complex<T> r = (A*t) * (t*A); // (1), Eq 2.2 - i hope, it generalizes the right way
    return r.real();                   // imag should be zero
  }

  /** Returns the probability for the given state A to be measured in "right" configuration. */
  static T getRightProbability(const rsQuantumSpin& A)
  {
    return getStateProbability(A, right());
  }

  /** Returns the probability for the given state A to be measured in "in" configuration. */
  static T getInProbability(const rsQuantumSpin& A)
  {
    return getStateProbability(A, in());
  }


  //rsQuantumSpin<T> r = right();


  // maybe have const static state members up,down,left,right,in,out and a general 
  // getStateComponent, getStateProbability functions and implelement
  // getUpProbability by using the general getStateProbability(A, targetState)

  // have a function to convert to bra - this is a complex conjugation of au, ad and also
  // turns the column vector into a row vector


  /** \name Operations */

  //void rotateUpComponent, rotateDownComponent, applyHadamard, etc.
  // measureUpState - should use p = getUpProbability(*this) and then put it into pure up state 
  // (with probability p) or down state (with probability 1-p) - the measurement destroys the 
  // superposition - the function probably get a pointer to a PRNG ...or maybe the class should 
  // have a PRNG pointer as member

  /** Applies a measurement operation to the state. This measurement will put the state vector 
  either into a pure "up" or pure "down" state and will return +1 in the former and -1 in the 
  latter case. Which one of the two it is is selected randomly (using our prng) according to the
  up-probability of our state. */
  T measureUpComponent()
  {
    //T Pu  = getUpProbability(*this); // optimize this!
    //T Pu  = au.real()*au.real() + au.imag()*au.imag(); // == conj(au) * au, (1) page 39
    T Pu = getSquaredNorm(au);
    T rnd = prng->getSample();
    if(rnd <= Pu) {       // should it be <= or < ?
      prepareUpState();
      return +1; }
    else {
      prepareDownState();
      return -1; }
  }

  T measureRightComponent()
  {
    //rsQuantumSpin<T> r = right();

    T Pr  = getRightProbability(*this);



  }


  // have similar measureRightComponent, measureInComponent functions and also for arbitrary 
  // observables described by a (Hermitian) 2x2 matrix of complex numbers (have a class 
  // "rsSpinOperator" for such matrices) - observables correspond to such matrices and the outcome
  // of any measurement will be one of the eigenvalues of the matrix - and if the state is an 
  // eigenvector, the measurement *will* be the corresponding eigenvalue (having an the eigenvector 
  // as state lets the probability for measuring the corresponding eigenvalue become one)
  // ...so that means, |u> and |d> are eigenvectors with eigenvalues +1 and -1 ...but what is the
  // matrix representing the observable? see (1) Eq 3.17 or 3.20 - the matrices seem to be the 
  // Pauli matrices (so they are not meant as operations to manipulate the state?)

  // todo: implement assignemnt-operator and copy-constructor - these should assign the prng
  // or maybe they should...yeah - i think, that's more convenient

protected:

  std::complex<T> au, ad;
    // our state consisting of the coefficients for up and down spin basis vectors

  rsNoiseGenerator<T>* prng = nullptr;
    // a pointer to a pseudo random number generator that is used in measurement operations which
    // destroy superposition, i.e. put the spin into a pure state - but which of the two possible
    // pure stats that is, is determined randomly - that's what this prng is used for




  static const T s;                // 1/sqrt(2)
  static const std::complex<T> i;  // imaginary unit
    // for convenience (we need these a lot)
};

template<class T>
const T rsQuantumSpin<T>::s = T(1) / sqrt(T(2));

template<class T>
const std::complex<T> rsQuantumSpin<T>::i = std::complex<T>(0, 1);


// define operators +,-

/** Computes the inner product between two states. Both states B,A are assumed to be ket
vectors - the operation of taking the inner product involves turning the first ket into a 
bra first (by complex conjugation of the au, ad coeffs) and then computing the sum of the 
products of corresponding elements. The important point is that you don't need to turn the
ket into a bra before using this - this is done internally by this operator. */
template<class T>
inline std::complex<T> operator*(const rsQuantumSpin<T>& B, const rsQuantumSpin<T>& A)
{
  return conj(B.getUpComponent())   * A.getUpComponent() 
       + conj(B.getDownComponent()) * A.getDownComponent(); // (1), page 20 ff
}
// 
// interchanging arguments leads to complex conjugation of the result
// maybe turn this into a function braKet(bra, ket)
// ..what abotu outer products? do we need such a thing?
// maybe have rsBra, rsKet classes (maybe as subclasses of some rsRowVector, rsColumnVector 
// classes)

/** Adds two kets. */
template<class T>
inline rsQuantumSpin<T> operator+(const rsQuantumSpin<T>& A, const rsQuantumSpin<T>& B)
{
  return rsQuantumSpin<T>(A.getUpComponent()   + B.getUpComponent(), 
                          A.getDownComponent() + B.getDownComponent() );
}

/** Subtracts two kets. */
template<class T>
inline rsQuantumSpin<T> operator-(const rsQuantumSpin<T>& A, const rsQuantumSpin<T>& B)
{
  return rsQuantumSpin<T>(A.getUpComponent()   - B.getUpComponent(), 
                          A.getDownComponent() - B.getDownComponent() );
}

/** Multiplies a scalar and a ket. */
template<class T>
inline rsQuantumSpin<T> operator*(const std::complex<T>& z, const rsQuantumSpin<T>& A)
{
  return rsQuantumSpin<T>(z * A.getUpComponent(), z * A.getDownComponent());
}

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
