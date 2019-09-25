#pragma once
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

...the design with teh pointer to the PRNG as member is a bit odd and may lead to bugs when the
pointer is unassigned (which may happen easily when you re-assign a a state variable)...come up 
with something better - maybe functions that need to create random numbers should get passed a 
pointer to the prng ...yes - that seems better - it also makes it more obvious from client code
which functions introduce random components

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

template<class T> class rsSpinOperator;

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

  void setState(const rsQuantumSpin<T>& newState)
  {
    au = newState.getUpComponent();
    ad = newState.getDownComponent();
  }


  
  void randomizeState() 
  {
    T Pu = prng->getSample();             // probability of "up"
    T Pd = T(1) - Pu;                     // probability of "down"
    T r  = T(1) / sqrt(Pu*Pu + Pd*Pd);    // normalizer
    Pu  *= r;
    Pd  *= r;
    T pu = T(2.0*PI) * prng->getSample(); // phase of au
    T pd = T(2.0*PI) * prng->getSample(); // phase of ad

    au = std::polar(Pu, pu);
    ad = std::polar(Pd, pd);

    //normalize(); // should already be normalized thanks to *= r
  }
  // needs nore tests - especially for the phase range
  // move to cpp file
  
  // maybe have an amount parameter between 0..1 - linearly interpolate between current state and
  // random new state - may be used to simulate decoherence

  /** Normalizes the state such that the total probability is unity - which it must be for a valid 
  state. */
  void normalize()
  {
    T r = sqrt(T(1) / getTotalProbability(*this));
    // or should we do 1/sqrt(t) instead of sqrt(1/t) - which one is better numerically?

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
  // can this be simplified? we could just call A.getUpComponent - actually, this function is 
  // redundant...

  static std::complex<T> getDownComponent(const rsQuantumSpin& A)
  { rsQuantumSpin<T> d; d.prepareDownState(); return d*A; }

  // make similar functions for left,right,in,out components




  /** Returns the probability to measure a target state t when a system is in state A. */
  static T getStateProbability(const rsQuantumSpin& A, const rsQuantumSpin& t)
  {
    std::complex<T> r = (A*t) * (t*A); // (1), Eq 3.11 (with lambda_i replaced by t)
    return r.real();                   // imag should be zero
  }

  /** Returns the probability for the given state A to be measured in "up" configuration. */
  static T getUpProbability(const rsQuantumSpin& A)
  {
    rsQuantumSpin<T> u;
    u.prepareUpState();
    std::complex<T> r = (A*u) * (u*A); // (1), Eq 2.2  ...can probably be optimized
    return r.real();                   // imag should be zero
  }

  /** Returns the probability for the given state A to be measured in "right" configuration. */
  static T getRightProbability(const rsQuantumSpin& A) { return getStateProbability(A, right()); }

  /** Returns the probability for the given state A to be measured in "in" configuration. */
  static T getInProbability(const rsQuantumSpin& A) { return getStateProbability(A, in());  }

  /** Tests, if the state A is close to this state withij the given tolerance. */
  bool isCloseTo(const rsQuantumSpin& A, T tol)
  {
    if(rsAbs(A.au-au) <= tol && rsAbs(A.ad-ad) <= tol)
      return true;
    return false;
  }



  // maybe have const static state members up,down,left,right,in,out and a general 
  // getStateComponent, getStateProbability functions and implelement
  // getUpProbability by using the general getStateProbability(A, targetState)

  // have a function to convert to bra - this is a complex conjugation of au, ad and also
  // turns the column vector into a row vector


  /** \name Spin Operators */

  //void rotateUpComponent, rotateDownComponent, applyHadamard, etc.
  // measureUpState - should use p = getUpProbability(*this) and then put it into pure up state 
  // (with probability p) or down state (with probability 1-p) - the measurement destroys the 
  // superposition - the function probably get a pointer to a PRNG ...or maybe the class should 
  // have a PRNG pointer as member


  /** \name Measurements */

  /** Measures the obervable variable that is associated with the given operator M. The result of
  this measurement will be one of the eigenvalues of M and after the measurement, the spin will be
  in a state given by the eigenvector that corresponds to the returned eigenvalue. The operator M 
  is supposed to be Hermitian (i.e. equal to itself transposed and conjugated): M = M^H which 
  implies its eigenvalues to be real numbers. It is also worth noting that, if the spin is in a 
  state of an eigenvector of M before the measurement, it will be guaranteed that the corresponding
  eigenvalue will result in the measurement (it's probability becomes one). This, together with the 
  fact that the act measurement will put the system in an eigenvector state of M, implies that 
  subsequent measurements of the same observable will always give the same result (assuming, of 
  course, that no manipulations of the state take place in between the measurements). */
  T measureObservable(const rsSpinOperator<T>& M); 
  // still buggy

  /** Applies a measurement operation to the state. This measurement will put the state vector 
  either into a pure "up" or pure "down" state and will return +1 in the former and -1 in the 
  latter case. Which one of the two it is is selected randomly (using our prng) according to the
  up-probability of our state. The operator/matrix that corresponds to that measurement is the 
  Pauli matrix sigma_z = [[1 0], [0,-1]]. */
  T measureUpComponent();

  /** Like measureUpComponent(), but for the "right" component represented by the Pauli matrix 
  sigma_x = [[0,1], [1,0]].  */
  T measureRightComponent(); 
  // not yet tested

  /** Like measureUpComponent(), but for the "in" component represented by the Pauli matrix 
  sigma_y = [[0,-i], [i,0]].  */

  T measureInComponent();    
  // not yet tested


protected:

  std::complex<T> au, ad;  // maybe rename to u,d
    // our state consisting of the coefficients for up and down spin basis vectors

  rsNoiseGenerator<T>* prng = nullptr;
    // a pointer to a pseudo random number generator that is used in measurement operations which
    // destroy superposition, i.e. put the spin into a pure state - but which of the two possible
    // pure stats that is, is determined randomly - that's what this prng is used for
  // get rid - pass the prng to measurement functions



  static const T s;                // 1/sqrt(2)
  static const std::complex<T> i;  // imaginary unit
    // for convenience (we need these a lot)


  friend class rsSpinOperator<T>;
};

template<class T>
const T rsQuantumSpin<T>::s = T(1) / sqrt(T(2));

template<class T>
const std::complex<T> rsQuantumSpin<T>::i = std::complex<T>(0, 1);


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

/** A class for representing linear operators on quantum spins (objects of class rsQuantumSpin). 
An operator is represented as complex valued 2x2 matrix M and applying the operator to a spin 
state v (represented by a complex valued row-vector, a.k.a. "ket") amounts to computing the 
matrix-vector product M*v. ...

There's a special class of operators associated with observable variables. Measuring a variable 
associated with such an operator will put the state into one of the eigenvectors of the operator 
and the result of the measurement will be one of its eigenvalues. Because physical measurements 
must be real numbers, an operator M corresponding to an observable must be Hermitian (i.e. 
M = M^H where M^H denotes the Hermitian transpose (= transpose and conjugate)). This ensures real
eigenvalues.....


todo: explain the unitarity stuff


*/

// move to RAPT::rsLinearAlagebra:

template<class T>
T rsEigenvalue2x2_1(T a, T b, T c, T d)
{
  return T(0.5) * (a + d - sqrt(a*a + T(4)*b*c - T(2)*a*d + d*d));
}

template<class T>
T rsEigenvalue2x2_2(T a, T b, T c, T d)
{
  return T(0.5) * (a + d + sqrt(a*a + T(4)*b*c - T(2)*a*d + d*d));
}

// move to RAPT::rsLinearAlgebra, make a function that computes bothe eigenvalues at once (we can
// re-use the value of the sqrt)
// Sage code to produce the formulas:
// var("a b c d")
// A = matrix([[a, b], [c, d]])
// A.eigenvalues()

template<class T>
void normalize(T& vx, T& vy)
{
  T rx = rsAbs(vx); rx *= rx;
  T ry = rsAbs(vy); ry *= ry;
  T s  = T(1) / sqrt(rx+ry);
  vx *= s;
  vy *= s;
}

template<class T>
void rsEigenvector2x2_1(T a, T b, T c, T d, T& vx, T& vy)
{
  if(b != T(0)) {
    vx = T(1);
    vy = T(0.5) * (a - d + sqrt(a*a + T(4)*b*c - T(2)*a*d + d*d)) / b; 
    normalize(vx, vy); }
  else {
    vx = T(0);
    vy = T(1); }
}
// ...needs tests

template<class T>
void rsEigenvector2x2_2(T a, T b, T c, T d, T& vx, T& vy)
{
  if(b != T(0)) {
    vx = T(1);
    vy = T(0.5) * (a - d - sqrt(a*a + T(4)*b*c - T(2)*a*d + d*d)) / b; 
    normalize(vx, vy); }
  else {
    if(a != d) {
      vx = T(1);
      vy = c/(a-d); 
      normalize(vx, vy); }
    else {
      vx = T(0);
      vy = T(1);  }} 
}

// the sqrt appears in all 4 formulas - what's its significance? maybe its worth to factor out and
// give it a name? maybe eigenSqrt...or has it to do with the determinant?


// var("a b c d")
// A = matrix([[a, b], [c, d]])
// A.eigenvectors_right()
// [(1/2*a + 1/2*d - 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2), [(1, -1/2*(a - d + sqrt(a^2 + 4*b*c - 2*a*d + d^2))/b)],  1),
//  (1/2*a + 1/2*d + 1/2*sqrt(a^2 + 4*b*c - 2*a*d + d^2), [(1, -1/2*(a - d - sqrt(a^2 + 4*b*c - 2*a*d + d^2))/b)],  1) ]
//
// special case when b=0 (leads to div-by-0 in formula above):
//   var("a b c d")
//   A = matrix([[a, 0], [c, d]])
//   A.eigenvectors_right()
//   [(d, [(0, 1)], 1), (a, [(1, c/(a - d))], 1)]
//
// needs further special case when a=d:
//   var("a b c d")
//   A = matrix([[a, 0], [c, a]])
//   A.eigenvectors_right()
//   [(a, [(0, 1)], 2)]

// these are the right eigenvectors - maybe have similar functions for the left eigenvectors


template<class T>
class rsSpinOperator // maybe rename to rsQuantumSpinOperator
{

public:

  static rsSpinOperator<T> pauliZ() { rsSpinOperator<T> z; z.setToPauliZ(); return z; }
  static rsSpinOperator<T> pauliX() { rsSpinOperator<T> x; x.setToPauliX(); return x; }
  static rsSpinOperator<T> pauliY() { rsSpinOperator<T> y; y.setToPauliZ(); return y; }


  /** \name Setup */

  /** Measurement operator for spin along the z-axis. Returns +1 for up, -1 for down. */
  void setToPauliZ() { a = T(1); b = T(0); c = T(0); d = T(-1); }

  /** Measurement operator for spin along the x-axis. Returns +1 for right, -1 for left. */
  void setToPauliX() { a = T(0); b = T(1); c = T(1); d = T(0);  }

  /** Measurement operator for spin along the y-axis. Returns +1 for in, -1 for out. */
  void setToPauliY() { a = T(0); b = -i;   c = i;    d = T(0);  }



  /** \name Inquiry */

  std::complex<T> getEigenvalue1() const { return rsEigenvalue2x2_1(a, b, c, d); }
  std::complex<T> getEigenvalue2() const { return rsEigenvalue2x2_2(a, b, c, d); }

  rsQuantumSpin<T> getEigenvector1() const
  {
    std::complex<T> vx, vy;
    rsEigenvector2x2_1(a, b, c, d, vx, vy);
    return rsQuantumSpin<T>(vx, vy);
  }

  rsQuantumSpin<T> getEigenvector2() const
  {
    std::complex<T> vx, vy;
    rsEigenvector2x2_2(a, b, c, d, vx, vy);
    return rsQuantumSpin<T>(vx, vy);
  }





  /** Access function (read/write) for the matrix elements. The indices i,j can both be 0 or 1. */
  //inline std::complex<T>& operator()(const int i, const int j) { return m[i][j]; }

  /** Applies this quantum spin operator to the given ket v and returns the resulting ket. */
  rsQuantumSpin<T> operator*(const rsQuantumSpin<T>& v) const
  {
    rsQuantumSpin<T> r;
    r.au = a * v.au  +  b * v.ad;
    r.ad = c * v.au  +  d * v.ad;
    return r;
  }



protected:

  std::complex<T> a, b, c, d; // matrix coefficients |a b|
                              //                     |c d|


  static const std::complex<T> i;  // imaginary unit
};

template<class T>
const std::complex<T> rsSpinOperator<T>::i = std::complex<T>(0, 1);