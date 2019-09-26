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

Unlike real quantum systems which are like a black box, we can look into the actual state which 
consists of two complex numbers. These two complex numbers are called probability amplitudes and 
their squared magnitudes represent the probabilities that the system will be found in an "up" or 
"down" state when the z-component of the spin is measured. The probability amplitudes of a quantum 
system actually behave totally deterministically - the random element only comes into play when you
do an actual measurement. In order to do so, you need to pass a pointer to a pseudo random number 
generator to the respective measurement function. You should make sure, that this generator is set
up to produce numbers between 0 and 1 with a uniform probability distribution.



...

To specify any state as a ket vector |A>, we express it as a linear combination of two (somewhat 
arbitrarily) choosen basis ket vectors |u> and |d> for "up" and "down" spin:

  |A> = au*|u> + ad*|d>


maybe rename to rsQuantumSpin (done) - what is used as the quantum "bit" is actually one of spin's 
components, such as the z-component (a pure "up" represents binary 1 and a pure "down" state
represents binary 0. Mixed states (with respect to the z-axis) represent a superposition of 0 and
1 (a pure left or right or in or out state is "mixed" with respect to the z axis)
done ..maybe make a simpler class rsQuantumBit that deals exclusively with the up/down direction
and doesn't consider the others at all - it should also not call the states "up" and "down" but
|1> and |0> respectively - it's more abstract and needs less "physical" features

References:
  (1) The Theoretical Minimum - Quantum Mechanics (Leonard Susskind, Art Friedman)  */

template<class T> class rsSpinOperator; // forward declaration

template<class T>
class rsQuantumSpin : public rsVector2D<std::complex<T>>
{

public:

  /** \name Construction */

  /** Default constructor. Creates a qubit in a pure "up" state. */
  rsQuantumSpin() { prepareUpState(); }

  /** Constructor to create a spin object with given up and down components. It does not verify, if
  these components specify a valid state. ...maybe do an assert... */
  rsQuantumSpin(const std::complex<T>& upComponent, const std::complex<T>& downComponent) 
  { 
    y = upComponent;
    x = downComponent;
  }

  /** Creates a spin object in pure "up" state. */
  static rsQuantumSpin<T> up()    { rsQuantumSpin<T> s; s.prepareUpState();    return s; }
  static rsQuantumSpin<T> down()  { rsQuantumSpin<T> s; s.prepareDownState();  return s; }
  static rsQuantumSpin<T> right() { rsQuantumSpin<T> s; s.prepareRightState(); return s; }
  static rsQuantumSpin<T> left()  { rsQuantumSpin<T> s; s.prepareLeftState();  return s; }
  static rsQuantumSpin<T> in()    { rsQuantumSpin<T> s; s.prepareInState();    return s; }
  static rsQuantumSpin<T> out()   { rsQuantumSpin<T> s; s.prepareOutState();   return s; }


  /** \name Setup */

  /** Puts the system into a pure "up" state. */
  void prepareUpState()    { y = 1; x =  0;   }
  void prepareDownState()  { y = 0; x =  1;   }
  void prepareRightState() { y = s; x =  s;   }  // (1) Eq 2.5
  void prepareLeftState()  { y = s; x = -s;   }  // (1) Eq 2.6
  void prepareInState()    { y = s; x =  s*i; }  // (1) Eq 2.10
  void prepareOutState()   { y = s; x = -s*i; }  // (1) Eq 2.10

  /** Assigns the coefficients (a.k.a. probability amplitudes) for "up" and "down" state to the 
  given values. It does not verify, if the numbers represent a valid state. */
  void setState(const std::complex<T>& newUpComponent, const std::complex<T>& newDownComponent)
  {
    y = newUpComponent;
    x = newDownComponent;
  }

  /** Sets this spin object into a given state copied from another spin object. */
  void setState(const rsQuantumSpin<T>& newState)
  {
    y = newState.getUpComponent();
    x = newState.getDownComponent();
  }

  /** Randomizes the state.... */
  void randomizeState(rsNoiseGenerator<T>* prng);
  // needs nore tests - especially for the phase range
  // move to cpp file
  // maybe have an amount parameter between 0..1 - linearly interpolate between current state and
  // random new state - may be used to simulate decoherence

  /** Normalizes the state such that the total probability is unity - which it must be for a valid 
  state. */
  void normalize()
  {
    T r = sqrt(T(1) / getTotalProbability(*this));  // or 1/sqrt(t) instead of sqrt(1/t) - which one is better numerically?
    y *= r;
    x *= r;
  }


  /** \name Inquiry */

  std::complex<T> getUpComponent()   const { return y; }
  std::complex<T> getDownComponent() const { return x; }
  // todo: getLeft/Right/In/Out Component - but these require more compilcated calculations
  // maybe rename all the "Component" functions to "Amplitude"

  /** Returns the squared norm (or magnitude, length, radius) of a complex number. */
  static T getSquaredNorm(const std::complex<T>& z)
  {
    return z.real()*z.real() + z.imag()*z.imag(); // == conj(z) * z, (1) page 39
  }
  // move to RAPT as rsSquaredNorm

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

  // make similar functions for left,right,in,out components and a general
  // getStateComponent


  /** Returns the probability to measure a target state t when a system is in state A. */
  static T getStateProbability(const rsQuantumSpin& A, const rsQuantumSpin& t)
  {
    std::complex<T> r = (A*t) * (t*A); // (1), Eq 3.11 (with lambda_i replaced by t)
    return r.real();                   // imag should be zero
  }



  /** Returns the probability for the given state A to be measured in "down" configuration. */
  static T getDownProbability(const rsQuantumSpin& A) { return getStateProbability(A, down()); }

  /** Returns the probability for the given state A to be measured in "up" configuration. */
  static T getUpProbability(const rsQuantumSpin& A)
  {
    rsQuantumSpin<T> u;
    u.prepareUpState();
    std::complex<T> r = (A*u) * (u*A); // (1), Eq 2.2  ...can probably be optimized
    return r.real();                   // imag should be zero
  }

  /** Returns the probability for the given state A to be measured in "left" configuration. */
  static T getLeftProbability(const rsQuantumSpin& A) { return getStateProbability(A, left()); }

  /** Returns the probability for the given state A to be measured in "right" configuration. */
  static T getRightProbability(const rsQuantumSpin& A) { return getStateProbability(A, right()); }


  /** Returns the probability for the given state A to be measured in "out" configuration. */
  static T getOutProbability(const rsQuantumSpin& A) { return getStateProbability(A, out());  }

  /** Returns the probability for the given state A to be measured in "in" configuration. */
  static T getInProbability(const rsQuantumSpin& A) { return getStateProbability(A, in());  }



  /** Tests, if the state A is close to the state of "this" with the given tolerance. */
  bool isCloseTo(const rsQuantumSpin& A, T tol)
  {
    if(rsAbs(A.y-y) <= tol && rsAbs(A.x-x) <= tol)
      return true;
    return false;
  }

  // have a function to convert to bra - this is a complex conjugation of v.y, ad and also
  // turns the column vector into a row vector


  /** \name Spin Operators */

  //void rotateUpComponent, rotateDownComponent, applyHadamard, etc.
  // measureUpState - should use p = getUpProbability(*this) and then put it into pure up state 
  // (with probability p) or down state (with probability 1-p) - the measurement destroys the 
  // superposition - the function probably get a pointer to a PRNG 

  /** \name Measurements */

  /** Measures the obervable variable that is associated with the given operator M. The result of
  this measurement will be one of the eigenvalues of M and after the measurement, the spin will be
  in a state given by the eigenvector that corresponds to the returned eigenvalue. The operator M 
  is supposed to be Hermitian (i.e. equal to itself transposed and conjugated): M = M^H which 
  implies its eigenvalues to be real numbers. It is also worth noting that, if the spin is in a 
  state of an eigenvector of M before the measurement, it will be guaranteed that the corresponding
  eigenvalue will result in the measurement (its probability becomes one). This, together with the 
  fact that the act measurement will put the system in an eigenvector state of M, implies that 
  subsequent measurements of the same observable will always give the same result (assuming, of 
  course, that no manipulations of the state take place in between the measurements). */
  T measureObservable(const rsSpinOperator<T>& M, rsNoiseGenerator<T>* prng); 

  /** Measures the spin along an arbitrary axis in 3D space given by the vector components
  nx, ny, nz. It is assumed that these components form a unit-length vector in 3-space. */
  T measureSpin(T nx, T ny, T nz, rsNoiseGenerator<T>* prng);
  // not yet tested

  /** Measures the z-component of the spin. This measurement will put the state vector either into
  a pure "up" or pure "down" state and will return +1 in the former and -1 in the latter case. 
  Which one of the two it is is selected randomly using the prng according to the up-probability 
  of our state. The operator/matrix that corresponds to that measurement is the 
  Pauli matrix sigma_z = [[1 0], [0,-1]]. */
  T measureSpinZ(rsNoiseGenerator<T>* prng);

  /** Similar to measureSpinZ, but for the x-component represented by the Pauli matrix 
  sigma_x = [[0,1], [1,0]]. */
  T measureSpinX(rsNoiseGenerator<T>* prng);

  /** Similar to measureSpinZ, but for the y-component represented by the Pauli matrix 
  sigma_y = [[0,-i], [i,0]]. */
  T measureSpinY(rsNoiseGenerator<T>* prng);



protected:

  //std::complex<T> au; //, ad;  // probability amplitudes for "up" and "down" - maybe rename to u,d
                           // coefficients for up and down spin basis vectors


  //rsVector2D<std::complex<T>> v; // for transitioning v.x == ad, v.y == au

  static const T s;                // 1/sqrt(2)
  static const std::complex<T> i;  // imaginary unit
    // for convenience (we need these a lot)

  friend class rsSpinOperator<T>;
};

template<class T> const T rsQuantumSpin<T>::s = T(1) / sqrt(T(2));
template<class T> const std::complex<T> rsQuantumSpin<T>::i = std::complex<T>(0, 1);

/** Computes the inner product between two states. Both states B,A are assumed to be ket
vectors - the operation of taking the inner product involves turning the first ket into a 
bra first (by complex conjugation of the au, ad coeffs) and then computing the sum of the 
products of corresponding elements. The important point is that you don't need to turn the
ket into a bra before using this - this is done internally by this operator. */
template<class T>
inline std::complex<T> operator*(const rsQuantumSpin<T>& B, const rsQuantumSpin<T>& A)
{
  return conj(B.getUpComponent())   * A.getUpComponent() 
       + conj(B.getDownComponent()) * A.getDownComponent(); // (1), page 30 ff
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
An operator is represented as complex valued 2x2 matrix M. There are two very distinct things that 
such operators are used for:

Firstly, an operator M can "act" on (or be applied to) a quantum state v. This means, that a new 
state is computed as the matrix-vector product w = M*v where v is the old state. Operators of that 
kind must be unitary matrices, because....

Secondly, operators may represent measurable or observable quantities. Measuring the value of an
observable associated with such an operator will put the state into one of the eigenvectors of the 
operator and the result of the measurement will be one of its eigenvalues. Because physical 
measurements must be real numbers, an operator M corresponding to an observable must be Hermitian 
(i.e. M = M^H where M^H denotes the Hermitian transpose (= transpose and conjugate)). This ensures 
real eigenvalues. 

Note that the act of setting the spin into an eigenstate of a measurement operator is *not*
the same thing as forming the matrix-vector product like it is done with the first kind of 
operator. Note also that it is only these measurements that involve setting the quantum state
into a randomly chosen one. Operations of the first kind act deterministically on the state 
consisting of the probability amplitudes.  

maybe factor out a baseclass rsMatrix2x2 - or maybe represent spins and operators generally as
rsVector2D and rsMatrix2x2 objects and all the functions and operatiions specific to quantum
stuff should be implemented procedurally - yes - that sounds like a better idea
*/

template<class T>
class rsSpinOperator : public rsMatrix2x2<std::complex<T>> // maybe rename to rsQuantumSpinOperator
{

public:

  /** Stadard constructor. You can pass the matrix elements. if you pass nothing, an identity 
  matrix will be created. */
  rsSpinOperator(std::complex<T> a = T(1), 
                 std::complex<T> b = T(0), 
                 std::complex<T> c = T(0), 
                 std::complex<T> d = T(1))
  {
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
  }

  static rsSpinOperator<T> pauliZ() { rsSpinOperator<T> z; z.setToPauliZ(); return z; }
  static rsSpinOperator<T> pauliX() { rsSpinOperator<T> x; x.setToPauliX(); return x; }
  static rsSpinOperator<T> pauliY() { rsSpinOperator<T> y; y.setToPauliZ(); return y; }
  // (1) says (on page 80 in the footnote) that these Pauli matrices together with the identity 
  // matrix are the quaternions - figure out what that means


  /** \name Setup */

  //void setValues(std::complex<T> a, std::complex<T> b, std::complex<T> c, std::complex<T> d)
  //{
  //  this->a = a; this->b = b; this->c = c; this->d = d;
  //}

  /** Measurement operator for spin along the z-axis. Returns +1 for up, -1 for down. */
  void setToPauliZ() { a = T(1); b = T(0); c = T(0); d = T(-1); }

  /** Measurement operator for spin along the x-axis. Returns +1 for right, -1 for left. */
  void setToPauliX() { a = T(0); b = T(1); c = T(1); d = T(0);  }

  /** Measurement operator for spin along the y-axis. Returns +1 for in, -1 for out. */
  void setToPauliY() { a = T(0); b = -i;   c = i;    d = T(0);  }



  /** \name Inquiry */

  //std::complex<T> getEigenvalue1() const { return rsEigenvalue2x2_1(a, b, c, d); }
  //std::complex<T> getEigenvalue2() const { return rsEigenvalue2x2_2(a, b, c, d); }

  /** Returns the first eigenvalue of this operator. */
  //std::complex<T> getEigenvalue1() const { return RAPT::rsLinearAlgebra::eigenvalue2x2_1(a, b, c, d); }

  /** Returns the second eigenvalue of this operator. */
  //std::complex<T> getEigenvalue2() const { return RAPT::rsLinearAlgebra::eigenvalue2x2_2(a, b, c, d); }

  /** Returns the first eigenvector of this operator. */
  rsQuantumSpin<T> getEigenvector1() const
  {
    std::complex<T> vx, vy;
    RAPT::rsLinearAlgebra::eigenvector2x2_1(a, b, c, d, &vx, &vy, true);
    return rsQuantumSpin<T>(vx, vy);
  }

  /** Returns the second eigenvector of this operator. */
  rsQuantumSpin<T> getEigenvector2() const
  {
    std::complex<T> vx, vy;
    RAPT::rsLinearAlgebra::eigenvector2x2_2(a, b, c, d, &vx, &vy, true);
    return rsQuantumSpin<T>(vx, vy);
  }

  /** Returns the expectation value for the observable M when the system is in state A. */
  static T getExpectedMeasurement(const rsSpinOperator<T>& M, const rsQuantumSpin<T>& A);


  // computes the commutator between two operators
  //static rsSpinOperator<T> commutator(const rsSpinOperator<T>& L, const rsSpinOperator<T>& M);




  /** Access function (read/write) for the matrix elements. The indices i,j can both be 0 or 1. */
  //inline std::complex<T>& operator()(const int i, const int j) { return m[i][j]; }








  /** Applies this quantum spin operator to the given ket v and returns the resulting ket. */
  rsQuantumSpin<T> operator*(const rsQuantumSpin<T>& v) const
  {
    rsQuantumSpin<T> r;
    r.y = a * v.y  +  b * v.x;
    r.x = c * v.y  +  d * v.x;
    return r;
  }



  //std::complex<T> a, b, c, d; // matrix coefficients |a b|
  //                            //                     |c d|
  // todo: either use rsMatrix2x2 or derive from rsMatrix2x2


  //rsMatrix2x2<std::complex<T>> M; // for transition/refactoring - later maybe derive from this class

protected:




  static const std::complex<T> i;  // imaginary unit
};

template<class T>
const std::complex<T> rsSpinOperator<T>::i = std::complex<T>(0, 1);

/** Multiplies a scalar and an operator. */
template<class T>
inline rsSpinOperator<T> operator*(const std::complex<T>& z, const rsSpinOperator<T>& A)
{
  return rsSpinOperator<T>(z * A.a, z * A.b, z*A.c, z*A.d);
}