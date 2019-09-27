#pragma once

/** Implements functions that deal with 2-dimensional complex vectors that represent the quantum 
state of a single spin (for example of an electron) and 2x2 complex matrices that represent 
operators on such state vectors. This is the most simple and prototypical quantum system and can 
also be used as quantum bit (qubit). The user can set the system into various predefined states and
has functions to manipulate the state. The state can also be measured in which case the system 
will - with a certain probability determined by the state - fall into one of two possible pure 
states corresponding to the measured variable. An example of such a measured observable is the spin
along the z-axis. The measured value will be either +1 or -1 (corresponding to "up" or "down") with
probabilities determined by the current state. After the measurement, however, this state will have 
been changed into a pure "up" or "down" state such that subsequent measurements of the z-component 
of the spin will always produce the same result again (with probability one).


States:

To specify any state as a ket vector |A>, we express it as a linear combination of two (somewhat 
arbitrarily) choosen basis ket vectors |u> = (1,0) and |d> = (0,1) for "up" and "down" spin:

|A> = au * |u> + ad * |d>    (1) Pg 38

where au and ad are the probability amplitudes to find the system und "up" or "down" state when
the z-component of the spin is measured. They can be computed for an arbitrary state A as:

au = <u|A>, ad = <d|A>   (1) Eq 2.1

where the inner product of two ket vectors A,B is defined as:

<A|B> = conj(A.x) * B.x + conj(A.y) * B.y.   (1) Eq 1.2 (with renaming)

which also implies that au = A.x and ad = A.y (this comes out due our choice |u> = (1,0), 
|d> = (0,1)). The numbers au, ad are called probability amplitudes and their squared magnitudes 
represent the probabilities that the system will be found in an "up" or "down" state when the 
z-component of the spin is measured.


Operators:

An operator is represented as complex valued 2x2 matrix M. There are two very distinct things that 
such operators are used for:

Firstly, an operator M can "act" on (or be applied to) a quantum state v. This means, that a new 
state is computed as the matrix-vector product w = M*v where v is the old state. Operators of that 
kind must be unitary matrices (i.e the inverse must be given by the conjugate transpose) because 
that's what the laws of quantum mechanics say (todo: add the *actual* explanation - why do they say
that?)

Secondly, operators may represent measurable or observable quantities. Measuring the value of an
observable associated with such an operator will put the state into one of the eigenvectors of the 
operator (randomly, with probabilities determined by the state) and the result of the measurement 
will be the corresponding eigenvalue. Because physical measurements must be real numbers, an 
operator M corresponding to an observable must be Hermitian (i.e. M = M^H where M^H denotes the 
Hermitian transpose (= transpose and conjugate)). This ensures real eigenvalues. 

Note that the act of setting the spin into an eigenstate of a measurement operator is *not*
the same thing as forming the matrix-vector product like it is done with the first kind of 
operator. Note also that it is only these measurements that involve setting the quantum state
into a randomly chosen one. Operations of the first kind act deterministically on the state.
The random element only comes into play when you do an actual measurement. In order to do so, you 
need to pass a pointer to a pseudo random number generator to the respective measurement function. 
You should make sure, that this generator is set up to produce numbers between 0 and 1 with a 
uniform probability distribution (TODO: this is a huge source of trouble - it's so easy to forget
and by default, my noise-generators produce values in -1...+1 and then one gets randomly wrong 
measurement results - use a specific kind of PRNG that can only produce numbers in 0..1)

References:
(1) The Theoretical Minimum - Quantum Mechanics (Leonard Susskind, Art Friedman) */

template<class T>
class rsQuantumSpin
{

public:

  typedef std::complex<T> Complex;
  typedef rsVector2D<std::complex<T>>  Vec;
  typedef rsMatrix2x2<std::complex<T>> Mat;
  typedef rsNoiseGenerator<T> PRNG;

  //-----------------------------------------------------------------------------------------------
  /** \name State setup */

  // pure state creation functions:
  static Vec up()    { Vec s; prepareUpState(s);    return s; }
  static Vec down()  { Vec s; prepareDownState(s);  return s; }
  static Vec right() { Vec s; prepareRightState(s); return s; }
  static Vec left()  { Vec s; prepareLeftState(s);  return s; }
  static Vec in()    { Vec s; prepareInState(s);    return s; }
  static Vec out()   { Vec s; prepareOutState(s);   return s; }

  static void prepareUpState(   Vec& A) { A.x = 1; A.y =  0;   }  // (1) Eq 2.11
  static void prepareDownState( Vec& A) { A.x = 0; A.y =  1;   }  // (1) Eq 2.12
  static void prepareRightState(Vec& A) { A.x = s; A.y =  s;   }  // (1) Eq 2.5
  static void prepareLeftState( Vec& A) { A.x = s; A.y = -s;   }  // (1) Eq 2.6
  static void prepareInState(   Vec& A) { A.x = s; A.y =  s*i; }  // (1) Eq 2.10
  static void prepareOutState(  Vec& A) { A.x = s; A.y = -s*i; }  // (1) Eq 2.10


  static void prepareState(Vec& A, T ur, T ui, T dr, T di)
  {
    A.x.real(ur); A.x.imag(ui); A.y.real(dr); A.y.imag(di);
  }

   /** Normalizes the state such that the total probability is unity - which it must be for a valid 
  state. */
  static void normalizeState(Vec& A);

  /** Randomizes the state....todo: allow to specify an amount between 0 and 1 - can be used to 
  inject noise to simulate quantum decoherence */
  static void randomizeState(Vec& v, PRNG* prng);


  //-----------------------------------------------------------------------------------------------
  /** \name Operator setup */

  /** Measurement operator for spin along the x-axis. Returns +1 for right, -1 for left. */
  static void setToPauliX(Mat& M) { M.a = T(0); M.b = T(1); M.c = T(1); M.d = T(0);  }

  /** Measurement operator for spin along the y-axis. Returns +1 for in, -1 for out. */
  static void setToPauliY(Mat& M) { M.a = T(0); M.b = -i;   M.c = i;    M.d = T(0);  }

  /** Measurement operator for spin along the z-axis. Returns +1 for up, -1 for down. */
  static void setToPauliZ(Mat& M) { M.a = T(1); M.b = T(0); M.c = T(0); M.d = T(-1); }
  // turn these into factory functions


  static Mat pauliX() { return Mat(T(0), T(1), T(1), T(0)); }
  static Mat pauliY() { return Mat(T(0),   -i,    i, T(0)); }
  static Mat pauliZ() { return Mat(T(1), T(0), T(0), T(1)); }
  // not yet tested

  // see here for more operators/"gates"
  // https://en.wikipedia.org/wiki/Quantum_logic_gate

  // for the multi-qubit gates, don't craete the matrices but instead directly operate on a 
  // number of input vectors like:
  // static void toffoli(Vec& v1, Vec& v2, Vec& v3);
  // maybe this should be done in another class rsQuantumGates


  /** Creates the Pauli vector which is the 3-vector of the 2x2 Pauli matrices. see (1) pg 83 or
  https://en.wikipedia.org/wiki/Pauli_matrices#Pauli_vector - i think, on wikipedia, the 
  (x,y,z)-hat vectors are simply the unit vectors in x,y,z directions of 3-space? */
  static rsVector3D<rsMatrix2x2<Complex>> pauliVector();
  // maybe move into a special section of factory functions

  // todo: void rotate(rotX, rotY), hadamard, cnot, etc


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Computes the inner product ("bracket") <A|B> of two "ket" vectors A,B. The left "ket" is 
  converted into a "bra" vector first (by complex conjugation and transposition). 
  interchanging arguments leads to complex conjugation of the result (verify)  */
  static std::complex<T> bracket(const Vec& A, const Vec& B) 
  { return conj(A.x) * B.x + conj(A.y) * B.y; }

  /** Returns the probability amplitude to measure an "up" configuration when z-spin is
  measured and the system is in state A. If we denote this amplitude by au, it is given by
  au = <u|A> = (u1, u2)^H * (A.x, A.y) = (1, 0)^H * (A.x, A.y)= A.x. 
  (see (1) Eq 2.1 with 2.11 for |u>) */
  static std::complex<T> getUpAmplitude(const Vec& A) { return A.x; }
  
  /** Returns the probability amplitude to measure a "down" configuration when z-spin is
  measured and the system is in state A. It's given by
  ad = <d|A> = (d1, d2)^H * (A.x, A.y) = (0, 1)^H * (A.x, A.y)= A.y  */
  static std::complex<T> getDownAmplitude(const Vec& A) { return A.y; }

  /** Returns the probability to measure a target state t when a system is in state A. */
  static T getStateProbability(const Vec& A, const Vec& t);

  /** Returns the probability for the given state A to be measured in "down" configuration. */
  static T getUpProbability(   const Vec& A) { return getStateProbability(A, up());    } 
  static T getDownProbability( const Vec& A) { return getStateProbability(A, down());  }
  static T getRightProbability(const Vec& A) { return getStateProbability(A, right()); }
  static T getLeftProbability( const Vec& A) { return getStateProbability(A, left());  }
  static T getInProbability(   const Vec& A) { return getStateProbability(A, in());    }
  static T getOutProbability(  const Vec& A) { return getStateProbability(A, out());   }

  /** Returns the total probability for given ket A, i.e. the probability to be in any state at 
  all - which must, of course, always return unity for a valid state. The function can be used for 
  sanity checks and/or to (re)normalize random states. */
  static T getTotalProbability(const Vec& A) { return rsAbsSquared(A.x) + rsAbsSquared(A.y); }

  /** Returns the expectation value for the observable M when the system is in state A. */
  static T getExpectedMeasurement(const Mat& M, const Vec& A);


  //-----------------------------------------------------------------------------------------------
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
  static T measureObservable(Vec& A, const Mat& M, rsNoiseGenerator<T>* prng); 

  /** Measures the spin along an arbitrary axis in 3D space given by the vector components
  nx, ny, nz. It is assumed that these components form a unit-length vector in 3-space
  (todo: relax that assumption - i.e. normalize internally, if necessarry).   */
  static T measureSpin(Vec& A, T nx, T ny, T nz, rsNoiseGenerator<T>* prng);
  // not yet tested

  /** Measures the z-component of the spin. This measurement will put the state vector either into
  a pure "up" or pure "down" state and will return +1 in the former and -1 in the latter case. 
  Which one of the two it is is selected randomly using the prng according to the up-probability 
  of our state. The operator/matrix that corresponds to that measurement is the 
  Pauli matrix sigma_z = [[1 0], [0,-1]]. */
  static T measureSpinZ(Vec& A, rsNoiseGenerator<T>* prng);

  /** Similar to measureSpinZ, but for the x-component represented by the Pauli matrix 
  sigma_x = [[0,1], [1,0]]. */
  static T measureSpinX(Vec& A, rsNoiseGenerator<T>* prng);

  /** Similar to measureSpinZ, but for the y-component represented by the Pauli matrix 
  sigma_y = [[0,-i], [i,0]]. */
  static T measureSpinY(Vec& A, rsNoiseGenerator<T>* prng);


  static const T s;                // 1/sqrt(2)
  static const std::complex<T> i;  // imaginary unit
};
template<class T> const T rsQuantumSpin<T>::s = T(1) / sqrt(T(2));
template<class T> const std::complex<T> rsQuantumSpin<T>::i = std::complex<T>(0, 1);



//=================================================================================================

/**


References:
 (1) Simulating Quantum Computers Using OpenCL: https://arxiv.org/pdf/1805.00988.pdf
     Source code: https://github.com/libtangle/qcgpu

*/

template<class T>
class rsQuantumComputer
{

public:

  typedef std::complex<T> Complex;
  typedef rsVector2D<std::complex<T>> Vec;
  typedef rsMatrix2x2<std::complex<T>> QGate;

  rsQuantumComputer() { allocateMemory(); }


  void applyGate(const QGate& g, int bitIndex);

  /** Returns the nth number where a given digit is cleared in the binary representation of the 
  number. */
  static int nth_cleared(int n, int target)  
  {
    int mask = (1 << target ) - 1;
    int not_mask = ~mask; 
    return (n & mask) | ((n & not_mask ) << 1);
  }



protected:

  void allocateMemory() { qbits.resize(numStates); }

  int numQBits  = 4;
  int numStates = 16; // = 2^numQBits
  std::vector<Vec> qbits;  // these are actually the states

};