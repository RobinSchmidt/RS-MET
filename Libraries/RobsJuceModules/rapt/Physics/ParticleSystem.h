#ifndef RAPT_PARTICLESYSTEM_H_INCLUDED
#define RAPT_PARTICLESYSTEM_H_INCLUDED

/** Class for representing particles that can move in 3D space. Each particle is characterized at 
each time-instant by its current position and velocity (time varying) and its mass and charge 
(constant). 

References:
-(1) The Feynman Lectures on Physics, Vol. 1, The New Millenium Edition
-(2) The Feynman Lectures on Physics, Vol. 2, The New Millenium Edition
-(3) Ein Jahr für die Physik (Thomsen), 2. Aufl.
-(4) Classical Mechanics - The Theoretical Minimum (Susskind)

*/

template<class T>
class rsParticle
{

public:

  /** Returns the kinetic energy of this particle */
  T getKineticEnergy() { return T(0.5) * mass * vel.getSquaredEuclideanNorm(); }

  /** Returns the momentum of this particle. */
  rsVector3D<T> getMomentum() { return mass * vel; }

  /** Returns the gravitational field at position vector p caused by this particle. Multiplying a 
  test mass (assumed to be at the given position) with the field gives the force on the test 
  mass. */
  rsVector3D<T> getGravitationalFieldAt(rsVector3D<T> p, T gravitationalConstant = 1);

  /** Returns the electric field at position vector p caused by this particle. */
  rsVector3D<T> getElecricFieldAt(rsVector3D<T> p, T electricConstant = 1);

  /** Returns the magnetic field at position vector p caused by this particle. */
  rsVector3D<T> getMagneticFieldAt(rsVector3D<T> p, T magneticConstant = 1);

  /** Returns the gravitational potential at position vector p due to this particle. The potential 
  is a scalar field whose negative gradient gives the gravitational (vector-) field. Multiplying 
  the potential with a test-mass (assumed to be at the given position) gives the potential energy 
  of the test-mass (...hmm...well...up to a factor of 2 and only in static situations). */
  T getGravitationalPotentialAt(rsVector3D<T> p, T gravitationalConstant = 1);

  /** Returns the electric potential at position vector p due to this particle. The potential is
  a scalar field whose negative gradient gives the electrical field. Multiplying the potential with
  a test-charge (assumed to be at the given position) gives the potential energy of the
  test-charge. */
  T getElectricPotentialAt(rsVector3D<T> p, T electricConstant = 1);

  /** Returns the magnetic potential at position vector p due to this particle. The magnetic 
  potential is a vector field whose curl gives the magnetic field. The vector potential is not 
  uniquely determined - you can add any vector field with zero curl and still get the same magnetic 
  field. Curl-free vector fields can be expressed as gradient of some scalar field, so the vector 
  potential is determined only up to a gradient of a scalar field. Here, we choose the gradient 
  field to be zero. */
  rsVector3D<T> getMagneticPotentialAt(rsVector3D<T> p, T magneticConstant = 1);


  /** Returns the gravitational potential energy of the passed particle p in the field of this 
  particle. */
  T getGravitationalEnergy(const rsParticle<T>& p, T gravitationalConstant = 1);

  /** Returns the electrical potential energy of the passed particle p in the field of this 
  particle. */
  T getElectricEnergy(const rsParticle<T>& p, T electricConstant = 1);

  /** Returns the magnetic potential energy of the passed particle p in the field of this 
  particle. */
  T getMagneticEnergy(const rsParticle<T>& p, T magneticConstant = 1);
  // todo: verify these formulas....they are supposed to hold up only in static situations (see 
  // (2), page 15.7 at the top)...maybe there are also formulas that work for dynamic situations? 
  // ..try to .figure out...see (2), page 15.15


  // todo: write function for magnetic potential (this is vector valued), figure out, if there's 
  // something like magnetic potential energy, have functions to compute potential energies of 
  // test-particles (passed as parameter)
  // see https://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Wiechert_potential
  // or (2), Eq. 15.24 for the magnetic vector potential and Eq. 15.20 for magnetic potential 
  // energy -> replace current density j in these formulas by charge * vel and the integrals can be 
  // ignored because we just look at one single particle (?)...stuff to figure out...

  rsVector3D<T> pos;  // position
  rsVector3D<T> vel;  // velocity
  T mass   = 1;       // todo: maybe distinguish inertial and gravitational mass
  T charge = 1;
  T size   = 1;       // to avoid force singularities
  // maybe have some kind of spin / angular-velocity (a vector)
};

//=================================================================================================

/** A class for simulating a system of particles that interact via gravitational, electric and
magnetic forces. 

todo:
-maybe include a frictional force
-maybe, when running the system in an iteration, regularly check the total energy and momentum
->this should remain constant (unless there is friction) - if it isn't, there might either be a 
bug or numerical errors. in the 2nd case, re-adjust the speeds such that the initial energy
is regained
*/

template<class T>
class rsParticleSystem
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsParticleSystem(int numParticles);

  /** Destructor. */
  //~rsParticleSystem(); 

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setNumParticles(int newNumParticles);

  /** Sets the stepsize for updating the state from one iteration to the next. */
  void setStepSize(T newStepSize) { stepSize = newStepSize; }

  /** Sets the constant used to scale all the gravitational forces. */
  void setGravitationalConstant(T newConstant) { cG = newConstant; }

  /** Sets the constant used to scale all the electric forces. */
  void setElectricConstant(T newConstant) { cE = newConstant; }

  /** Sets the constant used to scale all the magnetic forces. */
  void setMagneticConstant(T newConstant) { cM = newConstant; }

  /** Sets the exponent, by which the force law (asymptotically) depends inversely on the 
  distance. For example, with a value of 2, the force between two particles depends on the distance 
  d between them (asymptocically) proportionally to 1/d^2 - this is the physically correct 
  inverse-square law. The dependence is only asymptotic, because we also have the offsset set by
  setForceLawOffset - if this offset is zero, the dependence is not asymptotic but exact. */
  void setForceLawExponent(T newExponent) { p = newExponent+1; }
    // +1 because our formulas assume that the distance d appears in the numerator, because
    // we don't want to use the formula that uses the normalized distance vector (because the
    // normalization itself can cause a div-by-0)
    // todo: maybe express the exponent non-inverse, i.e. when the user wants an inverse-square
    // law, he should pass -2 instead of 2, then a value of 1 corresponds to Hooke's law
    // -> more convenient than the implicit inversion

  /** Sets an offset for the force-law which is mainly intended to avoid divisions-by-zero which
  would otherwise occur when the distance between two particles becomes zero. The value should be
  strictly positive. If you pass 0 (and the exponent is 2), you'll have the physically correct
  1/d^2 law - but this may lead to divisions by zero. */
  //void setForceLawOffset(T newOffset) { c = newOffset; }

  // makeTotalMomentumZero, makeCenterOfMassZero

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  T getKineticEnergy();

  T getGravitationalPotentialEnergy();

  T getElectricPotentialEnergy();

  T getMagneticPotentialEnergy();

  T getPotentialEnergy(); // this does not work yet

  T getTotalEnergy() { return getKineticEnergy() + getPotentialEnergy(); }

  rsVector3D<T> getTotalMomentum();

  //getCenterOfMass, getTotalMass
  //rsVector3D<T> getTotalAngularMomentum();

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Returns a scaler for the force based on the distance d. The complete distance-to-force law 
  will be the value returned by that function times the distance itself (this is, so we can use 
  formulas with unnormalized difference vectors - which have the distance in the numerator, 
  too). */
  T getForceScalerByDistance(T distance, T size1, T size2);

  /** Returns the force for a given distance. Mainly meant for making plots. */
  T getForceByDistance(T distance, T size1, T size2);

  /** Computes the force that particle p1 experiences due to the presence of particle p2. */
  rsVector3D<T> getForceBetween(const rsParticle<T>& p1, const rsParticle<T>& p2);

  /** Computes the forces that act on each particle and stores them in our "forces" member. */
  void updateForces();

  /** Updates the particle velocities according to the computed forces and the particle masses, 
  using F = m*a -> a = F/m -> v += a. */
  void updateVelocities();

  /** Updates the particle positions according to the current velocities by advancing their 
  positions one step into direction of the respective velocity. */
  void updatePositions();

  /** Updates the state of the system. */
  void updateState()
  {
    updateForces();
    updateVelocities();
    updatePositions();
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the positions and velocities of all particles to their initial values. */
  void reset();

  std::vector<rsVector3D<T>> initialPositions;
  std::vector<rsVector3D<T>> initialVelocities;
  std::vector<rsParticle<T>> particles;

protected:

  std::vector<rsVector3D<T>> forces;
  T stepSize = T(0.01);
  T cG = 1, cE = 1, cM = 1; // gravitational, electric, magnetic constants
  T p = 3;                  // force-by-distance-law exponent/power

};

#endif