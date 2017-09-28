#ifndef RAPT_PARTICLESYSTEM_H_INCLUDED
#define RAPT_PARTICLESYSTEM_H_INCLUDED

/** Class for representing particles that can move in 3D space. Each particle is characterized at 
each time-instant by its current position and velocity (time varying) and its mass and charge 
(constant). */

template<class T>
class rsParticle
{

public:

  /** Returns the kinetic energy of this particle */
  T getKineticEnergy() { return T(0.5) * mass * vel.getSquaredEuclideanNorm(); }

  /** Returns the momentum of this particle. */
  rsVector3D<T> getMomentum() { return mass * vel; }

  rsVector3D<T> pos;  // position
  rsVector3D<T> vel;  // velocity
  T mass = 1;
  T charge = 1;
};

//=================================================================================================

/** A class for simulating a system of particles that interact via gravitational, electric and
magnetic forces. 

References:
 -(1) The Feynman Lectures on Physics, Vol. 1, The New Millenium Edition

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

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  T getKineticEnergy();

  T getPotentialEnergy();

  T getTotalEnergy() { return getKineticEnergy() + getPotentialEnergy(); }

  rsVector3D<T> getTotalMomentum();

  //rsVector3D<T> getTotalAngularMomentum();

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

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
  T cG = 1, cE = 1, cM = 1;
};

#endif