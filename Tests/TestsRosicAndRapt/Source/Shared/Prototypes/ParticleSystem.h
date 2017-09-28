#ifndef ParticleSystem_h
#define ParticleSystem_h

#include "rapt/rapt.h"

/** This is just a skeleton/stub/idea for a particle system */



/** Class for representing 3-dimensional vectors. */

template<class T>
class rsVector3D
{

public:

  /** Constructor. Initializes coordinates with the passed values. */
  rsVector3D(T _x = 0, T _y = 0, T _z = 0) : x(_x), y(_y), z(_z) {}

  /** Returns the squared Euclidean norm of this vector. */
  T getSquaredEuclideanNorm() { return x*x + y*y + z*z; }

  /** Returns the Euclidean norm of this vector. */
  T getEuclideanNorm() { return sqrt(getSquaredEuclideanNorm()); }

  /** The 3 cartesian coordinate values. */
  T x, y, z;

};

/** Computes the cross-product between two 3D vectors v and w. Here, v is the left operand. That's 
important, because the cross-product is not commutative. Instead, we have (v*w) = -(w*v) where the
* symbol is used here to denote the cross-product. */
template<class T>
rsVector3D<T> cross(const rsVector3D<T>& v, const rsVector3D<T>& w)
{
  return rsVector3D<T>(v.y*w.z-v.z*w.y, v.z*w.x-v.x*w.z, v.x*w.y-v.y*w.x); 
}



/** A list of physical constants. Physical simulations may keep a pointer to an object of that type
and outlying code may modify this object. ...hmm...not sure if that's a good idea - maybe each 
simulation should keep the constants it actually needs as members.

References:
 -https://physics.nist.gov/cuu/Constants/Table/allascii.txt
*/

template<class T>
class rsPhysicalConstants
{

public:

  T boltzmann   = 1.38064852e-23;   // J K^-1, Boltzmann constant
  T gravitation = 6.67408e-11;      // m^3 kg^-1 s^-2, Newtonian constant of gravitation 
  T lightspeed  = 299792458;        // m s^-1, speed of light in vacuum
  T planck      = 6.626070040e-34;  // J s, Planck constant
  // ...more to come

};

/** Class for representing particles that can move in 3D space. Each particle is characterized at 
each time-instant by its current position and velocity (time varying) and its mass (constant). */

template<class T>
class rsParticle
{

public:

  rsVector3D<T> pos;  // position
  rsVector3D<T> vel;  // velocity
  T mass;
  T charge;

};


/**  */

template<class T>
class rsParticleSystem
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. */
  rsParticleSystem(size_t numParticles);

  /** Destructor. */
  //~rsParticleSystem(); 

  //-----------------------------------------------------------------------------------------------
  // \name setup:



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry:

  /*
  T getKineticEnergy();

  T getPotentialEnergy();

  T getTotalEnergy() { return getKineticEnergy() + getPotentialEnergy(); }

  T getTotalMomentum();

  T getTotalAngularMomentum();
  */

  //-----------------------------------------------------------------------------------------------
  // \name Processing:

  /** Computes the force between the two particle p1, p2. */
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
  // \name Misc:

  void reset();

protected:



  std::vector<rsParticle<T>> particles;
  std::vector<rsVector3D<T>> forces;
  T stepSize;


};

#endif
