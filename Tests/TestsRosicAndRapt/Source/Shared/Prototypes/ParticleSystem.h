#ifndef ParticleSystem_h
#define ParticleSystem_h

#include "rapt/rapt.h"

/** This is just a skeleton/stub/idea for a particle system */



/** Class for representing 3-dimensional vectors. */

template<class T>
class rsVector3D
{

public:

  T x, y, z;

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
  //T charge;

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



  //-----------------------------------------------------------------------------------------------
  // \name Processing:

  /** Computes the forces that act on each particle and stores them in our "forces" member. */
  void updateForces();

  /** Updates the velocities according to the computed forces. */
  void updateVelocities();

  /** Updates the positions according to the current velocities. */
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
