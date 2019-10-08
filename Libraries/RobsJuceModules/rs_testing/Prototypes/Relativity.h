#pragma once

/** A class to represent spacetime events from relativity theory. A spacetime event has 4 
coordinates: the 3 spatial coordinates x,y,z from regular 3D space and one time coordinate t. 

not yet tested
*/

template<class T>
struct rsSpaceTimeVector
{

  T t, x, y, z;

  typedef rsSpaceTimeVector<T> STV; // for convenience


  /** Constructor */
  rsSpaceTimeVector(T t, T x, T y, T z)
  { this->t = t; this->x = x; this->y = y; this->z = z; }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Computes the spacetime interval between two events "a" and "b". */
  static T getInterval(const STV& a, const STV& b)
  { T d = a - b; return sqrt(d.x*d.x + d.y*d.y + d.z*d.z - d.t*d.t); }
  // maybe factor out a function that returns the interval for a single event from the origin

  /** Computes the proper time between two events "a" and "b". */
  static T getProperTime(const STV& a, const STV& b)
  { return -getInterval(a, b); }
  // proper time and spacetime intervals are invariant with respect to the reference frame - they
  // are the same in the (x,t) and in the (x',t') system -> verify numerically!

  /** Computes the Lorentz factor gamma that occurs in converting between stationary and moving
  reference frames as a function of the relative speed v (expressed as fraction of the 
  lightspeed). */
  static T gamma(T v) { return T(1) / sqrt(T(1) - v*v); }
  // in common units: gamma = 1 / sqrt(1 - v^2/c^2)
  // todo: what's the vector version of that, i.e. the version that takes a velocity vector
  // maybe rename to lorentzFactor

  /** Converts an x-coordinate given in a stationary frame to the corresponding x' coordinate in 
  the moving frame at time t with relative velocity v between the frames. */
  static T xToMovingFrame(T t, T x, T v) { return gamma(v) * (x - v*t);  }  // 1.19
  // common units x' = (x-v*t) / sqrt(1 - v^2/c^2), 1.23

  static T tToMovingFrame(T t, T x, T v) { return gamma(v) * (t - v*x);  }  // 1.20
  // common units: t' = (t - v*x/c^2) / sqrt(1 - v^2/c^2), 1.24

  static T xToStationaryFrame(T tp, T xp, T v) { return gamma(v) * (xp + v*tp);  } // 1.21
   
  static T tToStationaryFrame(T tp, T xp, T v) { return gamma(v) * (tp + v*xp);  } // 1.22

  static bool isSeparationTimeLike(const STV& a, const STV& b)
  { return getInterval(a, b) < T(0); }
  // pg 58 ff

  static bool isSeparationSpaceLike(const STV& a, const STV& b)
  { return getInterval(a, b) > T(0); }

  static bool isSeparationLightLike(const STV& a, const STV& b)
  { return getInterval(a, b) == T(0); }


  /** Given a relative velocity v between a rest frame A and a moving frame B and another relative 
  velocity u between the moving frame B and another moving frame C, this function computes the 
  relative velocity betwen the rest frame A and the second moving frame C. This can not just be the 
  regular sum because that would potentially lead to velocities faster than light (take v=u=0.9c, 
  for example). */
  static T addVelocities(T u, T v) { return (u + v) / (T(1) + u*v); }  // 2.9
  // common units: w = (u + v) / (1 + u*v/c^2), 2.13



  //-----------------------------------------------------------------------------------------------


  STV operator-(const STV& b) const
  { return STV(t - b.t, x - b.x, y - b.y, z - b.z); }


};