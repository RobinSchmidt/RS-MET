#pragma once

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





  STV operator-(const STV& b) const
  { return STV(t - b.t, x - b.x, y - b.y, z - b.z); }


};