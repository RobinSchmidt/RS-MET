#pragma once

/** A class to represent spacetime events from relativity theory. A spacetime event has 4
coordinates: the 3 spatial coordinates x,y,z from regular 3D space and one time coordinate t. The
class implemented a bunch of formulas that are needed for converting time- and space coordinates
between reference frames that are moving with respect to one another. Because one can always
consider one of them stationary and only the other one to be moving, you'll find function names
referring to stationary and moving frames, but keep in mind that it's a matter of convention which
is which or even if you consider both moving.

not yet tested
*/

template<class T>
struct rsSpaceTimeVector
{

  T t, x, y, z;

  typedef rsSpaceTimeVector<T> STV; // for convenience


  /** Constructor */
  rsSpaceTimeVector(T t, T x, T y, T z)
  {
    this->t = t; this->x = x; this->y = y; this->z = z;
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Computes the spacetime interval between two events "a" and "b". */
  static T getInterval(const STV& a, const STV& b)
  {
    T d = a - b; return sqrt(d.x*d.x + d.y*d.y + d.z*d.z - d.t*d.t);
  }
  // maybe factor out a function that returns the interval for a single event from the origin

  /** Computes the proper time between two events "a" and "b". */
  static T getProperTime(const STV& a, const STV& b)
  {
    return -getInterval(a, b);
  }
  // proper time and spacetime intervals are invariant with respect to the reference frame - they
  // are the same in the (x,t) and in the (x',t') system -> verify numerically!

  /** Computes the Lorentz factor gamma that occurs in converting between stationary and moving
  reference frames as a function of the relative speed v (expressed as fraction of the
  lightspeed). */
  static T gamma(T v) { return T(1) / sqrt(T(1) - v*v); }
  // in common units: gamma = 1 / sqrt(1 - v^2/c^2)
  // todo: what's the vector version of that, i.e. the version that takes a velocity vector
  // maybe rename to lorentzFactor
  // maybe implement the approximations given in 3.15, 3.15 - or generally, in 3.11

  /* The Lorentz transformation is a linear transformation from the (x,t) coordinate system to the
  (x',t') coordinate system given by the matrix-vector-product:

  |x'| = | g   -vg| * |x|
  |t'|   |-vg   g |   |t|

  where g = 1/sqrt(1-v^2), vg = v*g. maybe we should put the t-coordinate first? ...maybe
  generalize to allow for separate vx, vy, vz velocities
  */

  static T gamma(T vx, T vy, T vz) { return T(1) / sqrt(T(1) - (vx*vx + vy*vy + vz*vz)); }


  /** Converts an x-coordinate given in a stationary frame to the corresponding x' coordinate in
  the moving frame at time t with relative velocity v (assumed to be in the x-direction) between
  the frames. */
  static T xToMovingFrame(T t, T x, T v) { return gamma(v) * (x - v*t); }  // 1.19
  // common units x' = (x-v*t) / sqrt(1 - v^2/c^2), 1.23

  static T tToMovingFrame(T t, T x, T v) { return gamma(v) * (t - v*x); }  // 1.20
  // common units: t' = (t - v*x/c^2) / sqrt(1 - v^2/c^2), 1.24

  static T xToStationaryFrame(T tp, T xp, T v) { return gamma(v) * (xp + v*tp); } // 1.21

  static T tToStationaryFrame(T tp, T xp, T v) { return gamma(v) * (tp + v*xp); } // 1.22

  static bool isSeparationTimeLike(const STV& a, const STV& b)
  {
    return getInterval(a, b) < T(0);
  }
// pg 58 ff

  static bool isSeparationSpaceLike(const STV& a, const STV& b)
  {
    return getInterval(a, b) > T(0);
  }

  static bool isSeparationLightLike(const STV& a, const STV& b)
  {
    return getInterval(a, b) == T(0);
  }


  /** Given a relative velocity v between a rest frame A and a moving frame B and another relative
  velocity u between the moving frame B and another moving frame C, this function computes the
  relative velocity betwen the rest frame A and the second moving frame C. This can not just be the
  regular sum because that would potentially lead to velocities faster than light (take v=u=0.9c,
  for example: v+u = 1.8c is not allowed). */
  static T addVelocities(T u, T v) { return (u + v) / (T(1) + u*v); }  // 2.9
  // common units: w = (u + v) / (1 + u*v/c^2), 2.13


  /** Given the 3 components of a regular 3D velocity vector v = (vx,vy,vz), this function computes
  the corresponding 4-velocity vector. */
  STV fourVelocity(T vx, T vy, T vz)
  {
    T g = gamma(vx, vy, vz);
    return STV(g, g*vx, g*vy, g*vz);  // 3.8, 3.9
  }
  // todo: verify 3.7: the proper time of a 4-velocity vector should always equal 1

  STV fourMomentum(T m, T vx, T vy, T vz) { return m * fourVelocity(vx, vy, vz); }
  // 3.31, 3.36

  T lagrangian(T m, T vx, T vy, T vz) { return -m * gamma(vx, vy, vz); } // 3.26


  T energy(T m, T vx, T vy, T vz)
  {
    return m * gamma(vx, vy, vz);
  } // 3.38 with c = 1
  // common units: E = m*c^2 / sqrt(1 - v^2/c^2), reduces to E = m*c^2 for v=0
  // maybe implement approximation Eq 3.39
  // what about 3.44 - what's the P^2 there?

  //-----------------------------------------------------------------------------------------------


  STV operator-(const STV& b) const
  {
    return STV(t - b.t, x - b.x, y - b.y, z - b.z);
  }


};

/** Multiplies a scalar and a spacetime vector. */
template<class T>
inline rsSpaceTimeVector<T> operator*(const T& s, const rsSpaceTimeVector<T>& v)
{
  return rsSpaceTimeVector<T>(s*v.t, s*v.x, s*v.y, s*v.z);
}

//=================================================================================================

template<class T>
class rsQuaternion
{

public:


  /** The quaternion is given by: q = a + b*i + c*j + d*k */
  T a, b, c, d;



  /** Constructor. Initializes elements such that the quaternion is given by:
  q = a + b*i + c*j + d*k */
  rsQuaternion(T a, T b, T c, T d) { setValues(a, b, c, d); }

  /** Standard constructor. Leaves elements uninitialized. */
  rsQuaternion() {}
  // try to get rid to enforce RAII



  /** \name Setup */

  void setValues(T a, T b, T c, T d) { this->a = a; this->b = b; this->c = c; this->d = d; }



  T getSquaredNorm() const { return a*a + b*b + c*c + d*d; }

  T getNorm() const { return sqrt(getSquaredNorm); }

  rsQuaternion<T> conj() const { return rsQuaternion<T>(a, -b, -c, -d); }

  rsQuaternion<T> reciprocal() const { return this->conj() / this->getSquaredNorm(); }
  // needs division by scalar to be defined

  /** Multiplies two quaternions. */
  rsQuaternion<T> operator*(const rsQuaternion<T>& q) const
  {
    return rsQuaternion<T>(
      a*q.a - b*q.b - c*q.c - d*q.d,
      a*q.b + b*q.a + c*q.d - d*q.c,
      a*q.c - b*q.d + c*q.a + d*q.b,
      a*q.d + b*q.c - c*q.b + d*q.a);
    // formula from: https://en.wikipedia.org/wiki/Quaternion#Hamilton_product
  }


  // todo: +,-,/


};
// this class should go somewhere into rapt/Math

// the idea is to construct an example electromagnetic vector potential A from q^2 where q is a
// quaternion and do experiments with that vector potential - i.e. derive electric and magnetic
// fields (E and B) from it and simulate the motion of a charged particle in that field - once with
// the regular maxwell equations and then also with the maxwell equations in tensor form. in the
// former case, use time t as independent parameter, in the latter, use proper time tau as
// independent parameter -> we may get a non-equidistant time-axis in the latter case
// the example vector potential is then:
// A(t,x,y,z) = |t^2 - x^2 - y^2 - z^2|
//              |t*x + x*t + y*z - z*y|
//              |t*y - x*z + y*t + z*x|
//              |t*z + x*y - y*x + z*t|
// -> verify and simplify...

