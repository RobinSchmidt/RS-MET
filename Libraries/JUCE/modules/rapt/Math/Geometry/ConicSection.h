#ifndef RAPT_CONIC_SECTION_H_INCLUDED
#define RAPT_CONIC_SECTION_H_INCLUDED

/** This is a class for dealing with conic sections, represented by the implicit equation: 
A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. */

template<class T>
class rsConicSection
{

public:

  /** Constructor. Creates a conic section object with the equation 
  A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. If you pass no parameters, it will by default create
  a unit circle. */
  rsConicSection(T A = 1, T B = 0, T C = 1, T D = 0, T E = 0, T F = -1);

  /** Given a line defined by the parameteric equations: x(t) = x + t*dx, y(t) = y + dy, this
  function computes the two values of t (t1, t2), where the line intersects this conic. If the line
  doesn't intersect the conic, t1 and t2 will be NaN (resulting from a (real) square root of a 
  negative number). */
  void lineIntersectionParameter(T x, T dx, T y, T dy, T* t1, T* t2);

  /** The coefficients in the equation. */
  T A, B, C, D, E, F;

};

#endif
