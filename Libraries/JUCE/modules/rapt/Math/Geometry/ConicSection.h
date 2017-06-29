#ifndef RAPT_CONIC_SECTION_H_INCLUDED
#define RAPT_CONIC_SECTION_H_INCLUDED

/** This is a class for dealing with conic sections, represented by the implicit equation: 
A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. */

template<class T>
class rsConicSection
{

public:

  /** Constructor. Creates a conic section object with the equation 
  A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. */
  rsConicSection(T A, T B, T C, T D, T E, T F);

  /** Given a line defined by the parameteric equations: x(t) = x + t*dx, y(t) = y + dy, this
  function computes the two values of t (t1, t2), where the line intersects this conic. If the line
  doesn't intersect the conic, t1 and t2 will be NaN (resulting from a (real) square root of a 
  negative number). */
  void lineIntersectionParameter(T x, T dx, T y, T dy, T* t1, T* t2);

protected:

  T A = 1, B = 0, C = 1, D = 0, E = 0, F = -1; // init as unit circle

};

#endif
