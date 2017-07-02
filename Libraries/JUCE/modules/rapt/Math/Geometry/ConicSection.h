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

  /** Given a line defined by the parameteric equations: x(t) = x0 + t*dx, y(t) = y0 + t*dy, this
  function computes the two values of t (t1, t2), where the line intersects this conic. If the line
  doesn't intersect the conic, the equation will formally yield two complex conjugate intersection
  points and the returned t1, t2 will be the real part of them. If there are two real intersection 
  points, they will be returned in ascending order. */
  void lineIntersectionParameter(T x0, T dx, T y0, T dy, T* t1, T* t2) const;

  /** Evaluates the conic equation A*x^2 + B*x*y + C*y^2 + D*x + E*y + F at the given coordinates 
  x,y. On the conic, this value equals zero, otherwise it is proportional to the signed distance of 
  the point from the conic. (question: what's the proportionality constant? can we normalize the 
  coeffs to make it 1? ...or wait...is it actually proportional or some other function of the 
  distance?...maybe derive a proper distance formula...)
  see here http://mathforum.org/library/drmath/view/52082.html
  */
  T evaluate(T x, T y) const;

  /** Given a point with coordinates x,y assumed to be on the conic, this function computes the
  coefficients for an implicit line equation A*x + B*y + C = 0 for the line that is tangent to the
  conic at the given point. */
  void getTangentCoeffs(T x, T y, T* A, T* B, T* C) const;

  //void getNormalVector
  //void getCurvature
  // see https://en.wikipedia.org/wiki/Implicit_curve#Tangent_and_normal_vector
  // for formulas


protected:

  /** The coefficients in the equation. */
  T A, B, C, D, E, F;

};

#endif
