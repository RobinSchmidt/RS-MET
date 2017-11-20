#ifndef RS_GEOMETRICFUNCTIONS_H
#define RS_GEOMETRICFUNCTIONS_H

namespace RSLib
{

  // functions related to geometric quantities like areas, circumfences, etc.

  /** Computes the area of a sector of an ellipse with major and minor half-semiaxes a and b. The
  angle parameter is the angle between the x-axis and a straight line from the origin outwards
  which itersects the ellipse. The returned area is the area of the elliptic segment with this
  angle. Note that this angle doe not equal the parameter t in the parametric equation for the
  ellipse but is related to it via t = atan(a*tan(angle)/b). */
  RSLib_API double rsEllipticSectorArea(double a, double b, double angle);

  // maybe the elliptic integrals and stuff would fit here as well?

}

#endif
