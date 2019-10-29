#ifndef RAPT_LINE2D_H_INCLUDED
#define RAPT_LINE2D_H_INCLUDED

/** This is a class for dealing with two-dimensional lines. 


\todo: 
-drag over code from RSLib
-use pointers for output parameters
-use const, where applicable
*/

template<class T>
class rsLine2D
{

public:

  /** Given the two end points of a line (x0,y0) and (x1,y1), this function computes the 
  coefficients for the implicit line equation A*x + B*y + C = 0. If normalize is set to true, the 
  coefficients are normalized such that A^2 + B^2 = 1. This implies that the right hand side of the 
  equation gives the (signed) distance of any point (x,y) from the line. If it's not normalized, 
  the right hand side will be proportional to the distance. (Not yet tested) */
  static void twoPointToImplicit(T x0, T y0, T x1, T y1, T* A, T* B, T* C, bool normalize = true);

  /** Given the two end points of a line (x1,y1) and (x2,y2), this function computes the 
  coefficients for the explicit line equation y = a*x + b. (not yet tested) */
  static inline void twoPointToExplicit(T x1, T y1, T x2, T y2, T* a, T* b)
  {
    *a = (y2-y1)/(x2-x1); // slope
    *b = y1 - *a * x1;    // offset
  }
    // todo: use pointers for output variables

  /** Given the coordinates of a point (x,y) and the parameters of an implicit line equation
  A*x + B*y + C = 0, this function computes the coordinates of a point that is reflected about this
  line and returns them in xr, yr. */
  static void reflectPointInLine(T x, T y, T A, T B, T C, T* xr, T* yr);

  /** Given left and right endpoint coordinates (xL,yL),(xR,yR), this function returns the zero
  crossing of the line that goes trough these points, i.e. the x-value where y=0. */
  static inline T zeroCrossing(T xL, T yL, T xR, T yR) { return xL - yL*(xR-xL)/(yR-yL); }
    // does it work also, if xL,yR is actually the right endpoint - check, i think, so, if so,
    // rename parameters to x1,y1,x2,y2


protected:

  //T x0 = 0, y0 = 0, x1 = 1, y1 = 1;

};

#endif
