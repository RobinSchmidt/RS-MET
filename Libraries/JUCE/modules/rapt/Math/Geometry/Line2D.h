#ifndef RAPT_LINE2D_H_INCLUDED
#define RAPT_LINE2D_H_INCLUDED

/** This is a class for dealing with two dimensional lines. */

template<class T>
class rsLine2D
{

public:

  /** Given the two end points of a line (x0,y0) and (x1,y1), this function computes the 
  coefficients for the implicit line equation A*x + B*y + C = 0. The coefficients are 
  normalized such that A^2 + B^2 = 1. This implies that the right hand side of the equation
  gives the distance of any point (x,y) from the line. (Not yet tested) */
  static void twoPointToImplicit(T x0, T y0, T x1, T y1, T& A, T& B, T& C);

protected:

  //T x0 = 0, y0 = 0, x1 = 1, y1 = 1;

};

#endif
