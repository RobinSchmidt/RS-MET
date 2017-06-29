#ifndef RAPT_ELLIPSE_H_INCLUDED
#define RAPT_ELLIPSE_H_INCLUDED

/** This is a class for dealing with ellipses. Ellipses are a special case of conic sections which 
have the general form A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. For ellipses, we have 
B^2 - 4*A*C < 0. The ellipse is set up in terms of user parameters defining a scale factor (with 
respect to the unit circle), an aspect ratio, determining the ration between the long and the short 
axis, a rotation angle, and a center point. Coordinates of points on the ellipse can be 
obtained by passing in a value of angle parameter p between 0 and 2*PI. */

template<class T>
class rsEllipse : protected rsConicSection<T>
{

public:

  rsEllipse(T scale = 1, T aspectRatio = 1, T angle = 0, T centerX = 0, T centerY = 0);

  // setup:

  /** Sets all the parameters that determine the ellipses shape, size, orientation and position
  at once. */
  void setParameters(T newScale = 1, T newAspectRatio = 1, T newAngle = 0, T newCenterX = 0, 
    T newCenterY = 0);

  //void setScaleFactor(T newFactor);
  //void setAspectRatio(T newRatio);
  //void setRotationAngle(T newAngle); // in radians
  //void setCenterX(T newX);
  //void setCenterY(T newY);

  // inquiry

  /** Computes a point on the ellipse corresponding to the given angle parameter. */
  void getPointOnEllipse(T angle, T* x, T* y) const;

  /** Returns the total area of the ellipse (not yet tested). */
  inline void getArea() const { return sqrt(scale) * T(PI); }

  /** Given coordinates of a point x,y assumed to be on the ellipse, this function returns the 
  implicit line equation coefficients of a line A*x + B*y + C = 0 that is tangent to the ellipse
  at this point. */
  void getTangentCoeffsAt(T x, T y, T* A, T* B, T* C);

  /** Updates the parametric and implicit equation coefficients according to user parameters. */
  void updateCoeffs();


protected:

  T scale = 1, ratio = 1, angle = 0, centerX = 0, centerY = 0;  // user parameters
  T Axc, Axs, Ayc, Ays;                                         // parametric equation coeffs

};

#endif
