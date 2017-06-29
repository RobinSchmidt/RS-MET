#ifndef RAPT_ELLIPSE_H_INCLUDED
#define RAPT_ELLIPSE_H_INCLUDED

/** This is a class for dealing with ellipses. Ellipses are a special case of conic sections which 
have the general form A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. For ellipses, we have 
B^2 - 4*A*C < 0. The ellipse is set up in terms of user parameters defining a scale factor (with 
respect to the unit circle), an aspect ratio, determining the ration between the long and the short 
axis, a rotation angle, and a center point. Coordinates of points on the ellipse can be 
obtained by passing in a value of angle parameter p between 0 and 2*PI. */

template<class T>
class rsEllipse : public rsConicSection<T>
{

public:


  // setup:

  void setScaleFactor(T newFactor);

  void setAspectRatio(T newRatio);

  void setRotationAngle(T newAngle); // in radians

  void setCenterX(T newX);

  void setCenterY(T newY);

  // inquiry

  /** Computes a point on the ellipse corresponding to the given angle parameter. */
  void getPointOnEllipse(T angle, T* x, T* y) const;

protected:

  T scale = 1, ratio = 1, angle = 0, centerX = 0, centerY = 0;

};

#endif
