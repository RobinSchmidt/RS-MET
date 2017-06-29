#ifndef ParticleBouncer_h
#define ParticleBouncer_h

#include "rosic/rosic.h"

/** A sound generator based on the idea of a particle that bounces around inside an enclosure. */

class ParticleBouncer
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  ParticleBouncer();   ///< Constructor.

  ~ParticleBouncer();  ///< Destructor.

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the initial position of the particle. */
  inline void setInitialPosition(double x, double y) { x0 = x; y0 = y; }

  /** Sets the aspect ratio for the enclosing ellipse, i.e. the ratio of its width to its 
  height. */
  inline void setEnclosureEllipseAspectRatio(double ratio) { a = sqrt(ratio); b = 1/a; }

  inline void setInitialIncrements(double newDeltaX, double newDeltaY)
  {
    dx0 = newDeltaX;
    dy0 = newDeltaY;
  }

  //void setAngle(double newAngle);
  //void setSpeed(double newSpeed);
  // use dx = speed * cos(angle), dy = speed * sin(angle)


  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns "a" parameter of enclosing ellipse. */
  inline double getEllipseA() const { return a; }

  /** Returns "b" parameter of enclosing ellipse. */
  inline double getEllipseB() const { return b; }

  /** Given a line-segment from (x,y) to (x+dx,y+dy) with the parameteric equations 
  x(t) = x + t*dx, y(t) = y + t*dy (t in 0..1) and reciprocals of the squares of
  ellipse parameters a2r = 1/a^2, b2r = 1/b^2 in an (implicit) ellipse equation 
  x^2/a^2 + y^2/b^2 - 1 = 0 (corresponding to the parametric equation x(t) = a * cos(2*PI*t),
  y(t) = b * sin(2*PI*t) (t in 0..1)), this function computes parameter t, where the line and the 
  ellipse intersect.  
  WARNING: it is assumed that the caller has already determined that such an intersection exists, 
  if there's none, you'll get a square-root of a negative number */
  double getLineEllipseIntersectionParameter(double x, double y, double dx, double dy,
    double a2r, double b2r);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  void getSampleFrame(double &x, double &y);

  //-----------------------------------------------------------------------------------------------
  // others:

  void reset();

protected:

  // parameters of enclosure ellipse:
  double a = 1, b = 1;          

  // particle coordinates:
  double x0 = 0, y0 = 0;        // initial
  double xc = 0, yc = 0;        // current

  // velocity components of particle (as increment per sample):
  double dx0 = 0.1, dy0 = 0.2;  // initial
  double dx  = 0.1, dy  = 0.2;  // current

};



#endif
