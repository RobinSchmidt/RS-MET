#ifndef ParticleBouncer_h
#define ParticleBouncer_h

#include "rosic/rosic.h"

/** A sound generator based on the idea of a particle that bounces around inside an enclosure. 
Maybe rename to EllipseBounncer */

class ParticleBouncer
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ParticleBouncer();   

  /** Destructor. */
  ~ParticleBouncer(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the initial position of the particle. */
  inline void setInitialPosition(double x, double y) { x0 = x; y0 = y; }

  /** Sets the aspect ratio for the enclosing ellipse, i.e. the ratio of its width to its 
  height. */
  inline void setEnclosureEllipseAspectRatio(double ratio) { a = sqrt(ratio); b = 1/a; }

  /** Sets the angle at which our particle is launched from its initial position (in degrees). */
  inline void setLaunchAngle(double newAngle) { angle = PI * newAngle / 180; }

  /** Sets the speed by which our particle moves around, i.e. the magnitude of the velocity, 
  which is a vector. */
  inline void setSpeed(double newSpeed) { speed = newSpeed; }
  // later replace this with setFrequency, setSampleRate...the speed is then proportional to
  // frequency/sampleRate

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
  y(t) = b * sin(2*PI*t) (t in 0..1)), this function computes the parameter value t, for which the 
  line and the ellipse intersect.  
  WARNING: It is assumed that the caller has already determined that such an intersection exists, 
  if there's none, you'll get a wrong result. Also, if there's one intersection, there's typically
  another one (unless the line just touches the ellipse) - the function returns the one with the
  larger value. */
  static double getLineEllipseIntersectionParameter(double x, double y, double dx, double dy,
    double a2r, double b2r);
  // todo: maybe make this a protected member function and do not pass a2r, b2r - instead make them
  // member variables

  /** Given the coordinates of a point (x,y) and the parameters of an implicit line equation
  A*x + B*y + C = 0, this function computes the coordinates of a point that is reflected about this
  line and returns them in xr, yr. */
  static void reflectPointInLine(double x, double y, double A, double B, double C,
    double *xr, double *yr);
  // todo: move this function into some geometry class (like rapt, Line2D)

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
  double dx  = 0.1, dy  = 0.2;  // current
  double speed = 0.2;
  double angle = 0.0;
};

#endif
