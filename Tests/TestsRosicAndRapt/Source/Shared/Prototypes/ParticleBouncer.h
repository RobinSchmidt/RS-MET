#ifndef ParticleBouncer_h
#define ParticleBouncer_h

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
  inline void setEnclosureEllipseAspectRatio(double ratio) { a = ratio; b = 1/a; }
  // check this

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns "a" parameter of enclosing ellipse. */
  inline double getEllipseA() const { return a; }

  /** Returns "b" parameter of enclosing ellipse. */
  inline double getEllipseB() const { return b; }

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  void getSampleFrame(double &x, double &y);

  //-----------------------------------------------------------------------------------------------
  // others:

  void reset();

protected:

  double a = 1, b = 1;        // parameters of enclosure ellipse
  double dx = 0.1, dy = 0.2;  // velocity components of particle (as increment per sample)
  double x0 = 0,   y0 = 0;    // initial particle coordinates
  double xc = 0,   yc = 0;    // current particle coordinates

};



#endif
