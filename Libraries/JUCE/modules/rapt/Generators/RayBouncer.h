#ifndef RAPT_RAYBOUNCER_H_INCLUDED
#define RAPT_RAYBOUNCER_H_INCLUDED


/** A sound generator based on the idea of a particle/ray that bounces around inside an elliptic 
enclosure. 
maybe rename this to RayBouncerCore

\todo: inlcude friction parameter -> replace "speed" by "initialSpeed" or "launchSpeed", re-use
speed variable for a decaying speed: speed -= friction*speed

*/

template<class T>
class rsRayBouncer
{

public:

  /** \name Construction/Destruction */

  /** Constructor.  */
  rsRayBouncer();


  /** \name Setup */

  /** Sets the initial position of the particle. */
  inline void setInitialPosition(T x, T y) { x0 = x; y0 = y; }

  /** Sets the angle at which our particle is launched from its initial position (in radians). */
  inline void setLaunchAngle(T newAngle) { angle = newAngle; }

  /** Sets the speed by which our particle moves around, i.e. the magnitude of the velocity, 
  which is a vector. */
  inline void setSpeed(T newSpeed) { speed = newSpeed; }
  // later replace this with setFrequency, setSampleRate...the speed is then proportional to
  // frequency/sampleRate - but maybe move this to some outside "driver" class
  // setting the speed variable is not enough - we also need to re-normalize the dx,dy

  /** Sets the maximum distance that the particle can travel before a reset is forced. This helps 
  to get more tonal sounds. */
  inline void setMaxDistance(T newMaxDistance) { maxDistance = newMaxDistance; }
   // maybe this should be set up in terms of cycles so the absolute maxDistance goes down when the
   // frequency goes up

  /** Sets the parameters of the enclosing ellipse. Note that we do not include a shear 
  transformation because that would be redundant, because a sheared ellipse is still an ellipse
  (affine transformations map conic sections into conic sections of the same type) */
  inline void setEllipseParameters(T newScale = 1, T newAspectRatio = 1, T newAngle = 0,
    T newCenterX = 0, T newCenterY = 0)
  {
    ellipse.setParameters(newScale, newAspectRatio, newAngle, newCenterX, newCenterY);
  }

  /** Sets the overall amount of bending effects on the velocity. This tends to drag the output to
  a more periodic waveform (corresponding to a harmonic spectrum). It can also create chaotic 
  behavior. */
  inline void setBend(    T newAmount) { bendAmount = newAmount; }
  inline void setBendX2Y( T newValue)  { xToY  = newValue; }
  inline void setBendY2X( T newValue)  { yToX  = newValue; }
  inline void setBendXX2X(T newValue)  { xxToX = newValue; }
  inline void setBendXX2Y(T newValue)  { xxToY = newValue; }
  inline void setBendXY2X(T newValue)  { xyToX = newValue; }
  inline void setBendXY2Y(T newValue)  { xyToY = newValue; }
  inline void setBendYY2X(T newValue)  { yyToX = newValue; }
  inline void setBendYY2Y(T newValue)  { yyToY = newValue; }



  /** \name Processing */

  /** Computes one x,y-pair of output values at a time. */
  void getSampleFrame(T &x, T &y);
  // use pointers for output variables

  /** Resets the internal state (position and velocity). */
  void reset();

  /** Resets position and velocity and advances the point a little bit into the direction of the
  velocity. Used for continuous response the changes of maxDistance. */
  void resetWithAdvance(T advance);


protected:

  /** \name Internal Functions */

  /** Computes the intersection point between the current line segment along 
  xc + t*dx, yc + t*dy and the ellipse. */
  inline void getLineEllipseIntersectionPoint(T* xi, T* yi);

  /** Given a point xt,yt assumed to be on the ellipse, this function reflects the point
  x,y about the tangent at that point. */
  inline void reflectInTangentAt(T xt, T yt, T* x, T *y);

  inline void ensurePointIsInEllipse(T xi, T yi);

  /** Updates dx, dy by taking the direction vector that points from the given line/ellipse 
  intersection point xi,yi to the current point as new direction and adjusting its length according 
  to the desired speed. */
  inline void updateVelocity(T xi, T yi);


  /** \name Data */

  rsEllipse<T> ellipse; // enclosing ellipse

  T tolerance = float(1.e-8);   // numerical tolerance for inside/outside checks

  // particle state:
  T x0 = 0, y0 = 0;      // initial position of particle
  T x , y, dx, dy;       // current position and velocity (as increment per sample):
  T distance = 0;        // total distance travelled

  // coefficients to scale products dx*dx, dx*dy, dy*dy to add them to the velocity vetor for
  // nonlinear effects:
  T xToY = 0, yToX = 0;
  T xxToX = 0, xyToX = 0, yyToX = 0;
  T xxToY = 0, xyToY = 0, yyToY = 0;
  T bendAmount = 0;  // global scaler for bending effects

  // user parameters: 
  T speed = T(0.2);          // speed (i.e. magnitude of velocity)
  T angle = T(0.0);          // launching angle
  T maxDistance = RS_INF(T); // maximum distance before a reset is forced



  // maybe introduce a maximum distance to travel and reset after that. per sample, we would update:
  // distance += speed;
  // if(distance > maxDistance)
  //   reset(distance-maxDistance)
  // the difference distance-maxDistance is used inside reset to advance the particle a little bit
  // into the initial direction to avoid jitter

  //friend class rsRayBouncerDriver;
};

//=================================================================================================

/** A class that encapsulates an rsRayBouncer object to let the user set it up in terms of musical
parameters (such as a pseudo frequency, etc.) and also some modulators that animate certain
parameters of the bouncer such as the parameters of the ellipse, an output transformation, etc. */

template<class T>
class rsRayBouncerDriver
{

public:

  /** \name Setup */

  /** Sets the frequency and sample rate in Hz */
  void setFrequencyAndSampleRate(T newFreq, T newRate);

  // ellipse parameters:
  inline void setEllipseSize(        T newSize)  { ellipseSize    = newSize;            }
  inline void setEllipseAspectRatio( T newRatio) { ellipseRatio   = newRatio;           }
  inline void setEllipseAngleDegrees(T newAngle) { ellipseAngle   = T(PI/180)*newAngle; } 
  inline void setEllipseCenterX(     T newX)     { ellipseCenterX = newX;               }
  inline void setEllipseCenterY(     T newY)     { ellipseCenterY = newY;               }

  // other parameters:
  void setStartX(T newX);
  void setStartY(T newY);
  inline void setLaunchAngleDegrees(T newAngle) { rayBouncer.setLaunchAngle(T(PI/180)*newAngle); }

  /** \name Processing */

  /** Computes one x,y-pair of output values at a time. */
  inline void getSampleFrame(T &x, T &y)
  {
    rayBouncer.setEllipseParameters(ellipseSize, ellipseRatio, ellipseAngle, ellipseCenterX,
      ellipseCenterY);
    rayBouncer.getSampleFrame(x, y);

    //// test - apply some high frequency texture derived form a nonlinear combination of the 
    //// outputs:
    //T k = 0.2;
    //T tx, ty;
    ////tx = x*(1-x*y);
    ////ty = y*(1-x*y);
    ////tx = y*(x+y);
    ////ty = x*(x+y);
    //tx = (1-x*y); tx *= tx; tx *= tx;
    //ty = tx;
    //x += k*tx;
    //y += k*ty;
  }

  /** Resets the internal state (position and velocity). */
  void reset();


  rsRayBouncer<T> rayBouncer;

protected:

  T ellipseSize = 1, ellipseRatio = 1, ellipseAngle = 0, ellipseCenterX = 0, ellipseCenterY = 0;
  T startX = 0, startY = 0;

  // rsSineOscillator lfoSize, lfoRatio, lfoX, lfoY;
  // rsNaiveSawOscillator lfoAngle;
  // maybe, instead of letting the center just x,y bob up and down, move the ellipse center in a 
  // circle or ellipse as well...or maybe do both...lots of movement in the ellipse should lead
  // to complex animated shape and sound
  // maybe when the ellipse size is animated, we should compensate for any frequency change that
  // is caused by this - or maybe compensate only partially according to a user parameter

};

//=================================================================================================

/** A simplified 1D version of the raybouncer concept. The particle just slides up and down on a 
one dimensional line and gets reflected whenever it hits the floor or the ceiling. 

todo:
-make it possible to change the shape by feedback (see elan's envelope for formulas)

*/

template<class T>
class rsRayBouncer1D
{

public:

  /** \name Setup */

  void setFloor(T newFloor) { floor = newFloor; }
  void setCeil( T newCeil)  { ceil = newCeil;   }  
  // todo: sanity checks (floor < ceil), maybe have something like abs(ceil-floor) > minDistance
   

  void setIncrement(T newIncrement) { inc = newIncrement; }

  void setShape(T newShape) { shape = newShape; }


  inline T getSample()
  {
    T out = x;

    // update:
    x  += dx;
    dx += shape*dx; //

    // reflection(s):
    while(x > ceil || x < floor)
    {
      if(x > ceil)
      {
        T over = x - ceil;
        x = ceil - over;
        dx = -inc;
      }
      if(x < floor)
      {
        T under = floor - x;
        x = floor + under;
        dx = inc;
      }
    }

    return out;
  }

  void reset() 
  { 
    x  = start; 
    dx = inc;
  }


protected:

  T floor = 0, ceil = 1;
  T sampleRate = 44100;
  T inc = T(0.01);
  T dx  = T(0.01);
  T x   = 0;
  T start = 0;
  T shape = 0;

};

#endif