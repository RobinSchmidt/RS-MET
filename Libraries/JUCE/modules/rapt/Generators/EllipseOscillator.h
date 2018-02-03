#ifndef RAPT_ELLIPSEOSCILLATOR_H_INCLUDED
#define RAPT_ELLIPSEOSCILLATOR_H_INCLUDED

/** Oscillator based on an ellipse in the xy-plane.

Credits:
 -Xoxos

more info:
 https://www.kvraudio.com/forum/viewtopic.php?p=6656238#p6656238
 https://gitlab.com/Hickler/Soundemote/issues/67
 https://github.com/RobinSchmidt/RS-MET/issues/72
play with parameters:
  https://www.desmos.com/calculator/7h9mknbv3q
i think, it works as follows:
-create x,y values on a circle (standard rotating phasor in the plane)
-convert x,y to values on an arbitrary ellipse
-(-project onto the x- and y-axis (i.e. take the x- and y-value))...really? 

the name is actually not really suitable - maybe rename to XoxosOscillator

...and make an EllipseOscillator as well

*/

template<class T>
class rsEllipseOscillator
{

public:

  inline void setA(T newA) { A = newA; }  // rename to setOffset

  inline void setB(T newB) { B = newB; }  // rename to setRotation

  inline void setC(T newC) { C = newC; }  // rename to setScale

  inline void setOmega(T newOmega) { w = newOmega; }

  inline void setRotationDegrees(T newRotation) { B = rsDegreeToRadiant(newRotation); }
   // no - this is not a rotation...maybe phase-shift?


  inline void updatePhase()
  {
    p += w;
    while(p > 2*PI)
      p -= 2*PI;

    // maybe do also a wraparound at 0 -> allow negative frequencies
  }

  // maybe have a function that returns x and y separately - maybe that's useful as stero signal?
  // if not, it is certainly helpful to figure out what the osc is doing - and then we can do 
  // mid/side (re)mixing

  inline void getSamplePair(T* x, T* y)
  {
    T c  = cos(p);  // x
    T s  = sin(p);  // y

    // precompute these:
    T sB = sin(B);
    T cB = cos(B);

    T Ac = A + c;    // x += a -> shift x
    T Cs = C * s;    // y *= c -> scale y


    T scl = 1 / sqrt(Ac*Ac + Cs*Cs);  // normalizer
    if(renorm != 1)
      scl = pow(scl, renorm);

    *x = scl*Ac*cB;  // phase(?) and
    *y = scl*Cs*sB;  // normalize

    updatePhase();
  }

  inline T getSample()
  {
    T x, y;
    getSamplePair(&x, &y);
    return 2*((1-mix)*x + mix*y); // precompute coeffs 2*(1-mix), 2*mix
    //return x + y;   // use mix*x + (1-mix)*y
  }

  inline void reset() { p = startPhase; }

protected:

  T startPhase = 0;  
  T p = 0;   // current phase
  T w = 0;   // normalized radian frequency
  // have independent phases and frequencies for x and y

  // waveshape parameters:
  T A = 0; 
  T B = 0;
  T C = 1;

  T renorm = 1;

  T mix = 0.5;
};

#endif