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
-(-project onto the x- and y-axis (i.e. take the x- and y-value))...really? */

template<class T>
class rsEllipseOscillator
{

public:

  inline void setA(T newA) { a = newA; }  // rename to setOffset

  inline void setB(T newB) { b = newB; }  // rename to setRotation

  inline void setC(T newC) { c = newC; }  // rename to setScale




  inline T getSample()
  {
    s  = sin(p);  // x
    c  = cos(p);  // y

    // precompute these:
    T sb = sin(b);
    T cb = cos(b);

    T ac = a + c;    // x += a -> shift y
    T cs = c * s;    // y *= c -> scale x

    T a = 1 / sqrt(ac*ac + cs*cs);  // normalizer
    T x = a*ac*cb;  // rotate and
    T y = a*cs*sb;  // normalize

    return x + y;   // use mix*x + (1-mix)*y
  }

protected:

  T p = 0;   // phase
  T w = 0;   // normalied radian frequency

  // waveshape parameters:
  T a = 0; 
  T b = 0;
  T c = 0;


};

#endif