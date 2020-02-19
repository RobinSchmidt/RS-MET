#ifndef RAPT_IMAGEGENERATOR_H_INCLUDED
#define RAPT_IMAGEGENERATOR_H_INCLUDED

/** A class for generating images by various algorithms. In particular, it may generate images that
visualize mathematical entities such as 2D functions z = f(x,y) or implicit or parameteric curves, 
so the algorithms here may serve as building blocks for creating mathematical plots. */

template<class TPix, class TVal>  // pixel and value types (for coordinates, heights, etc.)
class rsImageGenerator
{

public:




  /** Given plane coordinates x,y, this function computes a height above the plane that has the 
  shape of ridge (of height 1) that spirals around in a logarithmic spiral with the parametric 
  equation:
      x(t) = exp(a*t) * cos(sign * t + p); y(t) = exp(a*t) * sin(sign * t + p)
  For points that are on the spiral, the function will return zero and for points that are 
  "halfway" in between two "arcs", it will return 1. Halfway is to be understood in the logarithmic 
  sense - for example, if (1,0) and (2,0) are points on the spiral, the point (sqrt(2),0) would be 
  considered halfway between them. If the exponential growth parameter "a" is equal to 
  log(2)/(2*pi), the spiral will grow by a factor of 2 in each revolution. The "sign" parameter 
  should be +1 or -1 and determines the direction of the rotation. */
  TVal spiralRidge(TVal x, TVal y, TVal a = TVal(1), TVal p = TVal(0), TVal sign = TVal(1), 
    int profile = 0, TVal exponent = TVal(1));

  // move code for implicit function function plotting into this class


protected:

};


#endif