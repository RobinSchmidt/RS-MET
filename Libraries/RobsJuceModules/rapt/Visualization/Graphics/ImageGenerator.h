#ifndef RAPT_IMAGEGENERATOR_H_INCLUDED
#define RAPT_IMAGEGENERATOR_H_INCLUDED

/** A class for generating images by various algorithms. In particular, it may generate images that
visualize mathematical entities such as 2D functions z = f(x,y) or implicit or parameteric curves, 
so the algorithms here may serve as building blocks for creating mathematical plots. 

maybe rename to rsImagePlotter
*/

template<class TPix, class TVal>  // pixel and value types (for coordinates, heights, etc.)
class rsImageGenerator
{

public:

  void setRange(TVal minX, TVal maxX, TVal minY, TVal maxY)
  { xMin = minX; xMax = maxX; yMin = minY; yMax = maxY;  }


  /** Draws the curve defined by f(x,y) = c onto the image. It needs one solution x0,y0 for which
  f(x0,y0) = c holds as starting point. */
  void drawImplicitCurve(const std::function<TVal(TVal, TVal)>& f, TVal c, TVal x0, TVal y0,
    rsImage<TPix>& img, TPix color)
  { _drawImplicitCurve(f, c, x0, y0, img, color, false); }




  // todo: drawFunction (variants: y = f(x), z = f(x,y)), drawParametricCurve, drawCoordinateGrid



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






  rsImagePainter<TPix, TVal, TVal> painter;


protected:

  /** Internal function called by drawImplicitCurve - we need it, because we don't want the 
  "clockwise" parameter to be exposed to client code (it should be set to true only in the 
  recursive call of _draw... to itself). */
  void _drawImplicitCurve(const std::function<TVal(TVal, TVal)>& f, TVal c, TVal x0, TVal y0, 
    rsImage<TPix>& img, TPix color, bool clockwise); 


  TVal xMin, xMax, yMin, yMax;  // plotting range


};


#endif