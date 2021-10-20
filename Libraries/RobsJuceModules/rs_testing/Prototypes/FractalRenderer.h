#ifndef RS_FRACTALRENDERER_H
#define RS_FRACTALRENDERER_H

class rsFractalImageRenderer
{

public:

  using Vec2D = RAPT::rsVector2D<double>;
  using Color = rsFloat32x4;
 

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsFractalImageRenderer();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the iteration function to compute z[n+1] = f(z[n], p). The function should take as input
  the current value z[n] and from that produce a new value z[n+1]. The second argument p is some 
  optional parameter that the function can take. In Newton fractals, this plays no role and can be 
  ignored. In Mandelbrot fractals, it is used for the c in z[n+1] = z[n]^2 + c. */
  void setIterationFunction(const std::function<Vec2D(Vec2D z, Vec2D p)>& newFunc) 
  { iterationFunction = newFunc; }

  /** Sets the coloring function that determines a color (i.e. tuple of 4 floats) from a given 
  trajectory t. What the tuple represents (rgba, hsla, etc.) is up to client code. We don't 
  interpret it here, we just assign it according to the user-defined function. */
  void setColoringFunction(const std::function<Color(const std::vector<Vec2D>& t)> newFunc)
  { coloringFunction = newFunc; }
  // maybe the function should also take the initial value z0 and parameter p as arguments?



  // todo: setMaxNumIterations,
  // todo: set some presets with some well-known default fractals such as Mandelbrot, Nnewton


  //-----------------------------------------------------------------------------------------------
  /** \name Rendering */

  /** Given an array of locations, i.e. (x,y)-points in the plane, this function computes what the
  color should be at each such point according to the user defined iterationFunction, 
  coloringFunctions, etc. */
  std::vector<Color> render(const std::vector<Vec2D> locations, int maxNumIterations);
  // maybe use rsMatrix instead of std::vector - or just use raw pointers..hmm..no! rsMatrix
  // seems appropriate...or maybe have a function setRange(xMin, ...)
  // that seems more convenient



protected:



  std::function<Vec2D(Vec2D z, Vec2D p)> iterationFunction;


  std::function<Color(const std::vector<Vec2D>& t)> coloringFunction;
  /**< Used defined function that is supposed to take a trajectory of z-values as input and compute
  a color form that.  */

  std::vector<Vec2D> trajectory;
  /**< Records the trajectory of z for a given pixel. This is used as input to the coloring 
  algorithm. */

  double tol = 1.e-13;
  /**< Tolerance for convergence test.  */


};






#endif