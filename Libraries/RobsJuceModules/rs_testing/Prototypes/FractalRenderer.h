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

  /** Sets the function that determines the stopping criterion. It should take a trajectory as 
  input and return true, iff the iteration should stop. Examples of such criteria could be: 
  (1) the last element in the trajectory has a radius that exceeds some threshold, such that 
  divergence is  assured (used in the Mandelbrot fractal) or (2) the last element is very close to
  the second-to-last such that we can assume that we have converged to a fixed point (used in 
  Newton fractals). */
  void setStoppingCriterion(const std::function<bool(const std::vector<Vec2D>& t)> newFunc)
  { stoppingCriterion = newFunc; }

  /** Sets the range for the x- and y-coordinates. */
  void setCoordinateRange(double xMin, double xMax, double yMin, double yMax)
  { this->xMin = xMin; this->xMax = xMax; this->yMin = yMin; this->yMax = yMax; }

  /** Sets the size (width and height) of the image to be rendered in pixels. */
  void setImageSize(int width, int height) { w = width; h = height; }

  /** Sets an oversampling factor that is used in the rendering process. The final rendered image
  will then be downsampled using a local mean before the decimation. */
  void setOversampling(int newFactor) { oversample = newFactor; }
  // todo: maybe try using a median or some other rule

  /** Sets the maximum number of iterations to take until we stop. The iteration may stop earlier
  due to other stopping criteria such as convergence to a fixed point or divergence to infinity. */
  void setMaxNumIterations(int newMax) { maxIts = newMax; }

  // todo: 
  // -provide some factory functions for generating suitable iteration- and coloring functions
  // -provide a sort of setPreset(Preset newPreset) function that sets everything up for some 
  //  well-known default fractals such as Mandelbrot, Newton, etc. (these should use the factory
  //  functions)


  //-----------------------------------------------------------------------------------------------
  /** \name Rendering */

  /**  */
  rsImage<Color> render();
  // maybe use rsMatrix instead of std::vector - or just use raw pointers..hmm..no! rsMatrix
  // seems appropriate...or maybe have a function setRange(xMin, ...)
  // that seems more convenient



protected:


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  std::function<Vec2D(Vec2D z, Vec2D p)> iterationFunction;
  /**< User defined function f that is used to advance our iterate vector v[n] to the next time 
  step via v[n+1] = f(v[n]). */

  std::function<Color(const std::vector<Vec2D>& t)> coloringFunction;
  /**< User defined function that is supposed to take a trajectory of z-values as input and compute
  a color form that. */

  std::function<bool(const std::vector<Vec2D>& t)> stoppingCriterion;
  /**< user defined function for the stopping criterion ...tbc... */

  //std::vector<Vec2D> trajectory;
  /**< Records the trajectory of z for a given pixel. This is used as input to the coloring 
  algorithm. */

  double xMin = -1, xMax = +1, yMin = -1, yMax = +1;
  /**< Range for the coordinates  */

  int w = 480, h = 480;
  /** Image width and height in pixels. */

  int maxIts = 50;
  /** Maximum number of iterations. */

  int oversample = 1;
  /** Oversampling factor in the rendering. */

  //double tol = 1.e-13;
  /**< Tolerance for convergence test.  */

};

#endif