

rsFractalImageRenderer::rsFractalImageRenderer()
{
  // todo: 
  // -set it up in such a way as to render the Madelbrot fractal by calling setPreset with 
  //  the appropriate setting

  int dummy = 0;
}


rsImage<rsFloat32x4> rsFractalImageRenderer::renderRaw()
{
  rsAssert(oversample == 1, "Oversampling not yet implemented");

  int wo = w * oversample;
  int ho = h * oversample;
  rsImage<Color> img(wo, ho);         // image to render on, may be oversampled
  std::vector<Vec2D> t;               // trajectory of iterates
  t.reserve(maxIts);
  for(int j = 0; j < ho; j++)
  {
    for(int i = 0; i < wo; i++)
    {
      t.clear();
      double x = rsLinToLin((double)i, 0.0, double(wo-1), xMin, xMax); // optimize!
      double y = rsLinToLin((double)j, double(ho-1), 0.0, yMin, yMax); // dito
      Vec2D z(x, y);   // vector iterates
      Vec2D p(x, y);   // fixed parameter
      t.push_back(z);  // store 0-th iterate in trajectory
      for(int k = 0; k < maxIts-1; k++)
      {
        z = iterationFunction(z, p);
        if(stoppingCriterion(t))
          break;
        else
          t.push_back(z);
      }
      img(i, j) = coloringFunction(t);
    }
  }

  if(oversample == 1)
    return img;
  //else
  //  return rsImageProcessor::decimate(img, oversample, oversample, true); // true: use average


  return img;

  //// Decimate:
  //rsImage<Color> imgD(w, h);
  //return imgD;
  //// todo: factor out int rsImageProcessor::decimate such we can just do
  //// return rsImageProcessor::decimate(img, oversample, oversample);
}


/*

ToDo:

rsFractalImageRenderer:
-maybe we should have different variants of render() using different stopping criteria, e.g. 
 convergence, divergence ...or maybe that should be controlled by a setter
-Iteration Modes: 
 -NewtonComplex (as we do here)
 -NewtonVector (general 2D vector-field, using inverse Jacobian for Newton steps)
 -Holomorphic: general holomorphic iteration function, not necessarily based on Newton iteration
 -VectorField: general 2D -> 2D function
  ...actually they are all just general iteration functions - maybe use a general Vec2D -> Vec2D
  function in the iteration and provide factory methods that can create appropriate std::function
  objects that realize the special modes above
-Pixel Interpretation Modes:
 -initial value (as in Newton and Julia fractals)
 -parameter (as in Mandelbrot fractals)
 -maybe the iteration function should take a single 2D parameter (which may be ignored, if not
  needed)
-For each pixel, it should record the trajectory and then call a user-defined function with it
 which will determine the pixel "color" as float-quadruple, the interpretation of which is 
 relegated to later stages of the code
-It should be easy to vectorize, i.e. use with rsSimdVector - maybe we should use simd right 
 from the start
 -Maybe in this context, it's convenient, when it operates on a vector of (x,y)-coordinates 
  rather than on an image with on-the-fly computation of coordinates

-maybe allow the iteration function to depend explicitly on the iteration number, i.e.
   x[n+1] = f(x[n], n)
 this can be used to let the characteristics of the fractal change as function of the zoom level

-Provide some feedback about the state of the computation, maybe in % maybe with 2 decimal digits 
 after the point

-Implement a different kind of fractal renderer based on Lindenmayer systems - maybe 
 rsFractalVector/Drawing/LineRenderer 
 ...maybe the turtle can be extended to not only have a position and direction (i.e. velocity)
 but also a curvature (i.e. acceleration) and renders curved segments instead of straight lines

*/