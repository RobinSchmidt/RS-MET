


/*

ToDo:

rsFractalImageRenderer:
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


-Implement a different kind of fractal renderer based on Lindenmayer systems - maybe 
 rsFractalVector/Drawing/LineRenderer 
 ...maybe the turtle can be extended to not only have a position and direction (i.e. velocity)
 but also a curvature (i.e. acceleration) and renders curved segments instead of straight lines

*/