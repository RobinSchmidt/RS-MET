template<class T>
void rsLine2D<T>::twoPointToImplicit(T x0, T y0, T x1, T y1, T* A, T* B, T* C, bool normalize)
{
  *A = y0 - y1; // == -dy
  *B = x1 - x0; // ==  dx
  *C = -(x0 * *A + y0 * *B);
  if(normalize)
  {
    T s = 1 / sqrt(*A * *A + *B * *B);
    *A *= s;
    *B *= s;
    *C *= s;
  }
}

template<class T>
void rsLine2D<T>::reflectPointInLine(T x, T y, T A, T B, T C, T *xr, T *yr)
{
  T d = 2*(A*x + B*y + C) / (A*A + B*B);
  *xr = x - A*d;
  *yr = y - B*d;

  // formula taken from:
  // https://math.stackexchange.com/questions/1013230/how-to-find-coordinates-of-reflected-point
  // its also in Salomon's "Computer Graphics...", page 71, Eq. 3.12

  // maybe make a version that assumes A^2 + B^2 = 1 so we can get rid of the division
}

/*
Ideas:

-Provide functions to "construct" new points from a set of given points, as in the "synthetic
 geometry" of the ancient greeks. There, the allowed things to construct new points was to connect 
 any pair from a given set of points with a straight line (also, to continue the line as long as 
 desired) and taking distances between points with a compass and drawing circles centered at any of 
 the given points with radius being any of the distances that can be picked up somewhere. Any 
 intersections of these so formed lines and/or circles could be used as new points and added to the
 set (and the process could be applied recursively to the enlarged set).
-We could allow more constructions, like drawing an ellipse with its foci being two of the points 
 and the sum of the two strings being any of the distances that we can get. Or maybe parabolas and 
 hyperbolas with foci and directrix obtained from existing points/lines. Anything for which we can
 solve the intersection-equations in closed form. We could also find circles that pass through 3 
 given points (and maybe ellipses that pass through 4)?

*/