template<class Tx, class Ty>
void rsNumericDerivative(const Tx *x, const Ty *y, Ty *yd, int N, bool extrapolateEnds)
{
  rsAssert(y != yd, "cannot be used in place yet, y and yd have to be distinct");

  Tx dxl, dxr, dx;
  Ty a, b; 

  for(int n = 1; n < N-1; n++)
  {
    dxl   = x[n] - x[n-1];
    dxr   = x[n+1] - x[n];
    dx    = dxl + dxr;
    yd[n] = dxr * (y[n]-y[n-1])/(dxl*dx) + dxl * (y[n+1]-y[n])/(dxr*dx);
  }
  // todo: save the left weight and use it as right weight in the next iteration (save one division
  // per iteration)


  if( extrapolateEnds == true )
  {
    a = (yd[2] - yd[1]) / (x[2] - x[1]);
    b = yd[1] - a * x[1];
    yd[0] = a*x[0] + b;
    a = (yd[N-2] - yd[N-3]) / (x[N-2] - x[N-3]);
    b = yd[N-3] - a * x[N-3];
    yd[N-1] = a*x[N-1] + b;
      // maybe this can be simplified by using rsInterpolateLinear
  }
  else
  {
    yd[0]   = (y[1]   - y[0])   / (x[1]   - x[0]);
    yd[N-1] = (y[N-1] - y[N-2]) / (x[N-1] - x[N-2]);
  }
}
// todo:
// -make it possible to use the function in-place, i.e. y and yd point to the same memory
//  ->avoid using y[n-1] on the right hand side
// -before making those optimization, move the non-optimized version to prototypes and write a unit
//  test to compare optimized to non-optimized
// -make a simplified version, assuming equidistant abscissa values (maybe assume unit distance, 
//  derivatives for any other h can then be found by scaling
// -maybe allow different template types for x, y, yd so it can be used for complex or multivariate 
//  data as well. in the latter case, x,y would be vectors (of possibly different dimensionality) 
//  and the derivative yd would be the Jacobian matrix at each datapoint
// -maybe write a function that computes the numeric derivative at one particular datapoint and
//  ideally also higher order derivatives at that point - that's more convenient to use because
//  client code does not need to have buffers for all derivatives
// -try another approach: fit a polynomial of arbitrary order to a number of datapoints around
//  the n and return the derivative of the poynomial at that point (may this be equivalent to the
//  approach above when using 3 points for a quadratic polynomial?)



template<class Tx, class Ty>
void rsNumericIntegral(const Tx *x, const Ty *y, Ty *yi, int N, Ty c)
{
  Tx xo; 
  Ty yo, zo, tmp;
  xo = x[0]; yo = y[0]; zo = c; yi[0] = zo;    // "old" values (at index n-1)
  for(int n = 1; n < N; n++) {
    tmp = zo + (x[n]-xo)*(y[n]+yo)*Ty(0.5);    // compute integral by trapezoidal rule
    xo = x[n]; yo = y[n]; zo = tmp;            // update integrator state variables
    yi[n] = tmp;                               // write integral to output array
    //rsAssert(rsIsFiniteNumber(tmp));
  }
}
// todo:
// -implement a higher order method by making use of (numeric) derivative information to 
//  approximate the integral by cubic segments - this may use the yi array first for the numeric
//  derivative values (after numeric derivative is adapted for in-place use) and then overwrite 
//  them with the integral values)
// -for this, obtain natural cubic spline interpolation coeffs for all segments, integrate them to
//  quartic segments, obatin the segment integrals by evaluating the quartics at the starts/ends
//  and add them all up
// -make a simplified version that doesn't need an x-array (assume distance 1 between x-values)
// -implement path-integration - the path is defined by an array of vectors (taking the role of x)
//  and there should be a function value associated with each vector passed in another array 
//  (taking the role of y). in the above formula, the x[n]-xo term should be replaced by
//  norm(x[n]-xo) where "norm" should be the Euclidean norm (function values are multiplied by the 
//  lengths of the path-segments in the summation)
// -write N-dimensional integration functions that return the amount of N+1 space contained in
//  some hyperblock between x1, x2 (both of dimensionality N)
//  -maybe we somehow need a function that takes in a function of N variables and returns a 
//   function of N-1 variables (maybe using std::function)


template<class Tx, class Ty>
Ty rsNumericIntegrator<Tx, Ty>::integrate(const std::function<Ty(Tx)>& f, Tx a, Tx b)
{
  rsError("Not yet implemented");

  return Ty(0);
}

// Ideas: 
// -let the user set the sample evaluation points by passing a pointer to an array of Tx
// -alternatively, the user may set just a number and then the object auotmatically generates
//  the sample points
// -for this automatic sample poitn generation, the user may select between different algorithms,
//  by default, we just choose them equidistantly
// -use a (cubic) natural spline based on the datapoints and compute the integral as sum over the
//  integrals of the spline segments (for this, we need to split the spline generator such that it
//  can spit out arrays of polynomial coefficients like:
//  void getCubicNaturalSplineCoeffs(Tx* x, int N, Ty* a, Ty* b, Ty* c, Ty* d);
//  ...this may be also useful for rsInterpolatingFunction