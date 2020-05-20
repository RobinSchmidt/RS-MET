//-------------------------------------------------------------------------------------------------
// Differentiation:

template<class Ty>
//template<class F>
void rsNumericDifferentiator<Ty>::hessian(
  const std::function<Ty(Ty*)>& f, Ty* x, int N, Ty* pH, const Ty& h)
{
  // compute N diagonal elements:
  rsMatrixView<Ty> H(N, N, pH);
  Ty fc = f(x);
  for(int i = 0; i < N; i++) {
    Ty ti  = x[i];
    x[i]   = ti + h; Ty fp = f(x);
    x[i]   = ti - h; Ty fm = f(x);
    H(i,i) = (fm - Ty(2)*fc + fp) / (h*h);
    x[i]   = ti; }

  // compute (N^2-N)/2 off-diagonal elements:
  for(int i = 0; i < N; i++) {
    for(int j = i+1; j < N; j++) {
      Ty ti = x[i];
      Ty tj = x[j];
      x[i] = ti + h; x[j] = tj + h; Ty fpp = f(x);         // f_++
      x[i] = ti + h; x[j] = tj - h; Ty fpm = f(x);         // f_+-
      x[i] = ti - h; x[j] = tj + h; Ty fmp = f(x);         // f_-+
      x[i] = ti - h; x[j] = tj - h; Ty fmm = f(x);         // f_--
      H(i,j) = H(j,i) = (fpp + fmm - fpm - fmp) / (4*h*h); 
      x[i] = ti;
      x[j] = tj; }}

  // The formula for the diagonal elements is just the regular central difference for a 2nd 
  // derivative for one coordinate at a time. The formula for the off-diagonal elements was derived
  // by considering a bivariate function f(x,y) and computing its partial derivative with respect 
  // to x using a central difference:
  //   f_x ~= (f(x+h,y) - f(x-h,y)) / (2*h)
  // and then using a central difference with repect to y on f_x:
  //   f_xy ~= (f_x(x,y+h) - f_x(x,y-h)) / (2*h)
  // and then generalizing in the obvious way from the bivariate to the multivariate case by 
  // replacing x,y with i,j
}


template<class Ty>
template<class Tx>
void rsNumericDifferentiator<Ty>::derivative(
  const Tx *x, const Ty *y, Ty *yd, int N, bool extrapolateEnds)
{
  rsAssert(y != yd, "Cannot be used in place yet, y and yd have to be distinct");
  rsAssert(N >= 3, "Arrays need to be at least of length 3");
    // hmm - that's for the extrapolation - if it's not used, length = 2 would also work, i think

  Tx dxl, dxr, dx;
  Ty a, b; 

  for(int n = 1; n < N-1; n++) {
    dxl   = x[n] - x[n-1];
    dxr   = x[n+1] - x[n];
    dx    = dxl + dxr;
    yd[n] = dxr * (y[n]-y[n-1])/(dxl*dx) + dxl * (y[n+1]-y[n])/(dxr*dx);
  }
  // todo: save the left weight and use it as right weight in the next iteration (save one division
  // per iteration)

  if( extrapolateEnds == true ) { 
    a = (yd[2]   - yd[1]  ) / (x[2]   - x[1]  ); b = yd[1]   - a*x[1];   yd[0]   = a*x[0]   + b;
    a = (yd[N-2] - yd[N-3]) / (x[N-2] - x[N-3]); b = yd[N-3] - a*x[N-3]; yd[N-1] = a*x[N-1] + b;
  } else {
    yd[0]   = (y[1]   - y[0])   / (x[1]   - x[0]  );
    yd[N-1] = (y[N-1] - y[N-2]) / (x[N-1] - x[N-2]); }
  // maybe the first branch can be simplified by using rsInterpolateLinear, like
  //   yd[0]   = lerp(x[0],   x[1],   yd[1],   x[2],   yd[2]  );
  //   yd[N-1] = lerp(x[N-1], x[N-2], yd[N-2], x[N-3], yd[N-3]);

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
// -yet another approach: "invert" the trapezoidal integration algorithm, i.e. run it backwards in
//  order to get a numerical integration routine that is the inverse operation to trapezoidal 
//  integration
//  -it may return a value: the integration constant to be used, to get the original data back
// -make a numeric derivative routine that is the inverse of the trapezoidal integrator
//  rsDifferentiateTrapezoidal

template<class Ty>
void rsNumericDifferentiator<Ty>::stencilCoeffs(const Ty* x, int N, int d, Ty* c)
{
  rsAssert(d < N, "Stencil width must be greater than derivative order.");

  // establish matrix:
  Ty** A; rsMatrixTools::allocateMatrix(A, N, N);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      A[i][j] = pow(x[j], i);

  // establish right-hand-side vector:
  std::vector<Ty> rhs(N);
  rsFill(rhs, Ty(0));
  rhs[d] = rsFactorial(d);

  // compute coeffs by solving the linear system and clean up:
  rsLinearAlgebra::rsSolveLinearSystem(A, &c[0], &rhs[0], N);
  rsMatrixTools::deallocateMatrix(A, N, N);
  // In practice, the resulting coefficients have to be divided by h^d where h is the step-size and
  // d is the order of the derivative to be approximated. The stencil offsets in x are actually 
  // multipliers for some basic step-size h, i.e. a stencil -2,-1,0,1,2 means that we use values
  // f(x-2h),f(x-h),f(x),f(x+h),f(x+2h) to approximate the d-th derivative of f(x) at x

  // see: http://web.media.mit.edu/~crtaylor/calculator.html for explanation of the algorithm

  // todo: use rsMatrix, write unit test

  // One suboptimal thing is that we use floating point arithmetic where rational numbers could 
  // have been used instead to produce exact fractions...maybe a later refinement can do that - at
  // the end of the day, this computation is probably done offline anyway to obtain coeffs that 
  // will be hardcoded in some specific numerical derivative approximator...well - actually it's a
  // template and we could just isntantiate it with some rational number class - maybe try it in 
  // the testbed
}

/*

todo:
 -maybe it's better to not use x += h, x +- 2h, x +- 3h but instead 
  x +- h/3, x +- 2h/3, x +- h such that the total width of the stencil stays the same, 
  regardless of how many stencil-points we use. the goal is that the optimal choice of h 
  depends only on the problem, not on the number of stencil-points
 -try stencils where the points are not distributed equidistantly but maybe exponentially, for
  example: x +- h/4, x +- h/2, x +- h for a 5-point stencil
 -to avoid numerical error, it is desirable that the offsets (2, 2h, 3h, ..) are exactly 
  representable - also the final divisors should be exactly representable - maybe (inverse)
  powers of two are a good choice - at least, for the basic 3-point stencil x +- h, x +- 2h
  where the divisor is 1/(2h)
 -test it with some standard functions like exp, log, sin, cos, tan, 1/x, 1/(1+x^2) - maybe 
  polynomials (in which case we should get exact results, if the number of stencil points 
  matches the degree)
 -compute gradient of a multivariate function
  -this function should take a raw array (i.e. pointer) as input
 -derive and implement function to compute a Hessian matrix of a multivariate function

-figure out the accuracy experimentally - maybe this can be done by testing, how high 
 degree a polynomial can be such that we still get perfect results - i think, a 5-point stencil
 should/ be perfect for polynomials up to 5th degree (it's based on an interpolating polynomial
 of degree 5) 


 https://stackoverflow.com/questions/1174169/function-passed-as-template-argument
 https://stackoverflow.com/questions/10871100/pass-a-function-as-an-explicit-template-parameter
*/


//-------------------------------------------------------------------------------------------------
// Integration:

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
// -have a simpler version riemannSum which uses the midpoint of each interval as evaluation point 
//  - oh,  wait - this function is data-based and not based on a function that we may evaluate...but 
//  such an integration function should also be implemented - this one can then do Riemann sums or 
//  trapezoidal rule (and maybe higher order rules as well - Simpson, etc.)
// ...but we may still compute lower and upper Riemann sums


/*
template<class Tx, class Ty>
Ty rsNumericIntegrator<Tx, Ty>::integrate(const std::function<Ty(Tx)>& f, Tx a, Tx b)
{
  rsError("Not yet implemented");

  return Ty(0);
}
*/
// Ideas: 
// -let the user set the sample evaluation points by passing a pointer to an array of Tx
// -alternatively, the user may set just a number and then the object auotmatically generates
//  the sample points
// -for this automatic sample point generation, the user may select between different algorithms,
//  by default, we just choose them equidistantly but other choices may give better results, see: 
//  https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
//  https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature
//  https://en.wikipedia.org/wiki/Chebyshev_nodes
//  the best choice may depend on the problem
// -use a (cubic) natural spline based on the datapoints and compute the integral as sum over the
//  integrals of the spline segments (for this, we need to split the spline generator such that it
//  can spit out arrays of polynomial coefficients like:
//  void getCubicNaturalSplineCoeffs(Tx* x, int N, Ty* a, Ty* b, Ty* c, Ty* d);
//  ...this may be also useful for rsInterpolatingFunction

/*
In the Princeton Companion to Applied Mathematics is a formula for approximating the derivative of
analytic functions that produce real outputs for real inputs. It's
  f'(x) ~= Im( f(x + i*h) ) / h
It has an error of O(h^2) - opposed to O(h) when using f'(x) ~= f(x + i*h) / h. Why does this work?
Is it because the imaginary part of f(x) is zero?




Maybe rename to NumericAnalysis and include the interpolation stuff into this file as well 
because some interpolation stuff depends on numeric derivatives but some numeric derivatives/
integration stuff may depend on interpolation and if we templatize the functions, we need to 
take care that everything is defined before it gets used.


*/