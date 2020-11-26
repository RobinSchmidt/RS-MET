//-------------------------------------------------------------------------------------------------
// Differentiation:

template<class T>
template<class F>
void rsNumericDifferentiator<T>::hessian(const F& f, T* x, int N, T* pH, const T* h)
{
  // wrap raw pointer into matrix-view, evaluate f at x:
  rsMatrixView<T> H(N, N, pH);
  T fc = f(x);

  // compute N diagonal elements (2 evals for each):
  for(int i = 0; i < N; i++) {
    T ti   = x[i];
    x[i]   = ti + h[i]; T fp = f(x);
    x[i]   = ti - h[i]; T fm = f(x);
    H(i,i) = (fm - T(2)*fc + fp) / (h[i]*h[i]);
    x[i]   = ti; }

  // compute (N^2-N)/2 off-diagonal elements (4 evals for each):
  for(int i = 0; i < N; i++) {
    for(int j = i+1; j < N; j++) {
      T ti = x[i];
      T tj = x[j];
      x[i] = ti + h[i]; x[j] = tj + h[j]; T fpp = f(x);          // f_++
      x[i] = ti + h[i]; x[j] = tj - h[j]; T fpm = f(x);          // f_+-
      x[i] = ti - h[i]; x[j] = tj + h[j]; T fmp = f(x);          // f_-+
      x[i] = ti - h[i]; x[j] = tj - h[j]; T fmm = f(x);          // f_--
      H(i,j) = H(j,i) = (fpp + fmm - fpm - fmp) / (4*h[i]*h[j]);
      x[i] = ti;
      x[j] = tj; }}

  // # evaluations: 1 + 2*N + 4*(N^2-N)/2 = 2*N^2 + 1
  // maybe return that number - maybe that should all functions do - that information can be used
  // in higher-level algorithms such as those for numerical minimization when they use these 
  // functions from here - because in these, we are really interested in that number to assess 
  // their efficiency

  // The formula for the diagonal elements is just the regular central difference for a 2nd 
  // derivative for one coordinate at a time. The formula for the off-diagonal elements was derived
  // by considering a bivariate function f(x,y) and computing its partial derivative with respect 
  // to x using a central difference:
  //   f_x ~= (f(x+h,y) - f(x-h,y)) / (2*h)
  // and then using a central difference with respect to y on f_x:
  //   f_xy ~= (f_x(x,y+h) - f_x(x,y-h)) / (2*h)
  // and then generalizing in the obvious way from the bivariate to the multivariate case by 
  // replacing x,y with i,j
}
// ToDo:
// -can we make the formula for the off-diagonal elements more accurate by using fc = f(x)
//  -that would come at (almost) no cost because that value never needs to be re-evaluated
//  -maybe H(i,i) and H(j,j) could also be used
//  -maybe we could compute the A..F coeffs of a quadratic form (i.e. conic section) and
//   evaluate its partial derivatives? fpp,fpm,fmp,fmm,fc would determine 5 coeffs - maybe the
//   constant coeff F is irrelevant
// see:
// https://en.wikipedia.org/wiki/Hessian_matrix#Use_in_optimization
// https://en.wikipedia.org/wiki/Quasi-Newton_method

template<class T>
template<class Tx>
void rsNumericDifferentiator<T>::derivative(
  const Tx *x, const T *y, T *yd, int N, bool extrapolateEnds)
{
  rsAssert(y != yd, "Cannot be used in place yet, y and yd have to be distinct");
  rsAssert(N >= 3, "Arrays need to be at least of length 3");
    // hmm - that's for the extrapolation - if it's not used, length = 2 would also work, i think

  Tx dxl, dxr, dx;
  T  a, b; 

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


template<class T>
void rsNumericDifferentiator<T>::gradient2D(
  const rsGraph<rsVector2D<T>, T>& mesh, const T* u, int i, T* u_x, T* u_y)
{
  // Algorithm:
  // The algorithm is based on the fact that the directional derivative into the direction of an 
  // arbitrary vector a can be expressed as: D_a(u) = <g,a> where D_a(u) denotes the directional
  // derivative of the scalar function u(x,y) in the a-direction, g denotes the gradient of u and
  // <g,a> denotes the scalar product of vectors g and a. For the vertex i, we compute numerical 
  // estimates of the directional derivatives into the directions of all its neighbors and set up 
  // the above equation. When a vertex has more than 2 neighbors (which is the typical case), we 
  // have more equations than degrees of freedom (u_x, u_y), so the system is overdetermined and we
  // compute a weighted least squares solution where the weights are taken to be the values stored 
  // at the edges between the respective vertices. Meaningful weights could be either all ones 
  // (unweighted) or - probably better - something inversely proportional to (some power of) the 
  // distance between the vertices. The rationale behind this is that we expect the values obtained 
  // by the formula to be further off from the true values, the greater the distance between the 
  // vertices, but some experimentation for what kind of weighting gives the most accurate results 
  // is encouraged (for example, if Euclidean or Manhattan distance gives better results, etc.).
  // First experimental results indicate that w = pow(d, -p) seems to give most accurate results 
  // where d ist the (Euclidean) distance and p is the number of neighbors.

  using Vec2       = rsVector2D<T>;           // shorthand for convenience
  const Vec2& vi   = mesh.getVertexData(i);   // vertex i, at which we calculate the derivative
  int numNeighbors = mesh.getNumEdges(i);     // number of neighbors of vertex vi

  // If vi has no neighbors at all, we assign zeros to the partial derivatives:
  if(numNeighbors == 0) { *u_x = *u_y = T(0); return; }

  // If vi has only one neighbor, we have only one equation for our 2 degrees of freedom 
  // u_x, u_y, so there are infinitely many solutions. I'm not sure, if the minimum norm 
  // solution is the best thing to compute in such a case, but we should at least compute 
  // *something* that is *a* solution to the equation and the minimum norm solution is the only
  // meaningfully distinguished choice (right? or wrong?):
  if(numNeighbors == 1) {
    int j = mesh.getEdgeTarget(i, 0);
    const Vec2& vj = mesh.getVertexData(j);
    Vec2 dv = vj   - vi;                      // difference vector
    T    du = u[j] - u[i];                    // difference in function value
    rsLinearAlgebra::solveMinNorm(dv.x, dv.y, du, u_x, u_y);
    return; }

  // The typical case is that vi has >= 2 neighbors. In this case, we have either a critically
  // determined (numNeighbors == 2) or an overdetermined (numNeighbors > 2) system and we compute
  // a weighted least squares solution (which, in the case of a critically determined system, 
  // happens to be the exact solution...right?):
  static const T z = T(0);      // zero
  rsMatrix2x2<T> A(z,z,z,z);    // maybe rename to M = ATA (== A^T * A in most textbooks)
  Vec2 b(z,z), g;               // maybe rename b to Mb (== A^T * b in textbooks)
  for(int k = 0; k < numNeighbors; k++)       // loop over neighbors of vertex i
  {
    // Retrieve or compute intermediate variables:
    int j = mesh.getEdgeTarget(i, k);         // index of current neighbor of vi
    const Vec2& vj = mesh.getVertexData(j);   // current neighbor of vi
    Vec2 dv = vj - vi;                        // difference vector

    // Accumulate least-squares matrix and right-hand-side vector:
    T du = u[j] - u[i];                       // difference in function value
    T w  = mesh.getEdgeData(i, k);            // weight in weighted least squares
    A.a += w * dv.x * dv.x;
    A.b += w * dv.x * dv.y;
    A.d += w * dv.y * dv.y;
    b.x += w * dv.x * du;
    b.y += w * dv.y * du;
  }
  A.c = A.b;  // A.c is still zero - make A symmetric

  // Compute gradient that best explains the measured directional derivatives in the least 
  // squares sense and store it in outputs:
  rsMatrix2x2<T>::solveSave(A, g, b);  // g is the gradient vector that solves A*g = b
  *u_x = g.x; 
  *u_y = g.y;
}
// todo:
// -the "solveSave" call could be optimized - maybe we don't even need an explicit matrix and/or 
//  may make use of the symmetry of A (maybe a special solveSymmetric function could be used)
// -perhaps, solveSave should detect a zero determinant and switch between computing
//  a least-squares solution in case of an inconsistent RHS and a minimum norm solution in case
//  of a consistent RHS
// -do we need special treatment for 2 neighbours? I don't think so - the solution of the least 
//  squares problem should reduce to the exact solution if the number of equations equals the
//  number of unknowns
// -the matrix A can actually be precomputed (and perhaps stored as vertex data), likewise the
//  w*dv.x, w*dv.y coeffs used to establish the right-hand-side vector b - maybe the matrix 
//  elements a,b,c (or better: the elements of the inverse matrix) and the coeffs can be stored in
//  a data-structure rsMeshStencil2D at the nodes. Ultimately, we may just arrive at a scheme that
//  just forms a weighted sum of stencil values (i.e. the value at vertex i and its neighbors)

// See also: A Guide to Numerical Methods for Transport Equations (Dmitri Kuzmin):
// http://www.mathematik.uni-dortmund.de/~kuzmin/Transport.pdf  
// ...on page 20, it mentions fitting a (2D?) Taylor polynomial to the neighbors on 
// unstructured meshes. How does this relate to the directional derivative approach used here? Is 
// it equivalent? He says, fitting the Taylor polynomial is expensive and difficult to implement. 
// Using directional derivatives seems reasonably efficient and straightforward to implement. I 
// think, if the mesh is in fact rectangular, my approach reduces to using a central difference.

// Papers:
// https://www.researchgate.net/publication/254225242_Development_of_Irregular-Grid_Finite_Difference_Method_IFDM_for_Governing_Equations_in_Strong_Form
// https://www.semanticscholar.org/paper/Development-of-Irregular-Grid-Finite-Difference-(-)-GEORGE/0b8bcb2afdfee4d2fd8efc52f499a9d59f742f77
// http://www.fluidmal.uma.es/pdfs/JCOMP_2005.pdf
// https://www.researchgate.net/figure/Directional-derivatives-computed-using-one-sided-finite-differences-55-and-the_fig1_258919790

// Oh - the method is actually very similar to this:
// The finite difference method at arbitrary irregular grids and its application in applied 
// mechanics
// https://www.sciencedirect.com/science/article/abs/pii/0045794980901492
// the only difference is the use of directional derivatives in the derivation instead of an 
// interpolation polynomial

template<class T>
void rsNumericDifferentiator<T>::gradient2D(const rsGraph<rsVector2D<T>, T>& mesh, const T* u, 
  T* u_x, T* u_y)
{
  rsAssert(u_x != u && u_y != u, "rsNumericDifferentiator::gradient2D does not work in place");
  for(int i = 0; i < mesh.getNumVertices(); i++)
    gradient2D(mesh, u, i, &u_x[i], &u_y[i]);
}

template<class T>
void rsNumericDifferentiator<T>::laplacian2D(const rsGraph<rsVector2D<T>, T>& mesh, const T* u,
  T* L, T* w)
{
  int N = mesh.getNumVertices();
  T *t1 = &w[0], *t2 = &w[N], *t3 = &w[2*N];
  gradient2D(mesh, u,  L,  t1);   // L  = u_x,  t1 = u_y
  gradient2D(mesh, L,  t2, t3);   // t2 = u_xx, t3 = u_xy
  gradient2D(mesh, t1, t3, L );   // t3 = u_yx, L  = u_yy
  for(int i = 0; i < N; i++)
    L[i] += t2[i];                // L = u_xx + u_yy
}

template<class T>
void rsNumericDifferentiator<T>::laplacian2D_2(const rsGraph<rsVector2D<T>, T>& mesh, 
  const std::vector<T>& u, std::vector<T>& L)
{
  rsWarning("rsNumericDifferentiator::laplacian2D_2 is still under construction.");
  using Vec2 = rsVector2D<T>;
  int N = mesh.getNumVertices();
  rsAssert((int) u.size() == N);
  rsAssert((int) L.size() == N);
  for(int i = 0; i < N; i++) {                     // loop over all vertices
    Vec2 vi = mesh.getVertexData(i);               // current vertex location
    T uSum = T(0);                                 // weighted sum of neighbors
    T wSum = T(0);                                 // sum of weights
    for(int k = 0; k < mesh.getNumEdges(i); k++) { // loop over vi's neighbors
      int j = mesh.getEdgeTarget(i, k);            // index of current neighbor
      T w   = mesh.getEdgeData(  i, k);            // weight in weighted sum of neighbors
      Vec2 vj = mesh.getVertexData(j);             // location of current neighbor
      Vec2 dv = vj - vi;                           // difference vector
      T d2 = dv.x*dv.x + dv.y*dv.y;                // squared distance between vi and vj
      uSum += w*(u[j]-u[i])/d2;                    // accumulate sum of ...
      wSum += w;                                   // accumulate sum of weights
    }
    L[i] = T(4)*uSum/wSum;
  }
}
// this is still wrong - it works well only, if all the distances d2 are the same

template<class T>
void rsNumericDifferentiator<T>::stencilCoeffs(const T* x, int N, int d, T* c)
{
  rsAssert(d < N, "Stencil width must be greater than derivative order.");

  // establish matrix:
  T** A; rsMatrixTools::allocateMatrix(A, N, N);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      A[i][j] = pow(x[j], i);

  // establish right-hand-side vector:
  std::vector<T> rhs(N);
  rsFill(rhs, T(0));
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


 for the gradient:

 should the gradient be of type Tx or Ty? or maybe there should be just one type? what if
 y is a vector? then f would take an N-dim vector and produce a vector of possibly other 
 Consider z = f(x,y) where x,y are complex and z is real (for example f(x,y) = abs(x*y))
 Then, the gradient of would be given by grad(f) = (df/dx, df/dy) where 
    df/dx = lim_{h->0} (f(x+h) - f(h)) / h
 In this case, the stepsize h would also be complex, the numerator would be real (as a difference
 of real numbers) and the denominator would be complex, so the overall value df/dx would be 
 complex, which is Tx. On the other hand, if x,y are real and z is complex (for example, 
 f(x,y) = x + i*y), then the df/dx would also be complex (as before), but this time, this 
 corresponds to Ty. Maybe this makes really only sense, when we require Tx == Ty. Or maybe the 
 elements of the gradient should have their own type Tg, so it can be decided on instantiation, 
 which one it should be? Then, g[n] = (Tg(fp)-Tg(fm)) / (2*Tg(h));
 but we could also use:
    df/dx = lim_{dx->0} (f(x+dx) - f(h)) / |dx|
 this would be suitable, if x and dx are vectors - we take the norm in the denominator
 also: in gradient descent, we need the input vector x and the gradient vector to be of the same
 type -> maybe just use one type to start with, generalize when it becomes necessarry

 dimensionality - would this function then compute the Jacobian? i think, it would be natural, if
 it would -> try it using rsVector2D for Ty

 see also:
 https://github.com/wesselb/fdm

 for templatizing the function-type, see:
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