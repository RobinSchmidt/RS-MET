#ifndef RAPT_NUMERICCALCULUS_H
#define RAPT_NUMERICCALCULUS_H

//=================================================================================================

/** A class for computing numerical approximations to derivatives of functions. In some cases, the
functions are assumed to be given as a function object (aka functor), in other cases as an array of 
datapoints. The meaning of the 3 template parameters is:
  T:  general type for all values or, if input/output types of the function to be differentiated 
      are different, type of the output (ordinate type)
  Tx: type of the input to the function, if different from T (abscissa type)
  F:  type of the function object (if applicable)
T has to be given when the class template is instantiated, Tx and F may be given or inferred by 
the compiler later, when a particular member function is instantiated. This is because T can not 
be always inferred from the function call because T may only appear as output type in the function 
signature (see the function "derivative", for example). For example, Tx could be "double", T could 
be "rsVector2D<double>" and F could be "std::function<rsVector2D<double>(double)>". This would 
apply to functions that take real scalars (double) as input and produce 2D vectors (of real 
numbers) as output. Such functions define the parametric equation of 2D curves. They depend on a 
scalar parameter t (which is often interpreted as time) and produce a 2D vector as output for 
each t. In this case, the meaning of the first derivative would be the velocity vector and the 
second derivative would be the acceleration vector. Such derivatives are used a lot in the 
diffential geometry of curves (in the plane, but generalizes directly to higher dimensions), 
for example. 

If the situation is reversed, i.e. the input to the function is a vector and the output is a 
scalar, then we are dealing with a scalar field rather than a parametric curve, and the type of the 
derivative is actually also a vector (the gradient) and so is of the same type as the input, 
whereas in the parametric curve case, the derivative type was the same as the output...tbc
...also, if both input and output are vectors (vector fields), the type of the derivative is a 
matrix (the Jacobian)....tbc
.
 
 References:
   (1) http://web.media.mit.edu/~crtaylor/calculator.html  
   (2) https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences   */

template<class T>
class rsNumericDifferentiator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Functor derivatives

  /** Numeric approximation of the first derivative of function f at the value x with approximation
  step-size h. Uses a central difference with 2 evaluations of f and is 2nd order accurate in h. */
  template<class Tx, class F>
  static T derivative(const F& f, const Tx& x, const Tx& h)
  {
    return (f(x+h) - f(x-h)) / (Tx(2)*h);
  }
  // maybe rename to firstDerivative or firstDerivativeCentral
  // If, for example, Tx == T == real (float or double) and F is a function from real to real, this 
  // is the most basic numerical derivative computation. But as it stands, the function makes sense 
  // also when Tx is a scalar and T is a vector type. It would then compute the velocity vector of 
  // the parametric curve defined by f. But what if Tx is a vector and T is a scalar, i.e. f 
  // defines a scalar field - would the function make any sense in this case? I think, it probably 
  // should compute the gradient in this case - but would it? No, because the output is a scalar. 
  // What if input and output of f are vectors - the "derivative" would be the Jacobian matrix - in 
  // this case, this function makes no sense, i think...maybe we should have 3 template parameters: 
  // type of input of f, type of output of f, type of derivative...  maybe h should be of the 
  // derivative type then? no - we can't divide a vector by a vector - we could perhaps define 
  // division of a scalar by a vector as being element-wise...hmm...not so easy to design the most 
  // flexible and still convenient API...maybe cook up some applications that need such 
  // functionality to see it in action...
  // What if F is a function from complex to complex? I think, it would just work - but then h 
  // would also be complex number - which may not be a problem.


  /** Numeric approximation of the second derivative using 3 evaluations of f. 2nd order accurate 
  in h. */
  template<class Tx, class F>
  static T secondDerivative(const F& f, const Tx& x, const Tx& h)
  {
    return (f(x-h) - Tx(2)*f(x) + f(x+h)) / (h*h);
    // the Tx(2) may seem weird because it multiplies the result f(x) which is of type T, so one 
    // might expect T(2) here. However, if Tx is a scalar type and T is a vector type, this is
    // totally appropriate.
  }
  // Verify this and add to documentation:
  // I think, the accuracy is always equal to the number of evaluation points minus 1, for example 
  // here: 3 evaluation points -> 2nd order accurate in h. But: the formula for first derivative 
  // has only 2 evaluation points and is still 2nd order accurate. However, there is a concpetual 
  // 3rd evaluation point in the middle - but its coeff comes out as zero which is why the 
  // evaluation has been optimized away. So, it's the number of "conceptual" evaluation points 
  // minus 1. ...as said - i'm not totally sure about this - this needs to be verified...
  // ...seems like in the odd derivative formulas, the middle coeff always comes out as zero, 
  // which makes sense
  // make numerical experiments: plot approximation error of the various derivatives and formulas 
  // as function of h (maybe using a log-log scale)

  
  /** Numeric approximation of the third derivative using 4 evaluations of f. 3rd order accurate 
  in h (i think - verify - well, we have h^3 in the denominator, but that may not mean anything 
  because for the 5-point formula of the 1st derivative, we just have h^1 and it should be more 
  than 1st order accurate). */
  template<class Tx, class F>
  static T thirdDerivative(const F& f, const Tx& x, const Tx& h)
  {
    return (-f(x-Tx(2)*h) + Tx(2)*f(x-h) - Tx(2)*f(x+h) + f(x+Tx(2)*h)) / (Tx(2)*h*h*h);
  }
  // needs test
  // document the accuracy - maybe it's even 4th order accurate? with 5 "conceptual" evaluations?
  // coeffs found by:
  // http://web.media.mit.edu/~crtaylor/calculator.html
  // f_xxx = (-1*f[i-2]+2*f[i-1]+0*f[i+0]-2*f[i+1]+1*f[i+2])/(2*1.0*h**3)


  /** Computes 0th, 1st and 2nd derivative of f at x with stepsize h and stores them in f0, f1, f2 
  where we adopt the usual convention that the 0th derivative means just the function itself. This 
  is more efficient than using the separate functions. It uses 3 function evaluations whereas you 
  would need 6, if you would call the function itself (1 evaluation) and derivative (2 evaluations) 
  and secondDerivative (3 evaluations). The result is exactly the same - it just avoids to compute 
  some intermediate values values twice. */
  template<class Tx, class F>
  static void derivativesUpTo2(const F& f, const Tx& x, const Tx& h, T* f0, T* f1, T* f2)
  {
    T fp = f(x+h);  // "plus"
    T fm = f(x-h);  // "minus"
    *f0 = f(x);
    *f1 = (fp - fm) / (Tx(2)*h);
    *f2 = (fm - Tx(2)*(*f0) + fp) / (h*h);
  }
  // needs test

  /** Computes 0th, 1st, 2nd and 3rd derivative of f at x. Uses 5 evaluations of f with the 5-point 
  stencil: -2,-1,0,1,2. */
  template<class Tx, class F>
  static void derivativesUpTo3(const F& f, const Tx& x, const Tx& h, T* f0, T* f1, T* f2, T* f3)
  {
    // evaluate function at stencil points:
    T fm2 = f(x-2*h);  // "minus 2h", etc...
    T fm1 = f(x - h);
    T fc  = f(x    );  // centered ..get rid - assign to *f0 directly
    T fp1 = f(x + h);
    T fp2 = f(x+2*h);

    // form linear combinations to approximate derivatives:
    *f0 =                     fc;
    *f1 = ( fm2 -  8*fm1         +  8*fp1 - fp2) / (12*h);
    *f2 = (-fm2 + 16*fm1 - 30*fc + 16*fp1 - fp2) / (12*h*h);
    *f3 = (-fm2 +  2*fm1         -  2*fp1 + fp2) / (2*h*h*h);
  }
  // needs test
  // document the accuracy 
  // stencil: -2,-1,0,1,2
  // f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*h**1)
  // f_xx = (-1*f[i-2]+16*f[i-1]-30*f[i+0]+16*f[i+1]-1*f[i+2])/(12*1.0*h**2)
  // f_xxx = (-1*f[i-2]+2*f[i-1]+0*f[i+0]-2*f[i+1]+1*f[i+2])/(2*1.0*h**3)
  // maybe move to cpp file

  // todo: derivativesUpTo4 - this is as far as we may go with a 5-point stencil - for higher 
  // derivatives, we need more than 5 evaluation points


  /** Computes the partial derivative of the multivariate (N inputs) scalar function f with respect 
  to the n-th coordinate. It requires 2 evaluations of f. It's a bit inelegant that we can't 
  declare x as const, because we need to wiggle one of its elements in the internal computation. We 
  restore the exact same value before returning, so from the caller's I/O perspective, x can be 
  considered const. But such a perspective is valid only in situations, where the x-array is used 
  only in a single thread. The non-constness propagates out to any function that uses it, so it may 
  prevent some compiler optimizations elsewhere. We could be hacky and cast away the const but it's 
  probably a bad idea to do so. Preventing these optimizations is actually the point - not 
  preventing them could make the code non thread-safe. See:
  https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Res-casts-const
  https://visualstudiomagazine.com/articles/2016/04/26/dont-cast-away-const-in-cpp.aspx  */
  template<class F>
  static T partialDerivative(const F& f, T* x, int N, int n, const T h)
  {
    T t  = x[n];                  // temporary
    x[n] = t + h; T fp = f(x);    // fp = f(x0,x1,..,xn+h,..,x_M)
    x[n] = t - h; T fm = f(x);    // fm = f(x0,x1,..,xn-h,..,x_M)
    T f1 = (fp - fm) / (T(2)*h);  // df/dxn
    x[n] = t;                     // restore x[n]
    return f1;
  }

  /** Computes the partial derivative of the multivariate (N inputs) scalar function f with respect 
  to the n-th coordinate and also the diagonal element H(n,n) of the Hessian matrix, i.e. the 2nd 
  derivative with respect to coordinate n. It requires 3 evaluations of f. Here, like in 
  partialDerivative, the x-array is const from the caller's perspective but we can't declare it 
  const because we need to wiggle it temporarily internally. */
  template<class F>
  static void partialDerivativesUpTo2(const F& f, T* x, int N, int n, const T h, 
    T* f0, T* f1, T* f2)
  {
    *f0  = f(x);                             // f0 = f(x0,x1,..,x_M)    where M := N-1
    T t  = x[n];                             // temporary
    x[n] = t + h; T fp = f(x);               // fp = f(x0,x1,..,xn+h,..,x_M)
    x[n] = t - h; T fm = f(x);               // fm = f(x0,x1,..,xn-h,..,x_M)
    *f1  = (fp - fm) / (T(2)*h);             // df/dxn
    *f2  = (fm - T(2)*(*f0) + fp) / (h*h);   // d2f/dxn2
    x[n] = t;                                // restore input
  }

  /** Computes a central difference approximation of the gradient of a scalar-valued function 
  y = f(x) at the given N-dimensional position vector x by a central difference and stores the 
  result in g (also of length N). You also need to pass an array of approximation stepsizes along 
  the N directions in h. The function type F should take its input vector as a raw (pointer to an) 
  array of length N and return the scalar output. This function has 2*N evaluations of f. Note that 
  the input vector x is not const because we need to wiggle it internally, but at the end, the 
  original content of the x-vector will be restored (exactly - no roundoff error is involved) - so 
  it's quasi-const. */
  template<class F>
  static void gradient(const F& f, T* x, int N, T* g, const T* h)
  {
    rsAssert(x != g, "Can't be used in place");
    for(int n = 0; n < N; n++) 
      g[n] = partialDerivative(f, x, N, n, h[n]);
  }
  // move output variables to the end to make it consistent with the other functions like
  // derivativesUpTo2
  // The gradient gives the direction of steepest ascent of f and it's norm meausres, how fast f
  // would change when going along that direction
  // maybe make a function for the directional derivative - this is just the scalar product of
  // the gradient with a given vector (but that vector should be normalized to length 1, i think)

  /** Computes a numerical approximation of the Hessian matrix of the function f at the given 
  N-dimensional position vector x and writes the result int H which should be a pointer to the flat 
  data array of an a NxN matrix. The usage is similar to gradient. This function has 2*N^2 + 1 
  function evaluations of f. */
  template<class F>
  static void hessian(const F& f, T* x, int N, T* H, const T* h);
  // The Hessian matrix measures, how fast the gradient changes, when moving in a particular 
  // direction...explain this some more....

  /** Approximates the matrix-vector product H * v of the Hessian matrix H of f evaluated at 
  position x and an arbitrary given vector v and writes the result into Hv. This is the same as 
  v^T * H due to the symmetry of the Hessian matrix. The approximation is based on two gradient 
  evaluations at position vectors x + k*v and x - k*v for some scalar approximation stepsize k. The 
  array h is - as usual - the vector of stepsizes to numerically approximate the gradient itself. 
  The workspace must be of length 2*N. The two gradient evaluations each take 2*N evaluations of f, 
  so the function evaluates f 4*N times. As a side-note, when v is the i-th coordinate vector, i.e. 
  all zeros except v[i] = 1, then this will compute the i-th row (or column) of the Hessian matrix
  itself. But it's not algorithmically attractive to approximate the Hessian that way - compared
  to the implementation of the actual hessian function, it needs more evaluations of f, needs an
  additional memory workspace, needs the user to specify an additional parameter k and is 
  numerically less accurate (i think - i tried it with one specific choice of k in which case the
  estimate was less accurate).  For details of the general idea, see... */
  template<class F>
  static void hessianTimesVector(const F& f, T* x, T* v, int N, T* Hv, const T* h, T k, 
    T* workspace)
  {
    T *gp =  workspace;
    T *gm = &workspace[N];
    for(int n = 0; n < N; n++)  
      Hv[n] = x[n] + k*v[n];               // Hv temporarily used for x + k*v
    gradient(f, &Hv[0], N, &gp[0], h);     // gp: gradient at x + k*v
    for(int n = 0; n < N; n++)  
      Hv[n] = x[n] - k*v[n];               // Hv temporarily used for x - k*v
    gradient(f, &Hv[0], N, &gm[0], h);     // gm: gradient at x - k*v
    for(int n = 0; n < N; n++)
      Hv[n] = (gp[n] - gm[n]) / (2*k);     // Hv ~= (grad(x+k*v) - grad(x-k*v)) / (2*k)
  }

  // todo: compute numeric Jacobians


  //-----------------------------------------------------------------------------------------------
  // \name Data derivatives

  /** Cannot be used in-place yet: y and yd have to be distinct!
  Given an array of strictly monotonically increasing but not necessarily equidistant abscissa
  values in x and corresponding function values in y, this function fills the array yd with a
  numeric approximation of the derivative for each x value. All arrays are of length N. To
  compute the numeric derivative, we use a weighted average of the difference quotients left and
  right to the data point:

                y[n] - y[n-1]          y[n+1] - y[n]
  yd[n] = wL * --------------- + wR * ---------------
                x[n] - x[n-1]          x[n+1] - x[n]

  where the weights for the left and right difference quotients are determined by the distances
  dxL = x[n]-x[n-1] and dxR = x[n+1]-x[n] as wL = dxR/(dxL+dxR), wR = dxL/(dxL+dxR), such that
  the closer the x-axis value x[n-1] is to x[n], the more weight is given for the left quotient
  and vice versa. The weights add up to unity. If extrapolateEnds is true, the function will use
  linear extrapolation of the inner derivative values for the endpoints yd[0] and yd[N-1],
  otherwise it will use the (divided) forward difference at 0 and the backward difference at
  N-1. In a test with a sine function, the extrapolation gave more accurate results at the
  endpoints compared to simple differences, so it's probably better to use extrapolation. */
  template<class Tx>
  static void derivative(const Tx *x, const T *y, T *yd, int N, bool extrapolateEnds = true);
  // maybe rename to weightedCentralDifference

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Computes the stencil coefficients for a finite difference approximation of a derivative 
  according to: http://web.media.mit.edu/~crtaylor/calculator.html (todo: explain this here)
  Inputs: 
    x: array of normalized distances from the approximation point. "Normalized" means in this 
       context, that the step-size h is not yet included. So, if the x-array is given by: 
       x = [-2,-1,0,1,2], it means, we want to use the 5-point stencil: x-2h,x-h,x,x+h,x+2h.
    N: length of x
    d: derivative that should be approximated, e.g. d=2 for the 2nd derivative, must be less than N
  Output:
    c: The normalized coefficients, by which the function values f(x + k*h) must be multiplied 
       (k being one of the values from the x-array). Here, "normalized" means that the result must 
       be divided by h^d. The length of this array must also be N.  */
  //template<class T>
  static void stencilCoeffs(const T* x, int N, int d, T* c);


  //-----------------------------------------------------------------------------------------------
  // \name Convenience functions

  template<class F> 
  static std::vector<T> gradient(const F& f, std::vector<T>& x, const T* h)
  { std::vector<T> g(x.size()); gradient(f, &x[0], (int) x.size(), &g[0], h); return g; }
  // allocates


  template<class F>
  static void gradient(const F& f, T* x, int N, rsMatrix<T>& g, const T* h)
  { gradient(f, x, N, g.getDataPointer(), h); }
  // does not allocate - so this can be used in production as well

  template<class F>
  static void hessian(const F& f, T* x, int N, rsMatrix<T>& H, const T* h)
  { hessian(f, x, N, H.getDataPointer(), h); }
  // does not allocate - so this can be used in production as well



  /** Even more convenient convenience function...but allocates */
  template<class F>
  static rsMatrix<T> hessian(const F& f, T* x, int N, const T* h)
  { rsMatrix<T> H(N, N); hessian(f, x, N, H, h); return H; }
  // allocates

  template<class F>
  static rsMatrix<T> hessian(const F& f, std::vector<T>& x, const T* h)
  { return hessian(f, &x[0], (int) x.size(), h); }
  // allocates

  template<class F>
  static void hessianTimesVector(const F& f, T* x, T* v, int N, T* Hv, const T* h, T k)
  { std::vector<T> wrk(2*N); hessianTimesVector(f, x, v, N, Hv, h, k, &wrk[0]); }
  // allocates

  // maybe make convenience functions that take a scalar h - they should create temporary array and
  // set all values in it to h

};


//=================================================================================================

/** Computes the numerical integral of a function defined by data points, i.e. the function:
\f[ F(x) = \int_c^x f(t) dt \f] where the lower integration limit c can be passed as a parameter 
into the function. Usage is similar to rsNumericDerivative. The parameter c can also be seen as an 
integration constant and determines yi[0] that shifts the overall resulting function up or down 
along the y-axis. The algorithm uses the trapezoidal rule, i.e. it sums up the areas under the 
trapezoids defined by a piecewise linear interpolant that passes through the datapoints. */
template<class Tx, class Ty>
void rsNumericIntegral(const Tx *x, const Ty *y, Ty *yi, int N, Ty c = Ty(0));
// move to class rsNumericIntegrator and rename to trapezoidal,


/** just a stub, at the moment */

template<class Tx, class Ty>
class rsNumericIntegrator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setNumberOfSamplePoints(int newNumSamples) { numSamples = newNumSamples; }


  //-----------------------------------------------------------------------------------------------
  // \name Integration

  /** Computes the definite integral of f in the interval from a to b. */
  //Ty integrate(const std::function<Ty(Tx)>& f, Tx a, Tx b);


protected:

  int numSamples = 10;

};


#endif
