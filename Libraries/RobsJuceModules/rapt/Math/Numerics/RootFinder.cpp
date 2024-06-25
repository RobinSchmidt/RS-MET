
template<class T>
T rsRootFinder<T>::findRoot(const std::function<T(T)>& f, T y, T x0)
{
  T xL, xR;
  findBracket(f, &xL, &xR, y, x0);
  return bisection(f, xL, xR, y);   // ToDo: use better algo - maybe implement Brent or something
}

// Helper function:
template<class T>
inline bool isConvergedToRoot(T xL, T xR, T yM, T xTol, T yTol)
{
  rsAssert(xR >= xL, "Upper limit below lower limit in rsRootFinder)");

  if(abs(yM) <= yTol)    
    return true;
  // wait...this formula works only if the target-value yT for y is 0 - in general, we need to do: 
  // if(abs(yM-yT) <= yTol*abs(yT)) to make the y-tolerance relative, too


  T dx = xR-xL;             // for debug
  T lx = xTol*abs(xL+xR);

  if( (xR-xL) <= xTol*abs(xL+xR)*0.5) // x-tolerance is relative
    return true;
  // but this relative tolerance seems to be bad when the root or pole is exactly at 0

  return false;
}


template<class T>
T rsRootFinder<T>::bisection(const std::function<T(T)>& f, T xL, T xR, T y)
{
  static const int maxNumIterations = 100; // 100 should be enough for double-precision

  T tol = std::numeric_limits<T>::epsilon();
  // This is not always suitable! We really should let the caller pass this in

  T fL  = f(xL) - y;
  T xM, fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    xM = T(0.5)*(xL+xR);
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol))    // old
    //if(isConvergedToRoot(xL, xR, fM, tol, tol))
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM;          }
  }
  rsError("rsRootFinder::bisection failed to converge");
  return xM;
}

template<class T>
T rsRootFinder<T>::falsePosition(const std::function<T(T)>& f, T xL, T xR, T y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  T tol = std::numeric_limits<T>::epsilon();
  T fL  = f(xL) - y;
  T fR  = f(xR) - y;
  T xM, fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    xM = rsLine2D<T>::zeroCrossing(xL, fL, xR, fR); // = xL - (xR-xL) * fL / (fR-fL)
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) // old
    //if(isConvergedToRoot(xL, xR, fM, tol, tol))
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM; fR = fM; }
  }
  rsError("rsRootFinder::falsePosition failed to converge");
  return xM;
}

template<class T>
void rsRootFinder<T>::findBracket(const std::function<T(T)>& f, T* xL, T* xR, T y, T x0)
{
  // Init:
  T yL, yR, dL, dR;          // Declare internal variables.
  *xL = *xR = x0;            // Init bracket to degenerate interval [xL, xR] = [x0, x0].
  yL  =  yR = f(x0);         // Function values at the bracket points.
  dL  =  dR = T(1);          // Expansion deltas for left and right. They grow over time.

  // Helper functions to expand the current interval:
  auto expandRight = [&]() { *xR += dR; dR *= T(2); yR = f(*xR); };
  auto expandLeft  = [&]() { *xL -= dL; dL *= T(2); yL = f(*xL); };

  // Expand the interval to the left or right until yL <= y <= yR or yL >= y >= yR:
  while(true)
  {
    if(yL < yR)              // Function f goes upward inside current interval
    { 
      if(yR < y)             // yR is too small (and yL is even smaller)
        expandRight();       //   ..so let's expand right where yR gets bigger
      else if(yL > y)        // yL is too big (and yR is even bigger)
        expandLeft();        //   ..so let's expand left where yL gets smaller
      else
        break;               // We have reached yL <= y <= yR
    }
    else                     // yL >= yR so f goes down or is constant inside current interval.
    {
      if(yR > y)             // yR is too big (and yL is at least as big)
        expandRight();       //   ..so let's expand right where yR gets smaller
      else if(yL < y)        // yL is too small (and yR is at least as small)
        expandLeft();        //   ..so let's expand left where yL gets bigger
      else
        break;               // We have reached yL >= y >= yR
    }
  }

  // ToDo: 
  //
  // - Verify, if the <= and >= are used correctly in all the descriptions and code or if it 
  //   should be < and > etc., i.e. verify, if the edge-case behavior is working and documented 
  //   correctly.
  //
  // - Maybe allow the user to also select the initial values for dL, dR (both 1.0) and their 
  //   growth factors (both 2.0)
}




/*



ToDo: 
-Implement the modified false position method describen in Hamming's "Numerical Methods 
 for Scientists and Enginners (2nd Ed)", pages 65-67 - should converge better for functions with
 flat regions see also: https://en.wikipedia.org/wiki/False_position_method#The_Illinois_algorithm
 "When the new y-value has the same sign as the previous one, meaning that the data point before 
  the previous one will be retained, the Illinois version halves the y-value of the retained data
  point."

template<class T>
T rsRootFinder<T>::modifiedFalsePosition(std::function<T(T)>& f, T xL, T xR, T y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  T tol = std::numeric_limits<T>::epsilon();
  T fL  = f(xL) - y;
  T fR  = f(xR) - y;
  T xM, fM;
  for(int i = 1; i <= maxNumIterations; i++) {

    xM = rsLine2D<T>::zeroCrossing(xL, fL, xR, fR); // = xL - (xR-xL) * fL / (fR-fL)
    // not yet modified


    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) 
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM; fR = fM; }
  }
  rsError("rsRootFinder::modifiedFalsePosition failed to converge");
  return xM;
}


ToDo: 
-Maybe factor out common stuff...maybe calling a function getNewEstimate(xL, fL, xM, fM, xR, fR) 
 and call it as xM = getNewEstimate(xL, fL, ...) - such a structure can also accomodate for Brent
 algorithm later (fits (inverse) parabola through 3 points)
-Let the user pass the maximum number of iterations and (relative) tolerance. I think, we should 
 always work with relative rather than absolute tolerance because it just makes more sense.
-Implement Newton iteration with and without Schroeder's modification. The regular Newton 
 iteration uses the update rule:
   x[n+1] = x[n] + d[n]   where   d[n] = - f(x[n]) / f'(x[n])
 This rule has slow convergence when the root has a multiplicity m > 1. Schroeder's modification 
 fixes this problem by using the update rule:
   x[n+1] = x[n] + m*d[n]
 with d[n] defined as before. We could let the user pass the multiplicity m for cases where it is 
 known. If m is unknown, we could let the algorithm have a current estimate m[n] (initialized to 1) 
 and in each iteration, try 3 tentative update steps using m-1, m, m+1 and then actually use the
 step that brought us closest to zero and use as m[n+1] the value that gave rise to this best 
 result. We could also return our final estimate of m to the caller as estimate for the 
 multiplicity of the root. Maybe when m[n] = 1, don't try m-1. It would be useless because then
 we would do x[n+1] = x[n] + 0*d[n]. But maybe just do it anyway to keep the logic simple. It 
 isn't expensive. See Numerik (Meister, Sonar), page 249. When a too high m is used, it tends to
 create oscillations.
 Q: What other conditions besides root multiplicity can slow down convergence? I think, having a 
 flat region around the current root estimate (but not necessarily around the root itself) can 
 shoot the next estimate far away. IIRC, Numerical Recipies has some discussion about this (not 
 sure though). Can this issue be fixed by some modification, too?
-Also in the "Numerik" book by Meister/Sonar on page 278, there's the task to show that the 
 Schröder modification of Newton iteration applied to f(x) = (x+2)^2 * (x-1) * (x-7)^5 with a root
 of multiplicity 5 at x=7 yields an unstable iteration function: 
   x[n+1] = x[n] - 5 * f(x[n]) / f'(x[n])
 for any start value x[0] in the interval [6,8]. Try that experimentally. A general necessary 
 condition for any iteration to converge is that the derivative of the iteration function F (Phi in 
 the book) has a derivative with absolute value less than 1. I think, a sufficient condition is that
 |F(x) - F(y)| < |x-y| for all x,y in a neighborhood around the root. This is called the 
 contraction property (explained on page 274 and 257). The iteration function F is the function in
 x[n+1] = F(x[n]). In the case of Newton iteration, we have F(x[n]) = x[n] - f(x[n]) / f'(x[n]) but 
 the stability criterion applies to other kinds of fixed-point iterations as well. The book also 
 gives conditions for the edge case when |F'| = 1 in terms of 2nd and 3rd derivatives on pages 
 278/279. Maybe try these experimentally, too.
-Maybe implement an experiment that numerically computes the Feigenbaum constant (pg 279/280).
-On page 284, it is stated that the iteration x[n+1] = x[n] - (x[n]^p - q) / (p * x[x]^(p-1)) 
 converges to the the p-th root of q. Try than in an experiment. Is that Newton iteration applied 
 to the equation: x = q^(1/p)  ->  x^p = q  ->  x^p - q = 0? What would be a good starting value
 for the p-th root? Maybe right-shift the bits of the input by p positions or in case of floats, 
 divide the exponent by p? In decimal, we have sqrt(100) = 10, sqrt(10000) = 100, 
 sqrt(1000000) = 1000, etc. and the approximation would be very coarse but it woul be better in 
 binary, I think. Maybe implement such an algo for sqrt, cbrt for rsBigFloat.
-Implement the secant method. 
-Implement 2nd order methods without derivatives based on (inverse) parabolic interpolation. See:
 https://en.wikipedia.org/wiki/Brent%27s_method
 https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
 https://en.wikipedia.org/wiki/Muller%27s_method
 https://en.wikipedia.org/wiki/Ridders%27_method
 I think, inverse quadratic interpolation methods should find the root in one step when the 
 function f(x) is of the form: f(x) = a + b*sqrt(x+c). Here, c shifts the parabola left/right, a
 shifts it up/down and b determines the shape. Try this for example with these a,b,c:
   a = -6, b = 3, c = -2: has root at x = 6   https://www.desmos.com/calculator/n1wwd8cx3u
 Regular quadratic interpolation should find the root in one step for functions of the form:
 f(x) = a + b*(x+c)^2. Solving for x such that f(x) = 0 will require taking a square-root per 
 iteration. That could create a nonzero imaginary part. That may be a good or bad thing, depending 
 on whether a complex root may be expected or not. The Brent method should perhaps be the go-to 
 method for root-finding without derivatives
-Implement ternary search ("trisection"?) for finding minima or maxima. See:
 https://cp-algorithms.com/num_methods/ternary_search.html


*/
