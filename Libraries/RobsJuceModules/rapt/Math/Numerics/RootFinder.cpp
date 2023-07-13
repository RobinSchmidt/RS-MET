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
-Maybe factor out common stuf...maybe calling a function getNewEstimate(xL, fL, xM, fM, xR, fR) 
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
 isn't expensive. See Numerik (Meister, Sonar), page 249.
 Q: What other conditions besides root multiplicity can slow down convergence? I think, having a 
 flat region around the current root estimate (but not necessarily around the root itself) can 
 shoot the next estimate far away. IIRC, Numerical Recipies has some discussion about this (not 
 sure though). Can this issue be fixed by some modification, too?
   

*/
