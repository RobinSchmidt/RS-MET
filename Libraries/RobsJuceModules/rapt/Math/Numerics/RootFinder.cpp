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


// todo: implement the modified false position method describen in Hamming's "Numerical Methods 
// for Scientists and Enginners (2nd Ed)", pages 65-67 - should converge better for functions with
// flat regions
// see also: https://en.wikipedia.org/wiki/False_position_method#The_Illinois_algorithm
// "When the new y-value has the same sign as the previous one, meaning that the data point before the 
// previous one will be retained, the Illinois version halves the y-value of the retained data point."
/*
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
*/

// todo: factor out common stuf...maybe calling a function getNewEstimate(xL, fL, xM, fM, xR, fR) 
// and call it as xM = getNewEstimate(xL, fL, ...) - such a structure can also accomodate for brent
// algorithm later (fits (inverse) parabola through 3 points)

