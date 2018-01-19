template<class T>
T rsRootFinder<T>::bisection(std::function<T(T)>& f, T xL, T xR, T y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  T tol = std::numeric_limits<T>::epsilon();
  T fL  = f(xL) - y;
  T xM, fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    xM = T(0.5)*(xL+xR);
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) 
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM;          }
  }
  //rsError("rsRootFinder::bisection failed to converge");
  return xM;
}

template<class T>
T rsRootFinder<T>::falsePosition(std::function<T(T)>& f, T xL, T xR, T y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  T tol = std::numeric_limits<T>::epsilon();
  T fL  = f(xL) - y;
  T fR  = f(xR) - y;
  T xM, fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    //xM = xL - (xR-xL) * fL / (fR-fL);
    xM = rsLine2D<T>::zeroCrossing(xL, fL, xR, fR);
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) 
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM; fR = fM; }
  }
  //rsError("rsRootFinder::falsePosition failed to converge");
  return xM;
}
