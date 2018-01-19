/*
template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::bisection(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  int its = 0;                            // iteration counter
  Tx eps  = std::numeric_limits<Tx>::epsilon();
  Ty fL   = f(xL) - y;  
  Tx xM   = Tx(0.5)*(xL+xR);
  while(xR-xL > fabs(xM*eps)) {
    Ty fM = f(xM) - y;
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM;          }
    xM = Tx(0.5)*(xL+xR);
    if(its++ > maxNumIterations) {
      //rsError("rsRootFinder::bisection failed to converge");
      break; }
  }
  return xM;
}
*/

// old code above has been replaced with code that is structurally the same as falsePosition
// ...maybe something can be factored out
template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::bisection(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  Tx tol = std::numeric_limits<Tx>::epsilon();
  Ty fL  = f(xL) - y;
  Tx xM; Ty fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    xM = Tx(0.5)*(xL+xR);
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) 
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM;          }
  }
  //rsError("rsRootFinder::bisection failed to converge");
  return xM;
}


// move to somewhere else (maybe rsLine):
template<class T>
inline T lineZeroCrossing(T xL, T yL, T xR, T yR)
{
  return xL - yL*(xR-xL)/(yR-yL);

  //// old - code above is algebraically simplified:
  //T a = (yR-yL)/(xR-xL); // slope
  //T b = yL-a*xL;         // offset
  //return -b/a;
}

template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::falsePosition(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  Tx tol = std::numeric_limits<Tx>::epsilon();
  Ty fL  = f(xL) - y;
  Ty fR  = f(xR) - y;
  Tx xM; Ty fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    xM = lineZeroCrossing(xL, fL, xR, fR);
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) 
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM; fR = fM; }
  }
  //rsError("rsRootFinder::falsePosition failed to converge");
  return xM;
}
