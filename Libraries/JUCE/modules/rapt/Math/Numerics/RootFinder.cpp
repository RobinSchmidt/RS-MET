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

template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::falsePosition(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  static const int maxNumIterations = 60; // should be enough for double-precision
  Tx tol = std::numeric_limits<Tx>::epsilon();
  Ty fL  = f(xL) - y;
  Ty fR  = f(xR) - y;
  Tx xM; Ty fM;
  for(int i = 1; i <= maxNumIterations; i++) {
    xM = xL - (xR-xL) * Tx(fL/(fR-fL));
    fM = f(xM) - y;
    if(fM == 0 || xR-xL <= fabs(xM*tol)) 
      return xM; // done
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM; fR = fM; }
  }
  //rsError("rsRootFinder::falsePosition failed to converge");
  return xM;
}
