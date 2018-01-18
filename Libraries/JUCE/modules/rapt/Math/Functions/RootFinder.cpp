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
      //rsError("Bisection failed to converge");
      break; }
  }
  return xM;
}

// move to somewhere else:
template<class T>
inline T lineZeroCrossing(T xL, T yL, T xR, T yR)
{
  T a = (yR-yL)/(xR-xL); // slope
  T b = yL-a*xL;         // offset
  return -b/a;

  // maybe this can be simplified algebraically? ..like this:
  //return xL - yL*(xR-xL)/(yR-yL);
}

template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::falsePosition(std::function<Ty(Tx)>& func, Tx xLeft, Tx xRight, Ty y)
{

  return xLeft; // preliminary
}
