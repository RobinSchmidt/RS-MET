template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::bisection(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  T tol = T(0.000001);         // tolerance, make parameter..and/or use machine epsilon
  while( xR-xL > tol ) {
    T xM  = T(0.5) * (xL+xR);  // x in the middle betwen xL, xR
    T yM  = f(xM);             // function value at xM
    if( yM > y ) 
      xL = xM;
    else
      xR = xM; 
  }
  return T(0.5) * (xL+xR);
}
