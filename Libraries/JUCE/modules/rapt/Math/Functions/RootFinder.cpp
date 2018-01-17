template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::bisection(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  Tx eps = std::numeric_limits<Tx>::epsilon();
  Ty fL  = f(xL) - y;  
  Tx xM  = Tx(0.5)*(xL+xR);
  int its = 0;                // iteration counter - for development
  while(xR-xL > fabs(xM*eps)) {
    Ty fM = f(xM) - y;
    if(fL*fM > 0) { xL = xM; fL = fM; }
    else          { xR = xM;          }
    xM = Tx(0.5)*(xL+xR);
    its++; }
  return xM;
  // maybe have a maximum number of iterations after which to return and maybe trigger an assertion
  // when it is exceeded

  /*
  // too simple - works only for ascending functions:
  while( xR-xL > tol ) {
    Tx xM  = Tx(0.5) * (xL+xR);  // x in the middle betwen xL, xR
    Ty yM  = f(xM);              // function value at xM
    if( yM > y ) 
      xL = xM;
    else
      xR = xM; 
  }
  */
  //return Tx(0.5) * (xL+xR);
}
