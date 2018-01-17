template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::bisection(std::function<Ty(Tx)>& f, Tx xL, Tx xR, Ty y)
{
  //Tx tol = Tx(0.0000001);      // tolerance, make parameter..and/or use machine epsilon
  Tx tol = Tx(0.5) * std::numeric_limits<Tx>::epsilon();
  Ty fL = f(xL) - y;  
  int its = 0;                // iteration counter - for development
  while(xR-xL > tol) {
    Tx xM  = Tx(0.5)*(xL+xR);
    Ty fM = f(xM) - y;
    if(fL*fM > 0) 
    { 
      xL = xM; 
      fL = fM; 
    }
    else          
    { 
      xR = xM; 
    }
    its++;
  }
  return Tx(0.5)*(xL+xR);
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
