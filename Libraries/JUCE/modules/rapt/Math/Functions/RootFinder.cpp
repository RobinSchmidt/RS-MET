template<class Tx, class Ty>
Tx rsRootFinder<Tx, Ty>::bisection(std::function<Ty(Tx)>& f, Tx a, Tx b, Ty y)
{
  Tx tol = Tx(0.000001);      // tolerance, make parameter..and/or use machine epsilon
  Ty fa = f(a) - y;  
  //Ty fb = f(b) - y;
  int its = 0;                // iteration counter - for development
  while(b-a > tol) {
    Tx mid  = Tx(0.5)*(a+b);
    Ty fMid = f(mid) - y;
    if(fa*fMid > 0) 
    { 
      a  = mid; 
      fa = fMid; 
    }
    else          
    { 
      b  = mid; 
      //fb = fMid; // not needed?
    }
    its++;
  }
  return Tx(0.5)*(a+b);

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
