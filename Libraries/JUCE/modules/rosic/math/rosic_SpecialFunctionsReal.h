#ifndef rosic_SpecialFunctionsReal_h
#define rosic_SpecialFunctionsReal_h

//// rosic includes:
//
//#include "rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /** This file contains functions for real arguments that that can not be expressed as combinations of elementary functions, such as the 
  gamma-functions, Bessel-functions, Jacobian elliptic functions, etc. */

  // Bessel function J_n(x) */
  //INLINE double besselj(double x, double order);

  /** Returns the value of the modified Bessel-function I0. 
  ...quick and dirty implementation via power series - works well for small x. */
  double besselI0(double x);


  /** Returns the value of the Jacobian elliptic function cd(u,m) with argunment u and parameter m (m=k^2 with k being elliptic 
  modulus (?)) . */
  double cd(double u, double m);

  /** Returns the value of the Jacobian elliptic function cn(u,m) with argunment u and parameter m  m=k^2 with k being elliptic 
  modulus (?) ) . */
  double cn(double u, double m);

  /** Returns the value of the Jacobian elliptic function dn(u,m) with argunment u and parameter m (m=k^2 with k being elliptic 
  modulus (?) ) . */
  double dn(double u, double m);

  /** Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m), and dn(u|m) of parameter m between 0 and 1, and real argument u. 
  phi is the amplitdue of u. */
  int ellipticFunctions(double u, double m, double *sn, double *cn, double *dn, double* phi);

  /** Evaluates the complete elliptic integral of the first kind with elliptic modulus k. K will be assigned to the quarter period K(k) 
  and Kprime will be assigned to the quarter period K'(k) = K(k'), k' = sqrt(1-k^2). M is the number of Landen iterations which 
  determines the precision of the calculation. */
  void ellipticIntegral(double k, double *K, double *Kprime, int M=7);

  // gamma function
  //double gamma(double x);

  /** Returns the value of the Jacobian elliptic function sn(u,m) with argunment u and parameter m ( m=k^2 with k being elliptic 
  modulus (?) ) . */
  double sn(double u, double m);


  /** Finds a value x, such that f(x) = a*x + b*log(x) + c = 0. 
  \todo comment for which input combinations the function will work. 
  */
  double solveLinLogEquation(double a, double b, double c);

  /** Finds the solution x to the equation c*x + log(|x|) = y where c and y are given. 
  \todo test function for various input ranges - add comments for which inputs the function will 
  work.  */
  double solveLinLogEquationOld(double c, double y);




}

#endif 
