//#include "rosic_NumericalTests.h"
using namespace rotes;

#include "rosic/rosic.h"
using namespace rosic;

void rotes::testUnivariateScalarFunction()
{
  /*
  double a[4] = {4, 3, 2, 1};
  //double *pa = &a[0];

  Polynomial p(3, a);

  double f_0_1 = p.getValue(             1.0);
  double f_1_1 = p.approximateDerivative(1.0, 1, 1.e-2);
  double f_2_1 = p.approximateDerivative(1.0, 2, 1.e-2);
  double f_3_1 = p.approximateDerivative(1.0, 3, 1.e-2);
  double f_4_1 = p.approximateDerivative(1.0, 4, 1.e-2);
  */


  double a[5] = {-52, -23, 21, -7, 1};
  Polynomial p(4, a);
  double f1 = p.approximateDerivativeAt(0.0, 1, 1.e-2);
    
  double xr;
  //xr = p.findRootViaNewtonNonRobust(0.0);       //  8 iterations
  //xr = p.findRootViaChebychevNonRobust(0.0);    // 7 iterations

  //xr = p.findRootViaNewtonNonRobust(0.5);       // 11 iterations
  //xr = p.findRootViaChebychevNonRobust(0.5);    // 16 iterations
  //xr = p.findRootViaRidders(0.0, 10.0);    //  3 iterations


  //xr = p.findRootViaNewtonNonRobust(3.8);       // 5 iterations
  //xr = p.findRootViaChebychevNonRobust(3.8);    // 4 iterations
  //xr = p.findRootViaRidders(3.8, 4.3);          // 3 iterations


  // we use the polynomial p(x) = x^3 - x with roots at -1, 0, 1 (up, down, up)
  double a2[4] = {0, -1, 0, 1}; 
  Polynomial p2(3, a2);

  //xr = p2.findRootViaRidders(-3.0, -0.5);    // 6 iterations
  //xr = p2.findRootViaRidders(+0.5, +3.0);    // 6 iterations
  //xr = p2.findRootViaRidders(-0.9, +0.7);    // 6 iterations
  //xr = p2.findRootViaRidders(-0.7, +0.9);    // 6 iterations
  //xr = p2.findRootViaRidders(+0.0, +0.1);    // 0 iterations
  //xr = p2.findRootViaRidders(-0.1, +0.0);    // 0 iterations
  xr = p2.findRootViaRidders(-3.0, +6.0);    // 10 iterations (there 3 roots  between -3...+6)
  //xr = p2.findRootViaRidders(-6.0, +3.0);    // 9 iterations
  //xr = p2.findRootViaRidders(-1.1, +1.2);    // 7 iterations


  // try f(x) = 1 / (x-a)   -> a is a singularity instead of a root

  // todo: return a bool or better, merge with similar tests for rapt - requires to implement
  // the Ridders algo in rapt, too. The numerical analysis algos should all be templatized and go
  // to rapt anyway

  int dummy = 0;
}

/*
void rotes::testUnivariateRootFinder()
{
 
  int dummy = 0;
}
*/
