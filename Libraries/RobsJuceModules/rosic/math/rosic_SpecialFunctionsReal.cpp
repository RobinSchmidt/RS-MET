//#include "rosic_SpecialFunctionsReal.h"
//using namespace rosic;

#include "../_third_party/cephes_tweaked/double/mconf.h"
#include "../_third_party/cephes_tweaked/double/mtherr.c"
#include "../_third_party/cephes_tweaked/double/const.c"
#include "../_third_party/cephes_tweaked/double/sincos.c"
#include "../_third_party/cephes_tweaked/double/ellpj.c"

double rosic::besselI0(double x)
{
  // use the power series expansion (see Shenoi - Introduction to Digital Signal Processing and Filter Design, eq. 5.42, page 268):

  double x2 = 0.5*x;
  double powers[RAPT::rsNumInverseFactorials];
  powers[0] = 1.0;
  powers[1] = x2;
  for(int k=2; k<RAPT::rsNumInverseFactorials; k++)
    powers[k] = x2 * powers[k-1];

  double tmp[RAPT::rsNumInverseFactorials];
  tmp[0] = 0.0;
  for(int k=1; k<RAPT::rsNumInverseFactorials; k++)
  {
    tmp[k]  = RAPT::rsInverseFactorials[k] * powers[k];
    tmp[k] *= tmp[k];
  }

  double result = 1.0 + RAPT::rsArrayTools::sum(tmp, RAPT::rsNumInverseFactorials);

  return result;
}

double rosic::cd(double u, double m)
{
  double cdResult, snDummy, cnDummy, dnDummy, phiDummy;
  ellipticFunctions(u, m, &snDummy, &cnDummy, &dnDummy, &phiDummy);
  cdResult = cnDummy/dnDummy;
  return cdResult;
}

double rosic::cn(double u, double m)
{
  double snDummy, cnResult, dnDummy, phiDummy;
  ellipticFunctions(u, m, &snDummy, &cnResult, &dnDummy, &phiDummy);
  return cnResult;
}

double rosic::dn(double u, double m)
{
  double snDummy, cnDummy, dnResult, phiDummy;
  ellipticFunctions(u, m, &snDummy, &cnDummy, &dnResult, &phiDummy);
  return dnResult;
}

int rosic::ellipticFunctions(double u, double m, double *sn, double *cn, double *dn, double* phi)
{
  return ellpj(u, m, sn, cn, dn, phi);
}

void rosic::ellipticIntegral(double k, double *K, double *Kprime, int M)
{
  // this implementation is a port from the MatLab file ellipk by Sophocles Orfanidis

  double kmin = 1e-6; 
  double kmax = sqrt(1-kmin*kmin);

  double kp, L;
  double* v  = new double[M];
  double* vp = new double[M];

  int n;

  if( k == 1.0 )
   *K = INF;
  else if( k > kmax )
  {
    kp = sqrt(1.0-k*k);
    L  = -log(kp/4.0);
    *K = L + (L-1.0)*kp*kp/4.0; 
  }
  else
  {
    landen(k, M, v);
    for(n=0; n<M; n++)
      v[n] += 1.0;
    *K = RAPT::rsArrayTools::product(v, M) * 0.5*PI;
  }

  if( k == 0.0 )
   *Kprime = INF;
  else if( k < kmin )
  {
    L       = -log(k/4.0);
    *Kprime = L + (L-1.0)*k*k/4.0;
  }
  else
  {
    kp = sqrt(1.0-k*k);                             
    landen(kp, M, vp);
    for(n=0; n<M; n++)
      vp[n] += 1.0;
    *Kprime = RAPT::rsArrayTools::product(vp, M) * 0.5*PI;
  }

  delete[] v;
  delete[] vp;
}

double rosic::sn(double u, double m)
{
  double snResult, cnDummy, dnDummy, phiDummy;
  ellipticFunctions(u, m, &snResult, &cnDummy, &dnDummy, &phiDummy);
  return snResult;
}

double rosic::solveLinLogEquation(double a, double b, double c)
{
  rassert( a*b >= 0.0 ); 
    // "a" and "b" should have the same sign or one of them should be zero - otherwise there might 
    // not be a unique solution and the behavior is undefined

  int i      = 0;       // iteration counter
  int maxIt  = 1000;    // maximum number of iterations
  double tol = 1.e-13;  // tolerance

  // initial guess (check the math - maybe a better guess can be found):
  double x = 1.0;       
  if( fabs(a) < fabs(b) )
    x = exp(-c/b);
  else
  {
    x = fabs(c/a);  // try -c/a
  }

  // refinement via iteration:
  double xOld = x + 2*tol;  
  double f, fp;
  while( fabs(x-xOld) > tol && i <= maxIt ) // maybe use relative error: x*fabs(d) > 1.e-13
  {
    xOld = x;
    if( xOld <= 0.0 )
      x = exp(-(c+a*xOld)/b);       // fixed-point iteration step
    else
    {
      f  = a*xOld + b*log(xOld) + c;  // f(x)
      fp = a + b/xOld;                // f'(x)
      x  = xOld - f/fp;               // Newton step
    }
    i++;
  }

  rassert( i < maxIt ); // iteration did not converge
  return x;  
}

// does not work for negative c, for small c > 0 it behaves erratic
double rosic::solveLinLogEquationOld(double c, double y)
{
  // initial guess (check the math of the fabs - whether it is really the right thing to do):
  // check the function with lots of weird inputs
  double x;
  if( y > c )
  //if( y > 1/c )
  //if( y > 4 )
    x = fabs(y/c);           // x(y) is asymptotically linear in this regime
  else
  {
    //x = fabs(exp((y-c)/c));  // x(y) is asymptotically exponential in this regime
    x = fabs(exp(y-c));    // x(y) is asymptotically exponential in this regime
  }

  x = 1; // test


  int i     = 0;    // iteration counter
  int maxIt = 1000; // maximum number of iterations
  double tol = 1.e-13;  // tolerance


  

  /*
  // refine x via fixpoint iteration (version 1):
  double xOld = x + 2*tol;
  while( fabs(x-xOld) > tol && i <= maxIt ) // maybe use realtive error: x*fabs(d) > 1.e-13
  {
    xOld = x;
    x    = fabs((y-log(x))/c);
    i    = i+1;
  }
  */
  
 
  

  /*
  // refine x via fixpoint iteration (version 2):
  x = 0;
  double b = exp(y);
  double xOld = x + 2*tol;
  while( fabs(x-xOld) > tol && i <= maxIt ) // maybe use realtive error: x*fabs(d) > 1.e-13
  {
    xOld = x;
    x    = b*exp(-c*x);
    i    = i+1;
  }
  */
  
  
  
  


  
  // refine x via Newton iteration:
  double b = exp(y);
  double d = x;    // difference to previous iteration
  double xOld ;
  while( fabs(d) > tol && i <= maxIt ) // maybe use realtive error: x*fabs(d) > 1.e-13
  {
    xOld = x;
    double f  = c*x + log(x) - y;  // f(x)
    double fp = c + 1/x;                // f'(x)
    d = f/fp;

    //x = fabs(x-d);

    x = x-d;
    if( x < 0.0 )
    {
      //x = 0.5*xOld;
      //x = fabs((y-log(xOld))/c);
      x = b*exp(-c*xOld);  // a fixed point iteration step (instead of newton iteration)
    }

    i = i+1;
  }
  
  
  
  
  

  //rassert( i < maxIt ); // iteration did not converge

  return x;
}
