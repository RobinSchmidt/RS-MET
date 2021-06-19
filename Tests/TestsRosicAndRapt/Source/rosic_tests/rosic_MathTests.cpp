//#include "rosic_MathTests.h"
using namespace rotes;

//#include "../Shared/Plotting/rosic_Plotter.h"
using namespace rosic;

bool rotes::testComplexSqrt()
{
  bool result = true;
  int N = 500;  // number of tests

  RAPT::rsNoiseGenerator<double> ng;
  ng.setRange(-1000, +1000);
  double re, im;


  for(int n = 0; n < N; n++)
  {
    re = ng.getSample();
    im = ng.getSample();
    std::complex<double> z1 = std::complex<double>(re, im);
    rosic::Complex z2 = rosic::Complex(re, im);
    z1 = sqrt(z1);
    z2 = sqrtC(z2);

    // ...compare z1, z2

    z1 = std::complex<double>(0., im);
    z2 = rosic::Complex(0., im);
    z1 = sqrt(z1);
    z2 = sqrtC(z2);

    // ...compare z1, z2



  }

  // they seem to be the same..so what is wrong in the PoleZeroMapper?


  return result;
}

bool rotes::testCubicCoeffsTwoPointsAndDerivatives()
{
  bool result = true;

  double  x[2] = {-3, 2};
  double  y[2] = {-2, 5};
  double dy[2] = { 7, 3};
  double  a[4];

  RAPT::rsPolynomial<double>::cubicCoeffsTwoPointsAndDerivatives(a, x, y, dy);

  // check results:
  double yc, dyc;       // computed values
  double tol = 1.e-14;  // tolerance

  RAPT::rsPolynomial<double>::evaluateWithDerivative(x[0], a, 3, &yc, &dyc);
  result &= RAPT::rsIsCloseTo( yc,  y[0], tol);
  result &= RAPT::rsIsCloseTo(dyc, dy[0], tol);

  RAPT::rsPolynomial<double>::evaluateWithDerivative(x[1], a, 3, &yc, &dyc);
  result &= RAPT::rsIsCloseTo( yc,  y[1], tol);
  result &= RAPT::rsIsCloseTo(dyc, dy[1], tol);

  return result;
}

bool rotes::testPolynomialDiffAndInt()
{
  bool result = true;

  double a[6]  = {2, -1, 5, 7, -3, 2};
  double ad[5];
  double ai[7];

  RAPT::rsPolynomial<double>::derivative(a, ad, 5);
  result &= (ad[0] == -1);
  result &= (ad[1] == 10);
  result &= (ad[2] == 21);
  result &= (ad[3] == -12);
  result &= (ad[4] == 10);

  RAPT::rsPolynomial<double>::integral(a, ai, 5, 2.0);
  result &= (ai[0] ==  2.0);
  result &= (ai[1] ==  2.0/1.0);
  result &= (ai[2] == -1.0/2.0);
  result &= (ai[3] ==  5.0/3.0);
  result &= (ai[4] ==  7.0/4.0);
  result &= (ai[5] == -3.0/5.0);
  result &= (ai[6] ==  2.0/6.0);

  return result;
}

bool rotes::testPolynomialComposition()
{
  bool result = true;

  static const int na = 5;
  static const int nb = 4;
  static const int nc = na*nb;
  double a[na+1] = {2, -1, 5, 7, -3, 2}; // 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double b[nb+1] = {3,  1, 4, -5, 3};    // 3*x^4 - 5*x^3 + 4*x^2 + 1*x^1 - 3*x^0
  double c[nc+1];
  RAPT::rsPolynomial<double>::compose(a, na, b, nb, c);

  // check, if the composed c-polynomial returns the same result as applying the 2nd b-polynomial
  // to the result of the 1st a-polynomial:
  double x = -3.0; // input value
  double y1, y2;
  y1 = RAPT::rsPolynomial<double>::evaluate(x,  a, na);
  y1 = RAPT::rsPolynomial<double>::evaluate(y1, b, nb);
  y2 = RAPT::rsPolynomial<double>::evaluate(x,  c, nc);
  result &= (y1 == y2);

  return result;
}

bool rotes::testPolynomialWeightedSum()
{
  bool result = true;

  static const int pN = 5;
  static const int qN = 3;
  static const int rN = pN;              // == max(pN, qN);
  double p[pN+1] = {3, -1, 5, 7, -3, 2}; // p(x) = 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 3*x^0
  double q[qN+1] = {2, -3, 6, 2};        // q(x) =                 2*x^3 + 6*x^2 - 3*x^1 + 2*x^0
  double r[rN+1];

  // r(x) = 2*p(x) + 3*q(x) = 4*x^5 - 6*x^4 + 20*x^3 + 28*x^2 - 11*x^1 + 12*x^0
  RAPT::rsPolynomial<double>::weightedSum(p, pN, 2.0, q, qN, 3.0, r);
  result &= (r[0] ==  12);
  result &= (r[1] == -11);
  result &= (r[2] ==  28);
  result &= (r[3] ==  20);
  result &= (r[4] == - 6);
  result &= (r[5] ==   4);

  // exchange roles of function parameters (function takes the other branch, result should be the
  // same):
  RAPT::rsPolynomial<double>::weightedSum(q, qN, 3.0, p, pN, 2.0, r);
  result &= (r[0] ==  12);
  result &= (r[1] == -11);
  result &= (r[2] ==  28);
  result &= (r[3] ==  20);
  result &= (r[4] == - 6);
  result &= (r[5] ==   4);

  // use a truncated polynomial for p (such that p and q are of the same order):
  RAPT::rsArrayTools::fillWithZeros(r, rN+1);
  RAPT::rsPolynomial<double>::weightedSum(p, qN, 2.0, q, qN, 3.0, r);
  result &= (r[0] ==  12);
  result &= (r[1] == -11);
  result &= (r[2] ==  28);
  result &= (r[3] ==  20);
  result &= (r[4] ==   0);
  result &= (r[5] ==   0);

  return result;
}

bool rotes::testPolynomialIntegrationWithPolynomialLimits()
{
  bool result = true;

  static const int np = 5;
  static const int na = 3;
  static const int nb = 4;
  static const int nP = np+1;
  static const int nq = nP*nb;  // == nP*rmax(na, nb);

  double p[np+1] = {2, -1, 5, 7, -3, 2}; // 2*x^5 - 3*x^4 + 7*x^3 + 5*x^2 - 1*x^1 + 2*x^0
  double a[na+1] = {2, -3, 6, 2};        // 2*x^3 + 6*x^2 - 3*x^1 + 2*x^0
  double b[nb+1] = {3,  1, 4, -5, 3};    // 3*x^4 - 5*x^3 + 4*x^2 + 1*x^1 - 3*x^0
  double q[nq+1];

  double x = 1.5; // input value
  double y1, y2;

  // obtain coefficients of indefinite integral:
  double P[nP+1];
  RAPT::rsPolynomial<double>::integral(p, P, np);

  // compute integration limits for definite integral:
  double lowerLimit = RAPT::rsPolynomial<double>::evaluate(x, a, na);
  double upperLimit = RAPT::rsPolynomial<double>::evaluate(x, b, nb);

  // evaluate definite integral:
  y1 = RAPT::rsPolynomial<double>::evaluate(upperLimit, P, nP)
    - RAPT::rsPolynomial<double>::evaluate(lowerLimit, P, nP);


  RAPT::rsPolynomial<double>::integrateWithPolynomialLimits(p, np, a, na, b, nb, q);
  y2 = RAPT::rsPolynomial<double>::evaluate(x, q, nq);

  result &= RAPT::rsIsCloseTo(y2, y1, 1.e-13 * fabs(y1));

  return result;
}

bool rotes::testPolynomialRootFinder()
{
  // we use the polynomial p(x) = x^4 - 7x^3 + 21*x^2 - 23*x - 52 with roots at 2+3i, 2-3i, -1, 4 as test function
  double a1[5] = {-52, -23, 21, -7, 1};
  std::complex<double> r1[4];
  RAPT::rsPolynomial<double>::roots(a1, 4, r1);

  static const int maxN     = 20;
  static const int numTests = 1000;
  double range = 10.0;   // range for the real and imaginary parts of the roots
  //double tol   = 5.e-9; // tolerance
  double tol   = 5.e-8; // tolerance
  std::complex<double> a[maxN+1];     // polynomial coefficients
  std::complex<double> rTrue[maxN];   // true roots
  std::complex<double> rFound[maxN];  // roots that were found
  bool result = true;    // if this is still true at the end of the function, the test has passed
  RAPT::rsRandomUniform(-range, range, 0);  // set seed
  int i, j, k;
  for(i = 1; i <= numTests; i++)
  {
    // polynomial order for this test:
    int N = (int) RAPT::rsRandomUniform(1.0, maxN);

    // generate a bunch of random roots:
    for(k = 0; k < N; k++)
    {
      rTrue[k].real(random(-range, range));
      rTrue[k].imag(random(-range, range));
    }

    // obtain polynomial coeffs:
    RAPT::rsPolynomial<double>::rootsToCoeffs(rTrue, a, N);

    // find the roots:
    RAPT::rsPolynomial<double>::roots(a, N, rFound);

    // try to find a matching root in the found roots for each of the true roots:
    for(j = 0; j < N; j++)
    {
      bool matchFound = false;
      for(k = 0; k < N; k++)
      {
        if( abs(rFound[j]-rTrue[k]) < tol )
        {
          matchFound = true;
          break;
        }
      }
      rassert( matchFound == true );
      result &= matchFound;
    }
    rassert( result == true );
  }

  rassert( result == true );
  int dummy = 0;
  return result;
}



void rotes::testLinLogEquationSolver()
{
  static const int N = 100;
  double xMin = 0.01;
  double xMax = 100.0;
  double yMin = -50.0;
  double yMax = +50.0;

  double a = -10.0;
  double b = -2.0;
  int n;


  //double x[N], xLin[N], xLog[N], y[N], yLin[N], yLog[N], f[N];
  double x[N], y[N], yLin[N], yLog[N];


  // y as function of x:
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  for(n = 0; n < N; n++)
  {
    yLin[n] = a*x[n] + b;         // linear part
    yLog[n] = log(x[n]);          // logarithmic part
    y[n]    = yLin[n] + yLog[n];  // the sum of both parts (the objective function)
  }
  //Plotter::plotData(N, x, y, yLin, yLog);
  //Plotter::plotData(N, x, y);




  // x as function of y:
  RAPT::rsArrayTools::fillWithRangeLinear(y, N, yMin, yMax);
  for(n = 0; n < N; n++)
  {
    x[n] = rosic::solveLinLogEquation(a, b, -y[n]);
  }
  plotData(N, y, x);
}

void rotes::testLinLogEquationSolverOld()
{
  /*
  // automatic test:
  double x = 1.0;
  double c = 1.0;
  double y = c*x + log(x);
  double xc;  // computed value for x (from c and y)

  xc =  solveLinLogEquation(c, y);
  // check areAlmostEqual(x xc)

  // \todo: use a double-loop that loops through values of y and c in some range and check results
  */



  // interactive test:


  static const int N = 1000;
  double c = +0.001;
  double k = -0.5;
  double xMin = 0.01;
  double xMax = 10.0;
  int n;
  //double x[N], xLin[N], xLog[N], y[N], yLin[N], yLog[N], f[N];
  double x[N], y[N], yLin[N], yLog[N];

  // y as function of x:
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  for(n = 0; n < N; n++)
  {
    yLin[n] = c*x[n];          // linear part
    yLog[n] = log(x[n]);       // logarithmic part
    y[n]    = yLin[n] + yLog[n];  // the sum of both parts
  }
  plotData(N, x, y, yLin, yLog);

  /*
  // x as function of y:
  rosic::fillWithRangeLinear(y, N, -50, 50);
  for(n = 0; n < N; n++)
  {
    xLin[n] = y[n] / c;
    //xLog[n] = exp((y[n]-c)/c);
    xLog[n] = exp((y[n]-c));
    x[n]    = solveLinLogEquation(c, y[n]);
  }
  //Plotter::plotData(N, y, xLin, xLog);
  //Plotter::plotData(N, y, x, xLin, xLog);
  //Plotter::plotData(N, y, x, xLin);
  //Plotter::plotData(N, y, x, xLog);
  //Plotter::plotData(N, y, xLog);
  //Plotter::plotData(N, y, xLin);
  Plotter::plotData(N, y, x);
  */


  /*
  // g(x) = c*x + ln(x) + k
  rosic::fillWithRangeLinear(x, N, xMin, xMax);
  for(n = 0; n < N; n++)
  {
    yLin[n] = c*x[n] + k;              // linear part
    yLog[n] = log(x[n]);           // logarithmic part
    y[n] = yLin[n] + yLog[n] + k;  // the sum of both parts and the constant
  }
  Plotter::plotData(N, x, y, yLin, yLog);
  */

  int dummy = 0;
}

bool rotes::testLinearSystemSolver()
{
  bool ok = true;

  double x[3];
  double y[3]    = {-4, 15, 11};
  double A[3][3] = {{ 1, 2, -3},
                    {-5, 7,  2},
                    {-2, 2,  3}};
  double *pA[3];
  for(int i = 0; i < 3; i++)
    pA[i] = &A[i][0];

  solveLinearSystem(pA, x, y, 3);  // this is actually obsolete -> deprecate it

  double tol = 1.e-15;
  ok &= rsIsCloseTo(x[0], 1.0, tol);
  ok &= rsIsCloseTo(x[1], 2.0, tol);
  ok &= rsIsCloseTo(x[2], 3.0, tol);

  //printf("%s %.2f %.2f %.2f", "Result:", x[0], x[1], x[2]);
  //getchar();

  return ok;
}
