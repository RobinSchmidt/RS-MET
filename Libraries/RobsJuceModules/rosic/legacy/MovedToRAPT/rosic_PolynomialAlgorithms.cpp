/*
Complex rosic::evaluatePolynomialWithRoots(Complex s, Complex *r, int N)
{
  Complex result = 1.0;
  for(int i = 0; i < N; i++)
  {
    if( !r[i].isInfinite() )
      result *= (s - r[i]);
  }
  return result;
}
*/

// used in convergeToRootViaLaguerre:
double evaluatePolynomialWithTwoDerivativesAndError(Complex *a, int order, Complex z, Complex *P)
{
  P[0] = a[order];               // P(z)
  P[1] = Complex(0.0, 0.0);      // P'(z)
  P[2] = Complex(0.0, 0.0);      // P''(z)
  double err = P[0].getRadius(); // estimated roundoff error in evaluation of the polynomial
  double zA  = z.getRadius();    // absolute value of z
  for(int j = order-1; j >= 0; j--) 
  {
    P[2] = z * P[2] + P[1];
    P[1] = z * P[1] + P[0];
    P[0] = z * P[0] + a[j];
    err  = P[0].getRadius() + zA*err;
  }
  P[2] *= 2.0;
  return err;
}
/*
Complex rosic::convergeToRootViaLaguerre(Complex *a, int order, Complex initialGuess)
{
  const double eps = std::numeric_limits<double>::epsilon(); 

  static const int numFractions = 8; // number of fractions minus 1 (for breaking limit cycles)

  static const int itsBeforeFracStep = 10;  // number of iterations after which a fractional step 
                                            // is taken (to break limit cycles)

  static const int maxNumIterations = itsBeforeFracStep*numFractions;  

  // fractions for taking fractional update steps to break a limit cycles:
  static double fractions[numFractions+1] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0}; 
       
  Complex r = initialGuess; // the current estimate for the root
  for(int i = 1; i <= maxNumIterations; i++) 
  {
    Complex P[3];    // holds P, P', P'' 
    double  err = eps * evaluatePolynomialWithTwoDerivativesAndError(a, order, r, P);

    if( P[0].getRadius() <= err ) 
      return r;   
      // is this the "simplified stopping criterion due to Adams", referred to on page 373?
      // can we get rid of this? if so, we might also replace the above loop by 
      // evaluatePolynomialAndDerivativesAt

    // Laguerre's formulas:
    Complex G  = P[1]/P[0];                      // Eq. 9.5.6
    Complex H  = G*G - P[2]/P[0];                // Eq. 9.5.7
    Complex sq = sqrtC((order-1)*(order*H-G*G)); // the square-root Eq. 9.5.11
    Complex Gp = G + sq;                         // denominator in 9.5.11 with positive sign
    Complex Gm = G - sq;                         // denominator in 9.5.11 with negative sign

    // choose Gp or Gm according to which has larger magnitude (page 372, bottom), re-use Gp for 
    // the result:
    double GpA = Gp.getRadius();
    double GmA = Gm.getRadius();
    if( GpA < GmA ) 
    {
      Gp  = Gm;
      GpA = GmA;
    }

    // compute difference between old and new estimate for the root r (the 'a' variable in 
    // Eq. 9.5.8)
    Complex dr;                       
    if( GpA > 0.0 ) 
      dr = Complex(order, 0.0) / Gp;  // Eq. 9.5.11
    else
      dr = exp(log(1.0+r.getRadius())) * Complex(cos((double)i), sin((double)i)); // use sinCos()

    // compute new estimate for the root:
    Complex rNew = r - dr;   
    if( r == rNew )
      return r;  // converged

    // update our r-variable to the new estimate:
    if( i % itsBeforeFracStep != 0 ) 
      r = rNew;
    else 
      r = r - fractions[i/itsBeforeFracStep]*dr; // fractional step to break limit cycle
  }

  rassert(false);  // error - too many iterations taken, algorithm did not converge
  return 0.0;
}

void rosic::findPolynomialRoots(Complex *a, int order, Complex *roots)
{
  const double eps = 2.0e-14; // for float, it was 2.0e-6 - use template numeric_limit<T>

  // allocate memory for the coefficients of the deflated polynomial and initialize it as 
  // non-deflated polynomial:
  Complex *ad = new Complex[order+1];
  RAPT::rsArrayTools::copy(a, ad, order+1);

  // loop over the roots:
  for(int j = order; j >= 1; j--) 
  {
    // find a root of the deflated polynomial using the Laguerre-algorithm with 0 as initial guess:
    Complex r = convergeToRootViaLaguerre(ad, j, Complex(0.0, 0.0));    

    // polish the root by using the Laguerre method with the undeflated polynomial and the 
    // non-polished root as initial guess: 
    r = convergeToRootViaLaguerre(a, order, r);

    // maybe move into a member function Complex::zeroNegligibleImaginaryPart(double ratio); 
    // -> ratio = 2*eps:
    if( fabs(r.im) <= 2.0*eps*fabs(r.re) ) 
      r.im = 0.0;

    // store root in output array:
    roots[j-1] = r; 

    // deflate the deflated polynomial again by the monomial that corresponds to our most recently 
    // found root:
    Complex rem = ad[j];  // remainder - not used, needed for passing a dummy pointer
    //dividePolynomialByMonomialInPlace(ad, j, r, &rem);
    RAPT::rsPolynomial<Complex>::dividePolynomialByMonomialInPlace(ad, j, r, &rem);
  }

  rosic::sortComplexArrayByReIm(roots, order);
  delete[] ad;
}

void rosic::findPolynomialRoots(double *a, int order, Complex *roots)
{
  Complex *ac = new Complex[order+1];
  RAPT::rsArrayTools::convertBuffer(a, ac, order+1);
  findPolynomialRoots(ac, order, roots);
  delete[] ac;
}

rsArrayTools<Complex> rosic::getPolynomialCoefficientsFromRoots(rsArrayTools<Complex> roots)
{
  rsArrayTools<Complex> coeffs;

  coeffs.ensureAllocatedSize(roots.getNumElements()+1);
  coeffs.appendElement(1.0);

  if( roots.getNumElements() < 1 ) 
    return coeffs;

  for(int i=0; i<roots.getNumElements(); i++)
  {
    Complex z = roots[i];
    coeffs.appendElement(coeffs[i]);
    for(int j=i; j>=1; j--)
      coeffs[j] = coeffs[j-1] - z * coeffs[j];
    coeffs[0] = -z * coeffs[0];
  }

  return coeffs;
}

void rosic::rootsToCoeffs(Complex *r, Complex *a, int N)
{
  // use only the finite roots:
  Complex *rF = new Complex[N]; 
  int nF = copyFiniteValues(r, rF, N);

  RAPT::rsArrayTools::fillWithZeros(a, N+1);
  if( nF == 0 )
    a[0] = 1.0;
  else
  {
    a[0] = -rF[0];
    a[1] = 1.0;
    for(int M = 2; M <= nF; M++)
    {
      a[M] = a[M-1];
      Complex rM = rF[M-1];
      for(int n = M-1; n >= 1; n--)
        a[n] = a[n-1] - rM*a[n];
      a[0] = -rM*a[0];
    }
  }

  delete[] rF;
}

void rosic::rootsToCoeffs(Complex *r, double *a, int N)
{
  Complex *ac = new Complex[N+1];
  rootsToCoeffs(r, ac, N);
  for(int n = 0; n <= N; n++)
    a[n] = ac[n].re;
  delete[] ac;
}

double rosic::getRootOfLinearEquation(double a, double b)
{
  if( a == 0.0 )
  {
    DEBUG_BREAK;
    return 0.0;
  }
  else
    return -b/a;
}

rsArrayTools<Complex> rosic::getRootsOfQuadraticEquation(double a, double b, double c)
{
  // catch degenerate cases where the leading coefficient is zero:
  if( a == 0.0 )
  {
    rsArrayTools<Complex> roots(1);
    roots[0] = getRootOfLinearEquation(b, c);
    return roots;
  }

  rsArrayTools<Complex> roots(2);

  double D      = b*b - 4.0*a*c; // the discriminant
  double factor = 1.0 / (2.0*a); // a common factor that appears everywhere
  if( D > 0.0 )
  {
    // D > 0: two distinct real roots:
    double sqrt_D = sqrt(D);
    roots[0]      = factor * (-b+sqrt_D);
    roots[1]      = factor * (-b-sqrt_D);
  }
  else if( D == 0.0 )
  {
    // D == 0: a real root with multiplicity 2:
    roots[1] = roots[0] = Complex( -b * factor );
  }
  else
  {
    // D < 0: two complex conjugate roots:

    double imag = sqrt(-D) * factor;
    double real = -b       * factor;
    roots[0]    = Complex(real,  imag);
    roots[1]    = Complex(real, -imag);
  }

  return roots;
}

rsArrayTools<Complex> rosic::getRootsOfCubicEquation(double a, double b, double c, double d)
{
  // catch degenerate cases where the leading coefficient is zero:
  if( a == 0.0 )
    return getRootsOfQuadraticEquation(b, c, d);

  rsArrayTools<Complex> y(3);
  rsArrayTools<Complex> roots(3);

  // compute p,q as in the Bronstein page 40, Eq. 1.154c and the offset for the substitution
  // y = x + b/(3*a):
  double p = (3.0*a*c-b*b)/(9.0*a*a);
  double q = (b*b*b)/(27.0*a*a*a) - (b*c)/(6.0*a*a) + d/(2.0*a);

  double u, r, D, phi, ch, sh, re, im, tmp;

  if( p == 0.0 && q == 0.0 )
  {
    y[0] = y[1] = y[2] = 0.0;         // a triple real root at y=0
    // checked with y = (x-1)^3 = x^3-3*x^2+3*x-1
  }
  else if( p != 0.0 && q == 0.0 )
  {
    y[2] = 0.0;                       // a real root at y=0 and ...
    u    = -3.0*p;
    if( u > 0.0 )
    {
      tmp  =  sqrt(u);
      y[0] =  tmp;
      y[1] = -tmp;                    // ... two additional real roots or ...
      // checked with y = (x-1)*(x-2)*(x-3) = x^3-6*x^2+11*x-6
    }
    else // u < 0.0
    {
      tmp  =  sqrt(-u);
      y[0] =  Complex(0.0,  tmp);
      y[1] =  Complex(0.0, -tmp);     // ... two imaginary roots
      // checked with y = (x-4i)*(x+4i)*(x-0) = x^3+16*x
    }
  }
  else if( p == 0.0 && q != 0.0 )
  {
    u = -2.0*q;
    if( u > 0.0 )
    {
      tmp  = pow(u, 1.0/3.0);
      y[2] = tmp;                     // a real root at a positive y or ...
      phi  = (2.0/3.0)*PI;
      // checked with x^3+3*x^2+3*x
    }
    else // u < 0.0
    {
      tmp  = pow(-u, 1.0/3.0);
      y[2] = -tmp;                    // ... a real root at a negative y and ...
      phi  = PI/3.0;
      // checked with x^3+3*x^2+3*x+10
    }
    RAPT::rsSinCos(phi, &im, &re);
    re  *= tmp;
    im  *= tmp;
    y[0] = Complex(re,  im);
    y[1] = Complex(re, -im);           // ... two complex conjugate roots
  }
  else // both p and q are nonzero
  {
    r = RAPT::rsSign(q) * sqrt(fabs(p));
    if( p > 0.0 )
    {
      phi       = asinh( q/(r*r*r) );
      RAPT::rsSinhCosh(phi/3.0, &sh, &ch);
      y[0] = Complex(r*sh,  sqrt(3.0)*r*ch);
      y[1] = Complex(r*sh, -sqrt(3.0)*r*ch);
      y[2] = -2.0*r*sh;
      // checked with y = (x-i)*(x+i)*(x-1) = x^3-x^2+x-1
    }
    else // p < 0.0
    {
      D = q*q + p*p*p;
      if( D > 0.0 )
      {
        phi  = acosh( q/(r*r*r) );
        RAPT::rsSinhCosh(phi/3.0, &sh, &ch);
        y[0] = Complex(r*ch,  sqrt(3.0)*r*sh);
        y[1] = Complex(r*ch, -sqrt(3.0)*r*sh);
        y[2] = -2.0*r*ch;
        // checked with y = (x-i)*(x+i)*(x-3) = x^3-3*x^2+x-3
      }
      else // D <= 0.0
      {
        phi  = acos( q/(r*r*r) );
        y[0] =  2.0*r*cos(PI/3.0 + phi/3.0);
        y[1] =  2.0*r*cos(PI/3.0 - phi/3.0);
        y[2] = -2.0*r*cos(         phi/3.0);       // three distinct real roots
        // checked
      }
    }
  }

  // obtain the results for the original equation (back-substitution):
  double s = b/(3.0*a);
  roots[0] = y[0]-s;
  roots[1] = y[1]-s;
  roots[2] = y[2]-s;

  return roots;
}

double rosic::getCubicRootNear(double x, double a, double b, double c, double d,
                               double min, double max, int maxIterations)
{
  double f, df, xNew;

  f    = ((a*x+b)*x+c)*x+d;
  df   = (3.0*a*x+2.0*b)*x+c;
  xNew = x - f/df;
  int i = 1;
  while( xNew != x && i < maxIterations )
  {
    x    = xNew;
    f    = ((a*x+b)*x+c)*x+d;
    df   = (3.0*a*x+2.0*b)*x+c;
    xNew = x - f/df;
    i++;
  }

  return clip(xNew, min, max);
}

double rosic::getRootNear(double x, double *a, int order, double min, double max,
                          int maxIterations)
{
  // Newton/Raphson iteration:
  double f, df, xNew;
  RAPT::rsPolynomial<double>::evaluatePolynomialAndDerivativeAt(x, a, order, &f, &df);
  xNew  = x - f/df;
  int i = 1;
  while( xNew != x && i < maxIterations )
  {
    x    = xNew;
    RAPT::rsPolynomial<double>::evaluatePolynomialAndDerivativeAt(x, a, order, &f, &df);
    xNew = x - f/df;
    i++;
  }
  return clip(xNew, min, max);
}

void rosic::cubicCoeffsTwoPointsAndDerivatives(double *a, double *x, double *y, double *dy)
{
  // compute intermediate variables:
  double x0_2 = x[0]*x[0]; // x[0]^2
  double x0_3 = x0_2*x[0]; // x[0]^3
  double x1_2 = x[1]*x[1]; // x[1]^2
  double x1_3 = x1_2*x[1]; // x[1]^3
  double k1   = 3*x[0]*x1_2;
  double k2   = -3*x[1]*y[1];
  double k3   = dy[1]-dy[0];
  double s    = 1/(-x1_3+k1-3*x0_2*x[1]+x0_3);  // scaler

  a[0] =  s*(x0_2*(x1_2*k3+k2) + x0_3*(y[1]-x[1]*dy[1]) + x[0]*x1_3*dy[0] + y[0]*(-x1_3+k1));
  a[1] = -s*(x[0]*(x1_2*(2*dy[1]+dy[0])-6*x[1]*y[1]) - x0_3*dy[1] + x0_2*x[1]*(-dy[1]-2*dy[0])
             + x1_3*dy[0] + 6*x[0]*x[1]*y[0]);
  a[2] =  s*(x[0]*(x[1]*k3-3*y[1]) + x1_2*(dy[1]+2*dy[0]) + x0_2*(-dy[0]-2*dy[1]) + k2 
             + y[0]*(3*x[1]+3*x[0]));
  a[3] = -s*(x[1]*(dy[1]+dy[0]) + x[0]*(-dy[1]-dy[0]) - 2*y[1] + 2*y[0]);
}

void rosic::besselPolynomial(double *a, int order)
{
  int m, n;
  for(n=0; n<=order; n++)
    a[n] = 0.0;

  if( order == 0 )
  {
    a[0] = 1.0;
    return;
  }
  else if( order == 1 )
  {
    a[0] = 1.0;
    a[1] = 1.0;
    return;
  }
   
  // the general case is treated by recursion:
  a[0]       = 1.0;
  double *b1 = new double[order+1]; 
  double *b2 = new double[order+1]; 
  b2[0]      = 1.0;
  b2[1]      = 0.0;
  b1[0]      = 1.0;
  b1[1]      = 1.0;
  for(n=2; n<=order; n++)
  {
    double c = (double) (2*n-1);
    for(m=n; m>0; m--)
      a[m] = c*b1[m-1];
    a[0] = 0;
    for(m=0; m<n-1; m++)
      a[m] += b2[m];
    for(m=0; m<n; m++)
      b2[m] = b1[m];
    for(m=0; m<=n; m++)
      b1[m] = a[m];
  }
  delete[] b1; 
  delete[] b2; 
}

void rosic::legendrePolynomial(double *a, int order)
{
  if(order == 0) 
  {
    a[0] = 1.0;
    return;
  }
  if(order == 1) 
  {
    a[0] = 0.0;
    a[1] = 1.0;
    return;
  }

  a[0] = -0.5;
  a[1] =  0.0;
  a[2] =  1.5;
  if(order == 2) 
    return;

  int i, j;
  double *b1 = new double [order+1];
  double *b2 = new double [order+1];

  for(i = 0; i <= order; i++) 
  {
    b1[i] = b2[i] = 0.0;
  }
  b2[1] = 1.0;

  for(i = 3; i <= order; i++) 
  {
    for(j = 0; j <= i; j++) 
    {
      b1[j] = b2[j];
      b2[j] = a[j];
      a[j]  = 0.0;
    }
    for(j = i-2; j >= 0; j-=2) 
    {
      a[j] -= (i-1)*b1[j]/i;
    }
    for(j = i-1; j >= 0; j-=2) 
    {
      a[j+1] += (2*i-1)*b2[j]/i;
    }
  }
  delete [] b1;
  delete [] b2;
}

void rosic::maximumSlopeMonotonicPolynomial(double *w, int n)
{
  double *a,*p,*s,*v,c0,c1;
  int i,j,k;

  a = new double [n+1];
  p = new double [2*n+1];
  s = new double [2*n+1];
  v = new double [2*n+4];

  k = (n-1)/2;
  //  n = 2k + 1 for odd 'n' and n = 2k + 2 for even 'n':
  //  n: 1 2 3 4 5 6 ...
  //  k: 0 0 1 1 2 2 ...

  //  form vector of 'a' constants:
  if(n & 1)                   // odd
  {               
    for(i = 0; i <= k; i++) 
    {
      //a[i] = (2.0*i+1.0)/(M_SQRT2*(k+1.0));
      a[i] = (2.0*i+1.0)/(SQRT2*(k+1.0));
    }
  }                           // even
  else 
  {
    for(i = 0; i < k+1; i++) 
    {
      a[i] = 0.0;
    }
    if(k & 1) 
    {
      for(i = 1; i <= k; i+=2) 
      {
        a[i] = (2*i+1)/sqrt((double) ((k+1)*(k+2)));
      }
    }
    else 
    {
      for(i = 0; i <= k; i+=2) 
      {
        a[i] = (2*i+1)/sqrt((double) ((k+1)*(k+2)));
      }
    }
  }
  for(i = 0; i <= n; i++)
  {
    s[i] = 0.0;
    w[i] = 0.0;
  }

  // form s[] = sum of a[i]*P[i]
  s[0] = a[0];
  s[1] = a[1];
  for(i = 2; i <= k; i++) 
  {
    legendrePolynomial(p, i);
    for(j = 0; j <= i; j++) 
    {
      s[j] += a[i]*p[j];
    }
  }
  
  //  form v[] = square of s[]:
  for(i = 0; i <= 2*k+2; i++) 
  {
    v[i] = 0.0;
  }
  for(i = 0; i <= k; i++) 
  {
    for(j = 0; j <= k;j++) 
    {
      v[i+j] += s[i]*s[j];    
    }
  }
  
  // modify integrand for even 'n':
  v[2*k+1] = 0.0;
  if((n & 1) == 0) 
  {
    for(i = n; i >= 0; i--) 
    {
      v[i+1] += v[i];
    }
  }
  
  // form integral of v[]:
  for(i = n+1; i >= 0; i--) 
  {
    v[i+1] = v[i]/(double)(i+1.0);
  }
  v[0] = 0.0;

  // clear s[] for use in computing definite integral:
  for(i = 0; i < (n+2); i++)
  { 
    s[i] = 0.0;
  }
  s[0] = -1.0;
  s[1] =  2.0;

  // calculate definite integral:
  for(i = 1; i <= n; i++) 
  {
    if(i > 1) 
    {
      c0 = -s[0];
      for(j = 1; j < i+1; j++) 
      {
        c1 = -s[j] + 2.0*s[j-1];
        s[j-1] = c0;
        c0 = c1;
      }
      c1 = 2.0*s[i];
      s[i] = c0;
      s[i+1] = c1;
    }
    for(j = i; j > 0; j--) 
    {
      w[j] += (v[i]*s[j]);
    }
  }
  if((n & 1) == 0) 
    w[1] = 0.0;

  delete [] v;
  delete [] p;
  delete [] s;
  delete [] a;
}
*/