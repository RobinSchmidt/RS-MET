namespace RSLib
{

void rsLanden(double k, int M, double* v)
{
  // performs M iterations of the Landen transformation of an elliptic modulus k and returns the 
  // results in the array v which must be of length M+1
  int n;
  if(k == 0.0 || k == 1.0)
  {
    v[0] = k;
    for(n=1; n<M; n++)
      v[n] = 0.0;
    return;
  }
  else
  {
    for(n=0; n<M; n++)
    {
      k    = (k/(1.0+rsSqrt(1.0-k*k)));
      k   *= k;
      v[n] = k;
    }
  }
}

double rsIdentity(double x)
{
  return x;
}

  void rsEllipticIntegral(double k, double *K, double *Kprime, int M)
  {
    double kmin = 1e-6; 
    double kmax = rsSqrt(1-kmin*kmin);
    double kp, L;
    double* v  = new double[M];
    double* vp = new double[M];
    int n;

    if( k == 1.0 )
      *K = rsInfDouble;
    else if( k > kmax )
    {
      kp = rsSqrt(1.0-k*k);
      L  = -log(kp/4.0);
      *K = L + (L-1.0)*kp*kp/4.0; 
    }
    else
    {
      rsLanden(k, M, v);
      for(n=0; n<M; n++)
        v[n] += 1.0;
      *K = rsProduct(v, M) * 0.5*PI;
    }

    if( k == 0.0 )
      *Kprime = rsInfDouble;
    else if( k < kmin )
    {
      L       = -log(k/4.0);
      *Kprime = L + (L-1.0)*k*k/4.0;
    }
    else
    {
      kp = rsSqrt(1.0-k*k);                             
      rsLanden(kp, M, vp);
      for(n=0; n<M; n++)
        vp[n] += 1.0;
      *Kprime = rsProduct(vp, M) * 0.5*PI;
    }

    delete[] v;
    delete[] vp;
  }

  // Computes auxiliary functions f(x) and g(x) needed in the evaluation of Si(x) and Ci(x).
  // This approximation follows "Handbook of mathematical functions" by Abramowitz/Stegun, 
  // page 232.
  void rsSinCosIntegralAux(double x, double* f, double* g) // usable fo x >= 1
  {
    double x2 = x*x;    // x^2
    double x4 = x2*x2;  // x^4
    double x6 = x4*x2;  // x^6
    double x8 = x4*x4;  // x^8

    // Eq. 5.2.38:
    static const double 
      fa1 =  38.027264, fb1 =  40.021433,
      fa2 = 265.187033, fb2 = 322.624911,
      fa3 = 335.677320, fb3 = 570.236280,
      fa4 =  38.102495, fb4 = 157.105423;
    *f = (x8 + fa1*x6 + fa2*x4 + fa3*x2 + fa4) /
      (x*(x8 + fb1*x6 + fb2*x4 + fb3*x2 + fb4));

    // Eq. 5.2.39:
    static const double 
      ga1 =  42.242855, gb1 =   48.196927,
      ga2 = 302.757865, gb2 =  482.485984,
      ga3 = 352.018498, gb3 = 1114.978885,
      ga4 =  21.821899, gb4 =  449.690326;
    *g =  (x8 + ga1*x6 + ga2*x4 + ga3*x2 + ga4) /
      (x2*(x8 + gb1*x6 + gb2*x4 + gb3*x2 + gb4));
  }
  double rsSineIntegralViaAuxFunctions(double x) // usable for |x| > 1
  { 
    double xAbs = fabs(x);      
    double s, c, f, g;
    rsSinCos(xAbs, &s, &c);
    rsSinCosIntegralAux(xAbs, &f, &g);
    return rsSign(x) * (PI/2 - f*c - g*s); // Eq. 5.2.8
  }
  double rsSineIntegralViaPowerSeries(double x)  // recommended for |x| <= 3
  {
    double xAbs = fabs(x);  // |x|
    double sum  = 0.0;      // accumulator for result
    double mx2  = -x*x;     // -x^2
    int    k    = 1;        // k := 2*n+1
    double xk   = x;        // accumulator for x^k = x^(2*n+1)
    for(int n = 0; n < 15; n++)
    {
      double old = sum;
      sum += rsInverseFactorials[k] * xk / k;
      if( sum == old )
        break; // converged - maybe measure, if this really leads to efficiency improvements on 
               // the average -> throw arguments between 0 and 3 at the function, see if it performs 
               // better with or without the conditional early exit
               // also, if we know that |x| < 3, we might get away with fewer than 15 terms, 
               // in experiments, it never exceeded 13, when |x| < 3 ...but maybe we should shift
               // the crossover-point for the algorithm higher than 3
      xk *= mx2;
      k  += 2;
    }
    return sum;
  }
  double rsSineIntegralViaContinuedFractions(double x)  // recommended for |x| > 3
  {
    static const int maxIterations = 100;
    double xAbs = fabs(x);
    double a;
    rsComplexDbl b(1.0, xAbs);
    rsComplexDbl c = 1.0 / TINY;
    rsComplexDbl d = 1.0 / b;
    rsComplexDbl h = d;
    rsComplexDbl cd;
    for(int i = 1; i <= maxIterations; i++)
    {
      a   = -i*i;
      b  += 2.0;
      d   = 1.0 / (a*d + b);
      c   = b + a/c;
      cd  = c*d;             // factor
      h  *= cd;              // product
      if( fabs(cd.re-1.0) + fabs(cd.im) < EPS ) 
        break;
    }
    double si, co;
    rsSinCos(xAbs, &si, &co);
    h *= rsComplexDbl(co, -si); 
    return rsSign(x) * (0.5*PI + h.im);
  }
  double rsSineIntegral(double x)
  {
    double xAbs = fabs(x);
    if( xAbs <= 3.0 )
      return rsSineIntegralViaPowerSeries(x);
    else
    {
      return rsSineIntegralViaContinuedFractions(x); // exact but slow (really? -> measure)
      //return rsSineIntegralViaAuxFunctions(x);       // approximation
    }
  }

  double rsSinc(double x)
  {
    if( rsAbs(x) < EPS )
      return 1.0;
    return sin(x) / x;
  }

  double rsNormalizedSinc(double x)
  {
    return rsSinc(PI*x);
  }

}
