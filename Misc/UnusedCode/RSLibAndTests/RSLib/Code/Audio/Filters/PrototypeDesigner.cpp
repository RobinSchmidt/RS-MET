using namespace RSLib;

void RSLib::maximumSlopeMonotonicPolynomial(double *w, int n)
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
        a[i] = (2*i+1)/rsSqrt((double) ((k+1)*(k+2)));
      }
    }
    else 
    {
      for(i = 0; i <= k; i+=2) 
      {
        a[i] = (2*i+1)/rsSqrt((double) ((k+1)*(k+2)));
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


void updateLegendrePolynomial(double **P, double *P1, double *P2, int r)
{
  if( rsIsEven(r) )
  {
    rsLegendrePolynomialRecursion(P1, r, P2, P1);
    *P = P1;
  }
  else
  {
    rsLegendrePolynomialRecursion(P2, r, P1, P2);
    *P = P2;
  }
}
void RSLib::rsPapoulisPolynomial(double *v, int N)
{
  // temporary arrays for Legendre polynomials:
  double *P1 = new double[N/2+1];
  double *P2 = new double[N/2+1];
  double *P  = nullptr;  // pointer to the current P array

  // create integrand:
  int k, r;
  if( rsIsOdd(N) )
  {
    k = (N-1)/2;

    // create weighted sum of Legendre polynomials in v:
    rsFillWithZeros(v, k+1);
    for(r = 0; r <= k; r++)
    {
      updateLegendrePolynomial(&P, P1, P2, r);
      weightedSumOfPolynomials(v, r, 1.0, P, r, 2*r+1.0, v);
    }

    // square it:
    rsConvolve(v, k+1, v, k+1, v);
  }
  else
  {
    k = (N-2)/2;

    // generate Legendre polynomial of order k+1, store in P:
    for(r = 0; r <= k+1; r++)
      updateLegendrePolynomial(&P, P1, P2, r);

    // take the derivative, store in v:
    polyDerivative(P, v, k+1);

    // square it:
    rsConvolve(v, k+1, v, k+1, v);

    // multiply by (x+1):
    v[2*k+1] = 0;
    for(r = 2*k+1; r >= 1; r--)
      v[r] += v[r-1];
  }

  // integrate from -1 to 2*w^2-1:
  double a[1] = { -1 };  
  double b[3] = { -1, 0, 2};
  integratePolynomialWithPolynomialLimits(v, N-1, a, 0, b, 2, v);

  // scale, such that L^2(1) = 1:
  rsScale(v, 2*N+1, 1.0 / rsSum(v, 2*N+1));

  // clean up:
  delete[] P1;
  delete[] P2;
}

void rsHalpernU(double *a, int K)
{
  // Computes coefficients of the U-polynomials given in "Design and Analysis of Analog Filters", 
  // page 256-257 except for the scale factor in front.
  rsFillWithZeros(a, K+1);
  rsUint64 k, m;
  if( rsIsEven(K) )
  {
    k = K/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*m] = (-1)^(k-m) * (m+k)! / ( (k-m)! * (m!)^2 ):
      a[2*m] = (double) (rsProduct(k-m+1, m+k) / rsSquare(rsFactorial(m)));
      if( rsIsOdd(k-m) )
        a[2*m] *= -1;
    }
  }
  else
  {
    k = (K-1)/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*(k-m)+1] = (-1)^m * (2*k+1-m)! / (m! * (k+1-m)! * (k-m)!):
      a[2*(k-m)+1] = (double) (rsProduct(k-m+2, 2*k+1-m) / (rsFactorial(m)*rsFactorial(k-m)));
      if( rsIsOdd(m) )
        a[2*(k-m)+1] *= -1;
    }
  }
}
void RSLib::rsHalpernPolynomial(double *a, int N)
{  
  a[0] = 0;
  double *a1 = &a[1];           // index shift of one for multiplication by x in Eq. 8.19
  rsHalpernU(a1, N-1);          // create U-polynomial
  rsConvolve(a1, N, a1, N, a1); // square U-polynomial
  rsScale(a1, 2*N-1, 2*N);      // apply squared scale factor
  polyIntegral(a, a, 2*N-1);    // compute integral from 0 to w
}

void RSLib::rsGaussianPolynomial(double *a, int N, double wc)
{  
  rsFillWithZeros(a, 2*N+1);
  double g = log(2.0) / (wc*wc);  // gamma
  double s = 1;                   // scaler
  for(int k = 0; k <= N; k++)
  {
    a[2*k] = s;    // == g^k / k! == pow(g, k) / rsFactorial(k);
    s *= g/(k+1);
  }
}


//=================================================================================================
// class rsPrototypeDesigner

// construction/destruction:

rsPrototypeDesigner::rsPrototypeDesigner()
{
  L                     = 2;
  r                     = 0;
  N                     = L+r;
  approximationMethod   = BUTTERWORTH;
  prototypeMode         = LOWPASS_PROTOTYPE;
  numFinitePoles        = 2;
  numFiniteZeros        = 0;
  Ap                    = rsAmp2dB(sqrt(2.0));   // 3.01 dB passband ripple for lowpasses
  As                    = 60.0;                  // 60.0 dB stopband attenuation for lowpasses
  A                     = 0.0;                   // cut/boost in dB for shelvers
  A0                    = 0.0;                   // reference gain in dB for shelvers
  Rp                    = 0.95;                  // inner ripple as fraction of dB-peak-gain for shelv
  Rs                    = 0.05;                  // outer ripple as fraction of peak
  stateIsDirty          = true;                  // poles and zeros need to be evaluated
  updatePolesAndZeros();

  /*
  // for obtaining the scale-factors for the Bessel-filters:
  setApproximationMethod(BESSEL);
  setOrder(25);
  makeBesselLowpass();
  double p1 = p[0].getRadius();
  makeBesselLowpassFromTable();
  double p2 = p[0].getRadius();
  double factor = p2/p1;
  int dummy = 0;
  */
}

rsPrototypeDesigner::~rsPrototypeDesigner()
{

}

// setup:

void rsPrototypeDesigner::setOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder != N )
  {
    N = newOrder;
    if( rsIsOdd(N) )
    {
      r = 1;
      L = (N-1)/2;
    }
    else
    {
      r = 0;
      L = N/2;
    }
    stateIsDirty = true;
  }
}

void rsPrototypeDesigner::setApproximationMethod(int newApproximationMethod)
{
  if( newApproximationMethod < BUTTERWORTH || newApproximationMethod > PAPOULIS )
    rsError("Unknown approximation method");

  if( newApproximationMethod != approximationMethod )
  {
    approximationMethod = newApproximationMethod;
    stateIsDirty        = true;
  }
}

void rsPrototypeDesigner::setPrototypeMode(int newPrototypeMode)
{
  if( newPrototypeMode == LOWPASS_PROTOTYPE || newPrototypeMode == LOWSHELV_PROTOTYPE )
  {
    prototypeMode = newPrototypeMode;
    stateIsDirty  = true;
  }
  else
    rsError("Unknown mode");
}

void rsPrototypeDesigner::setPassbandRipple(double newPassbandRipple)
{
  if( newPassbandRipple >= 0.0 )
  {
    Ap           = newPassbandRipple;
    stateIsDirty = true;
  }
  else
    rsError("Ripple (in dB) must be >= 0");
}

void rsPrototypeDesigner::setStopbandRejection(double newStopbandRejection)
{
  if( newStopbandRejection >= 0.0 )
  {
    As           = newStopbandRejection;
    stateIsDirty = true;
  }
  else
    rsError("Ripple (in dB) must be >= 0");
}

void rsPrototypeDesigner::setGain(double newGain)
{
  if( newGain != A )
  {
    A            = newGain;
    stateIsDirty = true;
  }
}

void rsPrototypeDesigner::setReferenceGain(double newReferenceGain)
{
  if( newReferenceGain != A0 )
  {
    A0           = newReferenceGain;
    stateIsDirty = true;
  }
}

void rsPrototypeDesigner::setPassbandGainRatio(double newPassbandGainRatio)
{
  if( newPassbandGainRatio >= 1.0 || newPassbandGainRatio <= 0.0 || newPassbandGainRatio < Rs )
  {
    rsError("Bandwidth gain ratio makes no sense (inequation 51 is violated)");
    return;
  }
  if( newPassbandGainRatio != Rp  )
  {
    Rp           = newPassbandGainRatio;
    stateIsDirty = true;
  }
}

void rsPrototypeDesigner::setStopbandGainRatio(double newStopbandGainRatio)
{
  if( newStopbandGainRatio >= 1.0 || newStopbandGainRatio <= 0.0 || newStopbandGainRatio > Rp )
  {
    rsError("Stopband gain ratio makes no sense (inequation 51 is violated)");
    return;
  }
  if( newStopbandGainRatio != Rp  )
  {
    Rs           = newStopbandGainRatio;
    stateIsDirty = true;
  }
}

// static member functions:

void rsPrototypeDesigner::getNumBiquadsAndFirstOrderStages(int N, int &L, int &r)
{
  if( rsIsOdd(N) )
  {
    r = 1;
    L = (N-1)/2;
  }
  else
  {
    r = 0;
    L = N/2;
  }
}

double rsPrototypeDesigner::ellipdeg(int N, double k_1)
{
  int L;
  if( rsIsEven(N) )
    L = N/2;
  else
    L = (N-1)/2;

  double kmin = 1e-6;
  double k;

  if (k_1 < kmin)
    k = ellipdeg2(1.0 / (double) N, k_1);
  else
  {
    double kc = sqrt(1-k_1*k_1);			  // complement of k1

    double u_i;
    double prod = 1.0;
    for(int i=1; i<=L; i++)
    {
      u_i   = (double) (2*i-1) / (double) N;
      prod *= rsSnC(rsComplexDbl(u_i), kc).re; 
        // Can't we use the real-number version of the elliptic function here? Both arguments are
        // actually real.
    }
    prod      = prod*prod*prod*prod;
    double kp = pow(kc, (double) N) * prod; // complement of k
    k         = sqrt(1.0-kp*kp);
  }

  return k;
}

double rsPrototypeDesigner::ellipdeg1(int N, double k)
{
  int L;
  if( rsIsEven(N) )
    L = N/2;
  else
    L = (N-1)/2;

  double u_i;
  double prod = 1.0;
  for(int i=1; i<=L; i++)
  {
    u_i   = (double) (2*i-1) / (double) N;
    prod *= rsSnC(rsComplexDbl(u_i), k).re;
      // Can't we use the real-number version of the elliptic function here? Both arguments are
      // actually real.
  }
  prod      = prod*prod*prod*prod;
  double k1 = pow(k, (double) N) * prod;

  return k1;
}

double rsPrototypeDesigner::ellipdeg2(double N, double k)
{
  int M = 7;

  double K;
  double Kprime;
  rsEllipticIntegral(k, &K, &Kprime);

  double q  = exp(-PI*Kprime/K);
  double q1 = pow(q, N);

  int m;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for(m=1; m<=M; m++)
  {
    sum1 += pow(q1, (double) (m*(m+1)) );
    sum2 += pow(q1, (double) (m*m)     );
  }

  double tmp = (1.0+sum1)/(1.0+2.0*sum2);
  tmp       *= tmp;
  double k1  = 4.0 * sqrt(q1) * tmp;

  return k1;
}

double rsPrototypeDesigner::getRequiredButterworthOrder(double passbandFrequency, 
  double passbandRipple, double stopbandFrequency, double stopbandRipple)
{
  double Gp = pow(10.0, -passbandRipple/20.0);                      // (1),Eq.1
  double Gs = pow(10.0, -stopbandRipple/20.0);                      // (1),Eq.1
  double ep = sqrt(1.0 / (Gp*Gp) - 1.0);                            // (1),Eq.2
  double es = sqrt(1.0 / (Gs*Gs) - 1.0);                            // (1),Eq.2
  return log(es/ep) / log(stopbandFrequency/passbandFrequency);     // (1),Eq.9
}

double rsPrototypeDesigner::getRequiredChebychevOrder(double passbandFrequency, 
  double passbandRipple, double stopbandFrequency, double stopbandRipple)
{
  double Gp = pow(10.0, -passbandRipple/20.0);                           // (1),Eq.1
  double Gs = pow(10.0, -stopbandRipple/20.0);                           // (1),Eq.1
  double ep = sqrt(1.0 / (Gp*Gp) - 1.0);                                 // (1),Eq.2
  double es = sqrt(1.0 / (Gs*Gs) - 1.0);                                 // (1),Eq.2
  return rsAcosh(es/ep) / rsAcosh(stopbandFrequency/passbandFrequency);  // (1),Eq.9
}

double rsPrototypeDesigner::getRequiredEllipticOrder(double passbandFrequency, double passbandRipple,
  double stopbandFrequency, double stopbandRipple)
{
  double Gp = pow(10.0, -passbandRipple/20.0);                       // (1),Eq.1
  double Gs = pow(10.0, -stopbandRipple/20.0);                       // (1),Eq.1
  double ep = sqrt(1.0 / (Gp*Gp) - 1.0);                             // (1),Eq.2
  double es = sqrt(1.0 / (Gs*Gs) - 1.0);                             // (1),Eq.2
  double k  = passbandFrequency / stopbandFrequency;                 // (1),Eq.3
  double k1 = ep/es;                                                 // (1),Eq.3
  double  K, Kp, K1, K1p;
  rsEllipticIntegral(k,  &K,  &Kp);                                  // (1),Eq.19
  rsEllipticIntegral(k1, &K1, &K1p);
  return (K1p*K)/(K1*Kp);                                            // (1),Eq.34
}

void rsPrototypeDesigner::magSquaredNumAndDen(double *b, double *a, double *b2, double *a2, int N)
{
  double *am = new double[N+1];
  double *bm = new double[N+1];
  polyCoeffsForNegativeArgument(b, bm, N);  // coeffs of N(-s)
  polyCoeffsForNegativeArgument(a, am, N);  // coeffs of D(-s)
  multiplyPolynomials(b, N, bm, N, b2);     // coeffs of N(s)*N(-s)
  multiplyPolynomials(a, N, am, N, a2);     // coeffs of D(s)*D(-s)
  delete[] am;
  delete[] bm;
}

void rsPrototypeDesigner::shelvingMagSqrNumFromLowpassMagSqr(double *b2, double *a2, double k, 
  int N, double G0, double G, double *bShelf)
{
  rsWeightedSum(b2, a2, bShelf, 2*N+1, k*k*(G*G-G0*G0), G0*G0);

  /*
  // debug:
  double b2d[20], a2d[20], bsd[20];
  copyBuffer(b2,     b2d, 2*N+1);
  copyBuffer(a2,     a2d, 2*N+1);
  copyBuffer(bShelf, bsd, 2*N+1);
  int dummy = 0;
  */
}

// factor out shelvingMagSqrNumeratorFromLowpassMagSqr:
void rsPrototypeDesigner::shelvingMagSqrNumeratorFromLowpassTransfer(double *b, double *a, 
  double k, int N, double G0, double G, double *bShelf)
{
  double *a2 = new double[2*N+1];
  double *b2 = new double[2*N+1];

  // construct lowpass magnitude squared numerator and denominator N_LP(s)*N_LP(-s), 
  // D_LP(s)*D_LP(-s):
  magSquaredNumAndDen(b, a, b2, a2, N);

  // obtain coefficients for shelving filter's magnitude squared function numerator polynomial 
  // N_LS(s)*N_LS(-s):
  shelvingMagSqrNumFromLowpassMagSqr(b2, a2, k, N, G0, G, bShelf);

  delete[] a2;
  delete[] b2;
}

void rsPrototypeDesigner::scaleToMatchGainAtUnity(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double g)
{
  double  wc    = rsFilterAnalyzer::findAnalogFrequencyWithMagnitude(z, p, k, N, g, 1.0);
  double scaler = 1.0/wc;
  for(int n = 0; n < N; n++)
  {
    pNew[n] = scaler * p[n];
    zNew[n] = scaler * z[n];
  }
  int nz = rsGetNumFiniteValues(z, N);
  *kNew  = *k / pow(wc, N-nz);
}

void rsPrototypeDesigner::getInverseFilter(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N)
{
  rsComplexDbl *zTmp = new rsComplexDbl[N];
  rsCopyBuffer(z,    zTmp, N);
  rsCopyBuffer(p,    z,    N);
  rsCopyBuffer(zTmp, p,    N);
  *kNew = 1.0 / *k;
  delete[] zTmp;
}

int rsPrototypeDesigner::getLeftHalfPlaneRoots(double *a, rsComplexDbl *r, int N)
{
  rsComplexDbl *rTmp = new rsComplexDbl[N]; // maybe we can get rid of that temporary array
  findPolynomialRoots(a, N, rTmp);  
  int numLeftRoots = rsOnlyLeftHalfPlane(rTmp, r, N);

  /*
  // debug stuff:
  double aD[50];
  copyBuffer(a, aD, N+1);
  rsComplexDbl rL[50];
  fillWithValue(rL, 50, rsComplexDbl(rsInfDouble));
  numLeftRoots = onlyLeftHalfPlane(&rTmp[1], rL, N);
  */
  
  rsAssert(numLeftRoots == ceil(0.5*N)); // maybe take this out later
  delete[] rTmp;
  return numLeftRoots;
}

void rsPrototypeDesigner::getBesselLowpassZerosPolesAndGain(rsComplexDbl *z, rsComplexDbl *p, 
  double *k, int N)
{     
  // zeros are at infinity:
  rsFillWithValue(z, N, rsComplexDbl(rsInfDouble, 0.0));

  // find poles:
  double *a = new double[N+1];        // Bessel-Polynomial coefficients
  besselPolynomial(a, N);
  rsReverse(a, N+1);                  // we actually use a reverse Bessel polynomial

  findPolynomialRoots(a, N, p);  

  // set gain and scale poles to match Butterworth magnitude response asymptotically, if desired:
  bool matchButterworth = true; // maybe make this a parameter later
  if( matchButterworth == true )
  {
    double scaler = 1.0 / pow(a[0], 1.0/N);
    for(int n = 0; n < N; n++)
      p[n] *= scaler;
    *k = 1.0; 
  }
  else
    *k = a[0];

  delete[] a;
}

void rsPrototypeDesigner::getBesselLowShelfZerosPolesAndGain(rsComplexDbl *z, rsComplexDbl *p, double *k, int N, double G, double G0)
{
  // old version - needs root-finder twice:
  //getBesselLowpassZerosPolesAndGain(z, p, k, N);
  //PoleZeroMapper::sLowpassToLowshelf(z, p, k, z, p, k, N, G0, G);
  //return;

  // catch lowpass case:  
  if( G0 == 0.0 )
  {
    getBesselLowpassZerosPolesAndGain(z, p, k, N);
    *k *= G;
    return;
  } 

  // design boost filter and invert later, if a dip is desired:
  bool dip = false;
  if( G < G0 )
  {
    dip = true;
    G   = 1.0 / G;
    G0  = 1.0 / G0;
  }

  // construct lowpass denominator:
  double *a  = new double[N+1];   
  besselPolynomial(a, N);
  rsReverse(a, N+1);   // leaving this out leads to a modified Bessel filter response - maybe 
                       // experiment a bit, response looks good

  // find poles of the shelving filter:
  findPolynomialRoots(a, N, p);  

  // construct lowpass numerator:
  double *b = new double[N+1];        
  rsFillWithZeros(b, N+1);
  b[0] = a[0];

  // obtain magnitude-squared numerator polynomial for shelving filter:
  double *bS = new double[2*N+1];       
  shelvingMagSqrNumeratorFromLowpassTransfer(b, a, 1.0, N, G0, G, bS);

  // find left halfplane zeros (= zeros of the shelving filter):
  getLeftHalfPlaneRoots(bS, z, 2*N);  

  // set gain constant:
  *k = G0;          

  // now we have a shelving filter with correct low-frequency gain G and reference gain G0, but possibly still with wrong bandwidth gain 
  // GB at unity - now we adjust zeros/poles/gain to match GB:
  double GB = sqrt(G*G0);
  scaleToMatchGainAtUnity(z, p, k, z, p, k, N, GB);

  // invert filter in case of a dip:
  if( dip == true )
    getInverseFilter(z, p, k, z, p, k, N);

  // cleanup:
  delete[] a;
  delete[] b;
  delete[] bS;
}


void rsPrototypeDesigner::papoulisMagnitudeSquaredDenominator(double *a, int N)
{
  int n;
  rsFillWithZeros(a, 2*N+1);  // do we need this?

  // construct the polynomial L_N(w^2):
  maximumSlopeMonotonicPolynomial(a, N);  // does the same same as lopt(a, N); from C.R.Bond
    // use rsPapoulisPolynomial instead


  // flip sign of coeffs for odd powers (substitute w^2 with -s^2):
  for(n = 1; n <= N; n += 2)
    a[n] = -a[n];

  // convert polynomial in s^2 to the corresponding polynomial in s: 
  for(n = N; n >= 0; n--)
    a[2*n] = a[n];
  for(n = 1; n <= 2*N; n += 2)
    a[n] = 0.0;

  // add the constant 1 to the polynomial:
  a[0] += 1.0;
}

void rsPrototypeDesigner::getPapoulisLowpassZerosPolesAndGain(rsComplexDbl *z, rsComplexDbl *p, double *k, int N)
{ 
  // find poles:
  double *a2 = new double[2*N+1];     // coefficients of the magnitude-squared polynomial D(s)*D(-s)
  papoulisMagnitudeSquaredDenominator(a2, N);
  getLeftHalfPlaneRoots(a2, p, 2*N);  

  // zeros are at infinity:
  rsFillWithValue(z, N, rsComplexDbl(rsInfDouble, 0.0));

  // set gain at DC to unity:
  *k = sqrt(1.0/fabs(a2[2*N])); 

  delete[] a2;
}

void rsPrototypeDesigner::getPapoulisLowShelfZerosPolesAndGain(rsComplexDbl *z, rsComplexDbl *p, double *k, int N, double G, double G0)
{
  //getPapoulisLowpassZerosPolesAndGain(z, p, k, N);
  //PoleZeroMapper::sLowpassToLowshelf(z, p, k, z, p, k, N, G0, G);
  //return;


  // catch lowpass case:  
  if( G0 == 0.0 )
  {
    getPapoulisLowpassZerosPolesAndGain(z, p, k, N);
    *k *= G;
    return;
  } 

  // design boost filter and invert later, if a dip is desired:
  bool dip = false;
  if( G < G0 )
  {
    dip = true;
    G   = 1.0 / G;
    G0  = 1.0 / G0;
  }


  // factor out into a function getMagnitudeSquaredNumAndDen, the call this function here - this function can switch between Papoulis,
  // Gauss, etc. and fill the arrays accordingly, Bessel is a special case - it doesn't need to find poles of the mag-squared function,
  // there we can directly find the poles of the transfer function ...hmm...but maybe we can still factor out a function getPoles or 
  // something...

  // coefficients of the magnitude-squared polynomial D(s)*D(-s)
  double *a2 = new double[2*N+1];    
  papoulisMagnitudeSquaredDenominator(a2, N);
  getLeftHalfPlaneRoots(a2, p, 2*N);  

  // normalize denominator polynomial such that the leading coeff has unity as absolute value:
  double scaler = 1.0 / fabs(a2[2*N]);
  for(int n = 0; n <= 2*N; n++)
    a2[n] *= scaler;
    
  // construct lowpass numerator:
  double *b2 = new double[2*N+1];        
  rsFillWithZeros(b2, 2*N+1);
  b2[0] = 1.0;

  // end of "factor out" ...in general, we need to scale the b2-polynomial also by dividing through the leading coeff?


  // adjust lowpass DC gain via k:
  *k = sqrt(fabs(a2[0]));  // in general: sqrt(fabs(a2[0]/b2[0])) ?
  //*k = sign(a2[0]) * sqrt(fabs(a2[0]));

  // obtain magnitude-squared numerator polynomial for shelving filter:
  double *bS = new double[2*N+1];       
  shelvingMagSqrNumFromLowpassMagSqr(b2, a2, *k, N, G0, G, bS);

  // find left halfplane zeros (= zeros of the shelving filter):
  getLeftHalfPlaneRoots(bS, z, 2*N);  

  // set gain constant for shelving filter:
  *k = G0; 

  // adjust bandwidth:
  double GB = sqrt(G*G0);
  scaleToMatchGainAtUnity(z, p, k, z, p, k, N, GB);

  // invert filter in case of a dip:
  if( dip == true )
    getInverseFilter(z, p, k, z, p, k, N);

  delete[] a2;
  delete[] b2;
  delete[] bS;
}





void rsPrototypeDesigner::getEllipticLowpassZerosPolesAndGain(rsComplexDbl *z, rsComplexDbl *p, double *k, int N, double Gp, double Gs)
{
  int nz, L, r;
  if( rsIsEven(N) )
    nz = N;
  else
    nz = N-1;
  getNumBiquadsAndFirstOrderStages(N, L, r);

  // declare/assign/calculate some repeatedly needed variables:
  rsComplexDbl j(0.0, 1.0);                                    // imaginary unit
  double  ep  = sqrt(1.0/(Gp*Gp) - 1.0);                  // Eq. 2
  double  es  = sqrt(1.0/(Gs*Gs) - 1.0);                  // Eq. 2
  double  k1  = ep/es;                                    // Eq. 3
  double  kk  = ellipdeg(N, k1);                          // solve degree equation for k
  double  v_0 =  (-j*rsAsnC(j/ep, k1) / (double) N).re;   // from ellipap.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    //p[L+r-1] = -Omega_p/sinh(v_0*PI*0.5*kk);                   // Eq. 73
    p[N-1] = j * rsSnC(j*v_0, kk);                               // from ellipap.m - find Eq.
    z[N-1] = rsInfDouble;
  }

  // calculate the complex conjugate poles and zeros:
  double  u_i;
  rsComplexDbl zeta_i;
  for(int i=0; i<L; i++)
  {
    u_i      = (double) (2*(i+1)-1) / (double) N;                // Eq. 69
    zeta_i   = rsCdC(rsComplexDbl(u_i), kk);                     // from ellipap.m - find Eq.
    z[2*i]   = j / (kk*zeta_i);                                  // Eq. 62
    p[2*i]   = j*rsCdC((u_i-j*v_0), kk);
    z[2*i+1] = z[2*i].getConjugate();
    p[2*i+1] = p[2*i].getConjugate();
  }

  double H0 = pow(Gp, 1-r);    // preliminary - can be made simpler (without pow)
  rsComplexDbl n = rsProductOfFiniteFactors(p, N);
  rsComplexDbl d = rsProductOfFiniteFactors(z, N);
  *k        = H0 * (n/d).re;
}

/*
void rsPrototypeDesigner::getLowpassZerosPolesAndGain(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  int N, int approximationMethod)
{
  int dummy = 0;
}
*/

// inquiry:

rsComplexDbl rsPrototypeDesigner::getFilterResponseAt(rsComplexDbl s)
{
  rsComplexDbl num, den;
  rsComplexDbl tmp;
  int     Lz, Lp;

  // initialize the numerator and denominator:
  if( rsIsOdd(numFiniteZeros) )
  {
    num = -z[L+r-1].re;
    Lz  = (numFiniteZeros-1)/2;
  }
  else
  {
    num = 1.0;
    Lz  = numFiniteZeros/2;
  }
  if( rsIsOdd(numFinitePoles) )
  {
    den = -p[L+r-1].re;
    Lp  = (numFinitePoles-1)/2;
  }
  else
  {
    den = 1.0;
    Lp  = numFinitePoles/2;
  }

  // accumulate the product of the linear factors for the denominator (poles) and numerator 
  // (zeros):
  int i;
  for(i=0; i<Lz; i++)
    num *= ((s-z[i]) * (s-z[i].getConjugate()));
  for(i=0; i<Lp; i++)
    den *= ((s-p[i]) * (s-p[i].getConjugate()));

  // return the quotient of the calculated numerator and denominator as result:
  return num/den;
}

double rsPrototypeDesigner::getMagnitudeAt(double w)
{
  return getFilterResponseAt(rsComplexDbl(0.0, w)).getRadius();
}

double rsPrototypeDesigner::findFrequencyWithMagnitude(double magnitude, double wLow, double wHigh)
{
  // until we have something better, we search for the frequency at which the desired gain occurs 
  // by means of the bisection method:

  double wMid  = 0.5 * (wLow+wHigh);
  double mLow  = getMagnitudeAt(wLow);
  double mMid  = getMagnitudeAt(wMid);
  double mHigh = getMagnitudeAt(wHigh);

  while( wHigh-wLow > 0.0001 ) // introduce a threshold, maybe value should be lower
  {
    wMid  = 0.5 * (wLow+wHigh);
    mLow  = getMagnitudeAt(wLow);
    mMid  = getMagnitudeAt(wMid);
    mHigh = getMagnitudeAt(wHigh);

    if( mMid > magnitude )
      wLow = wMid;
    else
      wHigh = wMid;
  }

  return 0.5 * (wLow+wHigh); // preliminary
}

int rsPrototypeDesigner::getNumFinitePoles()
{
  return numFinitePoles;
}

int rsPrototypeDesigner::getNumFiniteZeros()
{
  return numFiniteZeros;
}

int rsPrototypeDesigner::getNumNonRedundantFinitePoles()
{
  if( rsIsEven(numFinitePoles) )
    return numFinitePoles/2;
  else
    return (numFinitePoles+1)/2;
}

int rsPrototypeDesigner::getNumNonRedundantFiniteZeros()
{
  if( rsIsEven(numFiniteZeros) )
    return numFiniteZeros/2;
  else
    return (numFiniteZeros+1)/2;
}

void rsPrototypeDesigner::getPolesAndZeros(rsComplexDbl *poles, rsComplexDbl *zeros)
{
  if( stateIsDirty == true )
    updatePolesAndZeros(); // re-calculate only if necesarry
  for(int i=0; i<(L+r); i++)
  {
    poles[i] = p[i];
    zeros[i] = z[i];
  }
}

bool rsPrototypeDesigner::hasCurrentMethodRippleParameter()

{
  if( prototypeMode == LOWPASS_PROTOTYPE )
  {
    if( (approximationMethod == ELLIPTIC) || (approximationMethod == CHEBYCHEV) )
      return true;
    else
      return false;
  }
  else
  {
    if( (approximationMethod == ELLIPTIC) || (approximationMethod == CHEBYCHEV)
      || (approximationMethod == INVERSE_CHEBYCHEV) )
      return true;
    else
      return false;
  }
}

bool rsPrototypeDesigner::hasCurrentMethodRejectionParameter()
{
  if( prototypeMode == LOWPASS_PROTOTYPE )
  {
    return (approximationMethod == ELLIPTIC) || (approximationMethod == INVERSE_CHEBYCHEV);
  }
  else
    return false;
}

bool rsPrototypeDesigner::needsSpecialHighShelvTransform()
{
  return (approximationMethod == BUTTERWORTH) 
    ||   (approximationMethod == CHEBYCHEV)    
    ||   (approximationMethod == INVERSE_CHEBYCHEV)
    ||   (approximationMethod == ELLIPTIC);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// pole/zero calculation:

void rsPrototypeDesigner::updatePolesAndZeros()
{
  if( stateIsDirty == true )
  {
    if( prototypeMode == LOWPASS_PROTOTYPE )
    {
      switch( approximationMethod )
      {
      case BUTTERWORTH:       makeButterworthLowpass();         break;
      case CHEBYCHEV:         makeChebychevLowpass();           break;
      case INVERSE_CHEBYCHEV: makeInverseChebychevLowpass();    break;
      case ELLIPTIC:          makeEllipticLowpass();            break;
      case BESSEL:            makeBesselLowShelv(  1.0, 0.0);   break; 
      case PAPOULIS:          makePapoulisLowShelv(1.0, 0.0);   break; 
      }
    }
    else if( prototypeMode == LOWSHELV_PROTOTYPE )
    {
      switch( approximationMethod )
      {
      case BUTTERWORTH:       makeButterworthLowShelv();                        break;
      case CHEBYCHEV:         makeChebychevLowShelv();                          break;
      case INVERSE_CHEBYCHEV: makeInverseChebychevLowShelv();                   break;
      case ELLIPTIC:          makeEllipticLowShelv();                           break;  
      case BESSEL:            makeBesselLowShelv(  rsDB2amp(A), rsDB2amp(A0));  break;
      case PAPOULIS:          makePapoulisLowShelv(rsDB2amp(A), rsDB2amp(A0));  break;
      }
    }
    stateIsDirty = false;
  }
}

void rsPrototypeDesigner::makeBypass()
{
  numFinitePoles = 0;
  numFiniteZeros = 0;
}

void rsPrototypeDesigner::makeButterworthLowpass()
{
  numFinitePoles = N;
  numFiniteZeros = 0;

  // intermediate variables:
  rsComplexDbl j(0.0, 1.0);                    // imaginary unit
  double  Gp     = sqrt(0.5);                  // -3.01 dB point gives cutoff frequency
  //double  Gp   = pow(10.0, -Ap/20.0);        // (1),Eq.1 - more general, cutoff gain specified
  double  ep     = sqrt(1.0/(Gp*Gp)-1.0);      // (1),Eq.2
  double  ep_pow = pow(ep, -1.0/(double) N);

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -ep_pow;                                    // Eq.70
    z[L+r-1] = rsInfDouble;                                // zero at infinity
  }
  // calculate the complex conjugate poles and zeros:
  double  u_i;
  for(int i=0; i<L; i++)
  {
    u_i  = (double) (2*(i+1)-1) / (double) N;              // Eq.69
    p[i] = ep_pow*j*rsExpC(j*u_i*PI*0.5);                  // Eq.70
    z[i] = rsInfDouble;                                    // zeros are at infinity
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeButterworthLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // catch some special cases:
  if( A0 == -rsInfDouble ) // lowpass-case
  {
    makeButterworthLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // intermediate variables:
  double G0   = rsDB2amp(A0);
  double G    = rsDB2amp(A);
  double GB   = sqrt(G0*G);                                // (2),Eq.52
  double ep   = sqrt( (G*G-GB*GB) / (GB*GB-G0*G0) );       // (2),Eq.12
  double g0   = pow(G0, 1.0 / (double) N);                 // (2),Eq.94
  double g    = pow(G,  1.0 / (double) N);                 // (2),Eq.94
  double wb   = 1.0;                                       // unit cutoff prototype
  double beta = wb * pow(ep, -1.0 / (double) N);           // (2),Eq.94

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -beta;                                      // (2),Eq.93
    z[L+r-1] = -g*beta/g0;                                 // (2),Eq.93
  }
  // calculate the complex conjugate poles and zeros:
  double phi, s, c;
  for(int i=0; i<L; i++)
  {
    phi     = (double) (2*(i+1)-1)*PI / (double) (2*N);    // (2),Eq.95
    s       = sin(phi);                                    // (2),Eq.95
    c       = cos(phi);                                    // (2),Eq.95
    z[i].re = -s*g*beta/g0;                                // (2),Eq.93
    z[i].im =  c*g*beta/g0;                                // (2),Eq.93
    p[i].re = -s*beta;                                     // (2),Eq.93
    p[i].im =  c*beta;                                     // (2),Eq.93
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeChebychevLowpass()
{
  numFinitePoles = N;
  numFiniteZeros = 0;

  // intermediate variables:
  double  Gp   = pow(10.0, -Ap/20.0);            // Eq. 1
  double  ep   = sqrt(1.0/(Gp*Gp) - 1.0);        // Eq. 2
  double  v_0  = rsAsinh(1.0/ep) / (N*PI*0.5);   // Eq. 72
  double  u_i;
  rsComplexDbl j(0.0, 1.0); // imaginary unit

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -sinh(v_0*PI*0.5);                    // Eq. 71
    z[L+r-1] = rsInfDouble;
  }
  // calculate the complex conjugate poles and zeros:
  for(int i=0; i<L; i++)
  {
    u_i  = (double) (2*(i+1)-1) / (double) N;        // Eq. 69
    p[i] = j*rsCosC((u_i-j*v_0)*PI*0.5);               // Eq. 71
    z[i] = rsInfDouble;                              // zeros at infinity
  }

  //gain = 1.0 / getFilterResponseAt(rsComplexDbl(0.0, 0.0)).getMagnitude();
  stateIsDirty = false;
}

void rsPrototypeDesigner::makeChebychevLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // calculate the linear gain-factors:
  double G0 = rsDB2amp(A0);
  double G  = rsDB2amp(A);

  // catch some special cases:
  if( A0 == -rsInfDouble ) // lowpass-case
  {
    makeChebychevLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // calculate intermediate variables:
  double Gp    = rsDB2amp(A0 + Rp*A);
  double ep    = sqrt( (G*G-Gp*Gp) / (Gp*Gp-G0*G0) );
  double g0    = pow(G0, 1.0 / (double) N);
  //double g     = pow(G,   1.0 / (double) N);
  double alpha = pow(1.0/ep + sqrt(1.0 + 1.0/(ep*ep)), 1.0/(double) N);
  double beta  = pow((G/ep + Gp*sqrt(1.0 + 1.0/(ep*ep)) ), 1.0/(double) N);
  double u     = log(beta/g0);
  double v     = log(alpha);
  double Gb    = sqrt(G0*G);
  double eb    = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );
  double wb    = 1.0 / cosh( rsAcosh(eb/ep) / (double)N ); // why 1/cosh(...) and not simply cosh?

  // calculate real pole and zero of the first order stage, if present and store them in the last array slots:
  if( r == 1 )
  {
    p[L+r-1] = -wb*sinh(v);
    z[L+r-1] = -wb*sinh(u);
  }

  // calculate the complex conjugate poles and zeros:
  double phi_i; //, s, c;
  rsComplexDbl j(0.0, 1.0); // imaginary unit
  for(int i=0; i<L; i++)
  {
    phi_i = (double) (2*(i+1)-1)*PI / (double) (2*N);
    z[i]  = j*wb*rsCosC(phi_i - j*u);
    p[i]  = j*wb*rsCosC(phi_i - j*v);
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeInverseChebychevLowpass()
{
  numFinitePoles = N;
  if( rsIsEven(N) )
    numFiniteZeros = N;
  else
    numFiniteZeros = N-1;

  // declare/assign/calculate some repeatedly needed variables:
  double  Gs = pow(10.0, -As/20.0);                      // Eq. 1
  double  es = sqrt(1.0/(Gs*Gs)-1.0);                    // Eq. 2
  double  v0 = rsAsinh(es) / (N*PI*0.5);                   // Eq. 74
  rsComplexDbl j(0.0, 1.0);                                   // imaginary unit

  double  wb = 1.0; // ...leads to a gain of Gs (stopband-gain) at unity (w=1), we rescale it here so as to have the -3 dB point at w=1:
  double  Gp = sqrt(0.5);
  double  ep = sqrt(1.0/(Gp*Gp)-1.0);
  wb         = cosh( rsAcosh(es/ep) / N );                      // (1),Eq.9

  // calculate the position of the real pole (if present):
  double  ui;
  if( r == 1 )
  {
    p[L+r-1] = -wb / sinh(v0*PI*0.5);                           // Eq.73 with k=1
    z[L+r-1] = rsInfDouble;
  }

  // calculate the complex conjugate poles and zeros:
  for(int i=0; i<L; i++)
  {
    ui   = (double) (2*(i+1)-1) / (double) N;                   // Eq.69
    z[i] = wb / (j*rsCosC(rsComplexDbl(ui*PI/2)));              // Eq.73 with k=1
    p[i] = wb / (j*rsCosC(rsComplexDbl((ui-j*v0)*PI*0.5)));     // Eq.73 with k=1
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeInverseChebychevLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // calculate the linear gain-factors:
  double G0 = rsDB2amp(A0);
  double G  = rsDB2amp(A);

  // catch some special cases:
  if( A0 == -rsInfDouble ) // lowpass-case
  {
    makeInverseChebychevLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // calculate intermediate variables (\todo check if the gains have reasonable values):
  //double Gs    = dB2amp(Rs*G + (1.0-Rs)*G0);
  double Gs    = rsDB2amp(A0 + Rs*A);
  double es    = sqrt( (G*G-Gs*Gs) / (Gs*Gs-G0*G0) );
  //double g0    = pow(G0, 1.0 / (double) N);
  double g     = pow(G,   1.0 / (double) N);
  double alpha = pow(es + sqrt(1.0 + es*es), 1.0 / (double) N);
  double beta  = pow((G0*es + Gs*sqrt(1.0 + es*es) ), 1.0/(double) N);
  double u     = log(beta/g);
  double v     = log(alpha);
  double Gb    = sqrt(G0*G);
  double eb    = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );
  double wb    = cosh( rsAcosh(es/eb) / (double) N );  // why not 1 / cosh(..)?

  // calculate real pole and zero of the first order stage, if present and store them in the last array slots:
  if( r == 1 )
  {
    z[L+r-1] = -wb/sinh(u);
    p[L+r-1] = -wb/sinh(v);
  }

  // calculate the complex conjugate poles and zeros:
  double  phi_i;
  rsComplexDbl j(0.0, 1.0); // imaginary unit
  for(int i=0; i<L; i++)
  {
    phi_i = (double) (2*(i+1)-1)*PI / (double) (2*N);
    z[i]  = wb / (j*rsCosC(phi_i-j*u));
    p[i]  = wb / (j*rsCosC(phi_i-j*v));
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeEllipticLowpass()
{
  numFinitePoles = N;
  if( rsIsEven(N) )
    numFiniteZeros = N;
  else
    numFiniteZeros = N-1;

  // declare/assign/calculate some repeatedly needed variables:
  rsComplexDbl j(0.0, 1.0);                                    // imaginary unit
  double  u_i;
  rsComplexDbl zeta_i;
  double  Gp  = pow(10.0, -Ap/20.0);                      // Eq. 1
  double  Gs  = pow(10.0, -As/20.0);                      // Eq. 1
  double  ep  = sqrt(1.0/(Gp*Gp) - 1.0);                  // Eq. 2
  double  es  = sqrt(1.0/(Gs*Gs) - 1.0);                  // Eq. 2
  double  k1  = ep/es;                                    // Eq. 3
  double  k   = ellipdeg(N, k1);                          // solve degree equation for k
  double  v_0 =  (-j*rsAsnC(j/ep, k1) / (double) N).re;     // from ellipap.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    //p[L+r-1] = -Omega_p/sinh(v_0*PI*0.5*k);                     // Eq. 73
    p[L+r-1] = j * rsSnC(j*v_0, k);                               // from ellipap.m
    z[L+r-1] = rsInfDouble;
  }
  // calculate the complex conjugate poles and zeros:
  for(int i=0; i<L; i++)
  {
    u_i    = (double) (2*(i+1)-1) / (double) N;                  // Eq. 69
    zeta_i = rsCdC(rsComplexDbl(u_i), k);                        // from ellipap.m
    z[i]   = j / (k*zeta_i);                                     // Eq. 62
    p[i]   = j*rsCdC((u_i-j*v_0), k);
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeEllipticLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // catch some special cases:
  if( A0 == -rsInfDouble ) // lowpass-case
  {
    makeEllipticLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // intermediate variables:
  double  G0  = rsDB2amp(A0);                          // reference amplitude
  double  G   = rsDB2amp(A);                           // boost/cut amplitude
  double  Gp  = rsDB2amp(A0 + Rp*A);                   // passband-amplitude (Rp near 1)
  double  Gs  = rsDB2amp(A0 + Rs*A);                   // stopband-amplitude (Rs near 0)
  double  Gb  = sqrt(G0*G);                            // (2),Eq.52 (gain at the bandedges)
  double  ep  = sqrt( (G*G-Gp*Gp) / (Gp*Gp-G0*G0) );   // (2),Eq.12
  double  es  = sqrt( (G*G-Gs*Gs) / (Gs*Gs-G0*G0) );   // (2),Eq.39
  double  eb  = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );   // (2),Eq.64
  double  k1  = ep/es;                                 // (2),Eq.39
  double  k   = ellipdeg(N, k1);                       // degree equation
  rsComplexDbl u = rsAcdC(rsComplexDbl(eb/ep), k1) / double(N);  // following text after (2),Eq.65
  double  wb  = 1.0 / rsCdC(u, k).re;                            // ...ditto
  rsComplexDbl j   = rsComplexDbl(0.0, 1.0);                     // imaginary unit
  rsComplexDbl ju0 = rsAsnC(j*G/(ep*G0), k1) / double(N);        // line 111 in hpeq.m
  rsComplexDbl jv0 = rsAsnC(j  / ep,     k1) / double(N);        // line 113 in hpeq.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = wb*(j*rsCdC(-1.0+jv0,k)).re;              // line 148 in hpeq.m
    z[L+r-1] = wb*(j*rsCdC(-1.0+ju0,k)).re;              // line 145 in hpeq.m
  }

  // calculate the complex conjugate poles and zeros:
  double ui;
  for(int i=0; i<L; i++)
  {
    ui   = (double) (2*(i+1)-1) / (double) N;          // (2),Eq.37
    p[i] = j*wb * rsCdC( (ui-jv0), k );                // line 179 in hpeq.m

    if( G0 == 0.0 && G != 0.0 )                        // lines 172-178 in hpeq.m
      z[i] = j*wb / (k*rsCdC(rsComplexDbl(ui),k));     // lowpass
    else if( G0 != 0.0 && G == 0.0 )
      z[i] = j*wb * rsCdC(rsComplexDbl(ui),k);         // highpass
    else
      z[i] = j*wb * rsCdC(ui-ju0,k);                   // low-shelv
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeBesselLowShelv(double G, double G0)
{
  rsFillWithZeros(p, maxNumNonRedundantPoles);
  rsFillWithZeros(z, maxNumNonRedundantPoles);
  numFinitePoles = N;
  if( G0 == 0.0 )
    numFiniteZeros = 0;
  else
    numFiniteZeros = N;

  rsComplexDbl zTmp[25];
  rsComplexDbl pTmp[25];
  double  kTmp;
  rsPrototypeDesigner::getBesselLowShelfZerosPolesAndGain(zTmp, pTmp, &kTmp, N, G, G0);

  // findPolynomialRoots returns the roots sorted by ascending real part. for a Bessel-polynomial, this ensures that the real pole, if 
  // present, is in pTmp[0] (it has the largest negative real part). this is importatnt for the next call:

  pickNonRedundantPolesAndZeros(zTmp, pTmp);
  stateIsDirty = false;
}


void rsPrototypeDesigner::makePapoulisLowShelv(double G, double G0)
{
  rsFillWithZeros(p, maxNumNonRedundantPoles);
  rsFillWithZeros(z, maxNumNonRedundantPoles);
  numFinitePoles = N;
  if( G0 == 0.0 )
    numFiniteZeros = 0;
  else
    numFiniteZeros = N;

  rsComplexDbl zTmp[25];
  rsComplexDbl pTmp[25];
  double  kTmp;
  rsPrototypeDesigner::getPapoulisLowShelfZerosPolesAndGain(zTmp, pTmp, &kTmp, N, G, G0);

  // findPolynomialRoots returns the roots sorted by ascending real part. for a Bessel-polynomial, this ensures that the real pole, if 
  // present, is in pTmp[0] (it has the largest negative real part). this is importatnt for the next call:

  pickNonRedundantPolesAndZeros(zTmp, pTmp);
  stateIsDirty = false;
}


void rsPrototypeDesigner::pickNonRedundantPolesAndZeros(rsComplexDbl *zTmp, rsComplexDbl *pTmp)
{
  rsZeroNegligibleImaginaryParts(pTmp, N, 1.e-11);
  rsZeroNegligibleImaginaryParts(zTmp, N, 1.e-11);
  rsOnlyUpperHalfPlane(pTmp, pTmp, N);
  rsOnlyUpperHalfPlane(zTmp, zTmp, N);
  rsCopyBuffer(pTmp, p, L+r);
  rsCopyBuffer(zTmp, z, L+r);

  // the caller is supposed to ensure that the real zero/pole, if present, is in zTmp[0], 
  // pTmp[0] - but we need it in the last positions z[L+r], p[L+r], so we reverse the arrays:
  rsReverse(p, L+r);
  rsReverse(z, L+r);
}
