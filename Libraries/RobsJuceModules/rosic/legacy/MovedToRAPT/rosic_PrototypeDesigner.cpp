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
  Ap                    = RAPT::rsAmpToDb(sqrt(2.0)); // 3.01 dB passband ripple for lowpasses
  As                    = 60.0;              // 60.0 dB stopband attenuation for lowpasses
  A                     = 0.0;               // cut/boost in dB for shelvers
  A0                    = 0.0;               // reference gain in dB for shelvers
  Rp                    = 0.95;              // inner ripple as fraction of dB-peak-gain for shelv
  Rs                    = 0.05;              // outer ripple as fraction of peak
  stateIsDirty          = true;              // poles and zeros need to be evaluated
  updatePolesAndZeros();
}

rsPrototypeDesigner::~rsPrototypeDesigner()
{

}

// parameter settings:

void rsPrototypeDesigner::setOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder != N )
  {
    N = newOrder;
    if( RAPT::rsIsOdd(N) )
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
    DEBUG_BREAK; // this is not one of the enumerated approximation methods

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
    DEBUG_BREAK; // this is not one of the enumerated modes
}

void rsPrototypeDesigner::setPassbandRipple(double newPassbandRipple)
{
  if( newPassbandRipple >= 0.0 )
  {
    Ap           = newPassbandRipple;
    stateIsDirty = true;
  }
  else
    DEBUG_BREAK; // ripple (in dB) must be >= 0
}

void rsPrototypeDesigner::setStopbandRejection(double newStopbandRejection)
{
  if( newStopbandRejection >= 0.0 )
  {
    As           = newStopbandRejection;
    stateIsDirty = true;
  }
  else
    DEBUG_BREAK; // ripple (in dB) must be >= 0
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
    DEBUG_BREAK; // this bandwidth gain ratio makes no sense (inequation 51 is violated)
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
    DEBUG_BREAK; // this stopband gain ratio makes no sense (inequation 51 is violated)
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
  if( RAPT::rsIsOdd(N) )
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
  if( RAPT::rsIsEven(N) )
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
    for(int i = 1; i <= L; i++)
    {
      u_i   = (double) (2*i-1) / (double) N;
      prod *= snC(u_i, kc).re;
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
  if( RAPT::rsIsEven(N) )
    L = N/2;
  else
    L = (N-1)/2;

  double u_i;
  double prod = 1.0;
  for(int i = 1; i <= L; i++)
  {
    u_i   = (double) (2*i-1) / (double) N;
    prod *= snC(u_i, k).re;
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
  ellipticIntegral(k, &K, &Kprime);

  double q  = exp(-PI*Kprime/K);
  double q1 = pow(q, N);

  int m;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for(m = 1; m <= M; m++)
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
  double Gp = pow(10.0, -passbandRipple/20.0);                       // (1),Eq.1
  double Gs = pow(10.0, -stopbandRipple/20.0);                       // (1),Eq.1
  double ep = sqrt(1.0 / (Gp*Gp) - 1.0);                             // (1),Eq.2
  double es = sqrt(1.0 / (Gs*Gs) - 1.0);                             // (1),Eq.2
  return acosh(es/ep) / acosh(stopbandFrequency/passbandFrequency);  // (1),Eq.9
}

double rsPrototypeDesigner::getRequiredEllipticOrder(double passbandFrequency, 
  double passbandRipple, double stopbandFrequency, double stopbandRipple)
{
  double Gp = pow(10.0, -passbandRipple/20.0);                       // (1),Eq.1
  double Gs = pow(10.0, -stopbandRipple/20.0);                       // (1),Eq.1
  double ep = sqrt(1.0 / (Gp*Gp) - 1.0);                             // (1),Eq.2
  double es = sqrt(1.0 / (Gs*Gs) - 1.0);                             // (1),Eq.2
  double k  = passbandFrequency / stopbandFrequency;                 // (1),Eq.3
  double k1 = ep/es;                                                 // (1),Eq.3
  double  K, Kp, K1, K1p;
  ellipticIntegral(k,  &K,  &Kp);                                    // (1),Eq.19
  ellipticIntegral(k1, &K1, &K1p);
  return (K1p*K)/(K1*Kp);                                            // (1),Eq.34
}

void rsPrototypeDesigner::magSquaredNumAndDen(double *b, double *a, double *b2, double *a2, int N)
{
  double *am = new double[N+1];
  double *bm = new double[N+1];
  RAPT::rsPolynomial<double>::polyCoeffsForNegativeArgument(b, bm, N);  // coeffs of N(-s)
  RAPT::rsPolynomial<double>::polyCoeffsForNegativeArgument(a, am, N);  // coeffs of D(-s)
  RAPT::rsPolynomial<double>::multiplyPolynomials(b, N, bm, N, b2);     // coeffs of N(s)*N(-s)
  RAPT::rsPolynomial<double>::multiplyPolynomials(a, N, am, N, a2);     // coeffs of D(s)*D(-s)
  delete[] am;
  delete[] bm;
}

void rsPrototypeDesigner::shelvingMagSqrNumFromLowpassMagSqr(double *b2, double *a2, double k, 
  int N, double G0, double G, double *bShelf)
{
  RAPT::rsArrayTools::weightedSum(b2, a2, bShelf, 2*N+1, k*k*(G*G-G0*G0), G0*G0);
}

// factor out shelvingMagSqrNumeratorFromLowpassMagSqr:
void rsPrototypeDesigner::shelvingMagSqrNumeratorFromLowpassTransfer(double *b, double *a, 
  double k, int N, double G0, double G, double *bShelf)
{
  double *a2 = new double[2*N+1];
  double *b2 = new double[2*N+1];

  // construct lowpass magnitude squared numerator and denominator 
  // N_LP(s)*N_LP(-s), D_LP(s)*D_LP(-s):
  magSquaredNumAndDen(b, a, b2, a2, N);

  // obtain coefficients for shelving filter's magnitude squared function numerator polynomial 
  // N_LS(s)*N_LS(-s):
  shelvingMagSqrNumFromLowpassMagSqr(b2, a2, k, N, G0, G, bShelf);

  delete[] a2;
  delete[] b2;
}

void rsPrototypeDesigner::scaleToMatchGainAtUnity(Complex *z, Complex *p, double *k, Complex *zNew, 
  Complex *pNew, double *kNew, int N, double g)
{
  double  wc    = rsFilterAnalyzer::findAnalogFrequencyWithMagnitude(z, p, k, N, g, 1.0);
  double scaler = 1.0/wc;
  for(int n = 0; n < N; n++)
  {
    pNew[n] = scaler * p[n];
    zNew[n] = scaler * z[n];
  }
  int nz = getNumFiniteValues(z, N);
  *kNew  = *k / pow(wc, N-nz);
}

void rsPrototypeDesigner::getInverseFilter(Complex* z, Complex* p, double* k, 
  Complex* zNew, Complex* pNew, double* kNew, int N)
{
  //rassert(false); // something seems wrong about this - we should write the inverted poles, zeros
  //                // and gain zeros into zNew, pNew, kNew

  Complex *zTmp = new Complex[N]; // to make it work, when the new arrays are equal to the old ones
  RAPT::rsArrayTools::copy(z,    zTmp, N);
  //copy(p,    z,    N);  // yes, this was wrong
  //copy(zTmp, p,    N);
  RAPT::rsArrayTools::copy(p,    zNew,    N);
  RAPT::rsArrayTools::copy(zTmp, pNew,    N);
  *kNew = 1.0 / *k;
  delete[] zTmp;
}

int rsPrototypeDesigner::getLeftHalfPlaneRoots(double *a, Complex *r, int N)
{
  Complex *rTmp = new Complex[N]; // maybe we can get rid of that temporary array
  findPolynomialRoots(a, N, rTmp);
  int numLeftRoots = onlyLeftHalfPlane(rTmp, r, N);

  rassert(numLeftRoots == ceil(0.5*N)); // maybe take this out later
  delete[] rTmp;
  return numLeftRoots;
}


void rsPrototypeDesigner::getBesselLowpassZerosPolesAndGain(Complex *z, Complex *p, double *k, 
  int N)
{
  // zeros are at infinity:
  RAPT::rsArrayTools::fillWithValue(z, N, Complex(INF, 0.0));

  // find poles:
  double *a = new double[N+1];        // Bessel-Polynomial coefficients
  //double a[20];
  RAPT::rsPolynomial<double>::besselPolynomial(a, N);
  RAPT::rsArrayTools::reverse(a, N+1);     // we actually use a reverse Bessel polynomial

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

void rsPrototypeDesigner::getBesselLowShelfZerosPolesAndGain(Complex *z, Complex *p, double *k, 
  int N, double G, double G0)
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
  RAPT::rsPolynomial<double>::besselPolynomial(a, N);
  RAPT::rsArrayTools::reverse(a, N+1);   // leaving this out leads to a modified Bessel filter response 
                                    // maybe experiment a bit, response looks good

  // find poles of the shelving filter:
  findPolynomialRoots(a, N, p);

  // construct lowpass numerator:
  double *b = new double[N+1];
  RAPT::rsArrayTools::fillWithZeros(b, N+1);
  b[0] = a[0];

  // obtain magnitude-squared numerator polynomial for shelving filter:
  double *bS = new double[2*N+1];
  shelvingMagSqrNumeratorFromLowpassTransfer(b, a, 1.0, N, G0, G, bS);

  // find left halfplane zeros (= zeros of the shelving filter):
  getLeftHalfPlaneRoots(bS, z, 2*N);

  // set gain constant:
  *k = G0;

  // now we have a shelving filter with correct low-frequency gain G and reference gain G0, but 
  // possibly still with wrong bandwidth gain GB at unity - now we adjust zeros/poles/gain to 
  // match GB:
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
  RAPT::rsArrayTools::fillWithZeros(a, 2*N+1);  // do we need this?

  // construct the polynomial L_N(w^2):
  RAPT::rsPolynomial<double>::maximumSlopeMonotonicPolynomial(a, N);  // does the same same as lopt(a, N); from C.R.Bond

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

void rsPrototypeDesigner::getPapoulisLowpassZerosPolesAndGain(Complex *z, Complex *p, double *k, 
  int N)
{
  // find poles:
  double *a2 = new double[2*N+1];     // coefficients of magnitude-squared polynomial D(s)*D(-s)
  papoulisMagnitudeSquaredDenominator(a2, N);
  getLeftHalfPlaneRoots(a2, p, 2*N);

  // zeros are at infinity:
  RAPT::rsArrayTools::fillWithValue(z, N, Complex(INF, 0.0));

  // set gain at DC to unity:
  *k = sqrt(1.0/fabs(a2[2*N]));

  delete[] a2;
}

void rsPrototypeDesigner::getPapoulisLowShelfZerosPolesAndGain(Complex *z, Complex *p, double *k, 
  int N, double G, double G0)
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

  // factor out into a function getMagnitudeSquaredNumAndDen, then call this function here - this 
  // function can switch between Papoulis, Gauss, etc. and fill the arrays accordingly, Bessel is 
  // a special case - it doesn't need to find poles of the mag-squared function, there we can 
  // directly find the poles of the transfer function ...hmm...but maybe we can still factor out a 
  // function getPoles or something...

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
  RAPT::rsArrayTools::fillWithZeros(b2, 2*N+1);
  b2[0] = 1.0;

  // end of "factor out" ...in general, we need to scale the b2-polynomial also by dividing through 
  // the leading coeff?

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

void rsPrototypeDesigner::getEllipticLowpassZerosPolesAndGain(Complex *z, Complex *p, double *k, 
  int N, double Gp, double Gs)
{
//  int nz;
//  if( isEven(N) )
//    nz = N;
//  else
//    nz = N-1;

  int L, r;
  getNumBiquadsAndFirstOrderStages(N, L, r);

  // declare/assign/calculate some repeatedly needed variables:
  Complex j(0.0, 1.0);                                    // imaginary unit
  double  ep  = sqrt(1.0/(Gp*Gp) - 1.0);                  // Eq. 2
  double  es  = sqrt(1.0/(Gs*Gs) - 1.0);                  // Eq. 2
  double  k1  = ep/es;                                    // Eq. 3
  double  kk  = ellipdeg(N, k1);                          // solve degree equation for k
  double  v_0 =  (-j*asnC(j/ep, k1) / (double) N).re;     // from ellipap.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    //p[L+r-1] = -Omega_p/sinh(v_0*PI*0.5*kk);                 // Eq. 73
    p[N-1] = j * snC(j*v_0, kk);                               // from ellipap.m - find Eq.
    z[N-1] = INF;
  }

  // calculate the complex conjugate poles and zeros:
  double  u_i;
  Complex zeta_i;
  for(int i = 0; i < L; i++)
  {
    u_i      = (double) (2*(i+1)-1) / (double) N;                // Eq. 69
    zeta_i   = cdC(u_i, kk);                                     // from ellipap.m - find Eq.
    z[2*i]   = j / (kk*zeta_i);                                  // Eq. 62
    p[2*i]   = j*cdC((u_i-j*v_0), kk);
    z[2*i+1] = z[2*i].getConjugate();
    p[2*i+1] = p[2*i].getConjugate();
  }

  double H0 = pow(Gp, 1-r);    // preliminary - can be made simpler (without pow)
  Complex n = productOfFiniteFactors(p, N);
  Complex d = productOfFiniteFactors(z, N);
  *k        = H0 * (n/d).re;
}

/*
void rsPrototypeDesigner::getLowpassZerosPolesAndGain(Complex *z, Complex *p, double *k, int N, 
int approximationMethod)
{
  int dummy = 0;
}
*/

// inquiry:

Complex rsPrototypeDesigner::getFilterResponseAt(Complex s)
{
  Complex num, den;
  Complex tmp;
  int     Lz, Lp;

  // initialize the numerator and denominator:
  if( RAPT::rsIsOdd(numFiniteZeros) ) {
    num = -z[L+r-1].re;
    Lz  = (numFiniteZeros-1)/2;  }
  else {
    num = 1.0;
    Lz  = numFiniteZeros/2; }
  if( RAPT::rsIsOdd(numFinitePoles) ) {
    den = -p[L+r-1].re;
    Lp  = (numFinitePoles-1)/2; }
  else {
    den = 1.0;
    Lp  = numFinitePoles/2; }

  // accumulate product of linear factors for denominator (poles) and numerator (zeros):
  int i;
  for(i = 0; i < Lz; i++) num *= ((s-z[i]) * (s-z[i].getConjugate()));
  for(i = 0; i < Lp; i++) den *= ((s-p[i]) * (s-p[i].getConjugate()));

  return num/den;
}

double rsPrototypeDesigner::getMagnitudeAt(double w)
{
  return getFilterResponseAt(Complex(0.0, w)).getRadius();
}

double rsPrototypeDesigner::findFrequencyWithMagnitude(double magnitude, double wLow, double wHigh)
{
  // until we have something better, we search for the frequency at which the desired gain occurs 
  // by means of the bisection method:

  double tolerance = 0.0001; // maybe make parameter
  double wMid, mMid;
  while( wHigh-wLow > tolerance )
  {
    wMid  = 0.5 * (wLow+wHigh);
    mMid  = getMagnitudeAt(wMid);
    if( mMid > magnitude )
      wLow = wMid;
    else
      wHigh = wMid;
  }
  return 0.5 * (wLow+wHigh);
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
  if( RAPT::rsIsEven(numFinitePoles) )
    return numFinitePoles/2;
  else
    return (numFinitePoles+1)/2;
}

int rsPrototypeDesigner::getNumNonRedundantFiniteZeros()
{
  if( RAPT::rsIsEven(numFiniteZeros) )
    return numFiniteZeros/2;
  else
    return (numFiniteZeros+1)/2;
}

void rsPrototypeDesigner::getPolesAndZeros(Complex *poles, Complex *zeros)
{
  if( stateIsDirty == true )
    updatePolesAndZeros(); // re-calculate only if necesarry
  for(int i = 0; i < (L+r); i++)
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
      case BUTTERWORTH:       makeButterworthLowShelv();                    break;
      case CHEBYCHEV:         makeChebychevLowShelv();                      break;
      case INVERSE_CHEBYCHEV: makeInverseChebychevLowShelv();               break;
      case ELLIPTIC:          makeEllipticLowShelv();                       break;
      case BESSEL:            makeBesselLowShelv(  RAPT::rsDbToAmp(A), RAPT::rsDbToAmp(A0));  break;
      case PAPOULIS:          makePapoulisLowShelv(RAPT::rsDbToAmp(A), RAPT::rsDbToAmp(A0));  break;
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
  Complex j(0.0, 1.0);                        // imaginary unit
  double  Gp     = sqrt(0.5);                 // use -3.01 dB point as cutoff frequency for Butterworths
  //double  Gp   = pow(10.0, -Ap/20.0);       // (1),Eq.1 - more general (cutoff gain can be specified), not used here
  double  ep     = sqrt(1.0/(Gp*Gp)-1.0);     // (1),Eq.2
  double  ep_pow = pow(ep, -1.0/(double) N);

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -ep_pow;                                    // Eq.70
    z[L+r-1] = INF;                                        // zero at infinity
  }
  // calculate the complex conjugate poles and zeros:
  double  u_i;
  for(int i = 0; i < L; i++)
  {
    u_i  = (double) (2*(i+1)-1) / (double) N;              // Eq.69
    p[i] = ep_pow*j*expC(j*u_i*PI*0.5);                    // Eq.70
    z[i] = INF;                                            // zeros are at infinity
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeButterworthLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // catch some special cases:
  if( A0 == NEG_INF ) // lowpass-case
  {
    makeButterworthLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // intermediate variables:
  double G0   = RAPT::rsDbToAmp(A0);
  double G    = RAPT::rsDbToAmp(A);
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
  for(int i = 0; i < L; i++)
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
  double  v_0  = asinh(1.0/ep) / (N*PI*0.5);     // Eq. 72
  double  u_i;
  Complex j(0.0, 1.0); // imaginary unit

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -sinh(v_0*PI*0.5);                    // Eq. 71
    z[L+r-1] = INF;
  }
  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    u_i  = (double) (2*(i+1)-1) / (double) N;        // Eq. 69
    p[i] = j*cosC((u_i-j*v_0)*PI*0.5);               // Eq. 71
    z[i] = INF;                                      // zeros at infinity
  }

  //gain = 1.0 / getFilterResponseAt(Complex(0.0, 0.0)).getMagnitude();
  stateIsDirty = false;
}

void rsPrototypeDesigner::makeChebychevLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // calculate the linear gain-factors:
  double G0 = RAPT::rsDbToAmp(A0);
  double G  = RAPT::rsDbToAmp(A);

  // catch some special cases:
  if( A0 == NEG_INF ) // lowpass-case
  {
    makeChebychevLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // calculate intermediate variables:
  double Gp    = RAPT::rsDbToAmp(A0 + Rp*A);
  double ep    = sqrt( (G*G-Gp*Gp) / (Gp*Gp-G0*G0) );
  double g0    = pow(G0, 1.0 / (double) N);
  //double g     = pow(G,   1.0 / (double) N);
  double alpha = pow(1.0/ep + sqrt(1.0 + 1.0/(ep*ep)), 1.0/(double) N);
  double beta  = pow((G/ep + Gp*sqrt(1.0 + 1.0/(ep*ep)) ), 1.0/(double) N);
  double u     = log(beta/g0);
  double v     = log(alpha);
  double Gb    = sqrt(G0*G);
  double eb    = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );
  double wb    = 1.0 / cosh( acosh(eb/ep) / (double)N ); // why 1/cosh(...) and not simply cosh?

  // calculate real pole and zero of the first order stage, if present and store them in the last array slots:
  if( r == 1 )
  {
    p[L+r-1] = -wb*sinh(v);
    z[L+r-1] = -wb*sinh(u);
  }

  // calculate the complex conjugate poles and zeros:
  double phi_i; //, s, c;
  Complex j(0.0, 1.0); // imaginary unit
  for(int i = 0; i < L; i++)
  {
    phi_i = (double) (2*(i+1)-1)*PI / (double) (2*N);
    z[i]  = j*wb*cosC(phi_i - j*u);
    p[i]  = j*wb*cosC(phi_i - j*v);
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeInverseChebychevLowpass()
{
  numFinitePoles = N;
  if( RAPT::rsIsEven(N) )
    numFiniteZeros = N;
  else
    numFiniteZeros = N-1;

  // declare/assign/calculate some repeatedly needed variables:
  double  Gs = pow(10.0, -As/20.0);                      // Eq. 1
  double  es = sqrt(1.0/(Gs*Gs)-1.0);                    // Eq. 2
  double  v0 = asinh(es) / (N*PI*0.5);                   // Eq. 74
  Complex j(0.0, 1.0);                                   // imaginary unit

  double  wb = 1.0; // ...leads to a gain of Gs (stopband-gain) at unity (w=1), we rescale it here
                    // so as to have the -3 dB point at w=1:
  double  Gp = sqrt(0.5);
  double  ep = sqrt(1.0/(Gp*Gp)-1.0);
  wb         = cosh( acosh(es/ep) / N );                        // (1),Eq.9

  // calculate the position of the real pole (if present):
  double  ui;
  if( r == 1 )
  {
    p[L+r-1] = -wb / sinh(v0*PI*0.5);                           // Eq.73 with k=1
    z[L+r-1] = INF;
  }

  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    ui   = (double) (2*(i+1)-1) / (double) N;                   // Eq.69
    z[i] = wb / (j*cosC(ui*PI/2));                              // Eq.73 with k=1
    p[i] = wb / (j*cosC((ui - j*v0)*PI/2));                     // Eq.73 with k=1
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeInverseChebychevLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // calculate the linear gain-factors:
  double G0 = RAPT::rsDbToAmp(A0);
  double G  = RAPT::rsDbToAmp(A);

  // catch some special cases:
  if( A0 == NEG_INF ) // lowpass-case
  {
    makeInverseChebychevLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // calculate intermediate variables (\todo check if the gains have reasonable values):
  //double Gs    = dB2amp(Rs*G + (1.0-Rs)*G0);
  double Gs    = RAPT::rsDbToAmp(A0 + Rs*A);
  double es    = sqrt( (G*G-Gs*Gs) / (Gs*Gs-G0*G0) );
  //double g0    = pow(G0, 1.0 / (double) N);
  double g     = pow(G,   1.0 / (double) N);
  double alpha = pow(es + sqrt(1.0 + es*es), 1.0 / (double) N);
  double beta  = pow((G0*es + Gs*sqrt(1.0 + es*es) ), 1.0/(double) N);
  double u     = log(beta/g);
  double v     = log(alpha);
  double Gb    = sqrt(G0*G);
  double eb    = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );
  double wb    = cosh( acosh(es/eb) / (double) N );  // why not 1 / cosh(..)?

  // calculate real pole and zero of the first order stage, if present and store them in the last
  // array slots:
  if( r == 1 )
  {
    z[L+r-1] = -wb/sinh(u);
    p[L+r-1] = -wb/sinh(v);
  }

  // calculate the complex conjugate poles and zeros:
  double  phi_i;
  Complex j(0.0, 1.0); // imaginary unit
  for(int i = 0; i < L; i++)
  {
    phi_i = (double) (2*(i+1)-1)*PI / (double) (2*N);
    z[i]  = wb / (j*cosC(phi_i-j*u));
    p[i]  = wb / (j*cosC(phi_i-j*v));
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeEllipticLowpass()
{
  numFinitePoles = N;
  if( RAPT::rsIsEven(N) )
    numFiniteZeros = N;
  else
    numFiniteZeros = N-1;

  // declare/assign/calculate some repeatedly needed variables:
  Complex j(0.0, 1.0);                                    // imaginary unit
  double  u_i;
  Complex zeta_i;
  double  Gp  = pow(10.0, -Ap/20.0);                      // Eq. 1
  double  Gs  = pow(10.0, -As/20.0);                      // Eq. 1
  double  ep  = sqrt(1.0/(Gp*Gp) - 1.0);                  // Eq. 2
  double  es  = sqrt(1.0/(Gs*Gs) - 1.0);                  // Eq. 2
  double  k1  = ep/es;                                    // Eq. 3
  double  k   = ellipdeg(N, k1);                          // solve degree equation for k
  double  v_0 =  (-j*asnC(j/ep, k1) / (double) N).re;     // from ellipap.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    //p[L+r-1] = -Omega_p/sinh(v_0*PI*0.5*k);                    // Eq. 73
    p[L+r-1] = j * snC(j*v_0, k);                                // from ellipap.m
    z[L+r-1] = INF;
  }
  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    u_i    = (double) (2*(i+1)-1) / (double) N;                  // Eq. 69
    zeta_i = cdC(u_i, k);                                        // from ellipap.m
    z[i]   = j / (k*zeta_i);                                     // Eq. 62
    p[i]   = j*cdC((u_i-j*v_0), k);
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeEllipticLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // catch some special cases:
  if( A0 == NEG_INF ) // lowpass-case
  {
    makeEllipticLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // intermediate variables:
  double  G0  = RAPT::rsDbToAmp(A0);                            // reference amplitude
  double  G   = RAPT::rsDbToAmp(A);                             // boost/cut amplitude
  double  Gp  = RAPT::rsDbToAmp(A0 + Rp*A);                     // passband-amplitude (Rp near 1)
  double  Gs  = RAPT::rsDbToAmp(A0 + Rs*A);                     // stopband-amplitude (Rs near 0)
  double  Gb  = sqrt(G0*G);                            // (2),Eq.52 (gain at the bandedges)
  double  ep  = sqrt( (G*G-Gp*Gp) / (Gp*Gp-G0*G0) );   // (2),Eq.12
  double  es  = sqrt( (G*G-Gs*Gs) / (Gs*Gs-G0*G0) );   // (2),Eq.39
  double  eb  = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );   // (2),Eq.64
  double  k1  = ep/es;                                 // (2),Eq.39
  double  k   = ellipdeg(N, k1);                       // degree equation
  Complex u   = acdC(eb/ep, k1) / N;                   // following text after (2),Eq.65
  double  wb  = 1.0 / cdC(u, k).re;                    // ...ditto
  Complex j   = Complex(0.0, 1.0);                     // imaginary unit
  Complex ju0 = asnC(j*G/(ep*G0), k1) / N;             // line 111 in hpeq.m
  Complex jv0 = asnC(j  / ep,     k1) / N;             // line 113 in hpeq.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = wb*(j*cdC(-1.0+jv0,k)).re;              // line 148 in hpeq.m
    z[L+r-1] = wb*(j*cdC(-1.0+ju0,k)).re;              // line 145 in hpeq.m
  }

  // calculate the complex conjugate poles and zeros:
  double ui;
  for(int i = 0; i < L; i++)
  {
    ui   = (double) (2*(i+1)-1) / (double) N;          // (2),Eq.37
    p[i] = j*wb * cdC( (ui-jv0), k );                  // line 179 in hpeq.m

    if( G0 == 0.0 && G != 0.0 )                        // lines 172-178 in hpeq.m
      z[i] = j*wb / (k*cdC(ui,k));   // lowpass
    else if( G0 != 0.0 && G == 0.0 )
      z[i] = j*wb * cdC(ui,k);       // highpass
    else
      z[i] = j*wb * cdC(ui-ju0,k);   // low-shelv
  }

  stateIsDirty = false;
}

void rsPrototypeDesigner::makeBesselLowShelv(double G, double G0)
{
  RAPT::rsArrayTools::fillWithZeros(p, maxNumNonRedundantPoles);
  RAPT::rsArrayTools::fillWithZeros(z, maxNumNonRedundantPoles);
  numFinitePoles = N;
  if( G0 == 0.0 )
    numFiniteZeros = 0;
  else
    numFiniteZeros = N;

  Complex zTmp[25];
  Complex pTmp[25];
  double  kTmp;
  rsPrototypeDesigner::getBesselLowShelfZerosPolesAndGain(zTmp, pTmp, &kTmp, N, G, G0);

  // findPolynomialRoots returns the roots sorted by ascending real part. for a Bessel-polynomial, 
  // this ensures that the real pole, if present, is in pTmp[0] (it has the largest negative real 
  // part). this is importatnt for the next call:

  pickNonRedundantPolesAndZeros(zTmp, pTmp);
  stateIsDirty = false;
}

void rsPrototypeDesigner::makePapoulisLowShelv(double G, double G0)
{
  RAPT::rsArrayTools::fillWithZeros(p, maxNumNonRedundantPoles);
  RAPT::rsArrayTools::fillWithZeros(z, maxNumNonRedundantPoles);
  numFinitePoles = N;
  if( G0 == 0.0 )
    numFiniteZeros = 0;
  else
    numFiniteZeros = N;

  Complex zTmp[25];
  Complex pTmp[25];
  double  kTmp;
  rsPrototypeDesigner::getPapoulisLowShelfZerosPolesAndGain(zTmp, pTmp, &kTmp, N, G, G0);

  // findPolynomialRoots returns the roots sorted by ascending real part. for a Bessel-polynomial, 
  // this ensures that the real pole, if present, is in pTmp[0] (it has the largest negative real 
  // part). this is importatnt for the next call:

  pickNonRedundantPolesAndZeros(zTmp, pTmp);
  stateIsDirty = false;
}


void rsPrototypeDesigner::pickNonRedundantPolesAndZeros(Complex *zTmp, Complex *pTmp)
{
  zeroNegligibleImaginaryParts(pTmp, N, 1.e-11);
  zeroNegligibleImaginaryParts(zTmp, N, 1.e-11);
  onlyUpperHalfPlane(pTmp, pTmp, N);
  onlyUpperHalfPlane(zTmp, zTmp, N);
  RAPT::rsArrayTools::copy(pTmp, p, L+r);
  RAPT::rsArrayTools::copy(zTmp, z, L+r);

  // the caller is supposed to ensure that the real zero/pole, if present, is in zTmp[0], pTmp[0] - 
  // but we need it in the last positions z[L+r], p[L+r], so we reverse the arrays:
  RAPT::rsArrayTools::reverse(p, L+r);
  RAPT::rsArrayTools::reverse(z, L+r);
}
