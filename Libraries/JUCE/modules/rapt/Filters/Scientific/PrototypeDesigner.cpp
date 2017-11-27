// helper functions (maybe make static member functions)

template<class T>
void updateLegendrePolynomial(T **P, T *P1, T *P2, int r)
{
  if( rsIsEven(r) )
  {
    rsPolynomial<T>::rsLegendrePolynomialRecursion(P1, r, P2, P1);
    *P = P1;
  }
  else
  {
    rsPolynomial<T>::rsLegendrePolynomialRecursion(P2, r, P1, P2);
    *P = P2;
  }
}

template<class T>
void rsHalpernU(T *a, int K)
{
  // Computes coefficients of the U-polynomials given in "Design and Analysis of Analog Filters", 
  // page 256-257 except for the scale factor in front.
  rsArray::fillWithZeros(a, K+1);
  rsUint64 k, m;
  if( rsIsEven(K) )
  {
    k = K/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*m] = (-1)^(k-m) * (m+k)! / ( (k-m)! * (m!)^2 ):
      a[2*m] = (T) (rsProduct(k-m+1, m+k) / rsSquare(rsFactorial(m)));
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
      a[2*(k-m)+1] = (T) (rsProduct(k-m+2, 2*k+1-m) / (rsFactorial(m)*rsFactorial(k-m)));
      if( rsIsOdd(m) )
        a[2*(k-m)+1] *= -1;
    }
  }
}



//=================================================================================================
// class rsPrototypeDesigner

// construction/destruction:

template<class T>
rsPrototypeDesigner<T>::rsPrototypeDesigner()
{
  L                   = 2;
  r                   = 0;
  N                   = L+r;
  approximationMethod = BUTTERWORTH;
  prototypeMode       = LOWPASS_PROTOTYPE;
  numFinitePoles      = 2;
  numFiniteZeros      = 0;
  Ap                  = (T)rsAmpToDb(sqrt(2)); // 3.01 dB passband ripple for lowpasses
  As                  = 60.0;     // 60.0 dB stopband attenuation for lowpasses
  A                   = 0.0;      // cut/boost in dB for shelvers
  A0                  = 0.0;      // reference gain in dB for shelvers
  Rp                  = (T)0.95;  // inner ripple as fraction of dB-peak-gain for shelv
  Rs                  = (T)0.05;  // outer ripple as fraction of peak
  stateIsDirty        = true;     // poles and zeros need to be evaluated
  updatePolesAndZeros();
}

template<class T>
rsPrototypeDesigner<T>::~rsPrototypeDesigner()
{

}

// parameter settings:

template<class T>
void rsPrototypeDesigner<T>::setOrder(int newOrder)
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

template<class T>
void rsPrototypeDesigner<T>::setApproximationMethod(int newApproximationMethod)
{
  if( newApproximationMethod < BUTTERWORTH || newApproximationMethod > PAPOULIS )
    RS_DEBUG_BREAK; // this is not one of the enumerated approximation methods

  if( newApproximationMethod != approximationMethod )
  {
    approximationMethod = newApproximationMethod;
    stateIsDirty        = true;
  }
}

template<class T>
void rsPrototypeDesigner<T>::setPrototypeMode(int newPrototypeMode)
{
  if( newPrototypeMode == LOWPASS_PROTOTYPE || newPrototypeMode == LOWSHELV_PROTOTYPE )
  {
    prototypeMode = newPrototypeMode;
    stateIsDirty  = true;
  }
  else
    RS_DEBUG_BREAK; // this is not one of the enumerated modes
}

template<class T>
void rsPrototypeDesigner<T>::setPassbandRipple(T newPassbandRipple)
{
  if( newPassbandRipple >= 0.0 )
  {
    Ap           = newPassbandRipple;
    stateIsDirty = true;
  }
  else
    RS_DEBUG_BREAK; // ripple (in dB) must be >= 0
}

template<class T>
void rsPrototypeDesigner<T>::setStopbandRejection(T newStopbandRejection)
{
  if( newStopbandRejection >= 0.0 )
  {
    As           = newStopbandRejection;
    stateIsDirty = true;
  }
  else
    RS_DEBUG_BREAK; // ripple (in dB) must be >= 0
}

template<class T>
void rsPrototypeDesigner<T>::setGain(T newGain)
{
  if( newGain != A )
  {
    A            = newGain;
    stateIsDirty = true;
  }
}

template<class T>
void rsPrototypeDesigner<T>::setReferenceGain(T newReferenceGain)
{
  if( newReferenceGain != A0 )
  {
    A0           = newReferenceGain;
    stateIsDirty = true;
  }
}

template<class T>
void rsPrototypeDesigner<T>::setPassbandGainRatio(T newPassbandGainRatio)
{
  if( newPassbandGainRatio >= 1.0 || newPassbandGainRatio <= 0.0 || newPassbandGainRatio < Rs )
  {
    RS_DEBUG_BREAK; // this bandwidth gain ratio makes no sense (inequation 51 is violated)
    return;
  }
  if( newPassbandGainRatio != Rp  )
  {
    Rp           = newPassbandGainRatio;
    stateIsDirty = true;
  }
}

template<class T>
void rsPrototypeDesigner<T>::setStopbandGainRatio(T newStopbandGainRatio)
{
  if( newStopbandGainRatio >= 1.0 || newStopbandGainRatio <= 0.0 || newStopbandGainRatio > Rp )
  {
    RS_DEBUG_BREAK; // this stopband gain ratio makes no sense (inequation 51 is violated)
    return;
  }
  if( newStopbandGainRatio != Rp  )
  {
    Rs           = newStopbandGainRatio;
    stateIsDirty = true;
  }
}

// static member functions:

template<class T>
void rsPrototypeDesigner<T>::getNumBiquadsAndFirstOrderStages(int N, int &L, int &r)
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

template<class T>
T rsPrototypeDesigner<T>::ellipdeg(int N, T k_1)
{
  int L;
  if( rsIsEven(N) )
    L = N/2;
  else
    L = (N-1)/2;

  T kmin = (T)1e-6;
  T k;

  if (k_1 < kmin)
    k = ellipdeg2(T(1) / T(N), k_1);
  else
  {
    T kc = sqrt(1-k_1*k_1);			  // complement of k1
    //T u_i;                      // old
    std::complex<T> u_i;          // use real argument later again
    T prod = 1.0;
    for(int i = 1; i <= L; i++)
    {
      u_i   = (T) (2*i-1) / (T) N;
      prod *= rsSnC(u_i, kc).real(); // use a real sn-function
    }
    prod = prod*prod*prod*prod;
    T kp = pow(kc, (T) N) * prod; // complement of k
    k    = sqrt(T(1)-kp*kp);
  }

  return k;
}

template<class T>
T rsPrototypeDesigner<T>::ellipdeg1(int N, T k)
{
  int L;
  if( rsIsEven(N) )
    L = N/2;
  else
    L = (N-1)/2;

  //T u_i;     
  std::complex<T> u_i;  // use real argument later again
  T prod = 1.0;
  for(int i = 1; i <= L; i++)
  {
    u_i   = (T) (2*i-1) / (T) N;
    prod *= rsSnC(u_i, k).real();
  }
  prod = prod*prod*prod*prod;
  T k1 = pow(k, (T) N) * prod;

  return k1;
}

template<class T>
T rsPrototypeDesigner<T>::ellipdeg2(T N, T k)
{
  int M = 7;

  T K;
  T Kprime;
  rsEllipticIntegral(k, &K, &Kprime);

  T q  = exp(-T(PI)*Kprime/K);
  T q1 = pow(q, N);

  int m;
  T sum1 = 0.0;
  T sum2 = 0.0;
  for(m = 1; m <= M; m++)
  {
    sum1 += pow(q1, (T) (m*(m+1)) );
    sum2 += pow(q1, (T) (m*m)     );
  }

  T tmp = (T(1)+sum1) / (T(1)+T(2)*sum2);
  tmp  *= tmp;
  T k1  = T(4) * sqrt(q1) * tmp;

  return k1;
}

template<class T>
T rsPrototypeDesigner<T>::getRequiredButterworthOrder(T passbandFrequency, T passbandRipple, 
  T stopbandFrequency, T stopbandRipple)
{
  T Gp = pow(T(10), -passbandRipple/T(20));                      // (1),Eq.1
  T Gs = pow(T(10), -stopbandRipple/T(20));                      // (1),Eq.1
  T ep = sqrt(T(1) / (Gp*Gp) - T(1));                            // (1),Eq.2
  T es = sqrt(T(1) / (Gs*Gs) - T(1));                            // (1),Eq.2
  return log(es/ep) / log(stopbandFrequency/passbandFrequency);  // (1),Eq.9
}

template<class T>
T rsPrototypeDesigner<T>::getRequiredChebychevOrder(T passbandFrequency, T passbandRipple, 
  T stopbandFrequency, T stopbandRipple)
{
  T Gp = pow(T(10), -passbandRipple/T(20));                         // (1),Eq.1
  T Gs = pow(T(10), -stopbandRipple/T(20));                         // (1),Eq.1
  T ep = sqrt(T(1) / (Gp*Gp) - T(1));                               // (1),Eq.2
  T es = sqrt(T(1) / (Gs*Gs) - T(1));                               // (1),Eq.2
  return acosh(es/ep) / acosh(stopbandFrequency/passbandFrequency); // (1),Eq.9
}

template<class T>
T rsPrototypeDesigner<T>::getRequiredEllipticOrder(T passbandFrequency, T passbandRipple, 
  T stopbandFrequency, T stopbandRipple)
{
  T Gp = pow(T(10), -passbandRipple/T(20));                     // (1),Eq.1
  T Gs = pow(T(10), -stopbandRipple/T(20));                     // (1),Eq.1
  T ep = sqrt(T(1) / (Gp*Gp) - T(1));                           // (1),Eq.2
  T es = sqrt(T(1) / (Gs*Gs) - T(1));                           // (1),Eq.2
  T k  = passbandFrequency / stopbandFrequency;                 // (1),Eq.3
  T k1 = ep/es;                                                 // (1),Eq.3
  T K, Kp, K1, K1p;
  rsEllipticIntegral(k,  &K,  &Kp);                             // (1),Eq.19
  rsEllipticIntegral(k1, &K1, &K1p);
  return (K1p*K)/(K1*Kp);                                       // (1),Eq.34
}

template<class T>
void rsPrototypeDesigner<T>::magSquaredNumAndDen(T* b, T* a, T* b2, T* a2, int N)
{
  T* am = new T[N+1];
  T* bm = new T[N+1];
  rsPolynomial<T>::polyCoeffsForNegativeArgument(b, bm, N);  // coeffs of N(-s)
  rsPolynomial<T>::polyCoeffsForNegativeArgument(a, am, N);  // coeffs of D(-s)
  rsPolynomial<T>::multiplyPolynomials(b, N, bm, N, b2);     // coeffs of N(s)*N(-s)
  rsPolynomial<T>::multiplyPolynomials(a, N, am, N, a2);     // coeffs of D(s)*D(-s)
  delete[] am;
  delete[] bm;
}

template<class T>
void rsPrototypeDesigner<T>::shelvingMagSqrNumFromLowpassMagSqr(T* b2, T* a2, T k, 
  int N, T G0, T G, T* bShelf)
{
  rsArray::weightedSum(b2, a2, bShelf, 2*N+1, k*k*(G*G-G0*G0), G0*G0);
}

// factor out shelvingMagSqrNumeratorFromLowpassMagSqr:
template<class T>
void rsPrototypeDesigner<T>::shelvingMagSqrNumeratorFromLowpassTransfer(T* b, T* a, T k, int N, 
  T G0, T G, T* bShelf)
{
  T* a2 = new T[2*N+1];
  T* b2 = new T[2*N+1];

  // construct lowpass magnitude squared numerator and denominator 
  // N_LP(s)*N_LP(-s), D_LP(s)*D_LP(-s):
  magSquaredNumAndDen(b, a, b2, a2, N);

  // obtain coefficients for shelving filter's magnitude squared function numerator polynomial 
  // N_LS(s)*N_LS(-s):
  shelvingMagSqrNumFromLowpassMagSqr(b2, a2, k, N, G0, G, bShelf);

  delete[] a2;
  delete[] b2;
}

template<class T>
void rsPrototypeDesigner<T>::scaleToMatchGainAtUnity(Complex* z, Complex* p, T* k, Complex* zNew, 
  Complex* pNew, T* kNew, int N, T g)
{
  T  wc    = rsFilterAnalyzer<T>::findAnalogFrequencyWithMagnitude(z, p, k, N, g, 1.0);
  T scaler = T(1)/wc;
  for(int n = 0; n < N; n++)
  {
    pNew[n] = scaler * p[n];
    zNew[n] = scaler * z[n];
  }
  int nz = rsGetNumFiniteValues(z, N);
  *kNew  = *k / pow(wc, N-nz);
}

template<class T>
void rsPrototypeDesigner<T>::getInverseFilter(Complex* z, Complex* p, T* k, Complex* zNew, 
  Complex* pNew, T* kNew, int N)
{
  //rassert(false); // something seems wrong about this - we should write the inverted poles, zeros
  //                // and gain zeros into zNew, pNew, kNew

  Complex *zTmp = new Complex[N]; // to make it work, when the new arrays are equal to the old ones
  rsArray::copyBuffer(z,    zTmp, N);
  rsArray::copyBuffer(p,    zNew, N);
  rsArray::copyBuffer(zTmp, pNew, N);
  *kNew = T(1) / *k;
  delete[] zTmp;
}

template<class T>
int rsPrototypeDesigner<T>::getLeftHalfPlaneRoots(T* a, Complex* r, int N)
{
  Complex *rTmp = new Complex[N]; // maybe we can get rid of that temporary array
  rsPolynomial<T>::findPolynomialRoots(a, N, rTmp);
  int numLeftRoots = rsOnlyLeftHalfPlane(rTmp, r, N);
  rsAssert(numLeftRoots == ceil(0.5*N)); // maybe take this out later
  delete[] rTmp;
  return numLeftRoots;
}

template<class T>
void rsPrototypeDesigner<T>::getBesselDenominatorCoeffs(T* a, int N)
{
  rsPolynomial<T>::besselPolynomial(a, N);
  rsArray::reverse(a, N+1);  // leaving this out leads to a modified Bessel filter response - maybe 
                             // experiment a bit, response looks good
}

template<class T>
void rsPrototypeDesigner<T>::getBesselLowShelfZerosPolesAndGain(Complex* z, Complex* p, T* k, 
  int N, T G, T G0)
{
  zpkFromTransferCoeffsLS(z, p, k, N, G, G0, &getBesselDenominatorCoeffs);
  return;
}

template<class T>
void rsPrototypeDesigner<T>::papoulisMagnitudeSquaredDenominator(T* a, int N)
{
  int n;
  rsArray::fillWithZeros(a, 2*N+1);  // do we need this?

  // construct the polynomial L_N(w^2):
  rsPolynomial<T>::maximumSlopeMonotonicPolynomial(a, N);  // does the same same as lopt(a, N); from C.R.Bond

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

//-----------------------------------------------
// new stuff:

template<class T>
void rsPrototypeDesigner<T>::papoulisPolynomial(T *v, int N)
{
  // temporary arrays for Legendre polynomials:
  T *P1 = new T[N/2+1];
  T *P2 = new T[N/2+1];
  T *P  = nullptr;  // pointer to the current P array

                    // create integrand:
  int k, r;
  if( rsIsOdd(N) )
  {
    k = (N-1)/2;

    // create weighted sum of Legendre polynomials in v:
    rsArray::fillWithZeros(v, k+1);
    for(r = 0; r <= k; r++)
    {
      updateLegendrePolynomial(&P, P1, P2, r);
      rsPolynomial<T>::weightedSumOfPolynomials(v, r, T(1), P, r, 2*r+T(1), v);
    }

    // square it:
    rsArray::convolve(v, k+1, v, k+1, v);
  }
  else
  {
    k = (N-2)/2;

    // generate Legendre polynomial of order k+1, store in P:
    for(r = 0; r <= k+1; r++)
      updateLegendrePolynomial(&P, P1, P2, r);

    // take the derivative, store in v:
    rsPolynomial<T>::polyDerivative(P, v, k+1);

    // square it:
    rsArray::convolve(v, k+1, v, k+1, v);

    // multiply by (x+1):
    v[2*k+1] = 0;
    for(r = 2*k+1; r >= 1; r--)
      v[r] += v[r-1];
  }

  // integrate from -1 to 2*w^2-1:
  T a[1] = { -1 };  
  T b[3] = { -1, 0, 2};
  rsPolynomial<T>::integratePolynomialWithPolynomialLimits(v, N-1, a, 0, b, 2, v);

  // scale, such that L^2(1) = 1:
  rsArray::scale(v, 2*N+1, 1.0 / rsArray::sum(v, 2*N+1));

  // clean up:
  delete[] P1;
  delete[] P2;
}

template<class T>
void rsPrototypeDesigner<T>::halpernPolynomial(T *a, int N)
{  
  a[0] = 0;
  T *a1 = &a[1];                // index shift of one for multiplication by x in Eq. 8.19
  rsHalpernU(a1, N-1);          // create U-polynomial
  rsArray::convolve(a1, N, a1, N, a1); // square U-polynomial
  rsArray::scale(a1, 2*N-1, 2*N);      // apply squared scale factor
  rsPolynomial<T>::polyIntegral(a, a, 2*N-1);    // compute integral from 0 to w
}

template<class T>
void rsPrototypeDesigner<T>::gaussianPolynomial(T *a, int N, T wc)
{  
  rsArray::fillWithZeros(a, 2*N+1);
  T g = log(T(2)) / (wc*wc);  // gamma
  T s = 1;                   // scaler
  for(int k = 0; k <= N; k++)
  {
    a[2*k] = s;    // == g^k / k! == pow(g, k) / rsFactorial(k);
    s *= g/(k+1);
  }
}

// end new
//-----------------------------------------------

template<class T>
void rsPrototypeDesigner<T>::getPapoulisLowpassZerosPolesAndGain(Complex* z, Complex* p, T* k, 
  int N)
{
  // find poles:
  T* a2 = new T[2*N+1];     // coefficients of magnitude-squared polynomial D(s)*D(-s)
  papoulisMagnitudeSquaredDenominator(a2, N);
  getLeftHalfPlaneRoots(a2, p, 2*N);

  // zeros are at infinity:
  rsArray::fillWithValue(z, N, Complex(RS_INF(T), 0.0));

  // set gain at DC to unity:
  *k = sqrt(T(1)/fabs(a2[2*N]));

  delete[] a2;
}

template<class T>
void rsPrototypeDesigner<T>::getPapoulisLowShelfZerosPolesAndGain(Complex* z, Complex* p, T* k, 
  int N, T G, T G0)
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
    G   = T(1) / G;
    G0  = T(1) / G0;
  }

  // factor out into a function getMagnitudeSquaredNumAndDen, then call this function here - this 
  // function can switch between Papoulis, Gauss, etc. and fill the arrays accordingly, Bessel is 
  // a special case - it doesn't need to find poles of the mag-squared function, there we can 
  // directly find the poles of the transfer function ...hmm...but maybe we can still factor out a 
  // function getPoles or something...

  // coefficients of the magnitude-squared polynomial D(s)*D(-s)
  T* a2 = new T[2*N+1];
  papoulisMagnitudeSquaredDenominator(a2, N);
  getLeftHalfPlaneRoots(a2, p, 2*N);

  // normalize denominator polynomial such that the leading coeff has unity as absolute value:
  T scaler = T(1) / fabs(a2[2*N]);
  for(int n = 0; n <= 2*N; n++)
    a2[n] *= scaler;

  // construct lowpass numerator:
  T* b2 = new T[2*N+1];
  rsArray::fillWithZeros(b2, 2*N+1);
  b2[0] = 1.0;

  // end of "factor out" ...in general, we need to scale the b2-polynomial also by dividing through 
  // the leading coeff?

  // adjust lowpass DC gain via k:
  *k = sqrt(fabs(a2[0]));  // in general: sqrt(fabs(a2[0]/b2[0])) ?
  //*k = sign(a2[0]) * sqrt(fabs(a2[0]));

  // obtain magnitude-squared numerator polynomial for shelving filter:
  T *bS = new T[2*N+1];
  shelvingMagSqrNumFromLowpassMagSqr(b2, a2, *k, N, G0, G, bS);

  // find left halfplane zeros (= zeros of the shelving filter):
  getLeftHalfPlaneRoots(bS, z, 2*N);

  // set gain constant for shelving filter:
  *k = G0;

  // adjust bandwidth:
  T GB = sqrt(G*G0);
  scaleToMatchGainAtUnity(z, p, k, z, p, k, N, GB);

  // invert filter in case of a dip:
  if( dip == true )
    getInverseFilter(z, p, k, z, p, k, N);

  delete[] a2;
  delete[] b2;
  delete[] bS;
}

//-------------------------------------------------------
// refactoring - not yet finsihed:

template<class T>
void rsPrototypeDesigner<T>::zpkFromTransferCoeffsLP(Complex* z, Complex* p, T* k, int N,
  void (*denominatorCoeffsFunction)(T* a, int N))
{
  // zeros are at infinity:
  rsArray::fillWithValue(z, N, Complex(RS_INF(T), 0.0));

  // find poles:
  T* a = new T[N+1];
  denominatorCoeffsFunction(a, N);
  rsPolynomial<T>::findPolynomialRoots(a, N, p);

  // set gain and scale poles to match Butterworth magnitude response asymptotically, if desired:
  bool matchButterworth = true; // maybe make this a parameter later
  if( matchButterworth == true )
  {
    T scaler = T(1) / pow(a[0], T(1)/N);
    for(int n = 0; n < N; n++)
      p[n] *= scaler;
    *k = T(1);
  }
  else
    *k = a[0];

  delete[] a;
}

template<class T>
void rsPrototypeDesigner<T>::zpkFromTransferCoeffsLS(Complex* z, Complex* p, T* k, int N, 
  T G, T G0, void (*denominatorCoeffsFunction)(T* a, int N))
{
  // catch lowpass case:
  if( G0 == 0.0 )
  {
    zpkFromTransferCoeffsLP(z, p, k, N, denominatorCoeffsFunction);
    *k *= G;
    return;
  }

  // design boost filter and invert later, if a dip is desired:
  bool dip = false;
  if(G < G0)
  {
    dip = true;
    G   = T(1) / G;
    G0  = T(1) / G0;
  }

  // construct lowpass denominator:
  T* a = new T[N+1];
  denominatorCoeffsFunction(a, N);

  // find poles of the shelving filter:
  rsPolynomial<T>::findPolynomialRoots(a, N, p);

  // construct lowpass numerator:
  T* b = new T[N+1];
  rsArray::fillWithZeros(b, N+1);
  b[0] = a[0];

  // obtain magnitude-squared numerator polynomial for shelving filter:
  T* bS = new T[2*N+1];
  shelvingMagSqrNumeratorFromLowpassTransfer(b, a, 1.0, N, G0, G, bS);

  // find left halfplane zeros (= zeros of the shelving filter):
  getLeftHalfPlaneRoots(bS, z, 2*N);

  // set gain constant:
  *k = G0;

  // now we have a shelving filter with correct low-frequency gain G and reference gain G0, but 
  // possibly still with wrong bandwidth gain GB at unity - now we adjust zeros/poles/gain to 
  // match GB:
  T GB = sqrt(G*G0);
  scaleToMatchGainAtUnity(z, p, k, z, p, k, N, GB);

  // invert filter in case of a dip:
  if( dip == true )
  getInverseFilter(z, p, k, z, p, k, N);

  // cleanup:
  delete[] a;
  delete[] b;
  delete[] bS;
}

// end refactoring
//-------------------------------------------------------


template<class T>
void rsPrototypeDesigner<T>::getEllipticLowpassZerosPolesAndGain(Complex* z, Complex* p, T* k, 
  int N, T Gp, T Gs)
{
//  int nz;
//  if( rsIsEven(N) )
//    nz = N;
//  else
//    nz = N-1;

  int L, r;
  getNumBiquadsAndFirstOrderStages(N, L, r);

  // declare/assign/calculate some repeatedly needed variables:
  Complex j(0.0, 1.0);                            // imaginary unit
  T  ep  = sqrt(T(1)/(Gp*Gp) - T(1));             // Eq. 2
  T  es  = sqrt(T(1)/(Gs*Gs) - T(1));             // Eq. 2
  T  k1  = ep/es;                                 // Eq. 3
  T  kk  = ellipdeg(N, k1);                       // solve degree equation for k
  T  v_0 =  (-j*rsAsnC(j/ep, k1) / T(N)).real(); // from ellipap.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    //p[L+r-1] = -Omega_p/sinh(v_0*PI*0.5*kk);                 // Eq. 73
    p[N-1] = j * rsSnC(j*v_0, kk);                               // from ellipap.m - find Eq.
    z[N-1] = RS_INF(T);
  }

  // calculate the complex conjugate poles and zeros:
  //T  u_i;
  Complex zeta_i;
  for(int i = 0; i < L; i++)
  {
    Complex u_i = (T) (2*(i+1)-1) / (T) N;             // Eq. 69
    zeta_i   = rsCdC(u_i, kk);                         // from ellipap.m - find Eq.
    z[2*i]   = j / (kk*zeta_i);                        // Eq. 62
    p[2*i]   = j*rsCdC((u_i-j*v_0), kk);
    z[2*i+1] = conj(z[2*i]); 
    p[2*i+1] = conj(p[2*i]);
  }

  T H0 = pow(Gp, 1-r);    // preliminary - can be made simpler (without pow)
  Complex n = rsProductOfFiniteFactors(p, N);
  Complex d = rsProductOfFiniteFactors(z, N);
  *k        = H0 * (n/d).real();
}

// inquiry:

template<class T>
std::complex<T> rsPrototypeDesigner<T>::getFilterResponseAt(Complex s)
{
  Complex num, den;
  Complex tmp;
  int     Lz, Lp;

  // initialize the numerator and denominator:
  if( rsIsOdd(numFiniteZeros) )
  {
    num = -z[L+r-1].real();
    Lz  = (numFiniteZeros-1)/2;
  }
  else
  {
    num = 1.0;
    Lz  = numFiniteZeros/2;
  }
  if( rsIsOdd(numFinitePoles) )
  {
    den = -p[L+r-1].real();
    Lp  = (numFinitePoles-1)/2;
  }
  else
  {
    den = 1.0;
    Lp  = numFinitePoles/2;
  }

  // accumulate product of linear factors for denominator (poles) and numerator (zeros):
  int i;
  for(i = 0; i < Lz; i++)
    num *= ((s-z[i]) * (s-conj(z[i])));
  for(i = 0; i < Lp; i++)
    den *= ((s-p[i]) * (s-conj(p[i])));

  return num/den;
}

template<class T>
T rsPrototypeDesigner<T>::getMagnitudeAt(T w)
{
  return abs(getFilterResponseAt(Complex(0.0, w)));
}

template<class T>
T rsPrototypeDesigner<T>::findFrequencyWithMagnitude(T magnitude, T wLow, T wHigh)
{
  // until we have something better, we search for the frequency at which the desired gain occurs 
  // by means of the bisection method:

  T tolerance = T(0.0001); // maybe make parameter
  T wMid, mMid;
  while( wHigh-wLow > tolerance )
  {
    wMid  = T(0.5) * (wLow+wHigh);
    mMid  = getMagnitudeAt(wMid);
    if( mMid > magnitude )
      wLow = wMid;
    else
      wHigh = wMid;
  }
  return T(0.5) * (wLow+wHigh);
}

template<class T>
int rsPrototypeDesigner<T>::getNumFinitePoles()
{
  updatePolesAndZeros();
  return numFinitePoles;
}

template<class T>
int rsPrototypeDesigner<T>::getNumFiniteZeros()
{
  updatePolesAndZeros();
  return numFiniteZeros;
}

template<class T>
int rsPrototypeDesigner<T>::getNumNonRedundantFinitePoles()
{
  updatePolesAndZeros();
  if( rsIsEven(numFinitePoles) )
    return numFinitePoles/2;
  else
    return (numFinitePoles+1)/2;
}

template<class T>
int rsPrototypeDesigner<T>::getNumNonRedundantFiniteZeros()
{
  updatePolesAndZeros();
  if( rsIsEven(numFiniteZeros) )
    return numFiniteZeros/2;
  else
    return (numFiniteZeros+1)/2;
}

template<class T>
void rsPrototypeDesigner<T>::getPolesAndZeros(Complex* poles, Complex* zeros)
{
  //if( stateIsDirty == true )// re-calculate only if necesarry ..is actually checked there internally
  updatePolesAndZeros(); 
  for(int i = 0; i < (L+r); i++)
  {
    poles[i] = p[i];
    zeros[i] = z[i];
  }
}

template<class T>
bool rsPrototypeDesigner<T>::hasCurrentMethodRippleParameter()

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

template<class T>
bool rsPrototypeDesigner<T>::hasCurrentMethodRejectionParameter()
{
  if( prototypeMode == LOWPASS_PROTOTYPE )
  {
    return (approximationMethod == ELLIPTIC) || (approximationMethod == INVERSE_CHEBYCHEV);
  }
  else
    return false;
}

template<class T>
bool rsPrototypeDesigner<T>::needsSpecialHighShelvTransform()
{
  return (approximationMethod == BUTTERWORTH)
    ||   (approximationMethod == CHEBYCHEV)
    ||   (approximationMethod == INVERSE_CHEBYCHEV)
    ||   (approximationMethod == ELLIPTIC);
}

// pole/zero calculation:

template<class T>
void rsPrototypeDesigner<T>::updatePolesAndZeros()
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
      case BUTTERWORTH:       makeButterworthLowShelv();                         break;
      case CHEBYCHEV:         makeChebychevLowShelv();                           break;
      case INVERSE_CHEBYCHEV: makeInverseChebychevLowShelv();                    break;
      case ELLIPTIC:          makeEllipticLowShelv();                            break;
      case BESSEL:            makeBesselLowShelv(  rsDbToAmp(A), rsDbToAmp(A0)); break;
      case PAPOULIS:          makePapoulisLowShelv(rsDbToAmp(A), rsDbToAmp(A0)); break;
      }
    }
    stateIsDirty = false;
  }
}

template<class T>
void rsPrototypeDesigner<T>::makeBypass()
{
  numFinitePoles = 0;
  numFiniteZeros = 0;
}

template<class T>
void rsPrototypeDesigner<T>::makeButterworthLowpass()
{
  numFinitePoles = N;
  numFiniteZeros = 0;

  // intermediate variables:
  Complex j(0.0, 1.0);                   // imaginary unit
  T  Gp     = sqrt(T(0.5));              // use -3.01 dB point as cutoff frequency for Butterworths
  //T  Gp   = pow(10.0, -Ap/20.0);       // (1),Eq.1 - more general (cutoff gain can be specified), not used here
  T  ep     = sqrt(T(1)/(Gp*Gp)-T(1));   // (1),Eq.2
  T  ep_pow = pow(ep, T(-1)/(T) N);

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -ep_pow;                                    // Eq.70
    z[L+r-1] = RS_INF(T);                                  // zero at infinity
  }
  // calculate the complex conjugate poles and zeros:
  //T  u_i;
  for(int i = 0; i < L; i++)
  {
    Complex u_i  = (T) (2*(i+1)-1) / (T) N;                // Eq.69
    //p[i] = ep_pow*j*rsExpC(j*u_i*T(PI*0.5));                  // old
    p[i] = ep_pow*j*exp(j*u_i*T(PI*0.5));                  // Eq.70
    z[i] = RS_INF(T);                                      // zeros are at infinity
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeButterworthLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // catch some special cases:
  if( A0 == -RS_INF(T) ) // lowpass-case
  {
    makeButterworthLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // intermediate variables:
  T G0   = rsDbToAmp(A0);
  T G    = rsDbToAmp(A);
  T GB   = sqrt(G0*G);                           // (2),Eq.52
  T ep   = sqrt( (G*G-GB*GB) / (GB*GB-G0*G0) );  // (2),Eq.12
  T g0   = pow(G0, T(1) / T(N));                 // (2),Eq.94
  T g    = pow(G,  T(1) / T(N));                 // (2),Eq.94
  T wb   = 1.0;                                  // unit cutoff prototype
  T beta = wb * pow(ep, T(-1) / (T) N);          // (2),Eq.94

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -beta;                                      // (2),Eq.93
    z[L+r-1] = -g*beta/g0;                                 // (2),Eq.93
  }
  // calculate the complex conjugate poles and zeros:
  T phi, s, c;
  for(int i = 0; i < L; i++)
  {
    phi = (T) (2*(i+1)-1)*T(PI) / T(2*N);    // (2),Eq.95
    s   = sin(phi);                          // (2),Eq.95
    c   = cos(phi);                          // (2),Eq.95
    z[i].real(-s*g*beta/g0);                 // (2),Eq.93
    z[i].imag( c*g*beta/g0);                 // (2),Eq.93
    p[i].real(-s*beta);                      // (2),Eq.93
    p[i].imag( c*beta);                      // (2),Eq.93
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeChebychevLowpass()
{
  numFinitePoles = N;
  numFiniteZeros = 0;

  // intermediate variables:
  T Gp   = pow(T(10), -Ap/T(20));            // Eq. 1
  T ep   = sqrt(T(1)/(Gp*Gp) - T(1));        // Eq. 2
  T v_0  = asinh(T(1)/ep) / T(N*PI*0.5);     // Eq. 72
  T u_i;
  Complex j(0.0, 1.0); // imaginary unit

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = -sinh(v_0*T(PI*0.5));         // Eq. 71
    z[L+r-1] = RS_INF(T);
  }
  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    u_i  = (T) (2*(i+1)-1) / (T) N;         // Eq. 69
    p[i] = j*rsCosC((u_i-j*v_0)*T(PI*0.5)); // Eq. 71
    z[i] = RS_INF(T);                       // zeros at infinity
  }

  //gain = 1.0 / getFilterResponseAt(Complex(0.0, 0.0)).getMagnitude();
  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeChebychevLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // calculate the linear gain-factors:
  T G0 = rsDbToAmp(A0);
  T G  = rsDbToAmp(A);

  // catch some special cases:
  if( A0 == -RS_INF(T) ) // lowpass-case
  {
    makeChebychevLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // calculate intermediate variables:
  T Gp    = rsDbToAmp(A0 + Rp*A);
  T ep    = sqrt( (G*G-Gp*Gp) / (Gp*Gp-G0*G0) );
  T g0    = pow(G0, T(1) / T(N));
  //T g     = pow(G,   1.0 / (T) N);
  T alpha = pow(T(1)/ep +  sqrt(T(1) + T(1)/(ep*ep)),   T(1)/(T) N);
  T beta  = pow((G/ep + Gp*sqrt(T(1) + T(1)/(ep*ep)) ), T(1)/(T) N);
  T u     = log(beta/g0);
  T v     = log(alpha);
  T Gb    = sqrt(G0*G);
  T eb    = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );
  T wb    = T(1) / cosh( acosh(eb/ep) / T(N) ); // why 1/cosh(...) and not simply cosh?

  // calculate real pole and zero of the first order stage, if present and store them in the last array slots:
  if( r == 1 )
  {
    p[L+r-1] = -wb*sinh(v);
    z[L+r-1] = -wb*sinh(u);
  }

  // calculate the complex conjugate poles and zeros:
  T phi_i; //, s, c;
  Complex j(0.0, 1.0); // imaginary unit
  for(int i = 0; i < L; i++)
  {
    phi_i = (T) (2*(i+1)-1)*T(PI) / T(2*N);
    z[i]  = j*wb*rsCosC(phi_i - j*u);
    p[i]  = j*wb*rsCosC(phi_i - j*v);
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeInverseChebychevLowpass()
{
  numFinitePoles = N;
  if( rsIsEven(N) )
    numFiniteZeros = N;
  else
    numFiniteZeros = N-1;

  // declare/assign/calculate some repeatedly needed variables:
  T  Gs = pow(T(10), -As/T(20));                    // Eq. 1
  T  es = sqrt(T(1) / (Gs*Gs)-T(1));                // Eq. 2
  T  v0 = asinh(es) / T(N*PI*0.5);                  // Eq. 74
  Complex j(0.0, 1.0);                              // imaginary unit

  T  wb = 1.0; // ...leads to a gain of Gs (stopband-gain) at unity (w=1), we rescale it here
                    // so as to have the -3 dB point at w=1:
  T  Gp = sqrt(T(0.5));
  T  ep = sqrt(T(1)/(Gp*Gp)-T(1));
  wb    = cosh( acosh(es/ep) / N );                             // (1),Eq.9

  // calculate the position of the real pole (if present):
  //T ui;
  if( r == 1 )
  {
    p[L+r-1] = -wb / sinh(v0*T(PI*0.5));                           // Eq.73 with k=1
    z[L+r-1] = RS_INF(T);
  }

  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    Complex ui = (T) (2*(i+1)-1) / (T) N;                            // Eq.69
    z[i] = wb / (j*rsCosC(ui*T(PI/2)));                              // Eq.73 with k=1
    p[i] = wb / (j*rsCosC((ui - j*v0)*T(PI/2)));                     // Eq.73 with k=1
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeInverseChebychevLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // calculate the linear gain-factors:
  T G0 = rsDbToAmp(A0);
  T G  = rsDbToAmp(A);

  // catch some special cases:
  if( A0 == -RS_INF(T) ) // lowpass-case
  {
    makeInverseChebychevLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // calculate intermediate variables (\todo check if the gains have reasonable values):
  //T Gs    = rsDbToAmp(Rs*G + (1.0-Rs)*G0);
  T Gs    = rsDbToAmp(A0 + Rs*A);
  T es    = sqrt( (G*G-Gs*Gs) / (Gs*Gs-G0*G0) );
  //T g0    = pow(G0, 1.0 / (T) N);
  T g     = pow(G,   T(1) / T(N));
  T alpha = pow(es + sqrt(T(1) + es*es), T(1) / T(N));
  T beta  = pow((G0*es + Gs*sqrt(T(1) + es*es) ), T(1)/T(N));
  T u     = log(beta/g);
  T v     = log(alpha);
  T Gb    = sqrt(G0*G);
  T eb    = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );
  T wb    = cosh( acosh(es/eb) / T(N) );  // why not 1 / cosh(..)?

  // calculate real pole and zero of the first order stage, if present and store them in the last
  // array slots:
  if( r == 1 )
  {
    z[L+r-1] = -wb/sinh(u);
    p[L+r-1] = -wb/sinh(v);
  }

  // calculate the complex conjugate poles and zeros:
  T  phi_i;
  Complex j(0.0, 1.0); // imaginary unit
  for(int i = 0; i < L; i++)
  {
    phi_i = (T) (2*(i+1)-1)*T(PI) / T(2*N);
    z[i]  = wb / (j*rsCosC(phi_i-j*u));
    p[i]  = wb / (j*rsCosC(phi_i-j*v));
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeEllipticLowpass()
{
  numFinitePoles = N;
  if( rsIsEven(N) )
    numFiniteZeros = N;
  else
    numFiniteZeros = N-1;

  // declare/assign/calculate some repeatedly needed variables:
  Complex j(0.0, 1.0);                                    // imaginary unit
  //T  u_i;
  Complex zeta_i;
  T  Gp  = pow(T(10), -Ap/T(20));               // Eq. 1
  T  Gs  = pow(T(10), -As/T(20));               // Eq. 1
  T  ep  = sqrt(T(1)/(Gp*Gp) - T(1));           // Eq. 2
  T  es  = sqrt(T(1)/(Gs*Gs) - T(1));           // Eq. 2
  T  k1  = ep/es;                               // Eq. 3
  T  k   = ellipdeg(N, k1);                     // solve degree equation for k
  T  v_0 =  (-j*rsAsnC(j/ep, k1)/(T)N).real();  // from ellipap.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    //p[L+r-1] = -Omega_p/sinh(v_0*PI*0.5*k);                    // Eq. 73
    p[L+r-1] = j * rsSnC(j*v_0, k);                                // from ellipap.m
    z[L+r-1] = RS_INF(T);
  }
  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    Complex u_i = (T) (2*(i+1)-1) / (T) N;                       // Eq. 69
    zeta_i = rsCdC(u_i, k);                                      // from ellipap.m
    z[i]   = j / (k*zeta_i);                                     // Eq. 62
    p[i]   = j*rsCdC((u_i-j*v_0), k);
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeEllipticLowShelv()
{
  numFinitePoles = N;
  numFiniteZeros = N;

  // catch some special cases:
  if( A0 == -RS_INF(T) ) // lowpass-case
  {
    makeEllipticLowpass();
    return;
  }
  else if( abs(A-A0) <  0.001 )
    makeBypass();

  // intermediate variables:
  T  G0  = rsDbToAmp(A0);                         // reference amplitude
  T  G   = rsDbToAmp(A);                          // boost/cut amplitude
  T  Gp  = rsDbToAmp(A0 + Rp*A);                  // passband-amplitude (Rp near 1)
  T  Gs  = rsDbToAmp(A0 + Rs*A);                  // stopband-amplitude (Rs near 0)
  T  Gb  = sqrt(G0*G);                            // (2),Eq.52 (gain at the bandedges)
  T  ep  = sqrt( (G*G-Gp*Gp) / (Gp*Gp-G0*G0) );   // (2),Eq.12
  T  es  = sqrt( (G*G-Gs*Gs) / (Gs*Gs-G0*G0) );   // (2),Eq.39
  T  eb  = sqrt( (G*G-Gb*Gb) / (Gb*Gb-G0*G0) );   // (2),Eq.64
  T  k1  = ep/es;                                 // (2),Eq.39
  T  k   = ellipdeg(N, k1);                       // degree equation
  //Complex u = rsAcdC(eb/ep, k1) / N;            // old
  Complex u = rsAcdC(Complex(eb/ep), k1) / T(N);  // following text after (2),Eq.65
  T  wb  = T(1) / rsCdC(u, k).real();             // ...ditto
  Complex j   = Complex(0.0, 1.0);                // imaginary unit
  Complex ju0 = rsAsnC(j*G/(ep*G0), k1) / T(N);   // line 111 in hpeq.m
  Complex jv0 = rsAsnC(j  / ep,     k1) / T(N);   // line 113 in hpeq.m

  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = wb*(j*rsCdC(T(-1)+jv0,k)).real();        // line 148 in hpeq.m
    z[L+r-1] = wb*(j*rsCdC(T(-1)+ju0,k)).real();        // line 145 in hpeq.m
  }

  // calculate the complex conjugate poles and zeros:
  //T ui;
  for(int i = 0; i < L; i++)
  {
    Complex ui = (T) (2*(i+1)-1) / (T) N;              // (2),Eq.37
    p[i] = j*wb * rsCdC( (ui-jv0), k );                // line 179 in hpeq.m

    if( G0 == 0.0 && G != 0.0 )                        // lines 172-178 in hpeq.m
      z[i] = j*wb / (k*rsCdC(ui,k));   // lowpass
    else if( G0 != 0.0 && G == 0.0 )
      z[i] = j*wb * rsCdC(ui,k);       // highpass
    else
      z[i] = j*wb * rsCdC(ui-ju0,k);   // low-shelv
  }

  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makeBesselLowShelv(T G, T G0)
{
  rsArray::fillWithZeros(p, maxNumNonRedundantPoles);
  rsArray::fillWithZeros(z, maxNumNonRedundantPoles);
  numFinitePoles = N;
  if( G0 == 0.0 )
    numFiniteZeros = 0;
  else
    numFiniteZeros = N;

  Complex zTmp[25];
  Complex pTmp[25];
  T  kTmp;
  rsPrototypeDesigner::getBesselLowShelfZerosPolesAndGain(zTmp, pTmp, &kTmp, N, G, G0);

  // findPolynomialRoots returns the roots sorted by ascending real part. for a Bessel-polynomial, 
  // this ensures that the real pole, if present, is in pTmp[0] (it has the largest negative real 
  // part). this is importatnt for the next call:

  pickNonRedundantPolesAndZeros(zTmp, pTmp);
  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::makePapoulisLowShelv(T G, T G0)
{
  rsArray::fillWithZeros(p, maxNumNonRedundantPoles);
  rsArray::fillWithZeros(z, maxNumNonRedundantPoles);
  numFinitePoles = N;
  if( G0 == 0.0 )
    numFiniteZeros = 0;
  else
    numFiniteZeros = N;

  Complex zTmp[25];
  Complex pTmp[25];
  T kTmp;
  rsPrototypeDesigner::getPapoulisLowShelfZerosPolesAndGain(zTmp, pTmp, &kTmp, N, G, G0);

  // findPolynomialRoots returns the roots sorted by ascending real part. for a Bessel-polynomial, 
  // this ensures that the real pole, if present, is in pTmp[0] (it has the largest negative real 
  // part). this is importatnt for the next call:

  pickNonRedundantPolesAndZeros(zTmp, pTmp);
  stateIsDirty = false;
}

template<class T>
void rsPrototypeDesigner<T>::pickNonRedundantPolesAndZeros(Complex *zTmp, Complex *pTmp)
{
  //const T thresh = T(1.e-11); // preliminary - old, worked for double
  T thresh = 10 * std::numeric_limits<T>::epsilon(); // ad-hoc - needs testing for double

  rsZeroNegligibleImaginaryParts(pTmp, N, thresh);
  rsZeroNegligibleImaginaryParts(zTmp, N, thresh);
  rsOnlyUpperHalfPlane(pTmp, pTmp, N);
  rsOnlyUpperHalfPlane(zTmp, zTmp, N);
  rsArray::copyBuffer(pTmp, p, L+r);
  rsArray::copyBuffer(zTmp, z, L+r);

  // the caller is supposed to ensure that the real zero/pole, if present, is in zTmp[0], pTmp[0] - 
  // but we need it in the last positions z[L+r], p[L+r], so we reverse the arrays:
  rsArray::reverse(p, L+r);
  rsArray::reverse(z, L+r);
}
