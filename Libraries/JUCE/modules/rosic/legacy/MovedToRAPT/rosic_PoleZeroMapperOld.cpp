#include <limits>

#include "rosic_PoleZeroMapper.h"
#include "../math/rosic_PolynomialAlgorithms.h"
#include "rosic_FilterAnalyzer.h"
using namespace rosic;


double findPrototypeCutoff(Complex *z, Complex *p, double *k, int N, double A, double initialGuess = 1.0)
{
  //if( A <= 0.0 || A >= 1.0 )
  //  rosic::error("A must be between 0 and 1 (exclusive) in findPrototypeCutoff");

  // determine the frequency interval in which the magnitude value A must fall (maybe wrap into a function):
  double wL = initialGuess;
  double wU = initialGuess;
  while( FilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wL) < A )
    wL = wL/2;
  while( FilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wU) > A )
    wU = wU*2;

  // find the frequency at which the magnitude is equal to A by bisection:  
  const double eps  = std::numeric_limits<double>::epsilon(); 
  int i       = 0;         // iteration counter
  int iMax    = 1000;      // maximum number of iterations
  double tol  = A*eps;     // amplitude tolerance - maybe use A*eps?
  double wM   = (wL+wU)/2; // midpoint
  double wOld = wM;
  double AM   = FilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wM);
  while( abs(AM-A) > tol && i < iMax ) 
  {
    if( AM > A )
      wL = wM;
    else if( AM < A )
      wU = wM;
    wOld = wM;
    wM   = (wL+wU)/2;
    AM = FilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wM);
    i++;

    if( abs(wM-wOld) < eps*min(wOld, wM) ) // 2nd convergence cirterion based on frequency difference - maybe use as only criterion
      break;  
  }

  if( i >= iMax )
    rosic::error("Slow convergence in findPrototypeCutoff");

  return wM;
}

void scaleToMatchGainAtCutoff(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, double *kNew, int N, double GC)
{
  double  wc    = findPrototypeCutoff(z, p, k, N, GC, 1.0);
  double scaler = 1.0/wc;
  for(int n = 0; n < N; n++)
  {
    pNew[n] = scaler * p[n];
    zNew[n] = scaler * z[n];
  }
  int nz = getNumFiniteValues(z, N);
  *kNew  = *k / pow(wc, N-nz);
}

// move to PrototypeDesigner, provide functions to design shelving prototypes directly there
void constructShelvingMagnitudeSquaredNumerator(Complex *b, Complex *a, double k, int N, double G0, double G, Complex *bShelf)
{
  Complex *am = new Complex[N+1];
  Complex *bm = new Complex[N+1];
  Complex *a2 = new Complex[2*N+1];
  Complex *b2 = new Complex[2*N+1];

  // construct N_LP(s)*N_LP(-s) and D_LP(s)*D_LP(-s):
  polyCoeffsForNegativeArgument(b, bm, N);  // coeffs of N_LP(-s)
  polyCoeffsForNegativeArgument(a, am, N);  // coeffs of D_LP(-s)
  multiplyPolynomials(b, N, bm, N, b2);     // coeffs of N_LP(s)*N_LP(-s)
  multiplyPolynomials(a, N, am, N, a2);     // coeffs of D_LP(s)*D_LP(-s)

  // obtain coefficients for shelving filter's magnitude squared function numerator polynomial N_LS(s)*N_LS(-s):
  weightedSum(b2, a2, bShelf, 2*N+1, Complex(k*k*(G*G-G0*G0), 0.0), Complex(G0*G0, 0.0));

  delete[] am;
  delete[] bm;
  delete[] a2;
  delete[] b2;
}

void PoleZeroMapper::sLowpassToLowshelf(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, double *kNew, 
                                        int N, double G0, double G)
{
  double kTmp;
  Complex *zTmp = new Complex[2*N];
  Complex *pTmp = new Complex[2*N];
  Complex *a    = new Complex[N+1];
  Complex *b    = new Complex[N+1];
  Complex *bS   = new Complex[2*N+1];  // shelving mgnitude squared numerator coeffs

  // we design a boost filter in any case - a dip filter will later be obtained by swapping poles and zeros:
  bool dip = false;
  if( G < G0 )
  {
    dip = true;
    G   = 1.0 / G;
    G0  = 1.0 / G0;
  }

  // scale poles/zeros/gain of the lowpass such that the resulting low-shelving filter has a gain of sqrt(G*G0) at unit frequency:
  double GB = sqrt(G0*G);                       // bandwidth gain  

  //double GC = sqrt((GB*GB-G0*G0)/(G*G-G0*G0));  // desired lowpass gain at unit frequency
  //scaleToMatchGainAtCutoff(z, p, k, zTmp, pTmp, &kTmp, N, GC);

  copyBuffer(z, zTmp, N);
  copyBuffer(p, pTmp, N);


  // obtain magnitude-squared polynomial for shelving filter:
  //rootsToCoeffs(zTmp, b, N, true);
  //rootsToCoeffs(pTmp, a, N, true);
  //constructShelvingMagnitudeSquaredNumerator(b, a, kTmp, N, G0, G, bS);
  rootsToCoeffs(z, b, N, true);
  rootsToCoeffs(p, a, N, true);
  constructShelvingMagnitudeSquaredNumerator(b, a, *k, N, G0, G, bS);


  // obtain zeros, poles and gain of the new shelving filter:
  findPolynomialRoots(bS, 2*N, &zTmp[-1]);
  int numLeftZeros = onlyLeftHalfPlane(zTmp, zTmp, 2*N); // maybe copy directly to zNew (->zNew as 2nd argument)


  if( dip == false )
    kTmp = sqrt(bS[2*N].getRadius());
  else
    kTmp = 1.0 / sqrt(bS[2*N].getRadius());

  scaleToMatchGainAtCutoff(zTmp, pTmp, &kTmp, zTmp, pTmp, &kTmp, N, GB);


  // if we make a dip-filter, poles and zeros exchange roles:
  if( dip == false )
  {
    copyBuffer(zTmp, zNew, N);
    copyBuffer(pTmp, pNew, N);
    *kNew = kTmp;
    
    //*kNew = sqrt(bS[2*N].getRadius());
  }
  else
  {
    copyBuffer(zTmp, pNew, N);
    copyBuffer(pTmp, zNew, N);
    *kNew = 1.0 / kTmp;

    //*kNew = 1.0 / sqrt(bS[2*N].getRadius());
  }

  delete[] zTmp;
  delete[] pTmp;
  delete[] a;
  delete[] b;
  delete[] bS;
}

void PoleZeroMapper::sLowpassToLowpass(Complex *z,    Complex *p,    double *k, 
                                       Complex *zNew, Complex *pNew, double *kNew, 
                                       int N, double wc)
{
  for(int n = 0; n < N; n++)
    zNew[n] = wc * z[n];
  for(int n = 0; n < N; n++)
    pNew[n] = wc * p[n];
  int nz = getNumFiniteValues(z, N); // number of finite zeros in prototype
  *kNew = *k * pow(wc, N-nz);
}



void PoleZeroMapper::sLowpassToHighpass(Complex *r, Complex *rNew, int N, double wc)
{
  for(int n = 0; n < N; n++)
  {
    if( r[n].isInfinite() )
      rNew[n] = Complex(0.0, 0.0);
    else
      rNew[n] = wc / r[n];
  }
}
void PoleZeroMapper::sLowpassToHighpass(Complex *z,    Complex *p,    double *k, 
                                        Complex *zNew, Complex *pNew, double *kNew, 
                                        int N, double wc)
{
  sLowpassToHighpass(z, zNew, N, wc);
  sLowpassToHighpass(p, pNew, N, wc);
  Complex n = productOfFiniteFactors(z, N);  
  Complex d = productOfFiniteFactors(p, N);
  *kNew     = *k * (n/d).re;  
  // \todo: include sign-factor: kNew = k * (-1)^(length(p) + length(z)) * prod(z) / prod(p); 
  // ..or maybe not? ..consult octave/matlab implementation
}

void PoleZeroMapper::sLowpassToBandpass(Complex *r, Complex *rNew, int N, double wc, double bw)
{
  Complex tmp1, tmp2;
  for(int n = 0; n < N; n++)
  {
    tmp1      = 0.5 * bw * r[n];
    tmp2      = -bw * r[n];

    // this is common with the bandreject case - factor out:
    tmp2      = 0.25 * tmp2*tmp2;
    tmp2      = sqrtC(tmp2 - wc*wc);
    rNew[n]   = tmp1 + tmp2;
    rNew[N+n] = tmp1 - tmp2;
  }
}
void PoleZeroMapper::sLowpassToBandpass(Complex *z,    Complex *p,    double *k, 
                                        Complex *zNew, Complex *pNew, double *kNew, 
                                        int N, double wl, double wu)
{
  double wc = sqrt(wl*wu);  // center (radian) frequency
  double bw = wu-wl;        // bandwidth
  sLowpassToBandpass(z, zNew, N, wc, bw);
  sLowpassToBandpass(p, pNew, N, wc, bw);
  for(int n = 0; n < N; n++)
  {
    if( z[n].isInfinite() )
    {
      zNew[n]   = Complex(0.0, 0.0);
      zNew[N+n] = Complex(INF, 0.0);
    }
  }
  int nz = getNumFiniteValues(z, N); // number of finite zeros in prototype
  *kNew  = *k * pow(bw, N-nz); 
}

void PoleZeroMapper::sLowpassToBandreject(Complex *r, Complex *rNew, int N, double wc, double bw)
{
  Complex tmp1, tmp2;
  for(int n = 0; n < N; n++)
  {
    tmp1      = 0.5 * bw / r[n];
    tmp2      = -bw / r[n];

    // this is common with the bandpass case - factor out:
    tmp2      = 0.25 * tmp2*tmp2;
    tmp2      = sqrtC(tmp2 - wc*wc);
    rNew[n]   = tmp1 + tmp2;
    rNew[N+n] = tmp1 - tmp2;
  }
}
void PoleZeroMapper::sLowpassToBandreject(Complex *z,    Complex *p,    double *k, 
                                          Complex *zNew, Complex *pNew, double *kNew, 
                                          int N, double wl, double wu)
{
  double wc = sqrt(wl*wu);  // center (radian) frequency
  double bw = wu-wl;        // bandwidth
  sLowpassToBandreject(z, zNew, N, wc, bw);
  sLowpassToBandreject(p, pNew, N, wc, bw);
  Complex sz = Complex(0.0, wc);
  for(int n = 0; n < N; n++)
  {
    if( z[n].isInfinite() )
    {
      zNew[n]   =  sz;
      zNew[N+n] = -sz;
    }
  }
  Complex n = productOfFiniteFactors(z, N);  
  Complex d = productOfFiniteFactors(p, N);
  *kNew     = *k * (n/d).re;  
  // \todo: include sign-factor - see m-file: kNew = k * (-1)^(length(p) + length(z)) * prod(z) / prod(p);
}



void PoleZeroMapper::prototypeToAnalogLowpass(Complex *poles, int numPoles, Complex *zeros, int numZeros, double *gain, 
                                                 double targetCutoff)
{
  int i;
  for(i=0; i<numPoles; i++)
    poles[i] = poles[i] * targetCutoff;
  for(i=0; i<numZeros; i++)
    zeros[i] = zeros[i] * targetCutoff;

  // \todo adjust gain
}

void PoleZeroMapper::sPlanePrototypeToLowpass(Complex *prototypePoles, Complex *prototypeZeros, Complex *targetPoles, 
                                                 Complex *targetZeros, int prototypeOrder, double targetCutoff)
{
  for(int i=0; i<prototypeOrder; i++)
  {
    targetPoles[i] = prototypePoles[i] * targetCutoff;
    targetZeros[i] = prototypeZeros[i] * targetCutoff;
  }
}

void PoleZeroMapper::prototypeToAnalogHighpass(Complex *poles, int numPoles, Complex *zeros, int numZeros, double *gain, 
                                                  double targetCutoff)
{
  int i;
  for(i=0; i<numPoles; i++)
    poles[i] = targetCutoff / poles[i];
  for(i=0; i<numZeros; i++)
    zeros[i] = targetCutoff / zeros[i];

  // there are numPoles-numZeros zeros at infinity in the lowpass-prototype and for for each such zero at infinity, we obtain a zero at 
  // s=0 in the highpass filter:
  for(i=numZeros; i<numPoles; i++)
    zeros[i] = Complex(0.0, 0.0);

  // \todo adjust gain
}

void PoleZeroMapper::sPlanePrototypeToHighpass(Complex *prototypePoles, Complex *prototypeZeros, Complex *targetPoles, 
                                                  Complex *targetZeros, int prototypeOrder, double targetCutoff)
{
  for(int i=0; i<prototypeOrder; i++)
  {
    targetPoles[i] = targetCutoff / prototypePoles[i];
    if( prototypeZeros[i].isInfinite() )
      targetZeros[i] = Complex(0.0, 0.0);
    else
      targetZeros[i] = targetCutoff / prototypeZeros[i];
  }
}

void PoleZeroMapper::prototypeToAnalogHighShelv(Complex *poles, int numPoles, Complex *zeros, int numZeros, double *gain, 
                                                   double targetCutoff)
{
  int i;
  for(i=0; i<numPoles; i++)
    poles[i] = zeros[i] * targetCutoff;
  for(i=0; i<numZeros; i++)
    zeros[i] = poles[i] * targetCutoff;

  // \todo adjust gain - huh= this is not supposed to work when overwriting the poles in the first array
}

void PoleZeroMapper::sPlanePrototypeToHighShelv(Complex *prototypePoles, Complex *prototypeZeros, Complex *targetPoles, 
                                                   Complex *targetZeros, int prototypeOrder, double targetCutoff)
{
  for(int i=0; i<prototypeOrder; i++)
  {
    targetPoles[i] = prototypeZeros[i] * targetCutoff;
    targetZeros[i] = prototypePoles[i] * targetCutoff;
  }
}

void PoleZeroMapper::prototypeToAnalogBandpass(Complex *poles, int numPoles, Complex *zeros, int numZeros, double *gain,
                                                  double targetLowCutoff, double targetHighCutoff)
{
  double wc = sqrt(targetLowCutoff*targetHighCutoff); // center (radian) frequency
  double bw = targetHighCutoff - targetLowCutoff;     // bandwidth
  Complex* tmpPoles = new Complex[2*numPoles];
  Complex* tmpZeros = new Complex[2*numPoles];

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the bandpass-filter:
  int     i;
  Complex tmp1, tmp2;
  for(i=0; i<numPoles; i++)
  {
    tmp1                     = 0.5 * bw * poles[i];
    tmp2                     = -bw * poles[i];
    tmp2                     = 0.25 * tmp2*tmp2;
    tmp2                     = sqrtC(tmp2 - wc*wc);
    //tmpPoles[2*i]   = tmp1 + tmp2;
    //tmpPoles[2*i+1] = tmp1 - tmp2;
    tmpPoles[i]              = tmp1 + tmp2;
    tmpPoles[2*numPoles-i-1] = tmp1 - tmp2;
  }
  for(i=0; i<numZeros; i++)
  {
    tmp1                     = 0.5 * bw * zeros[i];
    tmp2                     = -bw * zeros[i];
    tmp2                     = 0.25 * tmp2*tmp2;
    tmp2                     = sqrtC(tmp2 - wc*wc);
    //tmpZeros[2*i]   = tmp1 + tmp2;
    //tmpZeros[2*i+1] = tmp1 - tmp2;
    tmpZeros[i]              = tmp1 + tmp2;
    tmpZeros[2*numPoles-i-1] = tmp1 - tmp2;
  }

  // put additional zeros at 0.0 and infinity when there are less finite zeros than finite poles:
  //for(i=numZeros; i<(numPoles-numZeros); i++)
  for(i=numZeros; i<numPoles; i++)
  {
    //tmpZeros[2*i]   = 0.0;
    //tmpZeros[2*i+1] = INF;
    tmpZeros[i]              = 0.0;
    tmpZeros[2*numPoles-i-1] = INF;
  }

  //int         order = 2 * max(numPoles, numZeros);

  // copy the content of the temporary arrays into the output arrays:
  for(i=0; i<2*numPoles; i++)
  {
    zeros[i] = tmpZeros[i];
    poles[i] = tmpPoles[i];
  }

  // \todo: adjust gain

  // for debug:
  Complex pDbg[40];
  Complex zDbg[40];
  int k;
  for(k=0; k<2*numPoles; k++)
  {
    pDbg[k] = poles[k];
    zDbg[k] = zeros[k];
  }

  // free dynamically allocate memory:
  delete[] tmpPoles;
  delete[] tmpZeros;
}

void PoleZeroMapper::sPlanePrototypeToBandpass(Complex *prototypePoles, Complex *prototypeZeros, Complex *targetPoles, 
                                                  Complex *targetZeros, int prototypeOrder, double targetLowerCutoff, 
                                                  double targetUpperCutoff)
{
  double wc = sqrt(targetLowerCutoff*targetUpperCutoff); // center (radian) frequency
  double bw = targetUpperCutoff - targetLowerCutoff;     // bandwidth

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the bandpass-filter:
  Complex tmp1, tmp2;
  int k;
  for(k=0; k<prototypeOrder; k++)
  {
    // transform poles:
    tmp1 = 0.5 * bw * prototypePoles[k];
    tmp2 = -bw * prototypePoles[k];
    tmp2 = 0.25 * tmp2*tmp2;
    tmp2 = sqrtC(tmp2 - wc*wc);
    targetPoles[2*k]   = tmp1 + tmp2;
    targetPoles[2*k+1] = tmp1 - tmp2;

    // transform zeros:
    if( prototypeZeros[k].isInfinite() )
    {
      targetZeros[2*k]   = INF;
      targetZeros[2*k+1] = 0.0;
    }
    else
    {
      tmp1 = 0.5 * bw * prototypeZeros[k];
      tmp2 = -bw * prototypeZeros[k];
      tmp2 = 0.25 * tmp2*tmp2;
      tmp2 = sqrtC(tmp2 - wc*wc);
      targetZeros[2*k]   = tmp1 + tmp2;
      targetZeros[2*k+1] = tmp1 - tmp2;
    }
  }

  // re-arrange poles and zeros such that complex conjugate pairs (again) occupy successive slots:
  k = 1;
  while( k+1 < 2*prototypeOrder )
  {
    swap(targetPoles[k], targetPoles[k+1]);
    swap(targetZeros[k], targetZeros[k+1]);
    k += 4;
  }

  std::vector<Complex> pDbg, zDbg; // for debugging
  pDbg = toVector(targetPoles, 2*prototypeOrder);
  zDbg = toVector(targetZeros, 2*prototypeOrder);
  int dummy = 0;
}

void PoleZeroMapper::prototypeToAnalogBandstop(Complex *poles, int numPoles, Complex *zeros, int numZeros, double *gain, 
                                                  double targetLowCutoff, double targetHighCutoff)
{
  double wc = sqrt(targetLowCutoff*targetHighCutoff); // center (radian) frequency
  double bw = targetHighCutoff - targetLowCutoff;     // bandwidth

  int         order = 2 * rmax(numPoles, numZeros);
  Complex* tmpPoles = new Complex[order];
  Complex* tmpZeros = new Complex[order];

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the bandpass-filter:
  int     i;
  Complex tmp1, tmp2;
  for(i=0; i<numPoles; i++)
  {
    tmp1            = 0.5*bw / poles[i];
    tmp2            = -bw / poles[i];
    tmp2            = 0.25 * tmp2*tmp2;
    tmp2            = sqrtC(tmp2 - wc*wc);
    tmpPoles[2*i]   = tmp1 + tmp2;
    tmpPoles[2*i+1] = tmp1 - tmp2;
  }
  for(i=0; i<numZeros; i++)
  {
    tmp1            = 0.5*bw / zeros[i];
    tmp2            = -bw / zeros[i];
    tmp2            = 0.25 * tmp2*tmp2;
    tmp2            = sqrtC(tmp2 - wc*wc);
    tmpZeros[2*i]   = tmp1 + tmp2;
    tmpZeros[2*i+1] = tmp1 - tmp2;
  }

  // put additional zeros at s=j*wc and s=-j*wc when there are less finite zeros than finite poles:
  Complex j(0.0, 1.0); // imaginary unit
  for(i=numZeros; i<(numPoles-numZeros); i++)
  {
    tmpZeros[2*i]   =  j*wc;
    tmpZeros[2*i+1] = -j*wc;
  }

  // copy the content of the temporary arrays into the output arrays:
  for(i=0; i<order; i++)
  {
    zeros[i] = tmpZeros[i];
    poles[i] = tmpPoles[i];
  }

  // free dynamically allocate memory:
  delete[] tmpPoles;
  delete[] tmpZeros;

  // \todo: adjust gain
}

void PoleZeroMapper::sPlanePrototypeToBandreject(Complex *prototypePoles, Complex *prototypeZeros, Complex *targetPoles, 
                                                    Complex *targetZeros, int prototypeOrder, double targetLowerCutoff, 
                                                    double targetUpperCutoff)
{
  double wc = sqrt(targetLowerCutoff*targetUpperCutoff); // center (radian) frequency
  double bw = targetUpperCutoff - targetLowerCutoff;     // bandwidth

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the bandpass-filter:
  Complex tmp1, tmp2;  
  int k;
  for(k=0; k<prototypeOrder; k++)
  {
    // transform poles:
    tmp1 = 0.5 * bw / prototypePoles[k];
    tmp2 = -bw / prototypePoles[k];
    tmp2 = 0.25 * tmp2*tmp2;
    tmp2 = sqrtC(tmp2 - wc*wc);
    targetPoles[2*k]   = tmp1 + tmp2;
    targetPoles[2*k+1] = tmp1 - tmp2;

    // transform zeros:
    Complex j(0.0, 1.0); // imaginary unit
    if( prototypeZeros[k].isInfinite() )
    {
      targetZeros[2*k]   =  j*wc;
      targetZeros[2*k+1] = -j*wc;
    }
    else
    {
      tmp1 = 0.5 * bw / prototypeZeros[k];
      tmp2 = -bw / prototypeZeros[k];
      tmp2 = 0.25 * tmp2*tmp2;
      tmp2 = sqrtC(tmp2 - wc*wc);
      targetZeros[2*k]   = tmp1 + tmp2;
      targetZeros[2*k+1] = tmp1 - tmp2;
    }
  }

  // re-arrange poles and zeros such that complex conjugate pairs (again) occupy successive slots:
  k = 1;
  while( k+1 < 2*prototypeOrder )
  {
    swap(targetPoles[k], targetPoles[k+1]);
    int kp = (k-1)/4;
    if( !prototypeZeros[kp].isInfinite() )
      swap(targetZeros[k], targetZeros[k+1]); // don't want to swap the generated zeros at +-j*wc
    k += 4;
  }
}

void PoleZeroMapper::bilinearAnalogToDigital(Complex *poles, int numPoles, Complex *zeros, int numZeros, double sampleRate, 
                                                double *gain)
{
  int     i;
  Complex z;
  double  scaler = 0.5/sampleRate;

  for(i=0; i<numPoles; i++)
  {
    // we don't check against infinity (as we do with zeros) because infinite poles should actually never occur in practice
    z        = scaler * poles[i];
    poles[i] = (1.0+z)/(1.0-z);
  }

  for(i=0; i<numZeros; i++)
  {
    if( zeros[i].isInfinite() )
      zeros[i] = Complex(-1.0, 0.0);
    else
    {
      z        = scaler * zeros[i];
      zeros[i] = (1.0+z)/(1.0-z);
    }
  }
}

 


void PoleZeroMapper::zLowpassToLowpass(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, double *kNew, int N,
                                       double wc, double wt)
{
  // doo stufff....


  int dummy = 0;
}
