#ifndef RSLib_FilterAnalyzer_cpp
#define RSLib_FilterAnalyzer_cpp

namespace RSLib
{

  double onePoleMagnitudeAt(double b0, double b1, double a1, double w)
  {
    double c = 2*cos(w);
    double p = (b0*b0 + b1*b1 + b0*b1*c) / (1 + a1*a1 + a1*c);
    return rsSqrt(p);
  }

  double biquadMagnitudeAt(double b0, double b1, double b2, double a1, double a2, double w)
  {
    double c1  = 2*cos(w);
    double c2  = 2*cos(2*w);
    double num = b0*b0 + b1*b1 + b2*b2 + c1*(b0*b1 + b1*b2) + c2*b0*b2;
    double den = 1     + a1*a1 + a2*a2 + c1*(   a1 + a1*a2) + c2*   a2;
    double p   = num/den;
    return rsSqrt(p);
  }

  bool isBiquadStableAndMinimumPhase(double b0, double b1, double b2, double a1, double a2)
  {
    return areRootsOnOrInsideUnitCircle(b2, b1, b0) && areRootsOnOrInsideUnitCircle(a2, a1, 1);
    // reverse order, because coeffs multiply negative powers of z
  }

  double analogBiquadMagnitudeSquaredAt(double B0, double B1, double B2, double A0, 
    double A1, double A2, double w)
  {
    double w2  = w*w;   // w^2
    double w4  = w2*w2; // w^4
    double num = B0*B0 + (B1*B1-2*B0*B2)*w2 + B2*B2*w4;
    double den = A0*A0 + (A1*A1-2*A0*A2)*w2 + A2*A2*w4;
    return num/den;
  }

}

//-------------------------------------------------------------------------------------------------
// class rsFilterAnalyzer:

rsComplexDbl rsFilterAnalyzer::getAnalogFrequencyResponseAt(rsComplexDbl *z, rsComplexDbl *p, 
  double k, int N, double w)
{
  rsComplexDbl s = rsComplexDbl(0.0, w);
  rsComplexDbl n = evaluatePolynomialWithRoots(s, z, N); 
  rsComplexDbl d = evaluatePolynomialWithRoots(s, p, N); 
  return k * n/d;
}

double rsFilterAnalyzer::getAnalogMagnitudeResponseAt(rsComplexDbl *z, rsComplexDbl *p, double k, 
  int N, double w)
{
  return getAnalogFrequencyResponseAt(z, p, k, N, w).getRadius();
  // maybe it can be computed more efficiently?
}

void rsFilterAnalyzer::getAnalogMagnitudeResponse(rsComplexDbl *z, rsComplexDbl *p, double k, 
  int N, double *w, double *m, int numBins)
{
  for(int i = 0; i < numBins; i++)
    m[i] = getAnalogMagnitudeResponseAt(z, p, k, N, w[i]);
}

void rsFilterAnalyzer::getAnalogPhaseResponse(rsComplexDbl *z, rsComplexDbl *p, double k, int N, 
  double *w, double *phs, int numBins)
{
  for(int i = 0; i < numBins; i++)
    phs[i] = getAnalogFrequencyResponseAt(z, p, k, N, w[i]).getAngle();

  // \todo: unwrap the phase
}

double rsFilterAnalyzer::findAnalogFrequencyWithMagnitude(rsComplexDbl *z, rsComplexDbl *p, 
  double *k, int N, double A, double initialGuess)
{
  //if( A <= 0.0 || A >= 1.0 )
  //  rosic::error("A must be between 0 and 1 (exclusive) in findPrototypeCutoff");

  // determine the frequency interval in which the magnitude value A must fall (maybe wrap into a 
  // function):
  double wL = initialGuess;
  double wU = initialGuess;
  while( getAnalogMagnitudeResponseAt(z, p, *k, N, wL) < A )
    wL = wL/2;
  while( getAnalogMagnitudeResponseAt(z, p, *k, N, wU) > A )
    wU = wU*2;

  // find the frequency at which the magnitude is equal to A by bisection:  
  const double eps  = std::numeric_limits<double>::epsilon(); 
  int i       = 0;         // iteration counter
  int iMax    = 1000;      // maximum number of iterations
  double tol  = A*eps;     // amplitude tolerance - maybe use A*eps?
  double wM   = (wL+wU)/2; // midpoint
  double wOld = wM;        // to remember wM from previous iteration
  double AM   = getAnalogMagnitudeResponseAt(z, p, *k, N, wM);
  while( abs(AM-A) > tol && i < iMax ) 
  {
    if( AM > A )
      wL = wM;
    else if( AM < A )
      wU = wM;
    wOld = wM;
    wM   = (wL+wU)/2;
    AM   = getAnalogMagnitudeResponseAt(z, p, *k, N, wM);
    i++;

    // 2nd convergence cirterion based on frequency difference (maybe use as only criterion):
    if( abs(wM-wOld) < eps*std::min(wOld, wM) )                                                
      break;  
  }

  if( i >= iMax )
    rsError("Slow convergence in FilterAnalyzer::findAnalogFrequencyWithMagnitude");

  return wM;

  // todo: what we are doing here is bisection - maybe use class for root-finding without 
  // derivatives
}

double rsFilterAnalyzer::getBiquadMagnitudeAt(const double b0, const double b1, const double b2, 
  const double a1, const double a2, const double w)
{
  double c1  = cos(w);
  double c2  = 2.0*c1*c1-1.0;  // == cos(2.0*w) because cos(2*x) = 2*cos^2(x)-1
  double num = b0*b0 + b1*b1 + b2*b2 + 2.0 * (b0*b1 + b1*b2) * c1 + 2.0 * b0*b2*c2;
  double den = 1.0   + a1*a1 + a2*a2 + 2.0 * (   a1 + a1*a2) * c1 + 2.0 *    a2*c2;
  return sqrt(num/den);
}
    
void rsFilterAnalyzer::getBiquadMagnitudeResponse(const double b0, const double b1, 
  const double b2, const double a1, const double a2, double *w, double *mag, int numBins, 
  bool inDecibels)
{
  for(int k=0; k<numBins; k++)
    mag[k] = getBiquadMagnitudeAt(b0, b1, b2, a1, a2, w[k]);
  if(inDecibels)
    convertToDecibels(mag, numBins);
}

rsComplexDbl rsFilterAnalyzer::getBiquadTransferFunctionAt(const double b0, const double b1, 
  const double b2, const double a1, const double a2, const rsComplexDbl z)
{
  rsComplexDbl z1 = 1.0 / z;  // z^-1
  rsComplexDbl z2 = z1*z1;    // z^-2
  return (b0 + b1*z1 + b2*z2) / (1.0 + a1*z1 + a2*z2);
}

rsComplexDbl rsFilterAnalyzer::getBiquadCascadeTransferFunctionAt(double *b0, double *b1, 
  double *b2,double *a1, double *a2, int numBiquads, rsComplexDbl z)
{
  rsComplexDbl H = 1.0;
  for(int i=0; i<numBiquads; i++)
    H *= getBiquadTransferFunctionAt(b0[i], b1[i], b2[i], a1[i], a2[i], z);
  return H;
}

void rsFilterAnalyzer::getBiquadCascadeFrequencyResponse(double *b0, double *b1, double *b2, 
  double *a1, double *a2, int numBiquads, double *w, rsComplexDbl *H, int numBins, 
  int accumulationMode)
{
  if( accumulationMode == NO_ACCUMULATION )
  {
    rsFillWithValue(H, numBins, rsComplexDbl(1.0, 0.0));
    multiplyWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  }
  else if( accumulationMode == MULTIPLICATIVE_ACCUMULATION )
    multiplyWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  else if( accumulationMode == ADDITIVE_ACCUMULATION )
    addWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  else
    rsError("invalid accumulation mode");
}

void rsFilterAnalyzer::multiplyWithBiquadCascadeFrequencyResponse(double *b0, double *b1, 
  double *b2, double *a1, double *a2, int numBiquads, double *w, rsComplexDbl *H, int numBins)
{
  rsComplexDbl j(0.0, 1.0);
  for(int k=0; k<numBins; k++)
  {
    rsComplexDbl z = rsExpC(j*w[k]);
    H[k] *= getBiquadCascadeTransferFunctionAt(b0, b1, b2, a1, a2, numBiquads, z);
  }
}
   
void rsFilterAnalyzer::addWithBiquadCascadeFrequencyResponse(double *b0, double *b1, double *b2, 
  double *a1, double *a2, int numBiquads, double *w, rsComplexDbl *H, int numBins)
{
  rsComplexDbl j(0.0, 1.0);
  for(int k=0; k<numBins; k++)
  {
    rsComplexDbl z = rsExpC(j*w[k]);
    H[k] += getBiquadCascadeTransferFunctionAt(b0, b1, b2, a1, a2, numBiquads, z);
  }
}

void rsFilterAnalyzer::getBiquadCascadeMagnitudeResponse(double *b0, double *b1, double *b2, 
  double *a1, double *a2, int numBiquads, double *w, double *mag, int numBins, bool inDecibels, 
  bool accumulate)
{
  rsComplexDbl *H = new rsComplexDbl[numBins];
  double *tmp = new double[numBins];

  getBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  getMagnitudes(H, tmp, numBins);

  if(inDecibels)
  {  
    convertToDecibels(tmp, numBins, 0.0);    
    if( accumulate == false )
      rsCopyBuffer(tmp, mag, numBins);  
    else
      rsAdd(mag, tmp, mag, numBins);
    rsClipBuffer(mag, numBins, -200.0, rsInfDouble);  // avoid negative infinities
  }
  else
  {
    if( accumulate == false )
      rsCopyBuffer(tmp, mag, numBins);  
    else
      rsMultiply(mag, tmp, mag, numBins);
  }

  delete[] tmp;
  delete[] H;
}

void rsFilterAnalyzer::getBiquadCascadeMagnitudeResponse(double *b0, double *b1, double *b2,
  double *a1, double *a2, int numBiquads, double *frequencies, double sampleRate, 
  double *magnitudes, int numBins, bool inDecibels, bool accumulate)
{
  double *w = new double[numBins];
  double scaler = 2*PI / sampleRate;
  for(int k=0; k<numBins; k++)
    w[k] = scaler * frequencies[k];
  getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numBiquads, w, magnitudes, numBins, 
    inDecibels, accumulate);
  delete[] w;
}

double rsFilterAnalyzer::getBiquadPhaseResponseAt(double b0, double b1, double b2, double a1, 
  double a2, double w)
{
  rsComplexDbl z = rsExpC(rsComplexDbl(0.0, w));  // z = e^(j*w)
  rsComplexDbl H = getBiquadTransferFunctionAt(b0, b1, b2, a1, a2, z);  
  double  phase = H.getAngle();
  if( phase > 0.0 )
    phase -= 2.0*PI;
  return phase; 
}

double rsFilterAnalyzer::getBiquadCascadePhaseResponseAt(double *b0, double *b1, double *b2, 
  double *a1, double *a2, int numBiquads, double w)
{
  double phase = 0.0;
  for(int i = 0; i < numBiquads; i++)
    phase += getBiquadPhaseResponseAt(b0[i], b1[i], b2[i], a1[i], a2[i], w);
  return phase;
}

void rsFilterAnalyzer::getBiquadCascadePhaseResponse(double *b0, double *b1, double *b2, 
  double *a1, double *a2, int numBiquads, double *w, double *phases, int numBins, bool accumulate)
{
  if( accumulate == true )
  {
    for(int k = 0; k < numBins; k++)
      phases[k] += getBiquadCascadePhaseResponseAt(b0, b1, b2, a1, a2, numBiquads, w[k]);
  }
  else
  {
    for(int k = 0; k < numBins; k++)
      phases[k] = getBiquadCascadePhaseResponseAt(b0, b1, b2, a1, a2, numBiquads, w[k]);
  }
}

void rsFilterAnalyzer::getMagnitudes(rsComplexDbl *H, double *magnitudes, int length)
{
  for(int k=0; k<length; k++)
    magnitudes[k] = H[k].getRadius();
}

void rsFilterAnalyzer::getPhases(rsComplexDbl *H, double *phases, int length)
{
  for(int k=0; k<length; k++)
    phases[k] = H[k].getAngle();
}

void rsFilterAnalyzer::convertToDecibels(double *values, int length, double clipLowAmplitudeAt)
{
  for(int k=0; k<length; k++)
    values[k] = rsAmp2dBWithCheck(values[k], clipLowAmplitudeAt); 
}

void rsFilterAnalyzer::clampValuesAboveNyquist(double *frequencies, double *values, int length, 
  double sampleRate, double clampValue)
{
  for(int k=0; k<length; k++)
  {
    if( frequencies[k] > 0.5*sampleRate )
      values[k] = clampValue;
  }
}


#endif