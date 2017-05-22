//#include <limits>
//
//#include "rosic_FilterAnalyzer.h"
//#include "../math/rosic_PolynomialAlgorithms.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// filter analysis functions:

Complex FilterAnalyzer::getAnalogFrequencyResponseAt(Complex *z, Complex *p, double k, int N, double w)
{
  Complex s = Complex(0.0, w);
  Complex n = evaluatePolynomialWithRoots(s, z, N); 
  Complex d = evaluatePolynomialWithRoots(s, p, N); 
  return k * n/d;
}

double FilterAnalyzer::getAnalogMagnitudeResponseAt(Complex *z, Complex *p, double k, int N, double w)
{
  return getAnalogFrequencyResponseAt(z, p, k, N, w).getRadius();
  // maybe it can be computed more efficiently?
}

void FilterAnalyzer::getAnalogMagnitudeResponse(Complex *z, Complex *p, double k, int N, double *w, double *m, int numBins)
{
  for(int i = 0; i < numBins; i++)
    m[i] = getAnalogMagnitudeResponseAt(z, p, k, N, w[i]);
}

void FilterAnalyzer::getAnalogPhaseResponse(Complex *z, Complex *p, double k, int N, double *w, double *phs, int numBins)
{
  for(int i = 0; i < numBins; i++)
    phs[i] = getAnalogFrequencyResponseAt(z, p, k, N, w[i]).getAngle();

  // \todo: unwrap the phase
}

double FilterAnalyzer::findAnalogFrequencyWithMagnitude(Complex *z, Complex *p, double *k, int N, double A, double initialGuess)
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
  double wOld = wM;        // to remember wM from previous iteration
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

    if( abs(wM-wOld) < eps * min(wOld, wM) ) // 2nd convergence cirterion based on frequency difference - maybe use as only criterion
      break;  
    //if( abs(wM-wOld) < eps * std::min(wOld, wM) ) // 2nd convergence cirterion based on frequency difference - maybe use as only criterion
    //  break;  
  }

  if( i >= iMax )
    rosic::error("Slow convergence in FilterAnalyzer::findAnalogFrequencyWithMagnitude");

  return wM;
}













double FilterAnalyzer::getBiquadMagnitudeAt(const double b0, const double b1, const double b2, const double a1, const double a2, 
                                            const double w)
{
  double c1  = cos(w);
  double c2  = 2.0*c1*c1-1.0;  // == cos(2.0*w) by virtue of the identity cos(2*x) = 2*cos^2(x)-1
  double num = b0*b0 + b1*b1 + b2*b2 + 2.0 * (b0*b1 + b1*b2) * c1 + 2.0 * b0*b2*c2;
  double den = 1.0   + a1*a1 + a2*a2 + 2.0 * (   a1 + a1*a2) * c1 + 2.0 *    a2*c2;
  return sqrt(num/den);
}
    
void FilterAnalyzer::getBiquadMagnitudeResponse(const double b0, const double b1, const double b2, const double a1, const double a2,     
                                                double *w, double *mag, int numBins, bool inDecibels)
{
  for(int k=0; k<numBins; k++)
    mag[k] = getBiquadMagnitudeAt(b0, b1, b2, a1, a2, w[k]);
  if(inDecibels)
    convertToDecibels(mag, numBins);
}

Complex FilterAnalyzer::getBiquadTransferFunctionAt(const double b0, const double b1, const double b2, const double a1, const double a2, 
                                                    const Complex z)
{
  Complex z1 = 1 / z;  // z^-1
  Complex z2 = z1*z1;  // z^-2
  return (b0 + b1*z1 + b2*z2) / (1 + a1*z1 + a2*z2);
}

Complex FilterAnalyzer::getBiquadCascadeTransferFunctionAt(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads,                    
                                                           Complex z)
{
  Complex H = 1.0;
  for(int i=0; i<numBiquads; i++)
    H *= getBiquadTransferFunctionAt(b0[i], b1[i], b2[i], a1[i], a2[i], z);
  return H;
}

void FilterAnalyzer::getBiquadCascadeFrequencyResponse(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads,
                                                       double *w, Complex *H, int numBins, int accumulationMode)
{
  if( accumulationMode == NO_ACCUMULATION )
  {
    fillWithValue(H, numBins, Complex(1.0, 0.0));
    multiplyWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  }
  else if( accumulationMode == MULTIPLICATIVE_ACCUMULATION )
    multiplyWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  else if( accumulationMode == ADDITIVE_ACCUMULATION )
    addWithBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  else
    DEBUG_BREAK; // invalid accumulation mode
}

void FilterAnalyzer::multiplyWithBiquadCascadeFrequencyResponse(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads,
                                                                double *w, Complex *H, int numBins)
{
  Complex j(0.0, 1.0);
  for(int k=0; k<numBins; k++)
  {
    Complex z = expC(j*w[k]);
    H[k] *= getBiquadCascadeTransferFunctionAt(b0, b1, b2, a1, a2, numBiquads, z);
  }
}
   
void FilterAnalyzer::addWithBiquadCascadeFrequencyResponse(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads,                                                               
                                                           double *w, Complex *H, int numBins)
{
  Complex j(0.0, 1.0);
  for(int k=0; k<numBins; k++)
  {
    Complex z = expC(j*w[k]);
    H[k] += getBiquadCascadeTransferFunctionAt(b0, b1, b2, a1, a2, numBiquads, z);
  }
}

void FilterAnalyzer::getBiquadCascadeMagnitudeResponse(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads,
                                                       double *w, double *mag, int numBins, bool inDecibels, bool accumulate)
{
  Complex *H   = new Complex[numBins];
  double  *tmp = new double[numBins];

  getBiquadCascadeFrequencyResponse(b0, b1, b2, a1, a2, numBiquads, w, H, numBins);
  getMagnitudes(H, tmp, numBins);

  if(inDecibels)
  {  
    convertToDecibels(tmp, numBins, 0.0);    
    if( accumulate == false )
      copyBuffer(tmp, mag, numBins);  
    else
      add(mag, tmp, mag, numBins);
    clipBuffer(mag, numBins, -200.0, INF);  // avoid negative infinities
  }
  else
  {
    if( accumulate == false )
      copyBuffer(tmp, mag, numBins);  
    else
      multiply(mag, tmp, mag, numBins);
  }

  delete[] tmp;
  delete[] H;
}

void FilterAnalyzer::getBiquadCascadeMagnitudeResponse(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads, 
                                                       double *frequencies, double sampleRate, double *magnitudes, int numBins,                       
                                                       bool inDecibels, bool accumulate)
{
  double *w = new double[numBins];
  double scaler = 2*PI / sampleRate;
  for(int k=0; k<numBins; k++)
    w[k] = scaler * frequencies[k];
  getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numBiquads, w, magnitudes, numBins, inDecibels, accumulate);
  delete[] w;
}

double FilterAnalyzer::getBiquadPhaseResponseAt(double b0, double b1, double b2, double a1, double a2, double w)
{
  Complex z     = expC(Complex(0.0, w));  // z = e^(j*w)
  Complex H     = getBiquadTransferFunctionAt(b0, b1, b2, a1, a2, z);  
  double  phase = H.getAngle();
  if( phase > 0.0 )
    phase -= 2.0*PI;
  return phase; 
}

double FilterAnalyzer::getBiquadCascadePhaseResponseAt(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads, double w)
{
  double phase = 0.0;
  for(int i = 0; i < numBiquads; i++)
    phase += getBiquadPhaseResponseAt(b0[i], b1[i], b2[i], a1[i], a2[i], w);
  return phase;
}

void FilterAnalyzer::getBiquadCascadePhaseResponse(double *b0, double *b1, double *b2, double *a1, double *a2, int numBiquads, double *w, 
                                                   double *phases, int numBins, bool accumulate)
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

void FilterAnalyzer::getMagnitudes(Complex *H, double *magnitudes, int length)
{
  for(int k=0; k<length; k++)
    magnitudes[k] = H[k].getRadius();
}

void FilterAnalyzer::getPhases(Complex *H, double *phases, int length)
{
  for(int k=0; k<length; k++)
    phases[k] = H[k].getAngle();
}

void FilterAnalyzer::convertToDecibels(double *values, int length, double clipLowAmplitudeAt)
{
  for(int k=0; k<length; k++)
    values[k] = amp2dBWithCheck(values[k], clipLowAmplitudeAt); 
}

void FilterAnalyzer::clampValuesAboveNyquist(double *frequencies, double *values, int length, double sampleRate, double clampValue)
{
  for(int k=0; k<length; k++)
  {
    if( frequencies[k] > 0.5*sampleRate )
      values[k] = clampValue;
  }
}
