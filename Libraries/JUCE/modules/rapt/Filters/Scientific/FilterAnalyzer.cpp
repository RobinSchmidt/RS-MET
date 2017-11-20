template<class T>
std::complex<T> rsFilterAnalyzer<T>::getAnalogFrequencyResponseAt(Complex* z, Complex* p, T k, 
  int N, T w)
{
  Complex s = Complex(0.0, w);
  Complex n = rsPolynomial<T>::evaluatePolynomialWithRoots(s, z, N); 
  Complex d = rsPolynomial<T>::evaluatePolynomialWithRoots(s, p, N); 
  return k * n/d;
}

template<class T>
T rsFilterAnalyzer<T>::getAnalogMagnitudeResponseAt(Complex* z, Complex* p, T k, int N, T w)
{
  return abs(getAnalogFrequencyResponseAt(z, p, k, N, w));
  // maybe it can be computed more efficiently?
}

template<class T>
void rsFilterAnalyzer<T>::getAnalogMagnitudeResponse(Complex* z, Complex* p, T k, int N, T* w, 
  T* m, int numBins)
{
  for(int i = 0; i < numBins; i++)
    m[i] = getAnalogMagnitudeResponseAt(z, p, k, N, w[i]);
}

template<class T>
void rsFilterAnalyzer<T>::getAnalogPhaseResponse(Complex* z, Complex* p, T k, int N, T* w, T* phs, 
  int numBins)
{
  for(int i = 0; i < numBins; i++)
    phs[i] = getAnalogFrequencyResponseAt(z, p, k, N, w[i]).getAngle();

  // \todo: unwrap the phase
}

template<class T>
T rsFilterAnalyzer<T>::findAnalogFrequencyWithMagnitude(Complex* z, Complex* p, T* k, int N, T A, 
  T initialGuess)
{
  //if( A <= 0.0 || A >= 1.0 )
  //  rosic::error("A must be between 0 and 1 (exclusive) in findPrototypeCutoff");

  // determine the frequency interval in which the magnitude value A must fall (maybe wrap into a
  // function):
  T wL = initialGuess;
  T wU = initialGuess;
  while( rsFilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wL) < A )
    wL = wL/2;
  while( rsFilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wU) > A )
    wU = wU*2;

  // find the frequency at which the magnitude is equal to A by bisection:  
  const T eps  = std::numeric_limits<T>::epsilon(); 
  int i       = 0;         // iteration counter
  int iMax    = 1000;      // maximum number of iterations
  T tol  = A*eps;     // amplitude tolerance - maybe use A*eps?
  T wM   = (wL+wU)/2; // midpoint
  T wOld = wM;        // to remember wM from previous iteration
  T AM   = rsFilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wM);
  while( abs(AM-A) > tol && i < iMax ) 
  {
    if( AM > A )
      wL = wM;
    else if( AM < A )
      wU = wM;
    wOld = wM;
    wM   = (wL+wU)/2;
    AM = rsFilterAnalyzer::getAnalogMagnitudeResponseAt(z, p, *k, N, wM);
    i++;

    if( abs(wM-wOld) < eps * rsMin(wOld, wM) ) // 2nd convergence criterion based on frequency difference - maybe use as only criterion
      break;  
    //if( abs(wM-wOld) < eps * std::min(wOld, wM) ) // 2nd convergence cirterion based on frequency difference - maybe use as only criterion
    //  break;  
  }

  if(i >= iMax)
    RS_DEBUG_BREAK;
    //rsError("Slow convergence in FilterAnalyzer::findAnalogFrequencyWithMagnitude");

  return wM;
}

template<class T>
T rsFilterAnalyzer<T>::getBiquadMagnitudeAt(const T b0, const T b1, const T b2, const T a1, 
  const T a2, const T w)
{
  T c1  = cos(w);
  T c2  = 2.0*c1*c1-1.0;  // == cos(2.0*w) by virtue of the identity cos(2*x) = 2*cos^2(x)-1
  T num = b0*b0 + b1*b1 + b2*b2 + 2.0 * (b0*b1 + b1*b2) * c1 + 2.0 * b0*b2*c2;
  T den = 1.0   + a1*a1 + a2*a2 + 2.0 * (   a1 + a1*a2) * c1 + 2.0 *    a2*c2;
  return sqrt(num/den);
}

template<class T> 
void rsFilterAnalyzer<T>::getBiquadMagnitudeResponse(const T b0, const T b1, const T b2, 
  const T a1, const T a2, T* w, T* mag, int numBins, bool inDecibels)
{
  for(int k = 0; k < numBins; k++)
    mag[k] = getBiquadMagnitudeAt(b0, b1, b2, a1, a2, w[k]);
  if(inDecibels)
    convertToDecibels(mag, numBins);
}

template<class T>
std::complex<T> rsFilterAnalyzer<T>::getBiquadTransferFunctionAt(const T b0, const T b1, 
  const T b2, const T a1, const T a2, const Complex z)
{
  Complex z1 = 1 / z;  // z^-1
  Complex z2 = z1*z1;  // z^-2
  return (b0 + b1*z1 + b2*z2) / (1 + a1*z1 + a2*z2);
}

template<class T>
std::complex<T> rsFilterAnalyzer<T>::getBiquadCascadeTransferFunctionAt(T* b0, T* b1, T* b2, T* a1,
  T* a2, int numBiquads, Complex z)
{
  Complex H = 1.0;
  for(int i = 0; i < numBiquads; i++)
    H *= getBiquadTransferFunctionAt(b0[i], b1[i], b2[i], a1[i], a2[i], z);
  return H;
}

template<class T>
void rsFilterAnalyzer<T>::getBiquadCascadeFrequencyResponse(T* b0, T* b1, T* b2, T* a1, T* a2, 
  int numBiquads, T* w, Complex* H, int numBins, int accumulationMode)
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

template<class T>
void rsFilterAnalyzer<T>::multiplyWithBiquadCascadeFrequencyResponse(T* b0, T* b1, T* b2, T* a1, 
  T* a2, int numBiquads, T* w, Complex* H, int numBins)
{
  Complex j(0.0, 1.0);
  for(int k = 0; k < numBins; k++)
  {
    Complex z = expC(j*w[k]);
    H[k] *= getBiquadCascadeTransferFunctionAt(b0, b1, b2, a1, a2, numBiquads, z);
  }
}

template<class T> 
void rsFilterAnalyzer<T>::addWithBiquadCascadeFrequencyResponse(T* b0, T* b1, T* b2, T* a1, T* a2, 
  int numBiquads, T* w, Complex* H, int numBins)
{
  Complex j(0.0, 1.0);
  for(int k=0; k<numBins; k++)
  {
    Complex z = expC(j*w[k]);
    H[k] += getBiquadCascadeTransferFunctionAt(b0, b1, b2, a1, a2, numBiquads, z);
  }
}

template<class T>
void rsFilterAnalyzer<T>::getBiquadCascadeMagnitudeResponse(T* b0, T* b1, T* b2, T* a1, T* a2, 
  int numBiquads, T* w, T* mag, int numBins, bool inDecibels, bool accumulate)
{
  Complex *H   = new Complex[numBins];
  T  *tmp = new T[numBins];

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

template<class T>
void rsFilterAnalyzer<T>::getBiquadCascadeMagnitudeResponse(T* b0, T* b1, T* b2, T* a1, T* a2, 
  int numBiquads, T* frequencies, T sampleRate, T* magnitudes, int numBins, bool inDecibels, 
  bool accumulate)
{
  T *w = new T[numBins];
  T scaler = 2*PI / sampleRate;
  for(int k = 0; k < numBins; k++)
    w[k] = scaler * frequencies[k];
  getBiquadCascadeMagnitudeResponse(b0, b1, b2, a1, a2, numBiquads, w, magnitudes, numBins, 
    inDecibels, accumulate);
  delete[] w;
}

template<class T>
T rsFilterAnalyzer<T>::getBiquadPhaseResponseAt(T b0, T b1, T b2, T a1, T a2, T w)
{
  Complex z = expC(Complex(0.0, w));  // z = e^(j*w)
  Complex H = getBiquadTransferFunctionAt(b0, b1, b2, a1, a2, z);  
  T  phase = H.getAngle();
  if( phase > 0.0 )
    phase -= 2.0*PI;
  return phase; 
}

template<class T>
T rsFilterAnalyzer<T>::getBiquadCascadePhaseResponseAt(T* b0, T* b1, T* b2, T* a1, T* a2, 
  int numBiquads, T w)
{
  T phase = 0.0;
  for(int i = 0; i < numBiquads; i++)
    phase += getBiquadPhaseResponseAt(b0[i], b1[i], b2[i], a1[i], a2[i], w);
  return phase;
}

template<class T>
void rsFilterAnalyzer<T>::getBiquadCascadePhaseResponse(T* b0, T* b1, T* b2, T* a1, T *a2, 
  int numBiquads, T *w, T *phases, int numBins, bool accumulate)
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

template<class T>
void rsFilterAnalyzer<T>::getMagnitudes(Complex* H, T* magnitudes, int length)
{
  for(int k = 0; k < length; k++)
    magnitudes[k] = H[k].getRadius();
}

template<class T>
void rsFilterAnalyzer<T>::getPhases(Complex* H, T* phases, int length)
{
  for(int k = 0; k < length; k++)
    phases[k] = H[k].getAngle();
  // \todo unwrap(phases, length, 2*PI); maybe conditionally
}

template<class T>
void rsFilterAnalyzer<T>::convertToDecibels(T* values, int length, T clipLowAmplitudeAt)
{
  for(int k = 0; k < length; k++)
    values[k] = amp2dBWithCheck(values[k], clipLowAmplitudeAt); 
}

template<class T>
void rsFilterAnalyzer<T>::clampValuesAboveNyquist(T* frequencies, T* values, int length, 
  T sampleRate, T clampValue)
{
  for(int k = 0; k < length; k++)
  {
    if( frequencies[k] > 0.5*sampleRate )
      values[k] = clampValue;
  }
}
