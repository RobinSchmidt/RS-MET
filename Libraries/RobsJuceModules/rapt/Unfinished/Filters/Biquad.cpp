//-------------------------------------------------------------------------------------------------

template<class TCof>
rsBiquad<TCof>::rsBiquad()
{
  initializeCoefficients();             
}

template<class TCof>
rsBiquad<TCof>::~rsBiquad()
{

}

template<class TCof>
void rsBiquad<TCof>::setCoefficients(TCof newB0, TCof newB1, TCof newB2, TCof newA1, TCof newA2)
{
  b0 = newB0;
  b1 = newB1;
  b2 = newB2;
  a1 = newA1;
  a2 = newA2;
}

template<class TCof>
void rsBiquad<TCof>::initializeCoefficients()
{
  setCoefficients(1.0, 0.0, 0.0, 0.0, 0.0);
}

//-------------------------------------------------------------------------------------------------
// class rsBiquadDF1:

template<class TSig, class TCof>
rsBiquadDF1<TSig, TCof>::rsBiquadDF1()
{
  reset();                
}

template<class TSig, class TCof>
rsBiquadDF1<TSig, TCof>::~rsBiquadDF1()
{

}

template<class TSig, class TCof>
void rsBiquadDF1<TSig, TCof>::reset()
{
  x1 = x2 = y1 = y2 = 0.0;
}


//-------------------------------------------------------------------------------------------------
// class rsBiquadDesigner:

template<class T>
T rsBiquadDesigner::getBiquadMagnitudeAt(const T &b0, const T &b1, const T &b2, 
  const T &a1, const T &a2, const T &frequency, const T &sampleRate)
{
  T omega = 2.0 * PI * frequency / sampleRate; // optimize to one mul
  T c1  = cos(omega);
  T c2  = cos(2.0*omega);

  T a1m  = -a1; // we use the other sign-convention in the getSample-function
  T a2m  = -a2;

  T num = b0*b0 + b1*b1   + b2*b2   + 2.0*(b0*b1 + b1*b2)  *c1 + 2.0*b0*b2*c2;
  T den = 1.0   + a1m*a1m + a2m*a2m + 2.0*(  a1m + a1m*a2m)*c1 + 2.0*  a2m*c2;
  T mag = rsSqrt(num/den);

  return mag;
}

template<class T>
void rsBiquadDesigner::calculateFirstOrderLowpassCoeffsPrescribedNyquist(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &sampleRate, const T &frequency)
{
  T wc = 2.0*PI*frequency/sampleRate;        // normalized radian cutoff frequency
  T Wp = tan(wc/2.0);                        // pre-warped analog cutoff frequency
  T Wn = PI*sampleRate;                      // radian Nyquist frequency
  T Wc = 2.0*PI*frequency;                   // non-pre-warped analog cutoff frequency
  T k2 = 1.0 / (1.0 + ((Wn*Wn)/(Wc*Wc)) );   // gain of prototype at the Nyquist frequency

  // compute analog filter coefficients:
  T A2  = 1.0 / (Wp*Wp*(1.0-2.0*k2));        // A^2
  T B2  = k2*A2;                             // B^2
  T G02 = 1.0;                               // G0^2
  T A   = rsSqrt(A2);
  T B   = rsSqrt(B2);
  T G0  = rsSqrt(G02);                       // == 1.0 -> optimize out
  T rD  = 1.0 / (1.0+A);                     // reciprocal of denominator

  // compute digital filter coefficients:
  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1 -A); 
  a2 =  0.0;
}

template<class T>
void rsBiquadDesigner::calculateFirstOrderHighpassCoeffsPrescribedNyquist(T &b0, T &b1,
  T &b2, T &a1, T &a2, const T &sampleRate, const T &frequency)
{
  T wc = 2.0*PI*frequency/sampleRate;
  T Wp = tan(wc/2.0);
  T Wn = PI*sampleRate;
  T Wc = 2.0*PI*frequency;  

  T r   = (Wn*Wn) / (Wc*Wc);
  T k2  = r / (1+r);
  T A2  = -1 / (Wp*Wp*(1-2*k2));
  T B2  = k2*A2;
  T G02 = 0.0;

  T A  = rsSqrt(A2);
  T B  = rsSqrt(B2);
  T G0 = rsSqrt(G02);
  T rD = 1 / (1+A);  

  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1 -A); 
  a2 =  0.0;
}

template<class T>
void rsBiquadDesigner::calculateFirstOrderHighShelvCoeffsPrescribedNyQuist(T& b0, T& b1, 
  T& b2, T& a1, T& a2, const T& sampleRate, const T& frequency,                                                                          
  const T& gainFactor)
{
  // todo: catch special case when gainFactor == 1

  T wc = 2.0*PI*frequency/sampleRate;
  T Wp = tan(wc/2.0);
  T Wn = PI*sampleRate;
  T Wc = 2.0*PI*frequency;  

  // protoype coeffs (unmodified):
  T g0 = 1;
  T a  = 1 / (Wc*rsSqrt(gainFactor));
  T b  = gainFactor*a;

  // desired Nyquist frequency gain (squared):
  T k2 = (g0*g0 + b*b*Wn*Wn) / (1 + a*a*Wn*Wn);

  // modified analog filter coeffs:
  T G0 = g0;
  T A2 = (gainFactor-G0*G0) / (Wp*Wp*(k2-gainFactor));
  T B2 = k2*A2;
  T A  = rsSqrt(A2);
  T B  = rsSqrt(B2);

  // reciprocal of denominator:
  T rD =  1.0 / (1.0+A); 

  // digital filter coeffs:
  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1.0 -A);
  a2 =  0.0;
}

template<class T>
void rsBiquadDesigner::calculateFirstOrderLowShelvCoeffsPrescribedNyQuist(T& b0, T& b1,                                                                         
  T& b2, T& a1, T& a2, const T& sampleRate, const T& frequency, 
  const T& gainFactor)
{  
  // todo: catch special case when gainFactor == 1

  T wc = 2.0*PI*frequency/sampleRate;
  T Wp = tan(wc/2.0);
  T Wn = PI*sampleRate;
  T Wc = 2.0*PI*frequency;  

  // protoype coeffs (unmodified):
  T g0 = gainFactor;
  T a  = rsSqrt(gainFactor) / Wc;
  T b  = a;

  // desired Nyquist frequency gain (squared):
  T k2 = (g0*g0 + b*b*Wn*Wn) / (1 + a*a*Wn*Wn);

  // modified analog filter coeffs:
  T G0 = g0;
  T A2 = (gainFactor-G0*G0) / (Wp*Wp*(k2-gainFactor));
  T B2 = k2*A2;
  T A  = rsSqrt(A2);
  T B  = rsSqrt(B2);

  // reciprocal of denominator:
  T rD =  1.0 / (1.0+A); 

  // digital filter coeffs:
  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1.0-A);
  a2 =  0.0;
}
