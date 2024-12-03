#ifndef RAPT_STATEVARIABLEFILTER_H
#define RAPT_STATEVARIABLEFILTER_H

/** A zero delay feedback (ZDF) state variable filter (SVF). It offers all the frequency responses
from the RBJ biquad cookbook. The filter is parameterized in terms of the normalized radian 
frequency omega = 2*pi*frequency/sampleRate, the quality factor Q and for bell and shelving 
filters, the linear gain A. High values of Q generally mean "more resonance" or "narrower 
bandwidths". The filter is based on trapezoidal integration using TDF2 integrators. The filter
produces internally a lowpass, bandpass and highpass signal which are available at the same time.
You can produce these 3 signals using the getPartialOutputs() function. Alternatively, you can use
getSample() which produces a single signal that is a mix of these 3 signals using some mixing 
coefficients that are determined by the desired filter mode. Consider getSample() as the high level
API and getPartialOutputs() as a lower level API. Most of the time, client code will want to use 
getSample() but the 3 separate outputs are made available as well, just in case you want them. 

For more details about where all the formulas that we implement here come from, see:

  https://github.com/RobinSchmidt/RS-MET/blob/work/Notes/StateVariableFilter.txt

*/

template<class TSig, class TPar>       // Data types for signals and parameters
class rsStateVariableFilter
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setupMuted();
  void setupBypass();
  void setupLowpass(      TPar omega, TPar Q);
  void setupHighpass(     TPar omega, TPar Q);
  void setupBandpassSkirt(TPar omega, TPar Q);
  void setupBandpassPeak( TPar omega, TPar Q);
  void setupBandstop(     TPar omega, TPar Q);
  void setupAllpass(      TPar omega, TPar Q);
  void setupBell(         TPar omega, TPar Q, TPar A);
  void setupLowShelf(     TPar omega, TPar Q, TPar A);
  void setupHighShelf(    TPar omega, TPar Q, TPar A);

  /** Sets up the SVF coefficients such that the filter realizes the biquad transfer function 

    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2). 
  
  However, the formula that is implemented doesn't seem to work for all possible sets of biquad 
  coeffs. It works only when  (a1*a1 - a2*a2 - 2*a2 - 1) < 0, so it's not recommended to use it 
  blindly. It may fail. I have not yet figured out why that is and if something can be done to fix
  it. Maybe not. Maybe there's some inherent limitation in the SVF with regard to what biquads it 
  can simulate. Be careful! If it fails, it will trigger an rsAssert. */
  void setupFromBiquad(TPar b0, TPar b1, TPar b2, TPar a1, TPar a2);
  // Update:
  // Empirically, it seems to be the case that whenever a biquad is stable, then the formulas will
  // work, i.e. stability implies workability but not necessarily the other way around. That means,
  // there might be some unstable biquads that the SVF can also realize. Stability seems to be a
  // sufficient but not necessary condition for the formulas to work. So, for the filters that we 
  // usually care about, namely the stable ones, the formulas should be fine. Some more thorough 
  // research should be done on this, though. I once tried it with a billion random biquads. The 
  // unit test does only 1000 because it needs to be fast, but I once did the test with a billion
  // and it still passed.


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Produces the coefficients of an equivalent direct form biquad filter that implements the
  difference equation:

    y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

  and therefore has the transfer function:

    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)

  The biquad is equivalent in the sense that it has the same transfer function as this filter. */
  void convertToBiquad(TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2) const;

  /** Evaluates the magnitude response of this filter at a given normalized radian frequency 
  omega. This is useful for plotting it on a GUI. */
  TPar getMagnitudeAt(TPar omega) const;

  /** Evaluates the filter's z-domain transfer function H(z) value at the given value of z. */
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const;

  /** Produces the denominator coefficients of an equivalent direct form biquad filter. */
  void getBiquadDenominatorCoeffs(TPar* a1, TPar* a2) const;

  /** Produces the numerator coefficients of an equivalent direct form biquad filter. */
  void getBiquadNumeratorCoeffs(TPar* b0, TPar* b1, TPar* b2) const;

  /** Produces the numerator coefficients of the lowpass part of an equivalent direct form biquad 
  filter. */
  void getBiquadNumeratorCoeffsLP(TPar* b0, TPar* b1, TPar* b2) const;

  /** Produces the numerator coefficients of the bandpass part of an equivalent direct form biquad 
  filter. */
  void getBiquadNumeratorCoeffsBP(TPar* b0, TPar* b1, TPar* b2) const;

  /** Produces the numerator coefficients of the highpass part of an equivalent direct form biquad 
  filter. */
  void getBiquadNumeratorCoeffsHP(TPar* b0, TPar* b1, TPar* b2) const;


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Computes one sample at a time. Calls getPartialOutputs() and mixes the produced lowpass, 
  bandpass and highpass ouptuts according to the desired filter mode. */
  inline TSig getSample(TSig in);

  /** Returns the 3 outputs (lowpass, bandpass, highpass) of the core SVF in the output parameters 
  yL, yB, yH. */
  inline void getPartialOutputs(TSig in, TSig* yL, TSig* yB, TSig* yH);

  /** Resets the internal state. */
  void reset() { z1 = z2 = 0; }


protected:

  // State:
  TSig z1 = 0, z2 = 0;          // Integrator states. Maybe rename to u,v for consistency with derivation

  // Coeffs:
  TPar aL = 0, aB = 0, aH = 0;  // Mixing coeffs for lowpass, bandpass and highpass signals
  TPar g = 0;                   // Integrator gain
  TPar c = 0;                   // g + r (r is 2*R in Vadim's book, R is the damping coeff)
  TPar s = 1;                   // Scaler given by 1 / (1 + g*(g+r))

};

//-------------------------------------------------------------------------------------------------
// Implementation. Those functions that are typically called per sample are defined in the .h file 
// to facilitate inlining. The others are in the .cpp file.

// Setup:

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupMuted()
{
  // H(s) = 0

  g  = 0;
  c  = 0;
  s  = 0;
  aL = 0;
  aB = 0;
  aH = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupBypass()
{
  // H(s) = 1 

  g  = 0;
  c  = 0;
  s  = 1;
  aL = 0;
  aB = 0;
  aH = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupLowpass(TPar w, TPar Q)
{
  // H(s) = 1 / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 1;
  aB = 0;
  aH = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupHighpass(TPar w, TPar Q)
{
  // H(s) = s^2 / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 0;
  aB = 0;
  aH = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupBandpassSkirt(TPar w, TPar Q)
{
  // H(s) = s / (s^2 + s/Q + 1)   (constant skirt gain, peak gain = Q)

  TPar r = 1/Q;
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 0;
  aB = 1;
  aH = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupBandpassPeak(TPar w, TPar Q)
{
  // H(s) = (s/Q) / (s^2 + s/Q + 1)   (constant 0 dB peak gain)

  TPar r = 1/Q;
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 0;
  aB = r;
  aH = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupBandstop(TPar w, TPar Q)
{
  // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 1;
  aB = 0;
  aH = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupAllpass(TPar w, TPar Q)
{
  // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 1;
  aB = -r;
  aH = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupBell(TPar w, TPar Q, TPar A)
{
  // H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)

  TPar r = 1/(Q*A);
  g  = tan(0.5*w); 
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 1;
  aB = A*A*r;
  aH = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupLowShelf(TPar w, TPar Q, TPar A)
{
  // H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)

  TPar r = 1/Q;
  g  = tan(0.5*w) / sqrt(A);
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = A*A;
  aB = A*r;
  aH = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupHighShelf(TPar w, TPar Q, TPar A)
{
  // H(s) = A * (A*s^2 + (sqrt(A)/Q)*s + 1)/(s^2 + (sqrt(A)/Q)*s + A)

  TPar r = 1/Q;
  g  = tan(0.5*w) * sqrt(A);
  c  = g + r;
  s  = 1 / (1 + g*c);
  aL = 1;
  aB = A*r;
  aH = A*A;
}

// Processing:

template<class TSig, class TPar>
inline void rsStateVariableFilter<TSig, TPar>::getPartialOutputs(
  TSig in, TSig* yL, TSig* yB, TSig* yH)
{
  // Compute HP, BP and LP outputs:
  *yH = (in - c*z1 - z2) * s;            // == (in - (g+r)*z1 - z2) / (1 + g*(g+r))
  *yB = z1 + g * *yH; 
  *yL = z2 + g * *yB;

  // Update integrator state variables:
  z1 = 2 * *yB - z1;                     // Equivalent to: z1 += 2 * g * *yH
  z2 = 2 * *yL - z2;                     // Equivalent to: z2 += 2 * g * *yB
}

template<class TSig, class TPar>
inline TSig rsStateVariableFilter<TSig, TPar>::getSample(TSig in)
{
  TSig yL, yB, yH;
  getPartialOutputs(in, &yL, &yB, &yH);  // Produce LP, BP and HP signals
  return aH*yH + aB*yB + aL*yL;          // Mix them according to desired filter type
}

#endif