#ifndef RAPT_STATEVARIABLEFILTERMYSTRAN_H
#define RAPT_STATEVARIABLEFILTERMYSTRAN_H


/** A zero delay feedback (ZDF) state variable filter (SVF). It offers all the frequency responses
from the RBJ biquad cookbook. The filter is parameterized in terms of the normalized radian 
frequency omega = 2*pi*frequency/sampleRate, the quality factor Q and for bell and shelving 
filters, the linear gain A. High values of Q generally mean "more resonance" or "narrower 
bandwidths". The filter is based on trapezoidal integration using TDF2 integrators.

It implements this idea:

  https://www.kvraudio.com/forum/viewtopic.php?p=8992653#p8992653

See comments in the .cpp file for some more details.  */

template<class TSig, class TPar>       // Data types for signals and parameters
class rsStateVariableFilterMystran
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

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


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  // ToDo: move these to .cpp file

  //void convertToBiquad(TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2);


  /** Produces the coefficients of an equivalent direct form biquad filter that implements the
  difference equation:

    y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

  and therefore has the transfer function:

    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)

  The biquad coefficients are useful for evaluating the transfer function for a given z which in
  turn is useful to evaluate the frequency response, for example, for plotting it on a GUI. Maybe
  the transfer function can also be evaluated in other ways directly from our coefficients her. 
  However, the implementation of getTransferFunctionAt() makes use of this conversion 
  internally. */
  void convertToBiquad(TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2)
  {
    getBiquadNumeratorCoeffs(b0, b1, b2);
    getBiquadDenominatorCoeffs(  a1, a2);
  }
  // ToDo: implement also a setupFromBiquad method


  void getBiquadDenominatorCoeffs(TPar* a1, TPar* a2)
  {
    TPar s = scl;
    TPar c = gpr;
    *a1 =  2*(c*g + g*g)*s - 2;
    *a2 = -2*(c*g - g*g)*s + 1;
    // Simplify: factor out g, create variable for the common subexpression
  }

  void getBiquadNumeratorCoeffs(TPar* b0, TPar* b1, TPar* b2);


  void getBiquadNumeratorCoeffsLP(TPar* b0, TPar* b1, TPar* b2);
  /*
  {
    TPar s = scl;
    TPar c = gpr;
    *b0 =   s*g*g;
    *b1 = 2*s*g*g;
    *b2 =   s*g*g;
  }
  */

  void getBiquadNumeratorCoeffsBP(TPar* b0, TPar* b1, TPar* b2);
  /*
  {
    TPar s = scl;
    TPar c = gpr;
    *b0 =  g*s;
    *b1 =  0;
    *b2 = -g*s;
  }
  */

  void getBiquadNumeratorCoeffsHP(TPar* b0, TPar* b1, TPar* b2);
  /*
  {
    TPar s = scl;
    TPar c = gpr;
    *b0 =  s;
    *b1 = -2*s;
    *b2 =  s;
  }
  */

  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z);
  /*
  {
    TPar b0, b1, b2, a1, a2;
    convertToBiquad(&b0, &b1, &b2, &a1, &a2);
    rsComplex<TPar> d = TPar(1)/z, d2 = d*d;                              // d = z^-1, d2 = z^-2
    rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (TPar(1) + a1*d + a2*d2);
    return H;
  }
  */




  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Computes one sample at a time. Calls getOutputs() and mixes the produced lowpass, bandpass
  and highpass ouptuts according to the desired filter mode. */
  inline TSig getSample(TSig in);

  /** Returns the 3 outputs (lowpass, bandpass, highpass) of the core SVF. */
  inline void getOutputs(TSig in, TSig* yL, TSig* yB, TSig* yH);

  /** Resets the internal state. */
  void reset() { z1 = z2 = 0; }


protected:





  // State:
  TSig z1 = 0, z2 = 0;

  // Coeffs:
  TPar aL = 0, aB = 0, aH = 0;  // Mixing coeffs for lowpass, bandpass and highpass signals
  TPar g   = 0;                 // Integrator gain
  TPar gpr = 0;                 // g + r (r is 2*R in Vadim's book, R is the damping coeff)
  TPar scl = 1;                 // Scaler given by 1 / (1 + g*(g+r));

};

// Setup:

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBypass()
{
  // H(s) = 1 

  g   = 0; 
  gpr = 0;
  scl = 1;
  aL  = 0;
  aB  = 0;
  aH  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupLowpass(TPar w, TPar Q)
{
  // H(s) = 1 / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 1;
  aB  = 0;
  aH  = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupHighpass(TPar w, TPar Q)
{
  // H(s) = s^2 / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 0;
  aB  = 0;
  aH  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBandpassSkirt(TPar w, TPar Q)
{
  // H(s) = s / (s^2 + s/Q + 1)   (constant skirt gain, peak gain = Q)

  TPar r = 1/Q;
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 0;
  aB  = 1;
  aH  = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBandpassPeak(TPar w, TPar Q)
{
  // H(s) = (s/Q) / (s^2 + s/Q + 1)   (constant 0 dB peak gain)

  TPar r = 1/Q;
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 0;
  aB  = r;
  aH  = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBandstop(TPar w, TPar Q)
{
  // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 1;
  aB  = 0;
  aH  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupAllpass(TPar w, TPar Q)
{
  // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)

  TPar r = 1/Q;
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 1;
  aB  = -r;
  aH  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBell(TPar w, TPar Q, TPar A)
{
  // H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)

  TPar r = 1/(Q*A);
  g   = tan(0.5*w); 
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 1;
  aB  = A*A*r;
  aH  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupLowShelf(TPar w, TPar Q, TPar A)
{
  // H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)

  TPar r = 1/Q;
  g   = tan(0.5*w) / sqrt(A);
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = A*A;
  aB  = A*r;
  aH  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupHighShelf(TPar w, TPar Q, TPar A)
{
  // H(s) = A * (A*s^2 + (sqrt(A)/Q)*s + 1)/(s^2 + (sqrt(A)/Q)*s + A)

  TPar r = 1/Q;
  g   = tan(0.5*w) * sqrt(A);
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  aL  = 1;
  aB  = A*r;
  aH  = A*A;
}

// Processing:

template<class TSig, class TPar>
inline void rsStateVariableFilterMystran<TSig, TPar>::getOutputs(
  TSig in, TSig* yL, TSig* yB, TSig* yH)
{
  // Compute outputs:
  *yH = (in - gpr*z1 - z2) * scl;  // == (in - (g+r)*z1 - z2) / (1 + g*(g+r));
  *yB = z1 + g * *yH; 
  *yL = z2 + g * *yB;

  // State variable update:
  z1 = 2 * *yB - z1;               // Equivalent to: z1 += 2 * g * *yH
  z2 = 2 * *yL - z2;               // Equivalent to: z2 += 2 * g * *yB
}

template<class TSig, class TPar>
inline TSig rsStateVariableFilterMystran<TSig, TPar>::getSample(TSig in)
{
  TSig yL, yB, yH;
  getOutputs(in, &yL, &yB, &yH);   // Produce LP, BP and HP signals
  return aH*yH + aB*yB + aL*yL;    // Mix them according to desired filter type
}

#endif