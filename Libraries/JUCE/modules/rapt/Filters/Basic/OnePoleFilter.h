#ifndef RAPT_ONEPOLEFILTER_H_INCLUDED
#define RAPT_ONEPOLEFILTER_H_INCLUDED

/** This is an implementation of a simple one-pole filter unit.

\todo rename to FirstOrderFilter (it does not only have a pole but also a zero).
\todo make it possible to set up time constants in terms of dB/sec */

template<class TSig, class TPar>
class rsOnePoleFilter
{

public:

  /** This is an enumeration of the available filter modes. */
  enum modes
  {
    BYPASS = 0,
    LOWPASS_IIT,      // lowpass via impulse invariant transform
    HIGHPASS_MZT,     // highpass via matched-Z transform
    ALLPASS_BLT,      // allpass via bilinear transform
    LOWSHELV_NMM,     // low shelving via nyquist magnitude match
    HIGHSHELV_NMM,    // high shelving via nyquist magnitude match
    LOWSHELV_BLT,     // low shelving via bilinear transform
    HIGHSHELV_BLT,    // high shelving via bilinear transform
    LOWPASS_BLT,      // lowpass via bilinear transform ...needs test
    HIGHPASS_BLT      // highpass via bilinear transform ...needs test
  };
  // \todo (maybe): let the user choose between LP/HP versions obtained via bilinear trafo and 
  // impulse invariant trafo, magnitude-matched trafo
  // i.e. have options LOWPASS_IIT (impulse invariant transform), 
  // LOWPASS_BLT (bilinear tarnsform), LOWPASS_MZTI (matched z-transform improved), NFG (nyquist
  // frequency gain - a la orfanidis - maybe can also be called PMM for pointwise magnitude match)


  /** \name Construction/Destruction */

  /** Constructor. */
  rsOnePoleFilter();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(TPar newSampleRate);

  /** Chooses the filter mode. See the enumeration for available modes. */
  void setMode(int newMode);

  /** Sets the cutoff-frequency for this filter. */
  void setCutoff(TPar newCutoff);

  /** This will set the time constant 'tau' for the case, when lowpass mode is chosen. This is
  the time, it takes for the impulse response to die away to 1/e = 0.368... or equivalently, the
  time it takes for the step response to raise to 1-1/e = 0.632... */
  void setLowpassTimeConstant(TPar newTimeConstant) { setCutoff(1.0/(2*PI*newTimeConstant)); }

  /** Sets the gain factor for the shelving modes (this is not in decibels). */
  void setShelvingGain(TPar newGain);

  /** Sets the gain for the shelving modes in decibels. */
  void setShelvingGainInDecibels(TPar newGain);

  /** Sets the filter coefficients manually. */
  void setCoefficients(TPar newB0, TPar newB1, TPar newA1);

  /** Sets up the internal state variables for both channels. */
  void setInternalState(TSig newX1, TSig newY1);


  /** \name Inquiry */

  /** Returns the magnitude response of this filter at the given frqeuency (in Hz). */
  TPar getMagnitudeAt(TPar frequency);


  /** \name Audio Processing */

  /** Calculates a single filtered output-sample. */
  inline TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the internal buffers (for the \f$ x[n-1], y[n-1] \f$-samples) to zero. */
  void reset();

protected:

  /** \name Internal Functions */
  void calcCoeffs();  // calculates filter coefficients from filter parameters


  /** \name Data */

  // buffering:
  TSig x1, y1;

  // filter coefficients:
  TPar b0, b1; // feedforward coeffs
  TPar a1;     // feedback coeff

  // filter parameters:
  TPar cutoff;
  TPar shelvingGain;
  int    mode;

  TPar sampleRate;
  TPar sampleRateRec;  // reciprocal of the sampleRate

  // maybe factor out a baseclass that has only x1,y1,b0,b1,a1 as members

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
inline TSig rsOnePoleFilter<TSig, TPar>::getSample(TSig in)
{
  y1 = b0*in + b1*x1 + a1*y1;
  x1 = in;
  return y1;
}

#endif
