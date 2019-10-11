#ifndef rosic_DampingFilter_h
#define rosic_DampingFilter_h

namespace rosic
{


/** This class implements a filter which is meant be used as a frequency dependent gain inside a
feedback loop, so as to apply a frequency dependent damping inside the loop. It applies a global 
gain factor to a signal as well as a first order low-shelving and a first order high-shelving 
filter. In that respect, it is very much like the ToneControl class, but here we use a different 
definition of the corner-frequency. In the ToneControl class, the corner-frequency is defined to 
be the frequency at which the gain is the geometric mean between the gain the reference gain 
(which is unity) - here the reference gain does not need to be unity and the relative gain at the
corner-frequency can be specified arbitrarily. This facilitates the use of the filter inside the 
feedback-loop of a delay-line - here it may be desirable to define the corner-freq in terms of 
decay-time instead of in terms of gain. */

class DampingFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. */
  DampingFilter();

  /** Destructor. */
  ~DampingFilter();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Sets the samplerate. */
  void setSampleRate(double newSampleRate);

  /** Sets up a global gain factor (as factor, not in dB). */
  void setGlobalGainFactor(double newGlobalGainFactor);

  /** Sets the realative gain factor of the low-shelving filter  */
  void setLowGainFactor(double newLowGainFactor);

  /** Sets the corner frequency of the low-shelving filter */
  void setLowCrossoverFreq(double newLowCornerFreq);

  /** Sets the relative gain factor at which the crossover-frequency is measured for the
  low-shelving filter. See Orfanidis' paper about High Order Equalizers for details. */
  void setLowCrossoverGainFactor(double newLowCrossoverGainFactor);

  /** Sets the relative gain factor of the high-shelving filter. */
  void setHighGainFactor(double newHighGainFactor);

  /** Sets the corner frequency of the high-shelving filter. */
  void setHighCrossoverFreq(double newHighCornerFreq);

  /** Sets the relative gain factor at which the crossover-frequency is measured for the
  high-shelving filter. */
  void setHighCrossoverGainFactor(double newHighCrossoverGainFactor);

  //-----------------------------------------------------------------------------------------------
  // \name Processing:

  /** Calculates one output sample at a time. */
  INLINE double getSample(double in);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Resets the internal buffers of the filter to zero. */
  void reset();

  //===============================================================================================

protected:

  /** Calculates the coefficients for the filter. */
  void calculateCoefficients();

  /** Checks whether or not the combination of G and G_B is an allowed configuration - it must
  satisfy (G != 1) and (G_B != 1) and ( (1 < G_B < G) or (G < G_B < 1) ) in order to be
  allowed. */
  bool areGainsAllowed(double G, double G_B);

  // buffering:
  doubleA x1, x2;  // past input values
  doubleA y1, y2;  // past output values

  // overall 2-pole filter-coefficients:
  doubleA b0, b1, b2; // feedforward coeffs
  doubleA     a1, a2; // feedback coeffs

  // 1-pole filter-coefficients:
  //doubleA b0Ls, b1Ls, a1Ls; // low-shelf coeffs
  //doubleA b0Hs, b1Hs, a1Hs; // high-shelf coeffs

  // filter parameters:
  doubleA globalGainFactor;

  doubleA lowCrossoverFreq;
  doubleA lowCrossoverGainFactor;
  doubleA lowGainFactor;

  doubleA highCrossoverFreq;
  doubleA highCrossoverGainFactor;
  doubleA highGainFactor;

  doubleA sampleRate;
  doubleA sampleRateRec;  // reciprocal of the sampleRate

  //MutexLock mutex;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

INLINE double DampingFilter::getSample(double in)
{
  //mutex.lock();

  // calculate output-sample:
  double out = b0*in + b1*x1 + b2*x2 + a1*y1 + a2*y2; // + TINY;

  if(RAPT::rsIsNaN(out)) // do we still need this or is the issue fixed?
    DEBUG_BREAK;
  //if( _isnan(out) ) // we need to write our own, maybe isNan(double x ) { return x != x; }, see RSLib
  //  DEBUG_BREAK;

  // update buffer-variables:
  x2 = x1;
  x1 = in;
  y2 = y1;
  y1 = out;

  //mutex.unlock();

  return out;
}

} // end namespace rosic

#endif // rosic_DampingFilter_h
