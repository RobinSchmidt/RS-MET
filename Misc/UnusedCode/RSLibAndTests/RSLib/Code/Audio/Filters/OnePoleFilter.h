#ifndef RS_ONEPOLEFILTER_H
#define RS_ONEPOLEFILTER_H

namespace RSLib
{

  /**

  This is an implementation of a simple one-pole filter unit.

  \todo rename to FirstOrderFilter (it does not only have a pole but also a zero).
  \todo make it possible to set up time constants in terms of dB/sec

  */

  class RSLib_API rsOnePoleFilter
  {

  public:

    /** This is an enumeration of the available filter modes. */
    enum modes
    {
      BYPASS = 0,
      LOWPASS,         // rename to LOWPASS_IIT (impulse invariant transform)
      HIGHPASS,        // rename to HIGHPASS_MZT (matched Z transform)
      ALLPASS,
      LOWSHELV,        // rename to LOWSHELF_NMM "Nyquist Magnitude Match"
      HIGHSHELV,       // rename to LOWSHELF_NMM
      LOWSHELV_DAFX,   // rename to LOWSHELF_BLT
      HIGHSHELV_DAFX   // rename to HIGHSHELF_BLT
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
    void setSampleRate(double newSampleRate);

    /** Chooses the filter mode. See the enumeration for available modes. */
    void setMode(int newMode);

    /** Sets the cutoff-frequency for this filter. */
    void setCutoff(double newCutoff);

    /** This will set the time constant 'tau' for the case, when lowpass mode is chosen. This is 
    the time, it takes for the impulse response to die away to 1/e = 0.368... or equivalently, the
    time it takes for the step response to raise to 1-1/e = 0.632... */
    void setLowpassTimeConstant(double newTimeConstant) { setCutoff(1.0/(2*PI*newTimeConstant)); }

    /** Sets the gain factor for the shelving modes (this is not in decibels). */
    void setShelvingGain(double newGain);

    /** Sets the gain for the shelving modes in decibels. */
    void setShelvingGainInDecibels(double newGain);

    /** Sets the filter coefficients manually. */
    void setCoefficients(double newB0, double newB1, double newA1);

    /** Sets up the internal state variables for both channels. */
    void setInternalState(double newX1, double newY1);


    /** \name Inquiry */

    /** Returns the magnitude response of this filter at the given frqeuency (in Hz). */
    double getMagnitudeAt(double frequency);


    /** \name Audio Processing */

    /** Calculates a single filtered output-sample. */
    RS_INLINE double getSample(double in);


    /** \name Misc */

    /** Resets the internal buffers (for the \f$ x[n-1], y[n-1] \f$-samples) to zero. */
    void reset();

  protected:

    /** \name Internal Functions */
    void calcCoeffs();  // calculates filter coefficients from filter parameters


    /** \name Data */

    // buffering:
    double x1, y1;

    // filter coefficients:
    double b0, b1; // feedforward coeffs
    double a1;     // feedback coeff

    // filter parameters:
    double cutoff;
    double shelvingGain;
    int    mode;  

    double sampleRate; 
    double sampleRateRec;  // reciprocal of the sampleRate

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE double rsOnePoleFilter::getSample(double in)
  {
    y1 = b0*in + b1*x1 + a1*y1 + TINY;
    x1 = in;
    return y1;
  }

}

#endif
