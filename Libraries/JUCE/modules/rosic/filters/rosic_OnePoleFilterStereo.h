#ifndef rosic_OnePoleFilterStereo_h
#define rosic_OnePoleFilterStereo_h

// rosic-indcludes:
#include "../math/rosic_ComplexFunctions.h"

namespace rosic
{

  /**

  This is an implementation of a simple one-pole filter unit.

  */

  class OnePoleFilterStereo
  {

  public:

    /** This is an enumeration of the available filter modes. */
    enum modes
    {
      BYPASS = 0,
      LOWPASS,
      HIGHPASS,
      LOWSHELV,
      HIGHSHELV,
      ALLPASS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    OnePoleFilterStereo();   

    //---------------------------------------------------------------------------------------------
    // parameter settings:

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
    void setShelvingGainInDb(double newGainDb);

    /** Sets the filter coefficients manually. */
    void setCoeffs(double newB0, double newB1, double newA1);

    /** Sets up the internal state variables for both channels. */
    void setState(double newX1L, double newX1R, double newY1L, double newY1R);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the value of the filters transfer function at the complex value z. */
    Complex getTransferFunctionAt(Complex z);

    /** Returns the value of the magnitude response of the filter at the given frequency. */
    double getMagnitudeAt(double frequency);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample. */
    INLINE double getSample(double in);

    /** Calculates a single filtered output-sample. */
    INLINE void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal buffers (for the \f$ x[n-1], y[n-1] \f$-samples) to zero. */
    void resetBuffers();

    //=============================================================================================

  protected:

    // buffering:
    double x1L, x1R, y1L, y1R;

    // filter coefficients:
    double b0; // feedforward coeffs
    double b1;
    double a1; // feedback coeff

    // filter parameters:
    double cutoff;
    double shelvingGain;
    int    mode;  

    double sampleRate; 
    double sampleRateRec;  // reciprocal of the sampleRate

    // internal functions:
    void calcCoeffs();  // calculates filter coefficients from filter parameters

    // declare friends:
    friend class FunctionNodeLowpass6;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double OnePoleFilterStereo::getSample(double in)
  {
    // calculate the output sample:
    y1L = b0*in + b1*x1L + a1*y1L + TINY;

    // update the buffer variables:
    x1L = in;

    return y1L;
  }

  INLINE void OnePoleFilterStereo::getSampleFrameStereo(double *inL, double *inR, 
    double *outL, double *outR)
  {
    // calculate the output sample:
    y1L = b0*(*inL) + b1*x1L + a1*y1L + TINY;
    y1R = b0*(*inR) + b1*x1R + a1*y1R + TINY;

    // update the buffer variables:
    x1L = *inL;
    x1R = *inR;

    // store the output into its slots:
    *outL = y1L;
    *outR = y1R;
  }

} // end namespace rosic

#endif // rosic_OnePoleFilterStereo_h
