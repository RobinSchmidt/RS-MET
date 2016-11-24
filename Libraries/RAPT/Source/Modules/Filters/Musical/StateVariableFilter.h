#ifndef RAPT_STATEVARIABLEFILTER_H_INCLUDED
#define RAPT_STATEVARIABLEFILTER_H_INCLUDED

/** This is an implementation of a state variable filter using topology-preserving transform (TPT)
and zero-delay feedback (ZDF) technology. You can either use 3 outputs (lowpass, bandpass, 
highpass) from the SVF core by calling getOutputs() or let the filter itself form a linear 
combination of the these 3 to obtain a desired filter mode (in addition to the 3 modes above, there 
are also shelvers, a bell, etc.). */

class StateVariableFilter
{

public:

  /** Enumeration of the available filter modes. */
  enum modes
  {
    BYPASS = 0,
    LOWPASS,
    HIGHPASS,
    BANDPASS_SKIRT,  // constant skirt gain
    BANDPASS_PEAK,   // constant peak gain
    BANDREJECT,
    BELL,
    LOWSHELF,
    HIGHSHELF,
    ALLPASS,

    MORPH_LP_BP_HP,  // under construction

    NUM_MODES
  };


  /** \name Construction/Destruction */

  /** Constructor. */
  StateVariableFilter();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Chooses the filter mode. See the enumeration for available modes. */
  void setMode(int newMode);

  /** Sets the characteristic frequency (which is the cutoff-frequency for lowpass and highpass,
  the center-frequency for bandpass, bandreject and bell, and the halfgain-frequency for
  shelvers). */
  void setFrequency(double newFrequency);

  /** Sets the resonance gain (as linear gain factor) for low-, high- and (constant skirt gain)
  bandpass filters or the boost/cut gain for bell- and shelving filters. */
  void setGain(double newGain);

  /** Sets the bandwidth (in octaves) for (constant peak gain) bandpass filters and bell filters.
  In the case of shelving filters, this also determines the slope at the halfgain point.
  At B = (2*asinh(1/rsSqrt(2)))/log(2) = 1.899968626952992, the slope is as steep as it can be
  without overshooting. */
  void setBandwidth(double newBandwidth);

  /** Morphing parameter for the morphing filter types - this feature is under construction. */
  void setMorph(double newMorph);


  /** \name Inquiry */

  /** When you use getOutputs() directly to obtain the lowpass, bandpass and highpass signals
  of the core SVF, this function returns the value by which the bandpass signal has to be scaled
  in order to achieve lowpass + scaler*bandpass + highpass == original input. */
  inline double getBandpassScaler() const { return R2; }


  /** \name Audio Processing */

  /** Returns the 3 outputs (lowpass, bandpass, highpass) of the core SVF. */
  inline void getOutputs(double in, double &yL, double &yB, double &yH);

  /** Returns an appropriate linear combination of the 3 outputs of the core SVF in order to
  achieve various filter modes as selected via setMode(). */
  inline double getSample(double in);


  /** \name Misc */

  /** Resets the internal state buffers to zero. */
  void reset();

protected:

  /** \name Internal Functions */

  /** Calculates filter coefficients from the user parameters */
  void calcCoeffs();

  /** Computes damping coefficient R from desired bandwidth and the prewarped radian center
  frequency (for bandpass with constant peak gain, bandreject and allpass). */
  double bandwidthToR(double B);


  /** \name Data */

  // state:
  double s1, s2;

  // filter coefficients:
  double g;          // embedded integrator gain
  double R2;         // twice the damping coefficient (R2 == 2*R == 1/Q)
  double h;          // factor for feedback (== 1/(1+2*R*g+g*g))
  double cL, cB, cH; // coefficients for low-, band-, and highpass signals

  // parameters:
  double fs;    // sample-rate
  double fc;    // characteristic frequency
  double G;     // gain
  double B;     // bandwidth (in octaves)
  double m;     // morph parameter (0...1)
  int    mode;  // filter-mode

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

inline void StateVariableFilter::getOutputs(double in, double &yL, double &yB, double &yH)
{
  // compute highpass output via Eq. 5.1:
  yH = (in - R2*s1 - g*s1 - s2) * h;

  // compute bandpass output by applying 1st integrator to highpass output:
  yB = g*yH + s1;
  s1 = g*yH + yB; // state update in 1st integrator

  // compute lowpass output by applying 2nd integrator to bandpass output:
  yL = g*yB + s2;
  s2 = g*yB + yL; // state update in 2nd integrator

  // Remark: we have used two TDF2 integrators (Fig. 3.11) where one of them would be in code:
  // y = g*x + s; // output computation
  // s = g*x + y; // state update

  // as a cheap trick to introduce nonlinear behavior, we apply a nonlinearity to the states of 
  // the integrators (uncomment, if you want that):
  //s1 = tanh(s1);
  //s2 = tanh(s2);
}

inline double StateVariableFilter::getSample(double in)
{
  double yL, yB, yH;
  getOutputs(in, yL, yB, yH);
  return cL*yL + cB*yB + cH*yH;
}

#endif
