#ifndef RAPT_STATEVARIABLEFILTER_H_INCLUDED
#define RAPT_STATEVARIABLEFILTER_H_INCLUDED

/** This is an implementation of a state variable filter using topology-preserving transform (TPT)
and zero-delay feedback (ZDF) technology. You can either use 3 outputs (lowpass, bandpass, 
highpass) from the SVF core by calling getOutputs() or let the filter itself form a linear 
combination of the these 3 to obtain a desired filter mode (in addition to the 3 modes above, there 
are also shelvers, a bell, etc.). */

template<class TSig, class TPar> // signal, parameter types
class rsStateVariableFilter
{
  typedef const TSig& CRSig;
  typedef const TPar& CRPar;

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  rsStateVariableFilter();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(TPar newSampleRate);

  /** Enumeration of the available filter modes. */
  enum modes         // rename to Mode
  {
    BYPASS = 0,

    // Rober Bristow-Johnson's Cookbook designs:
    LOWPASS,         // rename to LOWPASS_RBJ or LowpassRBJ ...and the same also for HPF, BPF, etc.
    HIGHPASS,
    BANDPASS_SKIRT,  // constant skirt gain
    BANDPASS_PEAK,   // constant peak gain
    BANDREJECT,
    BELL,
    LOWSHELF,
    HIGHSHELF,
    ALLPASS,


    // UNDER CONSTRUCTION:
    // Martin Vicanek's hybrid MZT/manitude-matching designs
    LowpassMVS,      // The S stands for "simple", the cheaper design



    MORPH_LP_BP_HP,  // under construction

    NUM_MODES
  };


  /** Chooses the filter mode. See the enumeration for available modes. */
  void setMode(int newMode);

  /** Sets the characteristic frequency (which is the cutoff-frequency for lowpass and highpass,
  the center-frequency for bandpass, bandreject and bell, and the halfgain-frequency for
  shelvers). */
  void setFrequency(TPar newFrequency);

  /** Sets the resonance gain (as linear gain factor) for low-, high- and (constant skirt gain)
  bandpass filters or the boost/cut gain for bell- and shelving filters. */
  void setGain(TPar newGain);

  /** Sets the bandwidth (in octaves) for (constant peak gain) bandpass filters and bell filters.
  In the case of shelving filters, this also determines the slope at the halfgain point.
  At B = (2*asinh(1/rsSqrt(2)))/log(2) = 1.899968626952992, the slope is as steep as it can be
  without overshooting. */
  void setBandwidth(TPar newBandwidth);

  /** Morphing parameter for the morphing filter types - this feature is under construction. */
  void setMorph(TPar newMorph);

  /** Sets up the filter coefficients to simulate a biquad filter with given coeffs. */
  void setupFromBiquad(CRPar b0, CRPar b1, CRPar b2, CRPar a1, CRPar a2);
  // not yet tested


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** When you use getOutputs() directly to obtain the lowpass, bandpass and highpass signals
  of the core SVF, this function returns the value by which the bandpass signal has to be scaled
  in order to achieve lowpass + scaler*bandpass + highpass == original input. */
  inline TPar getBandpassScaler() const { return R2; }


  /** \name Audio Processing */

  /** Returns the 3 outputs (lowpass, bandpass, highpass) of the core SVF. */
  inline void getOutputs(TSig in, TSig &yL, TSig &yB, TSig &yH);
    // use pointer parameters for the outputs

  /** Returns an appropriate linear combination of the 3 outputs of the core SVF in order to
  achieve various filter modes as selected via setMode(). */
  inline TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the internal state buffers to zero. */
  void reset();

protected:

  /** \name Internal Functions */

  /** Calculates filter coefficients from the user parameters */
  void calcCoeffs();

  /** Computes damping coefficient R from desired bandwidth and the prewarped radian center
  frequency (for bandpass with constant peak gain, bandreject and allpass). */
  TPar bandwidthToR(TPar B);


  /** \name Data */

  // state:
  TSig s1, s2;

  // filter coefficients:
  TPar g;          // embedded integrator gain
  TPar R2;         // twice the damping coefficient (R2 == 2*R == 1/Q)
  TPar h;          // factor for feedback (== 1/(1+2*R*g+g*g))
  TPar cL, cB, cH; // coefficients for low-, band-, and highpass signals

  // parameters:
  TPar fs;         // sample-rate
  TPar fc;         // characteristic frequency
  TPar G;          // gain
  TPar B;          // bandwidth (in octaves)
  TPar m;          // morph parameter (0...1)
  int  mode;       // filter-mode

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
inline void rsStateVariableFilter<TSig, TPar>::getOutputs(TSig in, TSig &yL, TSig &yB, TSig &yH)
{
  // compute highpass output via Eq. 5.1:
  //yH = (in - R2*s1 - g*s1 - s2) * h;  // 3 mul, 3 sub
  yH = (in - (R2+g) * s1 - s2) * h;     // 2 mul, 2 sub, 1 add

  // compute bandpass output by applying 1st integrator to highpass output:
  yB = g*yH + s1;
  s1 = g*yH + yB; // state update in 1st integrator

  // compute lowpass output by applying 2nd integrator to bandpass output:
  yL = g*yB + s2;
  s2 = g*yB + yL; // state update in 2nd integrator

  // Notes:
  // we have used two TDF2 integrators (Fig. 3.11) where one of them would be in code:
  //   y = g*x + s; // output computation
  //   s = g*x + y; // state update
  // Implementation uses 6 mul, 5 add, 2 sub, 5 assign. A DF1 biquad uses 5 mul, 4 add, 5 assign so 
  // that's 1 mul less, 2 add/sub less. Also, an SVF needs the additional h coeff to cache a 
  // division result. So, DF biquads may still have a (small) advantage in terms of operation count 
  // and memory usage. ToDo: actually do some benchmarks of SVF vs DF1 vs DF2 vs TDF2 vs .... We may
  // also want to benchmark coefficient calculations then.
    
  // As a cheap trick to introduce nonlinear behavior, we apply a nonlinearity to the states of 
  // the integrators (uncomment, if you want that):
  //s1 = tanh(s1);
  //s2 = tanh(s2);
}

template<class TSig, class TPar>
inline TSig rsStateVariableFilter<TSig, TPar>::getSample(TSig in)
{
  TSig yL, yB, yH;
  getOutputs(in, yL, yB, yH);
  return cL*yL + cB*yB + cH*yH;
}

#endif
