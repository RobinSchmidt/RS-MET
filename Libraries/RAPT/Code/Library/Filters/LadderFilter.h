#ifndef RAPT_LADDERFILTER_H_INCLUDED
#define RAPT_LADDERFILTER_H_INCLUDED

//#include <complex>  // temporary - this should go away and be done in the RAPT.h ..or somewhere else
//using namespace std;

/**

This is an implementation of a ladder filter based on a chain of 4 one-pole lowpasses without a 
zero (i.e. leaky integrators) of the general form: y[n] = b*x[n] - a*y[n-1] where a is in the range 
(-1,0] and b = 1+a. It's a linear ladder with a unit delay in the feedback path. The coefficients 
of the lowpass stages (a,b) and the feedback gain are adjusted in order to compensate for the side 
effects of the unit delay, such that there's no shift of the resonant frequency and no instability 
artifacts for resonance-values < 1.

\todo: maybe provide a function setCutoffAndResonance such that when both are called at samplerate,
calcCoeffs is not called twice

\todo: Maybe introduce an automatic gain-scaling in dependence of the cutoff - especially for 
highpasses (but also bandpasses) it may be desirable to boost the signal, as the cutoff gets higher
to prevent the sound becoming too weak. But this should be done in a subclass (it complicates the
calculations a lot)

\todo: explore the possibilities of spreading the cutoff frequencies of the individual lowpass 
stages instead of letting them have all the same cutoff frequency

*/

template<class TSig, class TPar> // signal, parameter types
class LadderFilter
{

public:

  /** Enumeration of the available filter modes. */
  enum modes  // rename to LadderMode or just Mode, see  state handling here: 
              // https://www.juce.com/doc/tutorial_playing_sound_files
  {
    FLAT = 0,   
    LP_6,
    LP_12,
    LP_18,
    LP_24,
    HP_6,
    HP_12,
    HP_18,
    HP_24,
    BP_6_6,
    BP_6_12,
    BP_6_18,
    BP_12_6,
    BP_12_12,
    BP_18_6,

    //// experimental, from a KVR thread:
    //KVR_BP2,   // = BP_6_6 * 2 (2 times louder - fits the lowpass better)
    //KVR_BP4,   // = BP_12_12 * 2
    //KVR_NF2,   // notch/bandreject
    //KVR_NF4,   // notch/bandreject
    //KVR_PF2,   // 2nd order allpass? phase-shift?
    //KVR_PF4,   // 4th order allpass
    //// todo: do the math

    NUM_MODES
  };
  // add more modes: allpass1/2/3/4, notch(es), maybe peak if possible (perhaps requires gain 
  // parameter), shelf
  // maybe rename to Lowpass_6, etc., 


  /** \name Construction/Destruction */

  /** Constructor. */
  LadderFilter();


  /** \name Setup */

  /** Sets the cutoff frequency in Hz. */
  void setCutoff(TPar newCutoff);
    // maybe rename to setFrequency because we have different modes and the term cutoff applies
    // only to lowpass- and highpass filters (not to bandpass filters, for example)

  /** Sets the sample rate in Hz. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the resonance. This is the net amount of gain in the feedback loop - so the parameter is
  normalized to the range 0..1 where at 1, the stability/self-oscillation limit is reached. */
  void setResonance(TPar newResonance);

  /** Sets the coefficients by which the different signals, tapped off after successive 1-pole
  lowpass stages, are mixed together to form the final output. c0 multiplies the input, c1 the 
  output of the 1st lowpass stage, and so on. These mixing coefficients determine the general shape
  of the frequency response. For example, if c4 = 1 and all others are 0, we get a 24 dB/oct 
  lowpass - as in the classical Moog filter. */
  void setMixingCoefficients(TPar c0, TPar c1, TPar c2, TPar c3, TPar c4);

  /** Chooses the filter mode. See the enumeration for available modes. */
  void setMode(int newMode);


  /** \name Inquiry */

  /** Returns the gain factor that needs to be applied to the output signal to ensure unit gain
  at DC (in the lowpass case). Without such a compensation, the DC gain will typically drop with
  incresing resonance. This is only relevant, if you plan to use the getSampleNoGain functions in
  the subclasses, the regular getSample functions apply this compensation gain already.*/
  TPar getCompensationGain() { return g; }

  /** Writes the filter's internal state buffer into the passed array which must be of length 5. */
  void getState(TSig *state);

  /** Returns the filter's z-domain transfer function value at the given value of z. */
  complex<TPar> getTransferFunctionAt(complex<TPar> z);

  /** Returns the filter's magnitude response at the given frequency in Hz. */
  TPar getMagnitudeResponseAt(TPar frequency);

  // currently, the


  /** \name Audio Processing */

  /** Returns a single output sample without gain-compensation */
  inline TSig getSampleNoGain(TSig in)
  {
    //y[0] = in     - k*y[4];
    //y[0] = tanh(in - k*y[4]);   // test

    y[4] /= 1 + y[4]*y[4];     // nonlinearity applied to the feedback signal
    y[0]  = in - k*y[4];       // linear
    //y[0] /= 1 + y[0]*y[0];     // nonlineariry applied to input plus feedback signal
    y[1]  = b*y[0]  - a*y[1];
    y[2]  = b*y[1]  - a*y[2];
    y[3]  = b*y[2]  - a*y[3];
    y[4]  = b*y[3]  - a*y[4];
    return c[0]*y[0] + c[1]*y[1] + c[2]*y[2] + c[3]*y[3] + c[4]*y[4];

    // we should experiment with placing saturation at different points..
  }

  /** Returns a single output sample (with gain-compensation such that the DC-gain remains 
  unity, irrespective of the resonance) */
  inline TSig getSample(TSig in)
  {
    return g * getSampleNoGain(in);
  }

  /** Processes a buffer of given length. */
  inline void process(TSig *in, TSig *out, int length)
  {
    for(int n = 0; n < length; n++)
      out[n] = getSample(in[n]);
  }


  /** \name Misc */

  /** Resets the internal states to zero. To be overriden in subclasses. */
  void reset();


  /** \name Coefficient Computations */

  ///** Computes the desired compensation gain factor to get unit gain at DC with given filter- and
  //feedback coefficients a, b, k. */
  ////static TPar computeCompensationGain(TPar a, TPar b, TPar k);
  //static TPar computeCompensationGain(TPar k);

  /** Given some normalized net feedback loop gain fb (in the range 0..1 where 1 is the 
  self-oscillation/instability limit), cos(wc) and lowpass coefficients a, b, this function 
  computes the feedback factor k. */
  static TPar computeFeedbackFactor(TPar fb, TPar cosWc, TPar a, TPar b);

  /** Given a desired decay time for the resonance (defined as the time it takes to fall to the
  value 1/e = 0.36..) in seconds and a cutoff frequency in Hz, this function computes the desired
  normalized net feedback loop gain (in the range 0..1) to achieve this decay time. */
  static TPar resonanceDecayToFeedbackGain(TPar decay, TPar cutoff);

  /** Given a normalized radian cutoff frequency wc (in the range 0...pi) and a normalized overall
  feedback gain fb (in the range 0...1), this function computes the desired coefficients a, b for
  the one pole filter that realizes the difference equation: y[n] = b*x[n] - a*y[n-1], the
  feedback gain k by which the output of a chain of 4 such one-pole units should be fed back into
  the first unit and a compensation gain g that compensates for the loss of DC gain when turning up
  the feedback. */
  static void computeCoeffs(TPar wc, TPar fb, TPar *a, TPar *b, TPar *k, TPar *g);
  // todo: factor the function into a/b-computation, k-computation, g-computation - but leave 
  // this one as convenience function also

  /** Same as computeCoeffs(double wc, double fb, double *a, double *b, double *k, double *g) but 
  without the compensation gain computation.  */
  static void computeCoeffs(TPar wc, TPar fb, TPar *a, TPar *b, TPar *k);

  // make a static method for the output coefficients c[] 


protected:

  /** \name Internal Functions */
                        
  /** Calculates the one-pole coefficients, feedback-gain and output gain from the parameters. */
  virtual void calcCoeffs(); 

  /** \name Data */

  TPar c[5];        // mixing coeffs for stages 0..4
  TSig y[5];        // outputs of the stages 0..4
  TPar cutoff;      // cutoff frequency in Hz
  TPar resonance;   // resonance 0..1
  TPar sampleRate;  // samplerate in Hz
  int  mode;        // filter mode (see modes-enum)
  TPar k;           // feedback gain
  TPar g;           // output gain
  TPar a, b;        // leaky integrator coefficients for a stage: y[n] = b*x[n] - a*y[n-1]

  // maybe combine cutoff and samplerate into wc (normalized radian cutoff)
};

#endif
