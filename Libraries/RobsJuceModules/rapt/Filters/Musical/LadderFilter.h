#ifndef RAPT_LADDERFILTER_H_INCLUDED
#define RAPT_LADDERFILTER_H_INCLUDED

/** This is an implementation of a ladder filter based on a chain of 4 one-pole lowpasses without a 
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
stages instead of letting them have all the same cutoff frequency */

template<class TSig, class TPar> // signal, parameter types
class rsLadderFilter
{

public:
	//virtual ~LadderFilter() = default;

  // for convenience:
  typedef const TSig& CRSig;  // const reference to a signal value
  typedef const TPar& CRPar;  // const reference to a parameter value



  /** Enumeration of the available filter modes. */
  enum Mode
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

    NUM_MODES
  };
  // add more modes: allpass1/2/3/4, notch(es), maybe peak if possible (perhaps requires gain 
  // parameter), shelf
  // maybe rename to Lowpass_6, etc., 
  // see https://www.juce.com/doc/tutorial_playing_sound_files


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  rsLadderFilter();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the cutoff frequency in Hz. */
  void setCutoff(CRPar newCutoff);
    // maybe rename to setFrequency because we have different modes and the term cutoff applies
    // only to lowpass- and highpass filters (not to bandpass filters, for example)

  /** Sets the sample rate in Hz. */
  void setSampleRate(CRPar newSampleRate);

  /** Sets the resonance. This is the net amount of gain in the feedback loop - so the parameter is
  normalized to the range 0..1 where at 1, the stability/self-oscillation limit is reached. */
  void setResonance(CRPar newResonance);

  /** Sets the coefficients by which the different signals, tapped off after successive 1-pole
  lowpass stages, are mixed together to form the final output. c0 multiplies the input, c1 the 
  output of the 1st lowpass stage, and so on. These mixing coefficients determine the general shape
  of the frequency response. For example, if c4 = 1 and all others are 0, we get a 24 dB/oct 
  lowpass - as in the classical Moog filter. */
  void setMixingCoefficients(CRPar c0, CRPar c1, CRPar c2, CRPar c3, CRPar c4);

  /** Chooses the filter mode. See the Mode enum for the available modes. */
  void setMode(int newMode);

  /** Experimental feature - not yet ready for production! */
  void setBilinear(bool useBilinearTransform) 
  { bilinear = useBilinearTransform; updateCoefficients(); }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the gain factor that needs to be applied to the output signal to ensure unit gain
  at DC (in the lowpass case). Without such a compensation, the DC gain will typically drop with
  incresing resonance. This is only relevant, if you plan to use the getSampleNoGain functions in
  the subclasses, the regular getSample functions apply this compensation gain already.*/
  TPar getCompensationGain() { return g; }

  /** Writes the filter's internal state buffer into the passed array which must be of length 5. */
  void getState(TSig *state);

  /** Returns the filter's z-domain transfer function value at the given value of z. */
  std::complex<TPar> getTransferFunctionAt(const std::complex<TPar>& z); 
    // z needs to be a const-reference, too?

  /** Returns the filter's magnitude response at the given frequency in Hz. */
  TPar getMagnitudeResponseAt(CRPar frequency);
  // maybe rename to getMagnitudeAt

  /** Returns the transfer function rsRationalFunction object. */
  rsRationalFunction<TPar> getTransferFunction();

  /** Old implementation of getTransferFunction, using rsRationalFunction's arithmetic instead of
  just assigning the coeffs via analytically derived formulas (as the new one does). It's less 
  efficient and less precise that the new one, but nicely demonstrates how such a thing can be 
  done. It will also produce a function that is formally 8-pole, but fetaures pole/zero 
  cancellations. ToDo: move this eventually into the prototypes section. */
  rsRationalFunction<TPar> getTransferFunctionOld();

  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  /** Returns a single output sample without gain-compensation */
  //inline TSig getSampleNoGain(TSig in);
  inline TSig getSampleNoGain(CRSig in);

  /** Returns a single output sample (with gain-compensation such that the DC-gain remains 
  unity, irrespective of the resonance) */
  TSig getSample(CRSig in);

  /** Processes a buffer of given length. */
  void process(TSig *in, TSig *out, int length);

  /** Resets the internal states to zero. To be overriden in subclasses. */
  void reset();

  //-----------------------------------------------------------------------------------------------
  /** \name Coefficient Computations */

  /** Given some normalized net feedback loop gain fb (in the range 0..1 where 1 is the 
  self-oscillation/instability limit), cos(wc) and lowpass coefficient a, this function 
  computes the feedback factor k. */
  static TPar computeFeedbackFactor(CRPar fb, CRPar cosWc, CRPar a, bool bilinear);

  /** Given a desired decay time for the resonance (defined as the time it takes to fall to the
  value 1/e = 0.36..) in seconds and a cutoff frequency in Hz, this function computes the desired
  normalized net feedback loop gain (in the range 0..1) to achieve this decay time. */
  static TPar resonanceDecayToFeedbackGain(CRPar decay, CRPar cutoff);

  /** Given a normalized radian cutoff frequency wc (in the range 0...pi) and a normalized overall
  feedback gain fb (in the range 0...1), this function computes the desired coefficients a, b for
  the one pole filter that realizes the difference equation: y[n] = (1+a)*x[n] - a*y[n-1], the
  feedback gain k by which the output of a chain of 4 such one-pole units should be fed back into
  the first unit and a compensation gain g that compensates for the loss of DC gain when turning up
  the feedback. */
  static void computeCoeffs(CRPar wc, CRPar fb, CRPar s, TPar *a, TPar *k, TPar *g, bool bilinear);
  // todo: factor the function into a/b-computation, k-computation, g-computation - but leave 
  // this one as convenience function also

  /** Same as computeCoeffs(wc, fb, *a, *k, *g) but without the compensation gain computation. */
  static void computeCoeffs(CRPar wc, CRPar fb, TPar *a, TPar *k, bool bilinear);

  // make a static method for the output coefficients c[] and s

protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Internal Functions */
                        
  /** Calculates the one-pole coefficients, feedback-gain and output gain from the parameters. */
  //virtual void updateCoefficients();     // why virtual - can we get rid of this?
  void updateCoefficients(); 

  /** Computes the gain for a single filter stage in the bilinear case. In the no-zero case, 
  b = 1 + a. In the bilinear case, it's half that. The formula can be derived by defining: 
  c = cos(w), the magnitude squared is generally: (2*b^2 * (1+c)) / (1 + a^2 + 2*a*c) and setting 
  w = 0, such that c = 1 and then evaluating the DC gain as if b was 1 and taking the reciprocal.*/
  static TPar getBilinearB(TPar a) { return TPar(0.5) * (TPar(1) + a); }

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  // Algo parameters and state:
  TSig y[5];         // outputs of the stages 0..4 (filter state)
  TPar c[5];         // mixing coeffs for stages 0..4
  TPar a = TPar(0);  // coeff for a 1-pole stage (depends on cutoff)
  TPar k = TPar(0);  // feedback gain (depends on resonance)
  TPar g = TPar(1);  // compensation gain (depends on resonance)
  TPar s = TPar(1);  // scaler for k in compensation gain computation (depends on mode)

  // User parameters:
  TPar cutoff = TPar(1000);      // cutoff frequency in Hz
  TPar resonance = TPar(0);      // resonance 0..1
  TPar sampleRate = TPar(44100); // samplerate in Hz
  int  mode = Mode::FLAT;        // filter mode (see Mode enum)
  bool bilinear = false;         // flag to indicate using a zero at z=-1 per stage
  // todo: get rid of sample-rate, just maintain an omega = 2*PI*fc/fs

};

#endif
