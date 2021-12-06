#ifndef RAPT_LADDERFILTERS_H
#define RAPT_LADDERFILTERS_H

// rename to LadderVariations.h

// Various versions and extensions of the basic ladder filter. Currently experimental.

// \todo implement a rsLadderFilter baseclass, have subclasses rsLadderLinearZDF,
// rsLadderMystran, rsLadderLinear
//
// maybe have a version that provides a "spread" parameter that spreads the cutoff frequencies of
// the individual stages out..i guess, if it's implemented in such a way that the geometric mean
// of the cutoff frequencies equals the desired cutoff, the resonance won't move (but i'm not sure)
// perhaps, with such a spread parameter, peaking filters with adjustable badwidths could be
// possible?

//=================================================================================================

/** This is an implementation of a ladder filter based on a chain of 4 one-pole lowpasses without a
zero (i.e. leaky integrators) of the general form: y[n] = b*x[n] - a*y[n-1] where a is in the range
(-1,0] and b = 1+a. It's a linear ladder with a unit delay in the feedback path. The coefficients
of the lowpass stages (a,b) and the feedback gain are adjusted in order to compensate for the side
effects of the unit delay, such that there's no shift of the resonant frequency and no instability
artifacts for resonance-values < 1.

\todo merge this class with rsLadderFilter

*/

template<class TSig, class TPar>
class rsLadderFilter2
{

public:

  /** Enumeration of the available filter modes. */
  enum modes
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
  // parameter)


  /** \name Construction/Destruction */

  /** Constructor. */
  rsLadderFilter2();


  /** \name Setup */

  /** Sets the cutoff frequency in Hz. */
  void setCutoff(TPar newCutoff);

  /** Sets the sample rate in Hz. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the resonance. This is the net amount of gain in the feedback loop - so the parameter is
  normalized to the range 0..1 where at 1, the stability/self-oscillation limiti is reached. */
  void setResonance(TPar newResonance);

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


  /** \name Audio Processing */

  /** Returns a single output sample without gain-compensation */
  RS_INLINE TSig getSampleNoGain(TSig in);

  /** Returns a single output sample (with gain-compensation such that the DC-gain remains
  unity, irrespective of the resonance) */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the internal states to zero. To be overriden in subclasses. */
  void reset();


  /** \name Coeffiecient Computations */

  /** Computes the desired compensation gain factor to get unit gain at DC with given filter- and
  feedback coefficients a, b, k. */
  static TPar computeCompensationGain(TPar a, TPar b, TPar k);

  /** Given some normalized net feedback loop gain fb (in the range 0..1 where 1 is the
  self-oscillation/instability limit), cos(wc) and lowpass coefficients a, b, this function
  computes the feedback factor k. */
  static TPar computeFeedbackFactor(TPar fb, TPar cosWc, TPar a, TPar b);

  /** Given a desired decay time for the resonance (defined as the time it atkes to fall to the
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

  TPar c[5];           // mixing coeffs for stages 0..4
  TSig y[5];           // outputs of the stages 0..4
  TPar cutoff;         // cutoff frequency in Hz
  TPar resonance;      // resonance 0..1
  TPar sampleRate;     // samplerate in Hz
  int  mode;           // filter mode (see modes-enum)
  TPar k;              // feedback gain
  TPar g;              // ouput gain
  TPar a, b;           // leaky integrator coefficients for a stage: y[n] = b*x[n] - a*y[n-1]

  // maybe combine cutoff and samplerate into wc (normalized radian cutoff)
};

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilter2<TSig, TPar>::getSampleNoGain(TSig in)
{
  y[0] = in - k*y[4];
  y[1] = b*y[0] - a*y[1];
  y[2] = b*y[1] - a*y[2];
  y[3] = b*y[2] - a*y[3];
  y[4] = b*y[3] - a*y[4];
  return c[0]*y[0] + c[1]*y[1] + c[2]*y[2] + c[3]*y[3] + c[4]*y[4];
}

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilter2<TSig, TPar>::getSample(TSig in)
{
  return g * getSampleNoGain(in);
}

//=================================================================================================

/** This is a zero-delay-feedback (ZDF) version of the rsLadderFilter2 class. */

template<class TSig, class TPar>
class rsLadderFilterZDF : public rsLadderFilter2<TSig, TPar>
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsLadderFilterZDF();


  /** \name Audio Processing */

  /** Returns a single output sample without gain-compensation */
  RS_INLINE TSig getSampleNoGain(TSig in);

  /** Returns a single output sample (with gain-compensation such that the DC-gain remains
  unity, irrespective of the resonance) */
  RS_INLINE TSig getSample(TSig in);


protected:

  /** \name Internal Functions */

  /** Calculates the one-pole coefficients, feedback-gain and output gain from the parameters. */
  virtual void calcCoeffs();


  /** \name Data */

  TSig p[5];           // prediction coefficients for y[4]

};

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilterZDF<TSig, TPar>::getSampleNoGain(TSig in)
{
  // compute/predict output of final stage:
  this->y[4] = p[0]*in + p[1]*this->y[1] + p[2]*this->y[2] + p[3]*this->y[3] + p[4]*this->y[4];

  // update other filter states:
  this->y[0] =                   in - this->k * this->y[4];
  this->y[1] = this->b * this->y[0] - this->a * this->y[1];
  this->y[2] = this->b * this->y[1] - this->a * this->y[2];
  this->y[3] = this->b * this->y[2] - this->a * this->y[3];

  // form output (maybe factor out - this is the same as in the UDF case):
  return this->c[0] * this->y[0] +
         this->c[1] * this->y[1] +
         this->c[2] * this->y[2] +
         this->c[3] * this->y[3] +
         this->c[4] * this->y[4];

  // the ugly "this->" syntax is required by gcc for inherited members in template classes
}

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilterZDF<TSig, TPar>::getSample(TSig in)
{
  return this->g * getSampleNoGain(in);
}

//=================================================================================================

/**

This is a subclass of rsLadderFilterUDF that allows for saturation inside the feedback loop. The
user can set the lower and upper limiting value. In between the limits, an appropriately scaled and
shifted sigmoid function is applied. You can pass a function-poiter to the sigmoid to be used. It
should be normalized such that it limits to +-1 as the input goes to +-inf and have unit-slope at
the origin. When using a simple hardclipper that just clips at +-1, the function applied in the
feedback path is equivalent to y = clip(x, lowerLimit, upperLimit) - appropriate pre/post
scale/shift operations ensure this behavior. But you can also use tanh or other sigmoids that
satisfy the above normalization criteria.

ToDo:

maybe we need to override setCutoff, setDecay, etc. in order to apply a custom compensation
gain calculation here whenever one of these parameter changes - the goal is to always have unit
DC gain here, independently from the settings of decay/cutoff/drive etc. this also facilitates
the use in the resoshape filter.

rename to ...FeedbackClipped - or just ...Nonlinear

Currently, the total feedback can be seen a a multiplication of two factors, one from the
decay-setting "fbd" and one from the drive-setting that can be used to drive the filter into
self-oscillation "fbo", so we have: fb = fbd * fbo. Maybe other combinations, like fb = fbd + fbo,
fb = max(fbd, fbo) are worth to explore. The implications of the different formulas are:
fb = fbd * fbo:
 -fbo is ineffective when fbd = 0, i.e. decay is set to zero
 -when fbd depends on the cutoff (for example, for const-Q filters), the filter may cease to
  self-oscillate at low cutoffs
 -the feedback factor "k" can be separated into k_decay and k_drive (that's done in getSample)
fb = fbd + fbo:
 -when fbo >= 1, self-oscillation is ensured, regardless of fbd
 -when fbd > 0, the self-osc threshold for fbo is lowered
fb = max(fbd, fbo):
 -when fbd sweeps (for example as result of cutoff sweep in const-Q mode), the fb-sweep has a
  discontinuity, when fbd crosses fbo
 -it's always one value always that dominates, and within these two regimes, there are no
  interactions between the parameters
the formulas may be translated into each other, like
fbd * fbom = fbd + fboa -> fboa = fbd * fbom - fbd, fbom = (fbd + fboa) / fbd -  the latter formula
works only for fbd != 0
*/

template<class TSig, class TPar>
class rsLadderFilterFeedbackSaturated : public rsLadderFilter2<TSig, TPar>
{

public:

  /** Places where the feedback saturation function can be applied. */
  enum saturationPlaces
  {
    NOWHERE = 0,         // turns feedback saturation off
    PRE_FB_GAIN,         // before the "k" factor
    POST_FB_GAIN,        // after the "k" factor
    INPUT_AND_POST_GAIN, // to input and after k-factor
    POST_INPUT_ADD,      // after adding the input to the feedback signal
    POST_GAIN_AND_ADD,   // after "k" and again after adding input

    POST_1ST_STAGE,
    POST_2ND_STAGE,
    POST_3RD_STAGE,
    POST_4TH_STAGE,
    POST_EACH_STAGE,
    POST_EACH_AND_IN,

    IN_1ST_STAGE,
    IN_2ND_STAGE,
    IN_3RD_STAGE,
    IN_4TH_STAGE,
    IN_EACH_STAGE,

    NUM_SATPLACES,        //

    TEST1,                // experimental

  };
  // todo: figure out, if the pre-k/post-k modes are equivalent when adjusting the limits
  // appropriately -  i think, multiplying the limits by k should make them equivalent


  /** \name Construction/Destruction */

  /** Constructor. */
  rsLadderFilterFeedbackSaturated();


  /** \name Setup */

  /** Sets the lower clipping value for the feedback signal */
  void setLowerFeedbackLimit(TPar newLimit);

  /** Sets the upper clipping value for the feedback signal */
  void setUpperFeedbackLimit(TPar newLimit);

  /** Sets the drive factor for the feedback. */
  void setFeedbackDrive(TPar newDrive);

  /** Sets the gain of the saturating sigmoid function y = f(x) at x = 1. Allowed values are between
  0.5...1 (inclusive). */
  void setSaturationGainAt1(TPar newGain);

  /** Selects, at which point the feedback saturation is applied. @see saturationPlaces. */
  void setSaturationMode(int newMode);
    // maybe rename to setSaturationPlace


  /** \name Inquiry */

  /** Returns the compensation gain to be applied to the output of the filter to achieve unit gain
  at DC. This gain depends on the resonance decay and the feedback drive.
  \todo: move this function into the protected section, use it internally to compute the
  compensation gain. override setCutoff/etc. - anything that changes this gain */
  TPar getCompensationGain();

  TPar getLowerFeedbackLimit() { return loLimit; }
  TPar getUpperFeedbackLimit() { return hiLimit; }


  /** \name Audio Processing */

  RS_INLINE TSig getSampleNoGain(TSig in);
  RS_INLINE TSig getSample(TSig in);

  /** This function bypasses all the nonlinear additions introduced in this subclass and just
  accesses the inherited getSample function from the (linear) rsLadderFilter baseclass. We need
  this function for optimization purposes, when an object embedds a rsLadderFilterFeedbackSaturated
  object but wants to use only the linear features - this is the case in rsLadderresoShaped3 which
  inherits the nonlinear ladder filter from rsLadderresoShaped but needs only the linear features.
  It's a bit quirky, software-design wise...  */
  RS_INLINE TSig getSampleLinear(TSig in);

  /** Applies the saturation to the input */
  RS_INLINE TSig saturate(TSig x);

  // temporary - for test:
  rsOnePoleFilter<TSig, TPar> fbLpf;    // feedback lowpass - just for test
  TPar fbLpCutoff;        // relative cutoff for feedback lowpass

protected:

  /** Computes the scale/shift coefficients for the input and output of the normalized sigmoidal
  saturation function, so as to realize the desired range from loLimit to hiLimit. */
  void computeScaleAndShift();

  // user parameters:
  TPar loLimit, hiLimit;  // lower and upper limiting values
  TPar drive;             // saturation drive
  int mode;                 // decides, where the saturation is applied (temporary)

  // internal values:
  TPar scaleX, shiftX;    // scale and offset for the input for the normalized sigmoid
  TPar scaleY, shiftY;    // dito for normalized sigmoid's output

  // embedded objects:
  rsParametricSigmoid<TSig> sigmoid;
};

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilterFeedbackSaturated<TSig, TPar>::saturate(TSig x)
{
  return shiftY + scaleY * sigmoid.getValue(scaleX*x + shiftX);
}

// old:
//RS_INLINE double rsLadderFilterFeedbackSaturated::getSampleNoGain(double in)
//{
//  // decision, where the feedback saturation is applied:
//  switch(mode)
//  {
//  case NOWHERE:        y[0] = in - k*y[4];                 break; // no saturation
//  case PRE_FB_GAIN:    y[0] = in - k*saturate(drive*y[4]); break; // pre k
//  case POST_FB_GAIN:   y[0] = in - saturate(k*drive*y[4]); break; // post k
//  case POST_INPUT_ADD: y[0] = saturate(in - k*drive*y[4]); break; // post input add
//  // todo: have additional cases for post 1st stage, 2nd, 3rd, 4th
//  }
//
//  y[1] = b0*y[0] - a1*y[1];
//  y[2] = b0*y[1] - a1*y[2];
//  y[3] = b0*y[2] - a1*y[3];
//  y[4] = b0*y[3] - a1*y[4];
//
//  return c[0]*y[0] + c[1]*y[1] + c[2]*y[2] + c[3]*y[3] + c[4]*y[4];
//}

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilterFeedbackSaturated<TSig, TPar>::getSampleNoGain(TSig in)
{
  //// this is temporary test code:
  //fbLpf.setCutoff(fbLpCutoff * cutoff);
  //y[4] = fbLpf.getSample(y[4]);

  switch(mode)
  {
  case NOWHERE:                                  // no saturation
  {
    this->y[0] = in - this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case PRE_FB_GAIN:                              // pre k
  {
    this->y[0] = in - this->k * saturate(drive * this->y[4]);
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case POST_FB_GAIN:                             // post k
  {
    this->y[0] = in - saturate(this->k * drive * this->y[4]);
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case INPUT_AND_POST_GAIN:                   // input and post k
  {
    this->y[0] = saturate(in) - saturate(this->k * drive * this->y[4]);
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case POST_INPUT_ADD:                           // post input add
  {
    this->y[0] = saturate(in - this->k * drive * this->y[4]);
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case POST_GAIN_AND_ADD:
  {
    double tmp;
    tmp  = saturate(this->k * drive * this->y[4]);  // 1st saturation
    tmp  = saturate(in - tmp);      // 2nd saturation
                                    // todo: use different clipping limits for 1st and 2nd
    this->y[0] = tmp;
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;

  case POST_1ST_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = saturate(this->b * this->y[0] - this->a * this->y[1]);
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case POST_2ND_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = saturate(this->b * this->y[1] - this->a * this->y[2]);
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case POST_3RD_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = saturate(this->b * this->y[2] - this->a * this->y[3]);
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case POST_4TH_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = saturate(this->b * this->y[3] - this->a * this->y[4]);
  }
  break;
  case POST_EACH_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = saturate(this->b * this->y[0] - this->a * this->y[1]);
    this->y[2] = saturate(this->b * this->y[1] - this->a * this->y[2]);
    this->y[3] = saturate(this->b * this->y[2] - this->a * this->y[3]);
    this->y[4] = saturate(this->b * this->y[3] - this->a * this->y[4]);
  }
  break;

  case IN_1ST_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * saturate(this->y[1]);
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case IN_2ND_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * saturate(this->y[2]);
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case IN_3RD_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * saturate(this->y[3]);
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;
  case IN_4TH_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * saturate(this->y[4]);
  }
  break;
  case IN_EACH_STAGE:
  {
    this->y[0] = in - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * saturate(this->y[1]);
    this->y[2] = this->b * this->y[1] - this->a * saturate(this->y[2]);
    this->y[3] = this->b * this->y[2] - this->a * saturate(this->y[3]);
    this->y[4] = this->b * this->y[3] - this->a * saturate(this->y[4]);
  }
  break;
  case POST_EACH_AND_IN:
  {
    this->y[0] = saturate(in) - drive * this->k * this->y[4];
    this->y[1] = this->b * this->y[0] - this->a * saturate(this->y[1]);
    this->y[2] = this->b * this->y[1] - this->a * saturate(this->y[2]);
    this->y[3] = this->b * this->y[2] - this->a * saturate(this->y[3]);
    this->y[4] = this->b * this->y[3] - this->a * saturate(this->y[4]);
  }
  break;



  //case TEST1:
  //{
  //  y[0] = in - saturate(k*drive*y[4]);
  //  y[1] = b*y[0] - a*y[1];
  //  y[2] = b*y[1] - a*y[2];
  //  y[3] = b*y[2] - a*y[3];
  //  y[4] = b*y[3] - a*y[4];
  //}
  //break;
  //case TEST1:
  //{
  //  y[0] = in - k*y[4];
  //  y[1] = b*y[0] - a*y[1];
  //  y[2] = b*y[1] - a*y[2];
  //  y[3] = b*y[2] - a*y[3];
  //  y[4] = saturate(drive*(b*y[3] - a*y[4]));
  //}
  //break;
  case TEST1:
  {
    TSig C = -0.3;
    //y[0] = in - saturate(k*drive*y[4] + C*in) - saturate(C*in);
    this->y[0] = in - saturate(this->k * drive * this->y[4] + C * in);
    this->y[1] = this->b * this->y[0] - this->a * this->y[1];
    this->y[2] = this->b * this->y[1] - this->a * this->y[2];
    this->y[3] = this->b * this->y[2] - this->a * this->y[3];
    this->y[4] = this->b * this->y[3] - this->a * this->y[4];
  }
  break;

  default: return 0.0;
  }

  return this->c[0] * this->y[0] +
         this->c[1] * this->y[1] +
         this->c[2] * this->y[2] +
         this->c[3] * this->y[3] +
         this->c[4] * this->y[4];
}

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilterFeedbackSaturated<TSig, TPar>::getSample(TSig in)
{
  return this->g * getSampleNoGain(in);
    // maybe our gain here should be computed using a max(drive*feedback, 1.0) instead of the
    // regular normalized feedback
}

template<class TSig, class TPar>
RS_INLINE TSig rsLadderFilterFeedbackSaturated<TSig, TPar>::getSampleLinear(TSig in)
{
  return rsLadderFilter2<TSig, TPar>::getSample(in);
}

//=================================================================================================

/**

This is a variant of the ladder filter that splits off a purified resonance signal from the
(lowpass or otherwise) filtered part in order to separately process that pure resonance signal
further. The separation is achieved by running a resonant and a nonresonant filter and subtracting
the output of the nonresonant filter from the output of the resonant filter. The structure is:

in --------> NRF ------------>+----> out
      |             |         ^
      |            -|         |
      \----> RF ----+--> PP --/

NRF: nonresonant filter, RF: resonant filter, PP: post-processing

The main purpose of the separation of the resonance signal is the option to shape it further in
a postprocessing stage before adding it back to the nonresonant filter output. The post-processing
consists of passing the resonance signal through a resonator filter which applies a kind of attack
shape to the resonance, then adjusting the phase (by passing it through an allpass, giving a
quadrature (90° phase-shifted) version and forming a linear combination of the original resonance
and the quadrature component) and finally applying an overall gain to the resonance. Subclasses
may implement further post-processing steps, such as waveshaping, etc.


\todo: Adjust the cutoff frequency of the nonresonant filter in a way so as to minimize the
       bandwidth of the resonance signal - i.e. make it as pure as possible. I think, the cutoff of
       the lowpass should be increased a bit, but the amount of increase may depend on the cutoff
       frequency and/or resonance/decay-time. I suppose, that it should be a function of the
       relative bandwidth or Q. Some experimentation will be necessary. I think, minimizing the
       bandwidth is equivalent to minimizing the transient. The optimal frequency may also depend
       on the resonance gain (i.e. when it's 0, there should not be any shift of the cutoff
       frequency at all)

\todo: Switch to a (linear) ZDF ladder implementation.

\todo: Add a waveshaper to the resonance post-processing chain

Observations:
-the phase shift causes a dip in the spectrum below the resonance, it's most pronounced when
 the phase is around 135 degrees
-the attack may shift the dip - the higher the attack, the closer the dip moves to the resonance
 frequency

*/

template<class TSig, class TPar>
class rsLadderResoShaped
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. */
  rsLadderResoShaped();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the gain to be applied to the input as raw factor. */
  void setInputGain(TPar newGain);

  /** Sets the amount by which the filter input is leaked into the output. */
  void setInputLeak(TPar newLeak);

  /** Sets the cutoff-frequency. */
  void setCutoff(TPar newCutoff);

  /** Sets a decay-time for the resonance. This is the time, it takes for the resonance signal to
  decay down to 1/e. */
  void setResonanceDecay(TPar newDecay);

  /** Sets the amount of dependency of the resonance decay time from the resonance frequency. If
  it's set to zero, the decaytime will be independent from the resonant frequency, corresponding
  to a resonator filter with constant absolute bandwidth. If it's set to unity, the decaytime
  will be inversely proportional to the resonant frequency, corresponding to a resonator with
  constant relative bandwith (which implies constant Q as well). You may set it to values less
  than zero and larger than unity as well. In case of a nonzero dependency, we use a reference
  frequency of 1kHz - that means, at 1kHz (for resonance frequency), the decay time will be left
  unmodified. */
  void setDecayByFrequency(TPar newDecayByFreq);

  /** Sets the attack time for the resonance signal as multiplier for the (possibly frequency
  scaled) decay time of the ladder's resonance signal. Note that, for efficiency reasons, we
  define the attack time as the decay time constant of the two-pole filter that creates the
  attack. It's not exactly the time where the peak in the envelope occurs, because that would
  require a solution of an implicit equation and limit the factor to values < 1.
  \todo: maybe apply a mapping function, maybe look into findDecayScalerLess1 - the mapping should
  probably approximate that function. */
  void setResonanceAttack(TPar newAttack);

  /** Sets the gain of the resonance as raw multiplier. */
  void setResonanceGain(TPar newGain);

  /** Sets up a phase shift for the resonance (in radians) */
  void setResonancePhase(TPar newPhase);

  // feedback-saturation setters:
  void setFeedbackDrive(TPar newDrive);
  void setFeedbackLowerLimit(TPar newLimit);
  void setFeedbackUpperLimit(TPar newLimit);
  void setFeedbackSaturationGainAt1(TPar newGain);
  void setFeedbackSaturationPlace(int newPlace);


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes the two parts of the output signal: yf: the filtered signal, yr: the pure resonance
  signal, which when added together, give the final output. */
  virtual void getSignalParts(TSig in, TSig *yf, TSig *yr);

  /** Produces an output sample. */
  TSig getSample(TSig in);

  // todo: write a function that outputs also the intermediate signals - for experiments

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the internal states to zero. */
  void reset();


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Internal Functions */

  /** Sets up the scaled decay-time according to the parameters: decay, decayByFreq, cutoff and
  then calls setupAttackSmoother (because the attack time depends on the scaled decaytime). */
  void setupScaledDecay();

  /** Sets up the attack-smoothing filter. */
  void setupAttackSmoother();


  /** \name Embedded objects */

  rsLadderFilterFeedbackSaturated<TSig, TPar> resonant;
  rsLadderFilter2<TSig, TPar> nonResonant;
  rsTwoPoleFilter<TSig, TPar> attackSmoother;
  rsOnePoleFilter<TSig, TPar> allpass;
    // maybe later, instead of embedding a filter-object for the resonant and nonresonant signal,
    // we could derive this class from rsLadderFilterFeedbackSaturated
    // ...actually, eventually it may make more sense to implement the whole filter in one big
    // monolithic class. there might be some synergy-effects like not needing to compute certain
    // intermediate variables twice (or more) - values like cos(wc) are required in the coefficient
    // calculations of all embedded filter objects, so it makes sense to compute them only once.


  /** \name User Parameters */

  TPar sampleRate;            // sample rate
  TPar inGain;                // input gain factor
  TPar leak;                  // input leakage factor
  TPar cutoff;                // cutoff frequency
  TPar decay;                 // resonace decaytime in seconds
  TPar decayByFreq;           // dependency of decay from cutoff frequency
  TPar attack;                // resonance attack time (multiplier for (freq-scaled) decay time)
  TPar resGain;               // gain factor for resonance signal
  TPar phase;                 // phase-shift for resonance signal


  /** \name Internal data */

  double scaledDecay;  // decay, scaled by frequency

};

//=================================================================================================

/** Subclass of rsLadderResoShaped that adds nonlinear features such as waveshaping of the
resonance signal. */

template<class TSig, class TPar>
class rsLadderResoShaped2 : public rsLadderResoShaped<TSig, TPar>
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsLadderResoShaped2();


  /** \name Setup */

  /** Sets the saturation bypassed (or not). */
  //void setSaturationBypass(bool shouldBeBypassed);

  /** Selects a saturation mode. @see rsSidechainSaturator::modes. */
  void setSaturationMode(int newMode);

  ///** Sets the waveshaping function. You should pass a function-pointer to a function that takes
  //a "double" as argument and returns a "double" as return-value. */
  //void setSaturationFunction(double (*func)(double));

  ///** Sets the threshhold above which the nonlinear part of the saturation function begins. It
  //should be a value between 0..1. It can also be interpreted as a "hardness" parameter, where 0
  //amounts to the softest distortion and 1 amounts to full-blown hardclipping.  */
  //void setSaturationThreshold(double newThreshold);

  /** Sets the gain of the parametric sigmoid y = f(x) at x = 1. Allowed values are from 0.5...1. */
  void setSaturationGainAt1(TPar newGain);

  /** Sets the drive-multiplier wavershaper input. This is the factor by which the resonance signal
  gets multiplied before it goes into the waveshaper. */
  void setDrive(TPar newDrive);

  /** Sets the amount by which the pre-saturation gain (given by the drive factor) is to be undone
  post-saturation. If it's 0, there's no compensation, if it's 1, there's full compensation such
  that postGain = 1/preGain. Generally, we have postGain = (1/preGain)^compensation. */
  void setDriveCompensation(TPar newCompensation);

  /** Sets the offset to be added to the saturator sidechain input. */
  void setSaturationAddConstant(TPar newOffset);

  /** Sets the amount by which the input signal is added to the saturator sidechain input. */
  void setSaturationAddInput(TPar newAmount);

  /** Sets the amount by which the nonresonant filter output signal is added to the saturator
  sidechain input. */
  void setSaturationAddFiltered(TPar newAmount);

  /** Sets the sensitivity of teh gate as inverse factor by which the lower and upper saturation
  limits are scaled. 0: gating is effectively turned off, 1: limits equal to saturation limits -
  generally gateLimit = saturationlimit / gateSesitivity, which goes to infinity when the
  sensitivity goes to 0. */
  void setGateSensitivity(TPar newSensitivity);

  // 0: only input signal, 1: only lowpass signal
  void setGateMix(TPar newMix);

  // maybe we can also have a gate-amount parameter that crossfades between gated und non-gated
  // signal
  // REMOVE the gate stuff...

  /** \name Audio Processing */

  /** Computes the two parts of the output signal: yf: the filtered signal, yr: the pure resonance
  signal, which when added together, give the final output. */
  virtual void getSignalParts(TSig in, TSig *yf, TSig *yr); // override


protected:

  /** \name Embedded objects */

  rsSidechainSaturator<TSig, TPar> saturator;

  /** \name User Parameters */

  TPar drive;           // drive factor for saturator
  TPar addConst;        // constant offset to be added to saturator sidechain signal
  TPar addIn;           // multiplier for input signal for the sidechain
  TPar addFlt;          // multiplier for filter-output signal for the sidechain
  TPar driveComp;       // amount of gain compensation for drive
  TPar unDrive;         // post-saturation gain = (1/drive)^driveComp

  // remove:
  TPar gate;            // "amount" of gating of the resonance
  TPar gateMix;         // mix between input and lowpass signal in gate signal

};

//=================================================================================================

/**

Subclass of rsLadderResoShaped that adds the option of replacing the resonance waveform. This is
achieved by measuring instantaneous amplitude and phase of the resonance signal, that is, we model
the resonance signal as an enveloped sinusoid at the cutoff frequency: r = a * sin(p), where a
and p are instantaneous amplitude and phase which are measured from the resonance signal r and then
replace it by some other, arbitrary waveform: r' = a * wave(p). At the moment, the standard
waveforms: sine, triangle, square, saw-up, saw-down are provided but the structure is general
enough to support arbitrary waveforms.

todo:
-maybe try to do the measurement based on 2 succesive resonance samples using rsSineAmplitudeAndPhase
 instead of the allpass quadrature method - compare both methods, choose the better. ...but it think
 rsSineAmplitudeAndPhase will be much more expensive, so maybe it's better to stick with the allpass
 method.
-let the user choose an arbitrary waveform by setting wavetable - embedd a wavetable oscillator
 object to achieve this
 -anti-alias the wavetable oscillator by mip-mapping - we need to keep track of the (absolute value
  of) the phase difference between two successive samples to select the appropriately decimated
  wavetable
-maybe inherit from rsLadderResoShaped3 if we want to use waveshaping features from there

*/

template<class TSig, class TPar>
class rsResoReplacer : public rsLadderResoShaped<TSig, TPar>
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsResoReplacer();


  /** \name Setup */

  /** Selects a standard waveform for the resonance. 0: sine, 1: triangle, 2: square, 3: saw-up,
  4: saw-down. */
  void setResonanceWaveform(int newWaveform);

  /** Sets the multiplier by which the filter's cutoff frequency is to be multiplied to obtain the
  cutoff frequency for the filter that filters the resonance waveform. */
  void setResoCutoffMultiplier(TPar newMultiplier);

  /** Sets the amount by which the resonance waveform filter's cutoff frequency is modulated by the
  instantaneous amplitude of the resonance signal. */
  void setResoCutoffModulation(TPar newModulation);


  void setSampleRate(TPar newSampleRate);  // override


  /** \name Audio Processing */

  virtual void getSignalParts(TSig in, TSig *yf, TSig *yr);



  /** Given an input sample, this function computes the nonresonant filter output and the resonance
  signal reconstruction parameters (i.e. instantaneous amplitude and phase).*/
  virtual void getFltAmpPhs(TSig input, TSig *nonRes, TSig *resAmp, TSig *resPhase);


protected:





  /** Given an instantaneous amplitude and phase in the range 0..2*PI (or beyond 2*PI, but not below
  0), this function returns the intantaneous value of the reconstructed resonance waveform.  */
  virtual TSig getWaveForm(TSig amplitude, TSig phase);
    // todo: maybe change phase format (to 0..1) to make it more compatible with wavetable oscillator


  /** \name User Parameters */

  int waveform;  // preliminary - later use an embedded wavetable-osc that manages the waveform

  TPar resoCutoffMultiplier; // multiplier for the ...
  TPar resoCutoffModulation;

  // embedded objects:
  rsOnePoleFilter<TSig, TPar> resoFilter; // filters the resonance waveform

};


//=================================================================================================

/** Subclass of rsResoReplacer that extends it by adding an optional dynamic "phase-bumping" based
on the relation between the resonance-signal's phase and features of the input signal with the goal
to emulate growl-like effects that occur in ladder filters that involve a saturating feedback path.
*/

// todo: check, if the template type (TSig or TPar) is correct for all member variables

template<class TSig, class TPar>
class rsResoReplacerPhaseBumped : public rsResoReplacer<TSig, TPar>
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsResoReplacerPhaseBumped();


  /** \name Setup */


  void setSampleRate(TPar newSampleRate);  // override


  /** Sets the factor by which the bump function is multiplied. */
  void setBumpFactor(TPar newFactor);
    // rename to setChaosAmount

  /** Sets the cutoff frequencies for the two 2st order lowpass filters which are applied to the
  raw (spikey) bump function. */
  void setBumpCutoffs(TPar cutoff1, TPar cutoff2);

  /** Sets the factor by which the resonance signal is multiplied in the chaos-producing formula
  tmp *= in - c*yr; c is the chaos factor. Its value determines the nature of the chaos. */
  void setChaosParameter(TPar newParameter);


  /** Sets a limiting value for the instantaneous amplitude of the resonance waveform. We apply a
  soft-clipping function to the instantaneous amplitude in order simulate some aspect of the
  nonlinear ladder filter (namely, avoiding too large resonance amplitudes, when the cutoff
  frequency hits a partial of the input signal). The value should be a (strictly) positive number.
  We use the formula y = rsPositiveSigmoids::softClip(x/limit) * limit. */
  void setAmplitudeLimit(TPar newLimit);
    // move to baseclass

  /** Sets the range of input signals for which the resonance signal is (mostly) unaffected.
  Beyond this range, it gets progressively more attenuated to simulate the limiting effect in
  feedback saturated filters. */
  void setInputRange(TPar newRange);
    // move to baseclass

  /** Sets the amount by which the filter gets self-excited in order to emulate
  self-oscillation. In a linear filter, such as this, true self-oscillation is not possible
  because self-oscillation is based on an unstable filter which is re-stbilized by saturation
  (i.e. nonlinearity). */
  void setSelfExitation(TPar newExcitation);
    // this feature is experimental

  //void setResonanceFeedback(double newFeedback);

  /** Sets up a value that is added to the instantaneous amplitude measurement value before
  reconstructing the resonance signal. */
  void setAmplitudeOffset(TPar newOffset);



  /** \name Audio Processing */

  virtual void getSignalParts(TSig in, TSig *yf, TSig *yr);  // override


  // this is supposed to stay only during development, so we can look at the bump function
  TSig getPhaseBump()
  {
    return bump;
  }

protected:


  /** Function for updating the bump value. */
  void updateBump(TSig in, TSig yf, TSig yr);

  /** Returns a limited value for the instantaneous magnitude (from a raw, unlimited measurement
  value) according to our limiting settings here. */
  TSig limitMagnitude(TSig mag);

  /** Returns a multiplier for the resonance signal based on the input signal to simulate the
  bell-shaped envelope that is seen in self-oscillating feedback-saturated filters when there's a
  (sawtooth) input signal. */
  TSig inputScale(TSig in);

  /** Produces the self-excitation signal. */
  TSig getExcitation(TSig in);


  rsOnePoleFilter<TSig, TPar> inputDifferencer;
  rsOnePoleFilter<TSig, TPar> bumpSmoother1, bumpSmoother2;
    // if the concept turns out to be useful, replace the 2 bumpSmoother filters with a single
    // attack/decay filter that allows the user to set attack-time, decay-time and amplitude (i
    // think, i have such a filter somewhere in the legacy codebase)

  rsMovingAverage<TSig, TPar> bumpAverager;
  TPar averageSamples;           // number of samples for the moving averager
    // not used - remove

  TPar bumpCutoff1, bumpCutoff2; // later replace by attack/decay
  TPar bumpFactor;

  TSig bump;  // bump value for next sample
  TPar chaosParam;


  TPar ampLimit, ampLimitInv;
    // instantaneous amplitude limiting value and its reciprocal

  TPar inRange;
    // for the amplitude multiplier based on the input signal

  TPar selfExcite;  // amount of self-excitation
  TSig oldPhase;    // needed to compute self-exitation - rename to exitationPhase


  TPar ampOffset;   // a value that is added to the instananeous amplitude

};


//=================================================================================================

/** This is an implementation of a nonlinear model of a transistor ladder filter using
topology-preserving transform (TPT) and zero-delay feedback (ZDF) technology. The core of the
code is based on an implementation by Teemu Voipio published in the KVR-forum:
[...]

*/

template<class TSig, class TPar>
class rsLadderMystran : public rsLadderFilter2<TSig, TPar>
{

public:

  rsLadderMystran();
  RS_INLINE TSig getSample(TSig in);
  void reset();

  /** tanh(x)/x approximation, flatline at very high inputs so might not be safe for very large
  feedback gains [limit is 1/15 so very large means ~15 or +23dB] */
  RS_INLINE TSig tanhXdX(TSig x);

protected:

  TSig s[4]; // state variables
  TSig zi;   // previous input

};

//// LICENSE TERMS: Copyright 2012 Teemu Voipio
//
// You can use this however you like for pretty much any purpose,
// as long as you don't claim you wrote it. There is no warranty.
//
// Distribution of substantial portions of this code in source form
// must include this copyright notice and list of conditions.
//
template<class TSig, class TPar>
RS_INLINE TSig rsLadderMystran<TSig, TPar>::getSample(TSig in)
{
  // tuning and feedback
  TSig f = tan(PI * this->cutoff / this->sampleRate);  // move this calculation to setCutoff
  //TSig resonance = 0.5;                    // preliminary - make member, provide setter
  TSig r = (40.0/9.0) * this->resonance;         // maybe use 40/10 = 4 instead -> self-osc at r=1


  // some DC to get even harmonics (by RS):
  TSig dc = 0.1;
  in += dc;


  // input with half delay, for non-linearities:
  TSig ih = 0.5 * (in + zi); zi = in;

  // evaluate the non-linear gains:
  TSig t0 = tanhXdX(ih - r * s[3]);
  TSig t1 = tanhXdX(s[0]);
  TSig t2 = tanhXdX(s[1]);
  TSig t3 = tanhXdX(s[2]);
  TSig t4 = tanhXdX(s[3]);

  // g# the denominators for solutions of individual stages:
  TSig g0 = 1 / (1 + f*t1), g1 = 1 / (1 + f*t2);
  TSig g2 = 1 / (1 + f*t3), g3 = 1 / (1 + f*t4);

  // f# are just factored out of the feedback solution:
  TSig f3 = f*t3*g3, f2 = f*t2*g2*f3, f1 = f*t1*g1*f2, f0 = f*t0*g0*f1;

  // solve feedback:
  TSig y3 = (g3*s[3] + f3*g2*s[2] + f2*g1*s[1] + f1*g0*s[0] + f0*in) / (1 + r*f0);

  // then solve the remaining outputs (with the non-linear gains here):
  TSig xx = t0*(in - r*y3);
  TSig y0 = t1*g0*(s[0] + f*xx);
  TSig y1 = t2*g1*(s[1] + f*y0);
  TSig y2 = t3*g2*(s[2] + f*y1);

  // update state:
  s[0] += 2*f * (xx - y0);
  s[1] += 2*f * (y0 - y1);
  s[2] += 2*f * (y1 - y2);
  s[3] += 2*f * (y2 - t4*y3);

  return this->c[0]*xx + this->c[1]*y0 + this->c[2]*y1 + this->c[3]*y2 + this->c[4]*y3; // maybe use c0*in instead?
}

template<class TSig, class TPar>
RS_INLINE TSig rsLadderMystran<TSig, TPar>::tanhXdX(TSig x)
{
  TSig A = x*x;
  return ((A + 105)*A + 945) / ((15*A + 420)*A + 945);
}
// maybe make function static

#endif
