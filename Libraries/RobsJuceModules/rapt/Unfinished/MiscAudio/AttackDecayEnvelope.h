#pragma once  // use ifdef/define


//=================================================================================================

/** A filter that has an attack/decay shape as its impulse response. The user can adjust the attack
and decay time. It is based on the difference of two exponential decays with different time
constants. */

template<class T>
class rsAttackDecayFilter
{

public:


  rsAttackDecayFilter();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the attack time in samples, i.e. the number of samples it takes to reach the peak. */
  void setAttackSamples(T newAttack) { attackSamples = newAttack; coeffsDirty = true; }

  /** Sets the decay time constant in samples, i.e. the number of samples, it takes to decay to
  1/e for the more slowly decaying exponential. */
  void setDecaySamples(T newDecay) { decaySamples = newDecay; coeffsDirty = true; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the gain of this filter at DC. This value can be useful to know when you want to 
  create an envelope with sustain. You may then feed the reciprocal of that value as constant 
  input into the filter. */
  T getGainAtDC() const; // { return s*(cd-ca) / (T(1)+ca*cd-cd-ca); }

  /** Computes the reciprocal of the DC gain of this filter. This is the constant value, you want 
  to feed into the filter when you want to get a sustained output of 1. It's used for implementing 
  sustain in subclass rsAttackDecayEnvelope. */
  T getReciprocalGainAtDC() const; // { return (T(1)+ca*cd-cd-ca) / (s*(cd-ca)); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  T getSample(T in)
  {
    if(coeffsDirty)
      updateCoeffs();
    ya = in + ca * ya;
    yd = in + cd * yd;
    return s * (yd - ya);
  }

  /** Resets the internal state of both filters to zero. */
  void reset()
  {
    ya = yd = T(0);
  }

  /** Updates our filter coefficients according to user parameters. */
  void updateCoeffs();


protected:


  // Data:
  T ca, cd;  // coefficients for attack and decay filters
  T ya, yd;  // state of attack and decay filter
  T s;
  bool coeffsDirty = true;  // flag to indicate that coeffs need to be re-computed
  //std::atomic<bool> coeffsDirty = true;  // flag to indicate that coeffs need to be re-computed
  T attackSamples = T(20), decaySamples = 100; // sort of arbitrary
};
// maybe move into the same file where the modal filters are - it can be seen as a "modal" filter
// with zero frequency so it would fit in there

//=================================================================================================

/** An envelope generator based on rsAttackDecayFilter. It feeds the filter with a mix of a unit
impulse and a constant value, where the latter is responsible for a sustain phase */

template<class T>
class rsAttackDecayEnvelope : public rsAttackDecayFilter<T>
{

public:

  using Base = rsAttackDecayFilter<T>; // for conveniently calling basclass methods

  enum class AccuFormula // rename to AccumulationMode or RetriggerMode
  {
    none,        
    // Find better name - we do no scaling of the impulse at all. It has the effect that retriggers
    // pile up without any countermeasures. Maybe if we call this enum RetriggerMode, this mode 
    // could be called accumulative

    exact,       
    // Not yet usable - sometimes the iteration diverges.

    one_minus_yd
    // Find better name, maybe compByDec (for compensate by the decay feedback state)
    // Scales the impulse by (1-yd). Rationale: When the so scaled impulse is added back to the
    // yd value which is received in the feedback path, they add up to unity or at least almost. 
    // Maybe the yd value is one sample to old or new for the formula to be exact -> check that.
  };

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sustain level. This is the value that is added to the filter's input as long a note
  is being held. */
  void setSustain(T newSustain) { sustain = newSustain; }

  void setAccumulationMode(AccuFormula newMode)
  {
    accuFormula = newMode;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the constant value that should be fed into the filter during the sustain phase. */
  T getSustainInput() const { return sustain * Base::getReciprocalGainAtDC(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Event Handling */

  //void trigger();  // maybe let it take a "strength" parameter?

  void noteOn(int key, int vel)
  {
    currentNote = key; 
    // Maybe we should store the velocity, too? It coul be useful if we want to scale the sustain 
    // level by the velocity.

    T x = getAccuCompensatedImpulse();
    // Maybe x should be scaled according to key and vel? But maybe it should be the job of a 
    // higher level class to compute an overall strength which can be passed to this function. This
    // "strength" should include the key- and velocity scaling and directly multiply x here.

    Base::getSample(x);
    // We call geSample here to avoid the one sample delay due to the subtraction - the two
    // exponentials cancel each other at the very first sample
  }
  // Should we take into account the sustain here? The maximum excursion will increase, if we have
  // nonzero sustain. Can we compensate that by feeding some number other than 1? And if so, should 
  // we do it? Maybe that coupling is musically desirable?
  // Also: when anote-on is received before the envlope has decayed away, the overlapping impulse
  // responses will add up - can we scale the impulse in a way to compentsate this accumulative 
  // behavior (by a formula taking into account the previous output) - or better, continuously
  // blend between accumulative and non-accumulative behavior


  void noteOff(int key, int vel)
  {
    currentNote = -1;
  }
  // why do we expect key and vel parameters? ...maybe we can do something with the note-off 
  // velocity later, like adjusting the release time?


  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  T getSample()
  {
    //if(currentNote != -1)  return rsAttackDecayFilter<T>::getSample(sustain);
    if(currentNote != -1)  
      return Base::getSample(getSustainInput());
    else                   
      return Base::getSample(T(0));
  }
  // maybe we should have a conditional if(sustain == 0) to call simpler code, when sustain is not 
  // used - sustain calls getReciprocalGainAtDC which is moderately expensive

  void reset()
  {
    Base::reset();
    currentNote = -1;
  }


protected:


  /** Under construction. */
  T getExactAccuCompensation();

  /** Returns the desired input impulse height taking into account the current state of the filter.
  The goal is to compensate for the increasing height of the peaks when note-on events are received
  when the output has not yet decayed away. Without compensation, a quick succession of note-on 
  events will make the peak height grow. This function implements various formulas to be applied
  to the input impulse to counteract that effect. */
  T getAccuCompensatedImpulse();


  AccuFormula accuFormula = AccuFormula::none; // rename to retrigMode

  T sustain = T(0);
  int currentNote = -1;  // -1 is code for "none"

};

/*
todo:
-later extend it to allow for a two-stage decay (maybe in a subclass)
-make it polyphonic - allow to build a very basic subtractive synth from such envelopes, a simple
osc class (maybe TriSaw osc?) and the rsLadder filter
-the filter cutoff, osc-frequency and overall amplitude should be (polyphonically) modulated by
this simple envelope
-maybe also have some simple LFO class
*/









//=================================================================================================

// maybe this class should not be part of rapt but be moved to rosic - maybe also for the rsVoice
// and related classes

template<class T>
class rsAttackDecayEnvelopeVoice : public rsVoice
{

public:

  virtual void noteOn(int key, int vel) override;
  virtual void noteOff(int key, int vel) override;

protected:

  // template object containing shared state:
  rsAttackDecayEnvelope<T>* master;

  // per voice state:
  T ya, yd, ca, cd;  // state and coeffs for attack and decay

};

template<class T>
class rsTriSawVoice : public rsVoice
{

public:

  virtual void noteOn(int key, int vel) override;
  virtual void noteOff(int key, int vel) override;

protected:

  // template object containing shared state:
  rsTriSawOscillator<T>* master;

  // per voice state:
  T p = 0;     // current phase in 0..1
  T inc = 0;   // phase increment

};


template<class TSig, class TPar>
class rsLadderVoice : public rsVoice
{

public:

  virtual void noteOn(int key, int vel) override;
  virtual void noteOff(int key, int vel) override;

protected:

  // template object containing shared state:
  rsLadderFilter<TSig, TPar>* master;

  // per voice state:
  TSig y[5];   // outputs of the stages 0..4
  TPar a, b;   // leaky integrator coefficients for a stage: y[n] = b*x[n] - a*y[n-1]
  TPar k;      // feedback gain
  TPar g;      // output gain

};

// how can we handle to avoid duplicating the algorithmic code from the getSample functions of the
// underlying monophonic DSP classes ...and also all the coefficient calculation code? maybe by
// factoring the code out into (static) functions that do not operate on member data but instead
// get all their inputs and states as arguments - but that would uglify/complicate the monophonic
// implementations - something, i'd like to avoid - hmmm...but maybe, it's inavoidable?
// The rsLadderFilter class already has these static computeCoeffs functions - so it's already
// prepared for it - we would also need a static processSample function that works in a similar
// way




