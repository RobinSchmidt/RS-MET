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
  create an envelope with sustain - you may then feed the reciprocal of that value as constant 
  input into the filter. */
  T getGainAtDC() const 
  { 
    return s*(cd-ca) / (T(1)+ca*cd-cd-ca);
  }

  // todo: figure out the formula for the DC gain - should depend on ca,cd,s ..somthing like
  // s * (gd - ga) where gd, ga are the DC gains of the decay and attack filter...was it
  // gd = 1 / (1 + cd), ga = 1 / (1 + ca)? ...look it up
  // this value is needed for setting up a sustain...
  // todo: simplify to have only one division

  T getReciprocalGainAtDC() const
  {
    return (T(1)+ca*cd-cd-ca) / (s*(cd-ca));

    //return T(1) / getGainAtDC();
    // todo: this can be algebraically simplified such that we nee only 1 division instead of 3
    // -> do it
  }

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

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sustain level. This is the value that is added to the filter's input as long a note
  is being held. */
  void setSustain(T newSustain) { sustain = newSustain; }
  // doesn't work yet - we probably need to scale the sustain input according to the DC gain of the
  // filter


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  T getSustainInput() const { return sustain * Base::getReciprocalGainAtDC(); }
  // can the formula be optimized, i.e. algebraically simplified? 

  //-----------------------------------------------------------------------------------------------
  /** \name Event Handling */

  //void trigger();  // maybe let it take a "strength" parameter?

  void noteOn(int key, int vel)
  {
    currentNote = key;
    Base::getSample(T(1)); // maybe input should be scaled by vel?
    // We call geSample here to avoid the one sample delay due to the subtraction - the two
    // exponentials cancel each other at the very first sample
  }
  // Should we take into account the sustain here? the maximum excursion will increase, if we have
  // nonzero sustain. Can we compensate that by feeding some number other than 1? And if so, should 
  // we do it? Maybe that coupling is musically desirable?


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


protected:

  T sustain = T(0);
  int currentNote = -1;  // -1 is code for "none"

};

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




/*
todo:
-later extend it to allow for a two-stage decay (maybe in a subclass)
-make it polyphonic - allow to build a very basic subtractive synth from such envelopes, a simple
 osc class (maybe TriSaw osc?) and the rsLadder filter
 -the filter cutoff, osc-frequency and overall amplitude should be (polyphonically) modulated by
  this simple envelope
 -maybe also have some simple LFO class
*/
