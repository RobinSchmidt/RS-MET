#pragma once




/** A filter that has an attack/decay shape as its impulse response. The user can adjust the attack
and decay time. It is based on the difference of two exponential decays with different time 
constants. */

template<class T>
class rsAttackDecayFilter
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the attack time in samples, i.e. the number of samples it takes to reach the peak. */
  void setAttackSamples(T newAttack) { attackSamples = newAttack; coeffsDirty = true; }

  /** Sets the decay time constant in samples, i.e. the number of samples, it takes to decay to 
  1/e for the more slowly decaying exponen */
  void setDecaySamples(T newDecay) { decaySamples = newDecay; coeffsDirty = true; }





  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  T getSample(T in)
  {
    if(coeffsDirty)
      updateCoeffs();

    ya = in + ca * ya;
    yd = in + cd * yd;
    return yd - ya;
    // ....verify these formulas
  }

  /** Resets the internal state of both filters to zero. */
  void reset()
  {
    ya = yd = T(0);
  }

protected:


  /** Updates our filter coefficients according to user parameters. */
  void updateCoeffs();

  // Data:
  T ca, cd;  // coefficients for attack and decay filters
  T ya, yd;  // state of attack and decay filter
  bool coeffsDirty = true;  // flag to indicate that coeffs need to be re-computed
  //std::atomic<bool> coeffsDirty = true;  // flag to indicate that coeffs need to be re-computed
  T attackSamples = T(20), decaySamples = 100; // sort of arbitrary
};

//=================================================================================================

/** An envelope generator based on rsAttackDecayFilter. It feeds the filter with a mix of a unit 
impulse and a constant value, where the latter is responsible for a sustain phase */

template<class T>
class rsAttackDecayEnvelope : public rsAttackDecayFilter<T>
{

public:



  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sustain level. This is the value that is added to the filter's input as long a note
  is being held. */
  void setSustain(T newSustain) { sustain = newSustain; }

  //-----------------------------------------------------------------------------------------------
  /** \name Event Handling */

  //void trigger();  // maybe let it take a "strength" parameter?

  void noteOn(int key, int vel)
  {
    currentNote = key;
    impulse = T(1);   // maybe it should be scaled by vel?
    // ...and/or maybe we should call geSample once? (to avoid the one sample delay due to the 
    // subtraction - the two exponentials cancel each other at the very first sample)
  }

  void noteOff(int key, int vel)
  {
    currentNote = -1;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  T getSample()
  {
    T y(0);
    if(currentNote != -1)
      y = rsAttackDecayFilter::getSample(sustain+impulse);
    else
      y = rsAttackDecayFilter::getSample(impulse);
    impulse = T(0);
    return y;
  }





protected:


  T sustain = T(0);
  T impulse = T(0);
  int currentNote = -1;  // -1 is code for "none"

};

//=================================================================================================

template<class T>
class rsAttackDecayEnvelopePoly //: public rsVoice
{

public:


protected:

  rsAttackDecayEnvelope* master;

  T y, ca, cd;  // state and coeffs for attack and decay


};


/* 
todo: 
-implement  an attack/decay envelope based on a difference of exponentials
-use the formulas from the modal filter - there's code to compute the time-constants and weights 
 from attack/decay settings
-later extend it to allow for a two-stage decay (maybe in a subclass)
-make it polyphonic - allow to build a very basic subtractive synth from such envelopes, a simple 
 osc class (maybe TriSaw osc?) and the rsLadder filter
 -the filter cutoff, osc-frequency and overall amplitude should be (polyphonically) modulated by 
  this simple envelope
 -maybe also have some simple LFO class
*/