#pragma once

// template parameters: 
// T: type for signals and parameters, TOsc: class for the (blep-ready) osc, TBlep: class for blep
template<class T, class TOsc, class TBlep> 
// todo: have separate TSig, TPar template parameters

/** A class for producing waveforms like supersaw, supersquare, etc. using a blep-ready oscillator 
class as basis and applies a blep for anti-aliasing. 

todo: factor out a baseclass rsOscArray that is responsible for handling the array of increments,
i.e. contains all the members under "increment handling stuff - then we may later also have other 
subclasses, like a wavetable-based osc-array, etc. 

maybe avoid calling updateIncrements in the setters (or make it optional using a 2nd boolean 
parameter) - instead keep a std::atomic_bool incsUpToDate flag, set it false in the setters (and 
true in updateIncrements) check it in getSampleNaive and if it's false, call updateIncrements there
->good for simultaneous modulation of several parameters and thread-safe parameter changes (inc 
recalc will always be done in the audio thread and multiple parameter modulations per sample will 
not lead to multiple calls to updateIncrements (see rosic::TurtleSource/Sbowflake - there, i do it 
that way)  */ 

class rsBlepOscArray
{

public:


  rsBlepOscArray(RAPT::rsRatioGenerator<T>* ratioGenerator = nullptr);
  // get rid of passing the pointer to the ratio-generator

  /** Sets the reference phase increment which determines the center frequency of the osc stack 
  around which all other frequencies are arranged. */
  inline void setReferenceIncrement(T newIncrement) 
  { 
    refInc = newIncrement; 
    updateIncrements();
  }

  /** Sets the number of oscillators to use. */
  void setDensity(int newDensity)
  {
    rsAssert(newDensity <= getMaxDensity());
    numOscs = rsMin(newDensity, getMaxDensity());
    updateIncrements();
    // we need to update the increments here because the new number may be higher than the old and
    // the upper values may not yet contain valid data becuase on setReferenceIncrement, setDetune, 
    // etc. we only update the oscs up to numOscs in order to avoid computing increments that are 
    // not used
  }
  // rename to setDensity

  void setDetune(T newDetune) 
  { 
    detune = newDetune; 
    updateIncrements();
  }
  // todo: use an update strategy similar to TurtleSource - have an atomic bool that stores, if the
  // incs array is up to date, check in getSample if it is up to date and if not, update it - 
  // allows for efficient simultaneous modulation of frequency and detune from the mod-system
  // or maybe just have a function setIncrementAndDetune - to make it efficient in the mod-system,
  // a subclass (in rosic) shall be used that uses the bool - rapt is not the right place for such
  // infrastructe dependent decisions

  /** Sets up, how coherent the phases/positions of the phasors shall be after calling reset. A 
  value of 0 spreads the phases evenly in the phasor-interval 0..1, a value of 1 lets all saws 
  start coherently at the same phase zero, which will lead to the effect that there's a noticable
  attack transient at the start of the sound. */
  void setIntialPhaseCoherence(T newCoherence)
  {
    startPhaseDist = T(1) - newCoherence;
  }


  /** Sets the type of the frequency distribution that is used to arrange the individual saws 
  around the center frequency. */
  // why double? should be T!
  void setFrequencyDistribution(rsRatioGenerator<double>::RatioKind newDistribution)
  {
    rsAssert(ratioGenerator != nullptr); 
    ratioGenerator->setRatioKind(newDistribution);
    updateIncrements();
  }

  void setDistributionParameter(T newParam)
  {
    rsAssert(ratioGenerator != nullptr);
    ratioGenerator->setParameter1(newParam);
    updateIncrements();
  }

  /** Sets the maximum number of oscillators that can be used. May cause memory (re)allocation and
  should probably be called once at startup time. */
  void setMaxDensity(int newMaximum)
  {
    incs.resize(newMaximum);
    oscs.resize(newMaximum);
    numOscs = rsMin(numOscs, newMaximum);
    updateIncrements();
  }

  /** Returns the maximum number of oscillators that can be used. */
  int getMaxDensity() const { return (int) oscs.size(); } 


  inline T getSampleNaive()
  {
    T stepDelay = T(0);   // delay of step discontinuity
    T stepAmp   = T(0);   // amplitude of step discontinuity
    T out       = T(0);   // accumulator for (naive) output sample
    T amp = T(1) / sqrt(numOscs);  // todo: precompute
    for(int i = 0; i < numOscs; i++) {

      out += oscs[i].getSampleSaw(incs[i], &stepDelay, &stepAmp);
      // hmm..to make it work flexibly with other types of oscs, we need to call a generic 
      // getSample function - but then the osc would have to dispatch...based on what? maybe the
      // getSample function should take an int parameter?

      stepAmp *= amp;
      if(stepAmp != T(0))  
        blep.prepareForStep(stepDelay, stepAmp);
      // maybe try to do it without the branch and make performance test with both versions - it 
      // doesn't hurt to accumulate zero-valued signals into the corrector, so the code is valid
      // with or without the "if"
    }
    return amp * out;  
    // todo: scale the amplitude by 1/sqrt(numOscs) ..oh - but that factor should then also be 
    // applied to the stepAmp
  } 

  inline T getSample()
  {
    return blep.getSample(getSampleNaive());
  }

  void reset()
  {
    for(int i = 0; i < numOscs; i++) {
      T pos = startPhaseDist * T(i) / T(numOscs); 
      oscs[i].resetPhase(pos, incs[i]);
    }
    blep.reset();
  }

  // todo: have a reset function that allows to select an initial phase distribution...or maybe
  // just call it setPhases - it may take a parameter where 0 means all start at 0 and 1 means 
  // maximally sperad out (i.e. oscs[i].pos = i / numOscs

protected:

  /** Updates our array of phase increments according to the desired reference increment and 
  settings of detune, etc.. */
  void updateIncrements();

  // increment handling stuff:
  int numOscs = 1;        // current number of oscillators
  T refInc = T(0);        // reference increment
  T detune = T(0);
  T startPhaseDist = T(1); // distribution of start phases - 0: coherent, 1: maximally incoherent
  std::vector<T> incs;     // array of phase increments


  rsRatioGenerator<T>* ratioGenerator = nullptr; 
  // used to generate freq-ratios, shared among voices ..hmm - but maybe we should use a direct 
  // object and for the voices use the strategy outlined in Notes/Ideas.txt



  std::vector<TOsc> oscs; // oscillator array
  TBlep blep; 
  // a single blep is shared among all oscs...actually, it could even be shared among all voices
  // in a polyphonic situation - try to think of a way, how to do this - maybe by maintaining a 
  // pointer to the blep object instead of a direct object? ...we'll see...hmm...maybe that doesn't
  // make much sense, because in a synth, the osc-voices are not mixed before the filter - and 
  // applying the blep after the filter is invalid because the blep needs go through the filter, 
  // too - so it's probably best to keep things as is



  // todo: have an array of pan positions - if we do stereo later, we will need two blep objects 
  // - one for each channel
};


