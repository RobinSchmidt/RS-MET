#ifndef RAPT_PHASEDITHEROSCS_H_INCLUDED
#define RAPT_PHASEDITHEROSCS_H_INCLUDED

//=================================================================================================

/** A realtime oscillator that produces pitch-dithered waveforms. ...TBC... */

template<class T> 
class rsPitchDitherOsc
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Default constructor. It puts the object into a valid initial state by setting up a default
  period length of 100.0 samples and triggering the appropriate computations to set up our member
  variables that control the distribution of cycle lengths. */
  rsPitchDitherOsc();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the period, i.e. the desired length (in samples) of one cycle of the waveform. This is 
  a floating point value and it can be computed as  period = sampleRate / frequency  where 
  frequency is the desired oscillator frequency in Hz. This will immediately trigger a 
  recomputation of the probability distribution of the cycle lengths and update the currently used
  cycle length. */
  void setPeriod(T newPeriod) { setPeriodNoUpdate(newPeriod); updateCycleLength(); }

  /** Sets up a new period length just like setPeriod() does but without immediately updating the
  probability distribution and current cycle length. This results in the behavior that the new 
  period will not become effective immediately but only after finishing the currently running 
  cycle. */
  void setPeriodNoUpdate(T newPeriod)
  { calcCycleDistribution(newPeriod, &midLength, &probShort, &probMid); }

  /** Sets the seed for the pseudo random number generator. */
  void setRandomSeed(uint32_t newSeed) { seed = newSeed; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the average length of the cycles that are being produced. In general, if we have 3 
  integer cycle lengths given by L1,L2,L3 and cycles with these 3 lengths are produced with 
  probabilities p1,p2,p3 respectively, then the average cycle length P will be: 
  P = p1*L1 + p2*L2 + p3*L3. In our particular case here, we will have a given middle length L2 and
  the short and long lengths L1,L3 are then given as L1 = L2-1, L3 = L2+1. And, of course, for the 
  probabilities we will always have p1 + p2 + p3 = 1 such that p3 = 1 - (p1 + p2). */
  T getPeriod();
  // Needs tests. Maybe move explanation of the formula into implementation.

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns a sample of a phasor value, i.e. a value in the range 0..1 that can be used to create
  various waveforms. ...TBC... */
  inline T getSamplePhasor();
  // ToDo: Document, if the produced value is in the closed interval [0,1] or in the half-open 
  // interval [0,1). I think, it's the former and I also think that this might not really be the 
  // right thing to do - at least not for producing sine-waves. For producing saws, it may actually
  // e appropriate, though. I think, in general, whenever there is a jump discontinuity at the 
  // wrap-around point, we may want the closed interval and otherwise the half-open one.


  /** Returns a sample of an upward sawtooth wave */
  inline T getSampleSawUp()   { return T(-1) + T(2) * getSamplePhasor(); }

  /** Returns a sample of an downward sawtooth wave */
  inline T getSampleSawDown() { return T(+1) - T(2) * getSamplePhasor(); }

  /** Returns a sample of a pulse wave with given pulse-width. The default value of 0.5 produces a 
  square wave. */
  inline T getSamplePulse(T pw = T(0.5));

  /** Resets the internal state, i.e. the sample counter and the random generator. */
  inline void reset();

  //-----------------------------------------------------------------------------------------------
  // \name Helpers

  static void calcCycleDistribution(T period, T* midLength, T* probShort, T* probMid);


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internals

  /** Updates our cycleLength member by computing a new (pseudo) random cycle length to be used 
  for the next cycle. Called from getSample() after each cycle has been completed. */
  inline void updateCycleLength();

  //-----------------------------------------------------------------------------------------------
  // \name Data

  // Members that are accessed per sample:
  T phaseSlope;    // Phase increment per sample.
  T sampleCount;   // Is always in interval [0, cycleLength).
  T cycleLength;   // Is midLength or midLength + 1 or midLength - 1.

  // Members that are accessed per cycle:
  T midLength;     // The middle one of the 3 cycle lengths to be produced.
  T probShort;     // Probability to use midLength - 1.
  T probMid;       // Probability to use midLength.
  // Maybe rename cycleLength and midLength in lengthCurrent and lengthMid. Rationale: lengthMid 
  // would be more consistent with probShort and probMid (the "Mid" would be the suffix). Maybe
  // use lenCurrent (or lenNow) and lenMid. Maybe rename sawSlope to phaseSlope

  // Embedded DSP objects:
  rsRandomGenerator<T> prng;
  uint32_t seed = 0;

};

template<class T> 
inline T rsPitchDitherOsc<T>::getSamplePhasor()
{
  T p = phaseSlope * sampleCount;  // Compute output sample.
  sampleCount += T(1);             // Update counter. We will have produced 1 sample.
  if(sampleCount >= cycleLength)   // Is cycle finished?
  {                                // If so..
    sampleCount = T(0);            // ..Wrap around sample counter.
    updateCycleLength();           // ..Compute cycleLength and sawSlope for next cycle.
  }
  return p;                        // Return output sample.
}

template<class T> 
inline void rsPitchDitherOsc<T>::updateCycleLength()
{
  T r = prng.getSampleInUnitRange();         // Random number in interval [0,1).
  if(r < probShort)
    cycleLength = midLength - T(1);          // Next cycle is short.
  else if(r < probShort + probMid)
    cycleLength = midLength;                 // Next cycle is medium.
  else
    cycleLength = midLength + T(1);          // Next cycle is long.
  phaseSlope = T(1) / (cycleLength - T(1));  // Slope depends on cycle length.
}

template<class T> 
inline T rsPitchDitherOsc<T>::getSamplePulse(T pw) 
{ 
  T p = getSamplePhasor();
  if(p < pw)
    return T(-1);
  else
    return T(+1);

  // ToDo: Verify that this formula is what the user would expect. Maybe we should swap -1 and
  // +1? And/or maybe we should use if(p <= pw) rather than if(p < pw). Document these decsisions
  // and the reasons behind them. One reason to prefer to have the negative half-cycle first is
  // that this would be compatible with clipping a saw-up waveform and I think, the "up" variant
  // is the default expectation in case of a saw wave. Check what popular synthesizers do (Surge,
  // Serum, Diva, JP-8000, ...) and maybe do the same. Maybe to gigure out if < or <= is correct,
  // consider a square wave with an even integer cycle length. In such a case, we want the positive
  // and negative half-wave to have exactly the same number of samples. This may also depend on 
  // whether the phasor range is [0,1] or [0,1). 
}

template<class T> 
inline void rsPitchDitherOsc<T>::reset()
{
  sampleCount = T(0);
  prng.setState(seed);
  updateCycleLength();         // Important for correct initial cycleLength.  
}


#endif