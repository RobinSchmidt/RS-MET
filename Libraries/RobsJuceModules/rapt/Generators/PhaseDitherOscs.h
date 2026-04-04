#ifndef RAPT_PHASEDITHEROSCS_H_INCLUDED
#define RAPT_PHASEDITHEROSCS_H_INCLUDED

//=================================================================================================

/** A realtime oscillator that produces pitch-dithered waveforms. ...TBC...

ToDo: Explain the idea of pitch dithering

Warning: This class is not yet well tested and should be considered rather preliminary. Some 
details of the implementation may change. There may be bugs. */

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
  void setPeriod(T newPeriod, bool phasorRangeClosed) 
  { 
    setPeriodNoUpdate(newPeriod); 
    updateCycleLength(phasorRangeClosed);
  }
  // Maybe it should take a bool parameter "phasorRangeClosed"

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

  /** Returns the average length of the cycles that are being produced. */
  T getPeriod();
  // Needs tests.

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns a sample of a phasor value, i.e. a value in the range 0..1 that can be used to create
  various waveforms. ...TBC... */
  inline T getSamplePhasor(bool phasorRangeClosed);
  // ToDo: Document, if the produced value is in the closed interval [0,1] or in the half-open 
  // interval [0,1). I think, it's the former and I also think that this might not really be the 
  // right thing to do - at least not for producing sine-waves. For producing saws, it may actually
  // e appropriate, though. I think, in general, whenever there is a jump discontinuity at the 
  // wrap-around point, we may want the closed interval and otherwise the half-open one.


  /** Returns a sample of an upward sawtooth wave */
  inline T getSampleSawUp()   { return T(-1) + T(2) * getSamplePhasor(true); }

  /** Returns a sample of an downward sawtooth wave */
  inline T getSampleSawDown() { return T(+1) - T(2) * getSamplePhasor(true); }

  /** Returns a sample of a pulse wave with given pulse-width. The default value of 0.5 produces a 
  square wave. */
  inline T getSamplePulse(T pw = T(0.5));

  /** Resets the internal state, i.e. the sample counter and the random generator. */
  void reset(bool phasorRangeClosed);

  //-----------------------------------------------------------------------------------------------
  // \name Helpers

  static void calcCycleDistribution(T period, T* midLength, T* probShort, T* probMid);


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internals

  /** Updates our sampleCount member and takes care of appropriate wrap around with recomputation 
  of the new cycle length. */
  inline void updateSampleCount(bool phasorRangeClosed);

  /** Updates our cycleLength member by computing a new (pseudo) random cycle length to be used 
  for the next cycle. Called from updateSampleCount() after each cycle has been completed. */
  inline void updateCycleLength(bool phasorRangeClosed);

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
  // use lenCurrent (or lenNow) and lenMid. Maybe rename sawSlope to phaseSlope. Yes -  I think
  // lengthMid and lengthNow are the best choices.

  // Embedded DSP objects:
  rsRandomGenerator<T> prng;
  uint32_t seed = 0;

};

template<class T> 
inline T rsPitchDitherOsc<T>::getSamplePhasor(bool closed)
{
  T p = phaseSlope * sampleCount;  // Compute output sample.
  updateSampleCount(closed);       // Update sample counter. Possibly wraps around.
  return p;                        // Return output sample.
}

template<class T> 
inline void rsPitchDitherOsc<T>::updateSampleCount(bool closed)
{
  sampleCount += T(1);             // Update counter. We produce 1 sample at each update.
  if(sampleCount >= cycleLength)   // Is cycle finished?
  {                                // If so..
    sampleCount = T(0);            // ..Wrap around sample counter.
    updateCycleLength(closed);     // ..Compute cycleLength and phaseSlope for next cycle.
  }
}

template<class T> 
inline void rsPitchDitherOsc<T>::updateCycleLength(bool closed)
{
  T r = prng.getSampleInUnitRange();              // Random number in interval [0,1).
  if(r < probShort)
    cycleLength = midLength - T(1);               // Next cycle is short.
  else if(r < probShort + probMid)
    cycleLength = midLength;                      // Next cycle is medium.
  else
    cycleLength = midLength + T(1);               // Next cycle is long.
  phaseSlope = T(1) / (cycleLength - T(closed));  // Slope depends on cycle length.

  // Maybe as an optimization, pass the "closed" parameter not as bool but as type T so we can 
  // avoid the type conversion. Maybe we should assert that the value represents either T(0) or
  // T(1). Maybe add a function rsIsBoolean(T x) to the library that returns true iff x is 0 or 1
  // and use that function in a rsAssert here. Maybe have static const members phasorRangeClosed, 
  // phasorRangeHalfOpen of type T that are fixed to 0 and 1 such that the caller can uses these
  // as symbolic constants rather than having itself to make sure to only pass 0 or 1.
}

template<class T> 
inline T rsPitchDitherOsc<T>::getSamplePulse(T pw) 
{ 
  T p = getSamplePhasor(true);
  if(p < pw)
    return T(-1);
  else
    return T(+1);
}

template<class T> 
void rsPitchDitherOsc<T>::reset(bool closed)
{
  sampleCount = T(0);
  prng.setState(seed);
  updateCycleLength(closed);                 // Important for correct initial cycleLength.  
}


#endif