#ifndef RAPT_PITCHDITHEROSCS_H_INCLUDED
#define RAPT_PITCHDITHEROSCS_H_INCLUDED

//=================================================================================================

/** A realtime oscillator that produces pitch-dithered waveforms. ...TBC...

ToDo:

- Explain the idea of pitch dithering. Refer to the documents that I wrote up about the idea. They
  are currently in draft state, though. Document the phasorRangeClosed parameters that occur in 
  various places.

- Document what makes sense for the type T. I think, only scalar floating point types (i.e. 
  float, double, long double, etc.) are meaningful.

*/

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

  /** Sets the mean period (aka cycle length), i.e. the desired length (in samples) of one cycle of
  the waveform. The actually produced periods will be integers that straddle the desired mean
  period by following a suitable probability distribution. The mean period is a floating point value
  and it can be computed as  meanPeriod = sampleRate / frequency  where frequency is the desired
  oscillator frequency in Hz. This will immediately trigger a recomputation of the probability
  distribution of the cycle lengths and update the currently used cycle length. */
  void setMeanCycleLength(T newLength, bool phasorRangeClosed);

  /** Sets up a new period length just like setPeriod() does but without immediately updating the
  probability distribution and current cycle length. This results in the behavior that the new 
  period will not become effective immediately but only after finishing the currently running 
  cycle. */
  void setMeanCycleLengthNoUpdate(T newLength);

  /** Sets the seed for the pseudo random number generator. */
  void setRandomSeed(uint32_t newSeed) { seed = newSeed; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the average length of the cycles that are being produced. It is given by a weighted
  sum of the 3 integer cycle lengths that are produced where the weights are given by their 
  respective probabilities. The formula is cM = p1*c1 + p2*c2 + p3*c3 where cM is the mean cycle
  length, c1,c2,c3 are the actually produced integer cycle lengths and p1,p2,p3 are their
  associated probabilities. */
  T getMeanCycleLength();
  // Needs tests.

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Returns a sample of a phasor value, i.e. a value in the range 0..1 that can be used to create
  various waveforms. Via the parameter "phasorRangeClosed", callers can decide if they want the 
  phasor to be produced in the closed interval [0,1] or in the half-open interval [0,1). See the
  documentation of class rsWaveForms for some discussion in which circumtances one may want to opt
  for one or the other variant. */
  inline T getSamplePhasor(bool phasorRangeClosed);

  // Abbreviation for convenience in the functions below:
  using WF = rsWaveForms<T>;

  /** Returns a sample of an upward sawtooth wave */
  inline T getSampleSaw() { return WF::saw(getSamplePhasor(true)); }

  /** Returns a sample of a pulse wave with given pulse-width. The default value of 0.5 produces a 
  square wave. */
  inline T getSamplePulse(T pw = T(0.5)) { return WF::pulse(getSamplePhasor(true), pw); }

  /** Returns a sample of a sine wave. */
  inline T getSampleSine() { return WF::sine(getSamplePhasor(false)); }
  // Needs tests

  /** Resets the internal state, i.e. the sample counter and the random generator. */
  void reset(bool phasorRangeClosed);

  //-----------------------------------------------------------------------------------------------
  // \name Helpers

  /** Calculates the required probability distribution of the cycle lengths for the given desired
  mean cycle length given by the input parameter "meanCycleLength". In general, a cycle
  distribution for pitch dithering is determined by the 3 cycle lengths c1,c2,c3 to be produced
  along with their associated probabilities p1,p2,p3. However, in this setting here, it is always
  the case that c1 = c2 - 1, c3 = c2 + 1, p3 = 1 - (p1 + p2), so we have output parameters only for 
  c2 (= "lenMid"), p1 (= "probShort") and p2 (= "probMid"). Computing the rest, if really needed 
  (which it usually isn't), is up to the caller because this function is meant to be as efficient
  as possible because it's supposed to be called in a realtime context. The calculation implemented
  by this function is really the heart and soul of the pitch dithering idea that makes this 
  oscillator tick. That's why it has been made static and given output parameters rather than just
  operating directly on our member variables because we want to make the implementation re-usable
  by other oscillator code that also wants to implement pitch dithering. */
  static void calcCycleDistribution(T meanCycleLength, T* lenMid, T* probShort, T* probMid);


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internals

  /** Updates our sampleCount member and takes care of appropriate wrap around with recomputation 
  of the new cycle length if needed (which is the case when a cycle was just finished). */
  inline void updateSampleCount(bool phasorRangeClosed);

  /** Updates our cycleLength member by computing a new (pseudo) random cycle length to be used 
  for the next cycle. Called from updateSampleCount() after each cycle has been completed. */
  inline void updateCycleLength(bool phasorRangeClosed);

  //-----------------------------------------------------------------------------------------------
  // \name Data

  // Members that are accessed per sample:
  T phaseSlope;    // Increment per sample for phasor (which is in [0,1] or [0,1)).
  T sampleCount;   // Is integer and always in the interval [0, lenNow-1].
  T lenNow;        // Is lenMid or lenMid + 1 or lenMid - 1.

  // Members that are accessed per cycle:
  T lenMid;        // The middle one of the 3 integer cycle lengths to be produced.
  T probMid;       // Probability for lenMid.
  T probShort;     // Probability for lenMid - 1.
                   // Probability for lenMid + 1:  probLong = 1 - (probMid + probShort).

  // Embedded DSP objects:
  rsRandomGenerator<T> prng;
  uint32_t seed = 0;

};

//-------------------------------------------------------------------------------------------------
// Possibly inlined function implementations

template<class T> 
void rsPitchDitherOsc<T>::setMeanCycleLength(T newLength, bool closed)
{
  setMeanCycleLengthNoUpdate(newLength); 
  updateCycleLength(closed); 
}

template<class T> 
void rsPitchDitherOsc<T>::setMeanCycleLengthNoUpdate(T newLength)
{ 
  calcCycleDistribution(newLength, &lenMid, &probShort, &probMid); 
}

template<class T> 
inline T rsPitchDitherOsc<T>::getSamplePhasor(bool closed)
{
  T p = phaseSlope * sampleCount;      // Compute phasor output sample.
  updateSampleCount(closed);           // Increment with possible wraparound.
  return p;                            // Return phasor output sample.
}

template<class T> 
inline void rsPitchDitherOsc<T>::updateSampleCount(bool closed)
{
  sampleCount += T(1);                 // We produce 1 sample at each update.
  if(sampleCount >= lenNow)            // Is cycle finished?
  {                                    // If so..
    sampleCount = T(0);                // ..Wrap around sample counter.
    updateCycleLength(closed);         // ..Compute new lenNow and phaseSlope for next cycle.
  }
}

template<class T> 
inline void rsPitchDitherOsc<T>::updateCycleLength(bool closed)
{
  T r = prng.getSampleInUnitRange();   // Random number in interval [0,1).
  if(r < probShort)
    lenNow = lenMid - T(1);            // Next cycle is short.
  else if(r < probShort + probMid)
    lenNow = lenMid;                   // Next cycle is medium.
  else
    lenNow = lenMid + T(1);            // Next cycle is long.
  T maxCount = lenNow - T(closed);     // Maximum sample count until wrap around.
  phaseSlope = T(1) / maxCount;        // Phasor increment per sample.
}

template<class T> 
void rsPitchDitherOsc<T>::reset(bool closed)
{
  sampleCount = T(0);
  prng.setState(seed);
  updateCycleLength(closed);           // Important for correct initial lenNow.
}


#endif