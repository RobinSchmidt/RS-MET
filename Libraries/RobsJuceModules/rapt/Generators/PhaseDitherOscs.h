#ifndef RAPT_PHASEDITHEROSCS_H_INCLUDED
#define RAPT_PHASEDITHEROSCS_H_INCLUDED

//=================================================================================================

/** A class that factors out some functionality that is common to all generators that make use of 
pitch dithering. ...TBC... */

template<class T>
class rsPitchDitherHelpers
{

public:

  static void calcCycleDistribution(T period, T* midLength, T* probShort, T* probMid);
  // Maybe rename to cycleDistribEqualVariance (or ..EqVar) and maybe add functions to compute the
  // other distributions as well, although, I think, the other distributions are not really useful.
  // Not sure...maybe we shouldn't clutter the code with useless stuff. At least not the production
  // code. For research and prototype code, it's a different story. There, we may use the useless 
  // code to demonstrate in experiments that it is indeed useless.  


  // ToDo:
  //
  // - Add convenience functions to produce a whole signal vector of signals with various 
  //   waveforms. Maybe take the waveform as std::function or some callable template type F. 
  //
  // - Maybe create functions to produce various waveforms, including additively synthesized saw
  //   waves (maybe by using trig-recursions for an optimized implementation)
};

//=================================================================================================

/** A realtime oscillator that produces pitch-dithered sawtooth waves. */

template<class T> 
class rsPitchDitherSawOsc
{

public:

  using PDH = rsPitchDitherHelpers<T>;     // Shorthand for convenience

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Default constructor. It puts the object into a valid initial state by setting up a default
  period length of 100.0 samples and triggering the appropriate computations to set up our member
  variables that control the distribution of cycle lengths. */
  rsPitchDitherSawOsc()
  {
    setPeriod(T(100.0));                   // Triggers computations to set up members.
    reset();                               // Assigns sampleCount.
  }

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
  { PDH::calcCycleDistribution(newPeriod, &midLength, &probShort, &probMid); }

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
  T getPeriod()
  {
    T probLong = T(1) - (probShort + probMid);
    return probShort * (midLength - T(1)) + probMid * midLength + probLong * (midLength + T(1));
  }
  // Needs tests

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Produces one sample at a time. */
  inline T getSample()
  {
    T y = T(-1) + sawSlope * sampleCount;  // Compute output sample.
    sampleCount += T(1);                   // Update counter. We will have produced 1 sample.
    if(sampleCount >= cycleLength)         // Is cycle finished?
    {                                      // If so..
      sampleCount = T(0);                  // ..Wrap around sample counter.
      updateCycleLength();                 // ..Compute cycleLength and sawSlope for next cycle.
    }
    return y;                              // Return output sample.
  }

  /** Resets the internal state, i.e. the sample counter and the random generator. */
  void reset()
  {
    sampleCount = T(0);
    prng.setState(seed);
    updateCycleLength();                   // Important for correct initial cycleLength.
  }


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internals

  /** Updates our cycleLength member by computing a new (pseudo) random cycle length to be used 
  for the next cycle. Called from getSample() after each cycle has been completed. */
  inline void updateCycleLength()
  {
    T r = prng.getSampleInUnitRange();       // Random number in interval [0,1).
    if(r < probShort)
      cycleLength = midLength - T(1);        // Next cycle is short.
    else if(r < probShort + probMid)
      cycleLength = midLength;               // Next cycle is medium.
    else
      cycleLength = midLength + T(1);        // Next cycle is long.

    sawSlope = T(2) / (cycleLength - T(1));  // Slope depends on cycle length.
  }

  //-----------------------------------------------------------------------------------------------
  // \name Data

  // Members that are accessed per sample:
  T sawSlope;      // Increase of output value per sample.
  T sampleCount;   // Is always in interval [0, cycleLength).
  T cycleLength;   // Is midLength or midLength + 1 or midLength - 1.

  // Members that are accessed per cycle:
  T midLength;     // The middle one of the 3 cycle lengths to be produced.
  T probShort;     // Probability to use midLength - 1.
  T probMid;       // Probability to use midLength.

  // Embedded DSP objects:
  rsRandomGenerator<T> prng;
  uint32_t seed = 0;

  // Notes:
  //
  // - The members sampleCount, cycleLength and midLength are actually integer numbers but we let 
  //   them be of type T (which is typically float or double) anyway for optimization purposes. We 
  //   want to avoid int-to-float conversions because they are costly. The other members are indeed
  //   true floating point values that are not restricted to integers.
  // 
  // - There is no "probLong" member variable because that would be redundant. The probability to
  //   use the long cycle is always given by  1 - (probShort + probMid)  because probabilities must
  //   add up to 1.
  // 
  // - The member variables are deliberately not initialized in their declarations because some 
  //   initializations would require using some moderately complex formulas for a consistent, valid
  //   initial state. So we leave this member initialization to the constructor which calls some 
  //   functions to do the appropriate computations.
  //
  // - Maybe factor out all the stuff that has to do with the pitch-dithering into a separate class
  //   such that we can re-use the code for other types of pitch-dithering oscillators like, for 
  //   example, table lookup oscillators. Or maybe modify this class such that it can also produce
  //   a sawtooth in the range [0,1) that other oscillators can use as phasor. Maybe have a 
  //   function getPhase() and getSample() would just return 2*phase - 1. Check, if we currently 
  //   produce a saw in [-1,+1] or in [-1,+1) and document that. Maybe allow for different modes.
  //   Maybe the mode can be a compile-time parameter, i.e. a template parameter. Or maybe factor 
  //   out the stuff that is common to all variants into a baseclass and realize the different 
  //   variants as subclasses. The different subclasses need different implementations of 
  //   getSample() and updateCycleLength() but most of the code in these functions will be the same
  //   so maybe that should be factored out into functions. In getSample(), only the first line
  //   T y = T(-1) + sawSlope * sampleCount;  will be different. In the the 0..1 case, the T(-1)
  //   will be missing because we start at 0. In updateCycleLength(), only the last line 
  //   sawSlope = T(2) / (cycleLength - T(1));  will be different. In the 0..1 case, we need to 
  //   adapt the formula to T(1) / ... instead of T(2) / .... I think, if we want to produce closed
  //   intervals rather than half-open ones, we need to get rid of the " - T(1)" in the 
  //   denominator. ..but verify this! In general, it could make sense to have 4 variants with the 
  //   ranges [-1,+1], [-1,+1), [0,1], [0,1). Maybe the [0,1) version is the most important one. 
  //   This is the typical range for a phasor. This can be seen from what would happen if we would 
  //   use a phasor with range [0,1] with a sine wave produced as y = sin(2*PI*phasor). With the 
  //   closed interval, the 0 value would be repeated: Once it would occur at phasor = 0 and 
  //   secondly at phasor = 2*PI. That's clearly wrong, so [0,1) is the correct range for a phasor.
  //   To create a supersaw, it could actually be more efficient to just sum up the phasors and 
  //   then subtract sumOfAmplitudes once from the whole supersaw instead of subtracting 1 from
  //   each saw. So maybe it would be best to provide the two functions getPhasorSample() and 
  //   getSawSample() and the latter is just implemented as: "return 2 * getPhasorSample() - 1".
  //   Or maybe it should be called getSample(). But we could also rename this class to 
  //   rsPitchDitherOsc without limiting it to saw waves. in this case, getSawSample() would make 
  //   more sense. And then we could also have getSinSample(), getRectSample(), getTriSample(),
  //   getPulseSample(T pw), getTriSawSample(..), etc. If we do this, we may also rename sawSlope
  //   to phaseSlope. 
  // 
  // - Maybe that class can then also take the responsibility of rsPitchDitherHelpers which 
  //   currently just has this single static member function and I don't really think that it will
  //   need anything else (not sure, though). 
};


#endif