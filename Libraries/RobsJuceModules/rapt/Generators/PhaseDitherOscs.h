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
  // Not sure...maybe we shouldn't clutter the code with useless stuff. At leats not the production
  // code. For research and prototype code, it's a different story. There, we may use the useless 
  // code to demonstrate in experiments that it is indeed useless.  


  // ToDo:
  //
  // - Add convenience functions to produce a whole signal vector of signals with various 
  //   waveforms. Maybe take the waveform as std::function or some callable template type F. 
  //
  // - Maybe create functions to produce various wvaeforms, including additively syntehsized saw
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
    //prng.setRange(T(0), T(1));             // We use random numbers in the interval [0,1).
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
  void setPeriod(T newPeriod)
  {
    setPeriodNoUpdate(newPeriod);
    updateCycleLength();
  }

  /** Sets up a new period length just like setPeriod() does but without immediately updating the
  probability distribution and current cycle length. This results in the behavior that the new 
  period will not become effective immediately but only after finishing the currently running 
  cycle. */
  void setPeriodNoUpdate(T newPeriod)
  { 
    PDH::calcCycleDistribution(newPeriod, &midLength, &probShort, &probMid); 
  }

  /** Sets the seed for the pseudo random number generator. */
  void setRandomSeed(uint32_t newSeed)
  {
    seed = newSeed;
    //prng.setSeed(newSeed);
  }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the average length of the cycles that are being produced. If the 3 integer cycle 
  lengths are given by L1,L2,L3 and cycles with these 3 lengths are produced with probabilities 
  p1,p2,p3 respectively, then the average cycle length P will be: P = p1*L1 + p2*L2 + p3*L3. */
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
    //prng.reset();        // Old
    prng.setState(seed);   // New
    updateCycleLength();   // Also new. Not sure about that, though-
  }


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internals

  /** Updates our cycleLength member by computing a new (pseudo) random cycle length to be used 
  for the next cycle. Called from getSample() after each cycle has been completed. */
  inline void updateCycleLength()
  {
    //T r = prng.getSample();                  // Random number in interval [0,1).
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
  //rsNoiseGenerator<T> prng;

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
  // 
  // ToDo:
  //
  // - Replace the rsNoiseGenerator<T> member by a rsRandomGenerator<T>. This class does not exist 
  //   yet. It is supposed to factor out the integer random number generation from rsNoiseGenerator
  //   such that the object doesn't need to maintain the shift and scale members. It could provide
  //   a convience function for producing floating point outputs in the fixed range [0,1) but it 
  //   will not provide a user adjustable range. It could have a convenience function 
  //   getSampleFloat(T min, T max), though. But using that function in a hot loop is not 
  //   recommended for production code because it will need a costly division. Maybe it should have
  //   a normal getSample() function that produces floats in [0,1) and an additional getSampleRaw()
  //   or getSampleInt() function that produces the raw integer value.
  // 
  // - setRandomSeed() should probably eventually be replaced by a setRandomState function. Or 
  //   maybe complemented by it. We actually do want to have a setRandomSeed() function as well. 
  //   But when we switch to a PRNG that has no built in seed member, we will need a randomSeed 
  //   member variable here.
};


#endif