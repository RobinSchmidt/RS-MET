
template<class T>
rsPitchDitherOsc<T>::rsPitchDitherOsc()
{
  setMeanCycleLength(T(100.0), true);    // Triggers computations to set up members.
  reset(true);                           // Assigns sampleCount.
}

template<class T>  
T rsPitchDitherOsc<T>::getMeanCycleLength()
{
  T lenShort = lenMid - T(1);
  T lenLong  = lenMid + T(1);
  T probLong = T(1) - (probShort + probMid);
  return probShort*lenShort + probMid*lenMid + probLong*lenLong;

  // In general, if we have 3 integer cycle lengths given by c1,c2,c3 and cycles with these 3 
  // lengths are produced with probabilities p1,p2,p3 respectively, then the mean cycle length cM
  // will be: cM = p1*c1 + p2*c2 + p3*c3. In our particular case here, we will have a given middle 
  // length c2 and the short and long lengths c1,c3 are then given as c1 = c2-1, c3 = c2+1. And, of
  // course, for the probabilities we will always have p1+p2+p3 = 1 such that p3 = 1 - (p1 + p2).
}

template<class T>
void rsPitchDitherOsc<T>::calcCycleDistribution(T c, T* lenMid, T* probShort, T* probMid)
{
  // Compute cycle lengths c1,c2,c3:
  T ci = rsFloor(c);                   // Integer part of desired cycle length c
  T cf = c - ci;                       // Fractional part of it
  T c2 = ci;                           // Length of the middle length cycle
  if(cf >= T(0.5))                     // If fractional part of c is >= 0.5...
    c2 += T(1);                        // ..the mid length must be one sample longer
  T c1 = c2 - T(1);                    // Short cycles are one sample shorter than mid
  T c3 = c2 + T(1);                    // Long cycles are one sample longer than mid

  // Compute intermediates:
  T e1 = c1 - c;                       // Length error of short cycle
  T e2 = c2 - c;                       // Length error of mid cycle
  T e3 = c3 - c;                       // Length error of long cycle
  T v1 = e1 * e1;                      // Variance contribution from short cycles
  T v2 = e2 * e2;                      // Variance contribution from mid cycles
  T v3 = e3 * e3;                      // Variance contribution from long cycles
  T v  = T(0.25);                      // Target variance determined by the cf = 0.5 "worst case"
  T d1 = v - v1;                       // Deviation from target variance of short cycles
  T d2 = v - v2;                       // Deviation from target variance of mid cycles
  T d3 = v - v3;                       // Deviation from target variance of long cycles
  T s  = T(1) / (e3*(v1-v2) - e2*(v1-v3) + e1*(v2-v3));  // Common scaler for probabilities

  // Compute and assign outputs:
  *lenMid    = c2;
  *probShort = (d2*e3 - d3*e2) * s;    // Probability p1 to use short cycle with length c1
  *probMid   = (d3*e1 - d1*e3) * s;    // Probability p2 to use middle cycle with length c2
  //*probLong  = (d1*e2 - d2*e1) * s;  // That would be redundant. See below.
  
  // We don't have a probLong output parameter because that would be redundant. It would always be
  // given by 1 - (probShort + probMid). The derivation of these formulas can be found in the
  // textfile PitchDithering.txt in the research repo. ToDo: clean the derivation up and put it
  // into its own dedicated textfile here in the main repo! We actually already have now an .md
  // file but it's not yet finished. When it's done, reference it here.
}

//=================================================================================================
/*

Notes:

- The members sampleCount, lenNow and lenMid are actually integer numbers but we let 
  them be of type T (which is typically float or double) anyway for optimization purposes. We 
  want to avoid int-to-float conversions because they are costly. The other members are indeed
  true floating point values that are not restricted to integers.

- There is no "probLong" member variable because that would be redundant. The probability to
  use the long cycle is always given by  1 - (probShort + probMid)  because probabilities must
  add up to 1.

- The member variables are deliberately not initialized in their declarations because some 
  initializations would require using some moderately complex formulas for a consistent, valid
  initial state. So we leave this member initialization to the constructor which calls some 
  functions to do the appropriate computations.

- Partially done:
  Maybe factor out all the stuff that has to do with the pitch-dithering into a separate class
  such that we can re-use the code for other types of pitch-dithering oscillators like, for 
  example, table lookup oscillators. Or maybe modify this class such that it can also produce
  a sawtooth in the range [0,1) that other oscillators can use as phasor. Maybe have a 
  function getPhase() and getSample() would just return 2*phase - 1. Check, if we currently 
  produce a saw in [-1,+1] or in [-1,+1) and document that. Maybe allow for different modes.
  Maybe the mode can be a compile-time parameter, i.e. a template parameter. Or maybe factor 
  out the stuff that is common to all variants into a baseclass and realize the different 
  variants as subclasses. The different subclasses need different implementations of 
  getSample() and updateCycleLength() but most of the code in these functions will be the same
  so maybe that should be factored out into functions. In getSample(), only the first line
  T y = T(-1) + sawSlope * sampleCount;  will be different. In the the 0..1 case, the T(-1)
  will be missing because we start at 0. In updateCycleLength(), only the last line 
  sawSlope = T(2) / (cycleLength - T(1));  will be different. In the 0..1 case, we need to 
  adapt the formula to T(1) / ... instead of T(2) / .... I think, if we want to produce closed
  intervals rather than half-open ones, we need to get rid of the " - T(1)" in the 
  denominator. ..but verify this! In general, it could make sense to have 4 variants with the 
  ranges [-1,+1], [-1,+1), [0,1], [0,1). Maybe the [0,1) version is the most important one. 
  This is the typical range for a phasor. This can be seen from what would happen if we would 
  use a phasor with range [0,1] with a sine wave produced as y = sin(2*PI*phasor). With the 
  closed interval, the 0 value would be repeated: Once it would occur at phasor = 0 and 
  secondly at phasor = 2*PI. That's clearly wrong, so [0,1) is the correct range for a phasor.
  To create a supersaw, it could actually be more efficient to just sum up the phasors and 
  then subtract sumOfAmplitudes once from the whole supersaw instead of subtracting 1 from
  each saw. So maybe it would be best to provide the two functions getPhasorSample() and 
  getSawSample() and the latter is just implemented as: "return 2 * getPhasorSample() - 1".
  Or maybe it should be called getSample(). But we could also rename this class to 
  rsPitchDitherOsc without limiting it to saw waves. in this case, getSawSample() would make 
  more sense. And then we could also have getSinSample(), getRectSample(), getTriSample(),
  getPulseSample(T pw), getTriSawSample(..), etc. If we do this, we may also rename sawSlope
  to phaseSlope. ...BUT: I actually do thing that we produce the closed interval [0,1] here 
  and in the case of sawtooth waves, it is actually sort of appropriate because in the case
  of a jump discontinuity at the wrap around, it can make sense to return at the sample 
  instant zero one value and at the sample instant at the end of cycle the other value. 
  Hmmm...not sure what to do. Both variants have convincing arguments. Maybe we should 
  implement both and let the user choose? Maybe at compile time? The devil is in the detail!
  Maybe updateCycleLength could take a bool parameter closedInterval or something like that
  and we could give the use two versions of getSamplePhasor() like getSamplePhasorClosed(),
  getSamplePhasorHalfOpen(). Or maybe the "HalfOpen" version should go without qualification
  to indicate that this is the default. But will this lead to detuning? Is the current 
  implementation actually correctly tuned anyway? Maybe currently the cycles are one sample
  too short or too long? Verify this! ...done! Nope - it's alright. The period length is correct.

- Implement more waveforms: square, pulse, triangle, sine, trisaw, etc. Write into the 
  documentation that these standard waveforms can be used as examples for client code to 
  implement their own custom waveforms.

- Add a setPhase(T newPhase) function. It should set the sampleCounter to a phase 
  corresponding to the (rounded) newPhase value such that in the very next call to 
  getSamplePhasor(), we will get exactly that (rounded) newPhase value. We need to round 
  because our sampleCounter is an integer.

- Add convenience functions like setOmega(T newOmega), setFrequency(T newFreq, T sampleRate).
  setFrequency should perhaps just call setPeriod(sampleRate/newFreq)


ToDo:

- Add getSampleTriangle(), getSampleTriSaw(p, shape)

- Drag over the experiments and unit tests from the research repo into the main repo. Implement 
  some unit tests here.

- Implement a class rsPitchDitherSuperSawOsc. See comments in the experiments in the research repo
  for how to approach this. ..ok: we now have a class for that somewhere in the prototypes or the
  research repo.

- Maybe provide default values for the "phasorRangeClosed" parameter to make them optional. I 
  think, it should probably default to false because true is mostly to make saws look nicer but 
  they kinda also work with false whereas for a sine wave, using true will create audible 
  artifacts. I think, a half-open range is more common for a phasor. Maybe start a KVR thread like:
  Intervals for (normalized) osc-phases: Closed [0,1] or half-open [0,1)? It seems to me that the
  correct choice depends on the waveform to produce. Sines like the half-open version, saws like 
  the closed version. But that's awkward!

- To free the client code from this awkwardness, maybe the class rsPitchDitherOsc should be renamed
  to rsPitchDitherOscBase and we should provide a class rsPitchDitherOsc with an API containing a 
  function setWaveForm with a waveform from some enum and all the getSamplePhasor(), 
  getSampleSawUp(), etc. stuff should go away and we should just have the normal API with 
  getSample(), reset(), etc. Having to handle this half-open/closed phasor business is low level 
  stuff that client code should be able to tap in if needed but by default, that shouldn't be the 
  case. The enum should probably exist inside some class liek rsWaveForms such that it can be 
  re-used by other classes. Check the preliminary implementation of rsWaveForms in 
  rs_testing/Prototypes/Generators.h. This could be used as basis.

- Maybe as an optimization, pass the "phasorRangeClosed" parameter not as bool but as type T so we
  can avoid the type conversion in updateCycleLength(). Maybe we should assert that the value 
  represents either T(0) or T(1). Maybe add a function rsIsBoolean(T x) to the library that returns
  true iff x is 0 or 1 and use that function in a rsAssert here. Maybe have static const members 
  phasorRangeClosed, phasorRangeHalfOpen of type T that are fixed to 0 and 1 such that the caller 
  can use these as symbolic constants rather than having itself to make sure to only pass 0 or 1.

- Drag over the code for the pitch-dithered supersaw oscillator.

- When other oscillator classes that use that functions are added to the library, mention them 
  in the documentation of calcCycleDistribution(). For example, later we want to add the 
  pitch-dithered supersaw osc. We may also want to add pitchdithering to table lookup oscillators.

*/