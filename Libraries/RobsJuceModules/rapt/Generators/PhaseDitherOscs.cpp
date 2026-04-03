
template<class T>
rsPitchDitherOsc<T>::rsPitchDitherOsc()
{
  setPeriod(T(100.0));                   // Triggers computations to set up members.
  reset();                               // Assigns sampleCount.
}

template<class T>  
T rsPitchDitherOsc<T>::getPeriod()
{
  T probLong = T(1) - (probShort + probMid);
  return probShort * (midLength - T(1)) + probMid * midLength + probLong * (midLength + T(1));

  // In general, if we have 3 integer cycle lengths given by L1,L2,L3 and cycles with these 3 
  // lengths are produced with probabilities p1,p2,p3 respectively, then the average cycle length P
  // will be: P = p1*L1 + p2*L2 + p3*L3. In our particular case here, we will have a given middle 
  // length L2 and the short and long lengths L1,L3 are then given as L1 = L2-1, L3 = L2+1. And, of
  // course, for the probabilities we will always have p1 + p2 + p3 = 1 such that 
  // p3 = 1 - (p1 + p2).
}

template<class T>
void rsPitchDitherOsc<T>::calcCycleDistribution(
  T period, T* midLength, T* probShort, T* probMid)
{
  // Compute lengths:
  T floorLength = rsFloor(period);
  T fracLength  = period - floorLength;
  T L1, L2, L3;
  if(fracLength < T(0.5))
    L1 = floorLength - T(1);
  else
    L1 = floorLength;
  L2 = L1 + T(1);
  L3 = L2 + T(1);

  // Compute intermediates:
  T e1 = L1 - period;
  T e2 = L2 - period;
  T e3 = L3 - period;
  T m1 = e1*e1;
  T m2 = e2*e2;
  T m3 = e3*e3;
  T M  = T(0.25);
  T M1 = M - m1;
  T M2 = M - m2; 
  T M3 = M - m3;
  T S  = T(1) / (e3*(m1-m2) - e2*(m1-m3) + e1*(m2-m3));

  // Compute outputs:
  *midLength = L2;
  *probShort = (M2*e3 - M3*e2) * S;
  *probMid   = (M3*e1 - M1*e3) * S;
  //*probLong  = (M1*e2 - M2*e1) * S;  // Would be redundant. See Notes
  
  // Notes:
  // 
  // - We don't have a probLong parameter because that would be redundant. It would always be given
  //   by 1 - (probShort + probMid).
  //
  // - The derivation of these formulas can be found in the textfile PitchDithering.txt in the 
  //   research repo. ToDo: clean the derivation up and put it into its own dedicated textfile here
  //   in the main repo!
}

//=================================================================================================
/*

Notes:

- The members sampleCount, cycleLength and midLength are actually integer numbers but we let 
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

- Maybe factor out all the stuff that has to do with the pitch-dithering into a separate class
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
  too short or too long? Verify this!

- Maybe that class can then also take the responsibility of rsPitchDitherHelpers which 
  currently just has this single static member function and I don't really think that it will
  need anything else (not sure, though). 

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

- Drag over the experiments and unit tests from the research repo into the main repo

- Add unit test for rsNoiseGenerator, then factor out a class rsRandomGenerator that doesn't have
  members for seed, shift and scale, then use that class here instead of rsNoiseGenerator.

- Implement a class rsPitchDitherSuperSawOsc. See comments in the experiments in the reseatch repo
  for how to approach this.

- Clean up the derivations for the cycle distribution formulas and put them into a dedicated 
  textfile. They are currently in TempSketchPad.txt in the research repo and are rather messy.

- Add convenience functions to produce a whole signal vector of signals with various 
  waveforms. Maybe take the waveform as std::function or some callable template type F. 

- Maybe create functions to produce various waveforms, including additively synthesized saw
  waves (maybe by using trig-recursions for an optimized implementation)

*/