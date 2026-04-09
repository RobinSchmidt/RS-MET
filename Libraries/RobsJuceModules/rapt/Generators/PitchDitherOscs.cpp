
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
  T ci = rsFloor(c);                   // Integer part of desired cycle length c.
  T cf = c - ci;                       // Fractional part of it.
  T c2 = ci;                           // Length of the middle length cycle.
  if(cf >= T(0.5))                     // If fractional part of c is >= 0.5...
    c2 += T(1);                        // ..the mid length must be one sample longer.
  T c1 = c2 - T(1);                    // Short cycles are one sample shorter than mid.
  T c3 = c2 + T(1);                    // Long cycles are one sample longer than mid.

  // Compute intermediates:
  T e1 = c1 - c;                       // Length error of short cycle.
  T e2 = c2 - c;                       // Length error of mid cycle.
  T e3 = c3 - c;                       // Length error of long cycle.
  T v1 = e1 * e1;                      // Variance contribution from short cycles.
  T v2 = e2 * e2;                      // Variance contribution from mid cycles.
  T v3 = e3 * e3;                      // Variance contribution from long cycles.
  T v  = T(0.25);                      // Target variance determined by the cf = 0.5 "worst case".
  T d1 = v - v1;                       // Deviation from target variance of short cycles.
  T d2 = v - v2;                       // Deviation from target variance of mid cycles.
  T d3 = v - v3;                       // Deviation from target variance of long cycles.
  T s  = T(1) / (e3*(v1-v2) - e2*(v1-v3) + e1*(v2-v3));  // Common scaler for probabilities.

  // Compute and assign outputs:
  *lenMid    = c2;
  *probShort = (d2*e3 - d3*e2) * s;    // Probability p1 to use short cycle with length c1.
  *probMid   = (d3*e1 - d1*e3) * s;    // Probability p2 to use middle cycle with length c2.
  //*probLong  = (d1*e2 - d2*e1) * s;  // Would be redundant: p3 = 1 - (p1 + p2). See below.
  
  // We don't have a probLong output parameter because that would be redundant. It would always be
  // given by probLong = 1 - (probShort + probMid). We will also always have: 
  // lenShort = lenMid + 1, lenLong = lenMid + 1. The derivation of these formulas can be found in
  // the textfile PitchDithering.txt in the research repo. ToDo: clean the derivation up and put it
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

- Performance analysis: The per sample code when we call getSampleSaw() consists of the code in
  rsPitchDitherOsc::getSamplePhasor() and rsWaveForms::saw(). The former has 1 mul, 1 add, 1 if
  and the latter has 1 mul, 1 sub. So, in total, we have 2 mul, 1 add, 1 sub, 1 if per sample. To
  produce optimized code for a saw, it may be better to not first produce a phasor in 0..1 and 
  then converting it to -1..+1 but instead directly producing it in the range -1..+1. That would
  require to use  phaseSlope = T(2) / maxCount;  instead of  phaseSlope = T(1) / maxCount;  in 
  updateCycleLength(). But we really want the flexibility of having a phasor to be able to produce
  arbitrary waveforms so we accept this slight suboptimality. The per cycle calculations are:
  1 assign, 1 PRNG-evaluation, 3 add, 1 3-way branch (with 2 compares), 1 div. The nice thing is
  that when we use oversampling, the per-cycle calculations will be rarer when time is measured in
  terms of samples. The calling frequency will be constant when time is measured in absolute time.
  So only the per-sample calculations (2 mul, 1 add, 1 sub, 1 if) will have to be done at the 
  oversampled rate so pitch dithering should play really nice with oversampling. The calculations 
  in calcCycleDistribution() are only called when the user sets up a new frequency via 
  setMeanCycleLength() which could potentially be a thing that is getting called per sample when
  we have a pitch envelope or LFO or even frequency modulation. The somewhat high cost of pitch
  modulation could be mitigated by using setMeanCycleLengthNoUpdate() instead which would postpone
  the recalculation until the end of the currently running cycle. This seems to be a reasonable 
  thing to do.

- To create a supersaw, it could actually be more efficient to just sum up the phasors and 
  then subtract sumOfAmplitudes once from the whole supersaw instead of subtracting 1 from
  each saw. We would bypass the per-saw conversion from phasor to the actual saw and would instead
  just add up the phasors and then do the conversion for all of them at once. This is equivalent
  because everything is linear (Verify! Maybe it's not exactly linear but rather linear-with-shift
  aka affine - but that's something we can work with, too).

- Sound: At a sample rate of 44.1 kHz the low octaves sound pretty clean and the sound gets 
  progressively more noisy towards the higher octaves which is the behavior we expect. We can sense
  a clear pitch up to a fundamental around XXXX ...above that, it sounds more like high frequency
  noise without much tonality. With a highpass, it can make nice "mosquito" sounds. It's 
  interesting to feed it into a harsh waveshaper (hard-clip, fold, quantize, ...), cranking up the
  drive and and then applying an amplitude envelope. The great thing about this is that this 
  doesn't produce any new aliasing frequencies. They will again line up with the already existing 
  spectrum. We still get a perfectly pitch dithered waveform


ToDo:

- Implement more waveforms: square, pulse, triangle, sine, trisaw, etc. Write into the 
  documentation that these standard waveforms can be used as examples for client code to 
  implement their own custom waveforms.

- Add a setPhase(T newPhase) function. It should set the sampleCounter to a phase 
  corresponding to the (rounded) newPhase value such that in the very next call to 
  getSamplePhasor(), we will get exactly that (rounded) newPhase value. We need to round 
  because our sampleCounter is an integer.

- Add convenience functions like setOmega(T newOmega), setFrequency(T newFreq, T sampleRate).
  setFrequency should perhaps just call setPeriod(sampleRate/newFreq)

- Add getSampleTriangle(), getSampleTriSaw(p, shape)

- Maybe abbreviate calcCycleDistribution() as calcCycleDistrib().

- Drag over the experiments and unit tests from the research repo into the main repo. Implement 
  some unit tests here.

- Implement a class rsPitchDitherSuperSawOsc. See comments in the experiments in the research repo
  for how to approach this. ..ok: we now have a class for that somewhere in the prototypes or the
  research repo.

- Maybe provide default values for the "phasorRangeClosed" parameter to make them optional. I 
  think, it should probably default to false because true is mostly to make saws look nicer but 
  they kinda also work with false whereas for a sine wave, using true will create audible 
  artifacts. I think, a half-open range is more common for a phasor. Maybe have two functions:
  getSamplePhasor() and getSamplerPhasorClosed() where the former implementes the (standard, 
  default) half-open interval. Maybe start a KVR thread like: "Intervals for (normalized)
  osc-phases: Closed [0,1] or half-open [0,1)?" It seems to me that the correct choice depends on
  the waveform to produce. Sines like the half-open version, saws like the closed version. But
  that's awkward!

- To free the client code from this awkwardness, maybe the class rsPitchDitherOsc should be renamed
  to rsPitchDitherOscBase and we should provide a class rsPitchDitherOsc with an API containing a 
  function setWaveForm with a waveform from some enum and all the getSamplePhasor(), 
  getSampleSawUp(), etc. stuff should go away and we should just have the normal API with 
  getSample(), reset(), etc. Having to handle this half-open/closed phasor business is low level 
  stuff that client code should be able to tap in if needed but by default, that shouldn't be the 
  case. The enum should probably exist inside some class like rsWaveForms such that it can be 
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
  pitch-dithered supersaw osc. We may also want to add pitch dithering to table lookup oscillators.


Ideas:

- Combine pitch dithering with waveshaping. I think, when we waveshape a pitch-dithered waveform,
  we will not introduce aliasing. Instead, we will probably modify the noise spectrum. I think that
  because all that waveshaping does is to modify the instantaneous signal value. The output will 
  again be a pitch dithered signal with a new waveform. Verify that experimentally!

- Maybe try combining oversampled pitch-dithering with (non-oversampled) waveshaping. My guess is 
  that we may get aliasing but restrict it to subharmonics. For an oversampling factor of 2, I 
  would expect to see a suboctave, for a factor of 3 a subharmonic at 1/3 of the fundamental. Try
  that!

- Combine pitch dithering with hard-sync. I think, when the master osc is pitch dithered, we may 
  also get some results that don't show obvious aliasing. Try it!

- Currently, low notes sound cleaner than high notes. Maybe we could create a more complicated
  version of the idea in which we artificially dirtify the low notes by using a broader cycle 
  distribution that spans more that 3 integer cycle lengths. Maybe we should pick some high note as
  reference at which we use 3 lengths and then double the width of the distribution for each octave
  that we go down. If we use 3 lengths at 1 kHz, we could use 5 at 500 Hz, 9 at 250 Hz, 17 at 125
  Hz etc. The idea is that we use numbers of the from 2^k + 1 where the base case of 3 corresponds 
  to k = 1. The idea is that the effective width is actually 2 samples rather than 3 when we also 
  consider the probabilities as weights. At c = xxx.5, we do indeed only use 2 integer lengths 
  (both with probability 0.5). We could also say that we use 2 neighbors around the middle length.



*/