#ifndef RAPT_FAKERESONANCEFILTER_H
#define RAPT_FAKERESONANCEFILTER_H

/**

This is a filter that emulates a resonant lowpass filter by using a non-resonant lowpass in
parallel with a separate resonance path. The resonance path consists of a highpass filter in
series connection with a filter that has an enveloped sinusoid as its impulse response. This
envelope can be cotrolled by separate attack and decay parameters. The idea is that the highpass
converts edges in the input signal into impulse/spike like signals which in turn are converted
into enveloped sinusoids and added to the non-resonant lowpass output, thereby creating a
resonance-like effect. The advantage over more conventional techniques to produce filter
resonance, for example via feedback around a 4-pole ladder, is that the shape of the resonance
can be more directly controlled.

The structure is:

 in --------> LPF ----------->+----> out
       |                      ^
       |                      |
       \----> HPF ----> RF ---/

where: LPF: lowpass filter, HPF: highpass filter, RF: resonance filter

BUGS:
-when the delay is set to 0, the actual delay seems to be the length of the delayline - or
something

Ideas for improving the transient response::
use a resonant and a non-resonant lowpass filter and subtract their outputs to obtain the pure
resonance signal which can then be shaped further - or maybe this can just be used to analyze the
pure resonance signal - maybe then another filter can be put in parallell
that adds exactly what is missing in the resonance we create here....
so we do:
resonance = resonant - nonresonant (using regular filters), then we model
resonance = attackdecaysine + trnasientcorrection
where that "transientcorrection" would be an impulse response of another filter which would have
to figured out

things to try:
-attenuate/notch the resonance frequency from the lowpass path
-use an allpass for the resonance to shape  the transient

\todo maybe use a simplified highpass filter which is just a simple differencer ...but think
about the implications for what this means for using the filter at different samplerates.

\todo Replace the 4 individual 1st order lowpasses with a LadderFilter class - this will allow
for different filter modes in the upper path (not just lowpass but also highpass, bandpass, etc.)
as well as a secondary resonance. It will also clean up the code and make it more efficient.
Maybe provide a feedback parameter for that ladder filter, too - for easy comparing of the
feedback-resonance and the fake-resonance and to possibly further shape the resonance.

\todo maybe allow for negative delay times for the resonance path by using a delayline in the
lowpass path, too. whenever the user sets a negative delay, the lowpass path will be delayed
instead of the resonance path

\todo: for the waveshaping, maybe provide a control that blends between static and dynamic
waveshaping. Dynamic behavior is the normal case, where the output spectrum depends on input
amplitude. But we can have access to the amplitude envelope of the resonator (when we switch to
the spiraling-phasor implementation), so we may also provide a static waveshaping by dividing the
signal by its envelope before the waveshaper and multiplying the envelope back into it after the
waveshaper. The result should be an output spectrum that is (more or less) independent from input
amplitude. We could somehow "crossfade" continuously between both behaviors.

for the waveshaping of the resonance, let the user pass either a function-pointer or
functor, as explained here:
http://www.cprogramming.com/tutorial/functors-function-objects-in-c++.html

maybe introduce a dependence of the resonance amplitude from the decay time in such a way that
when it is set to unity, the amplitude is scaled with a factor that ensures the same overall
energy in the resonance signal, regardless of the decaytime (i think, the should be inversely
related to the decaytime then) ...or maybe let the amplitude directly depend on the frequency in
such a way that when freqToDecay and freqToAmp are both set to unity, the energy in the resonance
envelope is independent from frequency (decay scales inversely with frequency, amplitude
proportionally)

we could create a harmonic series like 2,4,8,16,... by using the 2nd chebychev polynomial
on the output of each respective previous polynomial ...but maybe successive squarers are better


\todo: optimizations:
-avoid the interpolation-method dispatch in the delayline
-get rid of the tempo-sync option in the rsFractionalDelayline class
-when going to stereo, use a SIMD-vector of 2 doubles as signal type - write a class rsFloat64x2
 and implement all the basic math operators between 2 values of that type and one value of the
 class-type and the other of type double ("broadcasting" operators)
-in the computation of the attack time constant, use a better initial guess for newton
 iteration (see comment there) */

template<class TSig, class TPar>
class rsFakeResonanceFilter
{

public:


  /** \name Construction/Destruction: */

  /** Constructor. */
  rsFakeResonanceFilter();


  /** \name Setup */

  /** Sets the sample rate in Hz. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the cutoff frequency of the lowpass in Hz. */
  void setLowpassCutoff(TPar newCutoff);

  /** Sets the cutoff frequency of the highpass in Hz. */
  void setHighpassCutoff(TPar newCutoff);

  /** Sets the pitch-shift (in semitones) of the resonance frequency with respect to the lowpass
  filter's cutoff frequency. If 0, the resonance will occur at the cutoff frequency, but with
  this parameter, it's possible to offset it from the cutoff. */
  void setResonanceShift(TPar newShift);

  /** Sets the resonance amplitude as a linear gain factor. */
  void setResonanceAmplitude(TPar newAmplitude);

  /** Sets the start phase of the resonance in radians. */
  void setResonancePhase(TPar newPhase);

  /** Sets the decay time of the resonance in seconds. */
  void setResonanceDecay(TPar newDecay);

  /** Sets the amount of dependency of the resonance decay time from the resonance frequency. If
  it's set to zero, the decaytime will be independent from the resonant frequency, corresponding
  to a resonator filter with constant absolute bandwidth. If it's set to unity, the decaytime
  will be inversely proportional to the resonant frequency, corresponding to a resonator with
  constant relative bandwith (which implies constant Q as well). You may set it to values less
  than zero and larger than unity as well. In case of a nonzero dependency, we use a reference
  frequency of 1kHz - that means, at 1kHz (for resonance frequency), the decay time will be left
  unmodified. */
  void setDecayByFrequency(TPar newDecayByFreq);

  /** Sets the attack time of the resonance as a scale factor k to be applied to the decay time,
  i.e. the absolute attack time is taken to be this factor times the decay time:
  attack = k * decay. For technical resons, the value must satisfy: 0 < k < 1. */
  void setResonanceAttack(TPar newAttack);

  /** Sets the delaytime for the resonance path with respect to the filter path. The unit is a
  multiplier for an optimal delay time that is given by the reciprocal of the lowpass cutoff
  frequency. Choosing the delay in this way gives the best alignment of the resonance start with
  lowpassed signal edges. */
  void setResonanceDelay(TPar newDelay);


  /** \name Audio Processing */

  /** Calculates a single filtered output-sample. */
  RS_INLINE TSig getSample(TSig in);

  /** Computes the lowpass signal (upper path) */
  RS_INLINE TSig getLowpassOutput(TSig in);

  /** Computes the resonance signal (lower path) */
  RS_INLINE TSig getResonanceOutput(TSig in);



  /** \name Misc */

  /** Resets the internal buffers. */
  void reset();


protected:


  /** Sets up the resonance filter from the user parameters. */
  void setupResonance();

  /** Sets up the delay in the resonance path from the filter's cutoff frequency and the delay
  parameter. */
  void setupResonanceDelay();


  // processors in filter path:
  rsOnePoleFilter<TSig, TPar> lpf1, lpf2, lpf3, lpf4;
    // Four 1st order lowpasses in series for the lowpass path
    // \todo: Maybe replace this series of 4 individual 1st order filters by an actual 
    // implementation of the Moog-ladder structure.


  // processors in resonance path:
  rsFractionalDelayLine<TSig, TPar>   dl;     // delayline
  rsOnePoleFilter<TSig, TPar>         hpf;    // highpass
  rsModalFilterWithAttack<TSig, TPar> rf;     // resonator


  /** \name Internal Functions */


  /** \name Data */

  TPar fs;      // sample rate in Hz
  TPar cutoff;

  // maybe use more descriptive names:
  TPar shift;        // resonance pitch-shift with respect to cutoff
  TPar gr;           // resonance gain as raw factor
  TPar pr;           // resonance start phase in radians
  TPar dr;           // resonance decay in seconds
  TPar decayByFreq;  // scaling of decay-time by resonance frequency
  TPar ar;           // resonance attack as factor 0 < k < 1 to be applied to dr to arrive at
                     // the absolute attack time
  TPar delay;        // resonance delay

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsFakeResonanceFilter<TSig, TPar>::getSample(TSig in)
{
  TSig yl = getLowpassOutput(in);
  TSig yr = getResonanceOutput(in);

  // at this point, we may later further shape the resonance signal's amplitude envelope, for
  // example by "ducking" it with a signal that is derived from the input (and/or lowpass output)
  // in order to mimic saturation effects in analog filters...

  return yl + yr;
  // maybe use a weighted sum - introduce a "mix" parameter that mixes between lowpass and 
  // resonance output
}

template<class TSig, class TPar>
RS_INLINE TSig rsFakeResonanceFilter<TSig, TPar>::getLowpassOutput(TSig in)
{
  TSig yl;
  yl = lpf1.getSample(in);
  yl = lpf2.getSample(yl);
  yl = lpf3.getSample(yl);
  yl = lpf4.getSample(yl);
  return yl;
}

template<class TSig, class TPar>
RS_INLINE TSig rsFakeResonanceFilter<TSig, TPar>::getResonanceOutput(TSig in)
{
  TSig yr;
  yr = hpf.getSample(in);
  yr = dl.getSample(yr);
  yr = rf.getSample(yr);
  return yr;

  // todo: try to change the order of resonator/delayline/highpass filter - as long as everything
  // is LTI, this should not make any difference, but as soon as nonlinearities are introduced
  // (saturation, sweep, etc.), and/or we use modulation, there will be a difference
  // i tend to think: 
  // delayline -> resonator -> highpass 
  // might be best: the delayline gets typically low-frequency rich signal which is good for the
  // interpolator and the post-resonator highpass may remove low-frequency aliasing components 
  // which are typically the most objectionable ones
  // but no - we can't put the highpass after the resonater if the resonator should introduce
  // waveshaping stuff.. so maybe: delay -> high -> reso is best

  // for the waveshaper: maybe use Chebychev polynomials to control the harmonics directly
}

#endif
