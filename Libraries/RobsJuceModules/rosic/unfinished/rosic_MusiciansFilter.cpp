

void rsResoSplitFilter::processFrame(float inL, float inR,
  float* fltOutL, float* fltOutR, float* resOutL, float* resOutR)
{
  rsFloat32x4 x(inL, inR, inL, inR);
  rsFloat32x4 y = ladder.getSample(x);
  *fltOutL = y[2];
  *fltOutR = y[3];
  *resOutL = y[0] - *fltOutL;
  *resOutR = y[1] - *fltOutR;
}

//=================================================================================================

rsResoWaveFilter::rsResoWaveFilter()
{


}

void rsResoWaveFilter::setCutoff(float newCutoff)
{


}

void rsResoWaveFilter::setResonance(float r)
{
  rsFloat32x4 reso(r, r, 0.f, 0.f);

}

//=================================================================================================

/*
Notes:

-the name "MusiciansFilter" should contrast it with "EngineersFilter"
-this filter is made for musical purposes, mainly for use as filter in the context of subtractive 
 synthesis
-in contrast to EngineersFilter, the following aspects are important:
 -it should respond nicely to modulation
 -it should have a good, "musical" response to changes to control parameters
  -> the control parameters should have a perceptually linear scaling behavior 
  -> the control parameters should be perceptually decoupled from one another
 -it may show some weird, nonlinear behaviors - ideally, the amount of these behaviors should be 
  continuously adjustable
-the following aspects are less important:
 -we don't care much about technical parameters like filter slope, ripple, etc.
 -we don't need strict linearity and cleanliness
 
-some perceptual dimensions for which we may want to give a slider to the user could be:
 growly, noisy, bubbly, silky, screamy, gnarly, warm, etc.
 
-some features of the filter's resonance that are important to character are: 
 -the amplitude of the resonance
 -the bandwidth of the resonance, i.e. the ringing time
 -the micro-envelope of the resonance waveform (i.e. it's amplitude envelope)
 -the micro-envelope of the frequency of the resonance 
 -the shapes of these micro-envelopes
 -micro-envelopes can be derived from the (filtered) input
 -tendency to lock into modes of the input vs doing its own thing (i.e. resonate at its own 
  resonance frequency with regard to the input signal), i think, this is related to resonance 
  bandwidth - post filtering it more strongly makes it narrower and less prone to lock into input 
  modes (->verify!)
 -waveshape of the resonance. the basic ladder produces a sinusoidal waveshape
 -start-phase of the waveshape
 -the (pseudo)-randomness of that start-phase in different cycles of the input signal
 -how the resonance scales with frequency - constant-Q vs constant bandwidth, etc.

-The basic setup of the filter is as follows:
 -We use two ladder filters in parallel, one with and one without feedback
 -We subtract the signal without resonance for the signal with resonance to obtain the pure 
  resonance signal.
 -We post-process the resonance signal in various ways
 -We produce the output as sum of the non-resonant signal and the post-processed resonance signal.
 
The post-processing should include the following
 -Adjustment of overall amplitude (reso-gain)
 -Narrowing the bandwidth of the resonance signal (which also amounts to increasing its ringing 
  time) by passing it through a resonator of bandpass
 -This resonator/bandpass could be realized using the state-vector filter or the (nonlinear) modal
  filter
 -Estimate instantaneous amplitude and phase of the resonance (...should this be done on the raw or
  narrowed resonance? If the latter, we should make sure that the nonlinearities do not destroy 
  the sinusoidal shapes - at least not too much)
  -> The instantaneous phase will be a sawtooth like signal. Maybe we can anti-alias it with a 
     PolyBlep (we would need to detect the phase wrap-around of the reso-wave between samples)
  

 
 

 
 
Some hypotheses/speculations, how the perceptual parameters could correlate to signal features:
-the "silky/smooth/creamy" dimension could be related to the tendency of the resonance to lock into 
 the modes of the input signal - does it happen at all, how quickly, etc. The opposite would be,
 if the filter would just do its own thing (resonate at its own resonance freq, regardless of 
 whether an input mode is near it). The less the filter tends toward locking into input modes, the
 smoother a frequency sweep wil sound - the resonance does not "hop" from one freq to another but
 transitions smoothly. This can be achieved by passing the resonance through another resonator 
 filter - the narrower that filter is, the less the resonance is affected by the input modes (i 
 think)
-"bubbly" could be related to the frequency envelope of the resonance
-"growly" could be related to a pseudo-random change of the resonance phase between cycles

User parameters could be:
-Frequency, Feedback, Reso-Gain, Reso-Bandwidth, waveshape, (micro)env-amount, phase-randomization 
 (via some nonlinear feedback algo), waveshape as function of input (i.e. somehow derived form the
 filtered input)
 
-The following signals are readily available: overall input, the output taps of the non-resonating 
 ladder, the output taps of the resonating ladder, the raw resonance signal, the post-processed 
 resonance signal (actually, at all stages of the post-processing)
-The following signals could be made available with relative ease: a quadrature component of the 
 resonance (via an allpass), the instantaneous amplitude of the raw resonance (uses the quadrature 
 component and the original reso - square both, add, take sqrt - or maybe use the squared env 
 directly -> cheaper)
-The following signals could be made available with some more processing: the input amplitude 
 (using an envelope follower), ...
->these signals are the ingredients that we can work with to achieve a desired behavior of the 
  resonance

More ideas:
-Maybe one could use the input signal itself for the resonance waveform? Maybe always buffer the most 
 recent N samples where N = fc/fs and use that as wavetable?
 
ToDo:
-maybe call the filter ResoWave
-maybe allow (maybe in a subclass) to use single-cycle waveforms for the resonance
-maybe allow to use a wavetable where the table index is selected by some feature of the input signal
 like overall amplitude
-make an efficient version of the Ladder filter using SIMD to compute the resonant and non-resonant
 part simultaneously - maybe use rsFloat32x4, for stereo...oh - that means, this class should go to 
 rosic
-move the experimental filters from rapt/Unfinished/Filters into prototypes





-maybe factor out an rsResoSplitFilter that just as the ladder and does the splitting between 
 resonant and nonresonant
-this could actually be dragged to rapt and turned into a template having one template parameter
 T that can be float or double, we ould then use rsSimdVector<T, 4>

*/