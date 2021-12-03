


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
 -the (pseudo)-randomness of that star-phase in different cycles of the input signal

-the basic setup of the filter is as follows:

 
 

 
 
Some hypotheses, how the perceptual parameters could correlate to signal features:
-the "silky/smooth/creamy" dimension could be related to the tendency of the resonance to lock into 
 the modes of the input signal - does it happen at all, how quickly, etc. The opposite would be,
 if the filter would just do its own thing (resonate at its own resonance freq, regardless of 
 whether an input mode is near it). The less the filter tends toward locking into input modes, the
 smoother a frequency sweep wil sound - the resonance does not "hop" from one freq to another but
 transitions smoothly. This can be achieved by passing the resonance through another resonator 
 filter - the narrower that filter is, the less the resonance is affected by the input modes (i 
 think)
-


*/