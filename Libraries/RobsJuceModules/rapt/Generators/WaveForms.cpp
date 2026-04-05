
//=================================================================================================
/*

Notes:

- To keep the API lean, we deliberately do not provide waveforms that can be obtained by trivial 
  transformations of the existing ones such as sawDown (which is just a negative saw) or cosine
  (which is just a phase-shifted sine). These shall be created by the calling code from the given
  ones.


ToDo:

- Add the waveforms: triangle, triSaw, ...
 
- pulseNoDc(): Should subtract the DC that would otherwise be there. I think, we need to 
  add or subtract (pw - 0.5) and maybe we need to scale it. We should then have a unit test 
  that verifies that the DC is indeed zero for all pulse-widths in 0..1.

- Add functions to produce additively synthesized (brickwall lowpassed) saw and pulse waves (maybe
  by using trig-recursions for the sines of the various frequencies for an optimized 
  implementation). I think, there is some prototype code in in the rs_testing module that can 
  render sawtooth and pulse-waves using truncated Fourier series. That code could be helpful for 
  this. But thes code there is naively implemented, i.e. doesn't use trig-recursions, so we may 
  keep that as prototype implementation and use it in unit tests to render reference signals for 
  the optimized implementations.

- See also other places where we have implemented similar functionality, for example:

    rsBlepReadyOscBase  in  rapt/Unfinished/MiscAudio/BlepBlampOscs
    rsTriSawOscillator  in  rapt/Generators/VariousOscillators.h
    rsTriSaw            in  rapt/AudioBasics/AudioFunctions.h
    rsPulseWave, ...    in  rapt/Math/Functions/RealFunctions.h

  This code should ideally all be consolidated in this new class rsWaveForms. There may be 
  more places. I thnk, somewhere in rosic we have waveform rendering functions. I think, there's
  also an enum for enumerating the waveforms. If the new class also gets such an enum, it should be
  compatible with the old one. Then deprecate the old functions. Maybe do an AI assisted search for
  all places where we deal with waveforms and try to consolidate as much as possible of this code
  here.

- Maybe add an enum class to enumerate the different waveforms. But maybe that should belong to
  rosic

- Add to the class documentation a paragraph that explains for some example waveforms (maybe use 
  saw and sine) why some waveforms prefer a closed interval [0,1] and others require an open 
  interval [0,1). I was deliberate about the wording "prefer" in the first case. I think, 
  they do not strictly "require" it. They would still kinda work with [0,1) but using [0,1] is 
  just a little bit nicer and more perfect. I think, the general rule is: when there is a jump
  discontinuity at the wrap-around point, the closed interval is preferred. If the waveform 
  smoothly joins back to itself, we want the half-open interval because otherwise, a value gets 
  repeated. It would occur once at the end of the old cycle and then once again at the start of the
  new cycle - and this is very wrong. In the other case where we use an interval of [0,1) for waves
  that prefer [0,1], we would just modify the height of the jump and perhaps introduce some DC both
  of which is much less problematic, I think.

- Create functions for realtime processing that take a phasor and covert it into a waveform
  and also create functions that produce or manipulate the whole waveform. Maynipulations 
  could be things like reversal, negation, circular shift, fractalization, etc. Maybe for 
  things like fractalization, it could make sense to compute in double precision even when the
  type is float. Fractalization should mix in waveforms of octaves above the original. These 
  higher waveforms could themselves be manipulated with shifts, negation, reversal, etc. using 
  alternating patterns. We should get self-similar waveshapes out of this.

- Implement functions like reverse, shift, etc. also for phasors. I think, reverse would just
  be 1-p, shift would be (p+s) % 1 where % 1 would be fmod - but maybe a special variant that
  leaves 1 as is i.e. does not if(x >= 1) return x-1  but rather if(x > 1) return x-1. Maybe 
  we could also "fractalize" a phasor value? Try it!

- Verify that the formula for pulse is what the user would expect. Maybe we should swap -1 and
  +1? And/or maybe we should use if(p <= pw) rather than if(p < pw). Document these decsisions
  and the reasons behind them. One reason to prefer to have the negative half-cycle first is
  that this would be compatible with clipping a saw-up waveform and I think, the "up" variant
  is the default expectation in case of a saw wave. Check what popular synthesizers do (Surge,
  Serum, Diva, JP-8000, ...) and maybe do the same. Maybe to figure out if < or <= is correct,
  consider a square wave with an even integer cycle length. In such a case, we want the
  positive and negative half-wave to have exactly the same number of samples. This may also 
  depend on whether the phasor range is [0,1] or [0,1). 

- Maybe the functions to manipulate prototype waveforms that we have in Straightliner's 
  oscillators could also be moved here.

- Add waveform rendering functions like renderSaw(T* buffer, int length). Or maybe make another
  class rsWaveFormRenderer for these rendering functions. Maybe the class rsWaveForms should be
  renamed to rsPhasorToWaveConverter. Splitting the functionalities of converting phasors and 
  rendering full waveforms into buffers into separate classes may potentially reduce code size 
  when some parts of the functionalities is needed only for float and some others for both float
  and double because we can then have our template instantiations at a finer granularity so we can
  be more picky about what functionality gets instantiated for which data type.

- Let's assume that client code reverses the waveform by passing 1-p instead of p. Will that play
  nicely with phasors in the interval [0,1)? The reversal would convert the interval to (0,1]. It 
  should probably be unproblematic, I guess. With the [0,1] interval, it's immediately obvious that
  there's no problem because it will get converted into the exact same interval [0,1] by the
  reversal. Create an experiment that tests this (perhaps using saw and sine as examples for both
  cases) and if it's all ok, add to the documentation that client code can reverse the waveforms in
  this way.

*/