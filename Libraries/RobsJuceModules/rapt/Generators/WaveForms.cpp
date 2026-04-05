
//=================================================================================================
/*

ToDo:

- Maybe add an enum class to enumerate the different waveforms. But maybe that should belong to
  rosic
 
- Document the functions sawUp(), sawDown(), etc. The documentation should tell if the phasor is
  expected in the interval [0,1) or [0,1). I think, the default should be [0,1).

- Maybe sort the phasor conversion functions into two categories: Those that require [0,1) and 
  those that prefer [0,1]. I was deliberate about the wording "prefer" in the 2nd case. I think, 
  they do not strictly "require" it. They would still kinda work with [0,1) but using [0,1] is 
  just a little bit nicer and more perfect.

- triangle, sin, triSaw, 
 
- pulseNoDc(): Should subtract the DC that would otherwise be there. I think, we need to 
  add or subtract (pw - 0.5) and maybe we need to scale it. We should then have a unit test 
  that verifies that the DC is indeed zero for all pulse-widths in 0..1.

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

- Add functions to produce various waveforms, including additively synthesized (brickwall 
  lowpassed) saw and pulse waves (maybe by using trig-recursions for the sines of the various
  frequencies for an optimized implementation)

- See also other places where we have implemented similar functionality, for example:
    rsBlepReadyOscBase  in  rapt/Unfinished/MiscAudio/BlepBlampOscs
    rsTriSawOscillator  in  rapt/Generators/VariousOscillators.h
    rsTriSaw            in  rapt/AudioBasics/AudioFunctions.h
    rsPulseWave, ...    in  rapt/Math/Functions/RealFunctions.h
  This code should ideally all be consolidated in this new class rsWaveForms. There may be 
  more places. I think in RealFunctions may also be some waveform-producing functions. Also,
  somewhere in rosic in the waveform rendering functions. I think, there's also an enum for 
  enumerating the waveforms. If the new class also gets such an enum, it should be compatible
  with the old one. Then deprecate the old functions.

- Maybe the functions to manipulate prototype waveforms that we have in Straightliner's 
  oscillators could also be moved here.



*/