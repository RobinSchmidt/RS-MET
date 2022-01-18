In this file, I collect some ideas for possible new features of the sampler engine.

Ideas for new opcodes
=====================

-reverse:  
Toggle reverse playback in the RegionPlayer. The user may to also want to use the offset opcode to 
determine where to start. If not defined, it starts at the last sample by default. Maybe 
loop_start and loop_end should also reverse roles in reverse mode.

-normalize:  
Apply an additional scale factor (directly at the RegionPlayer) that normalizes the sample on the 
fly according to different criteria. Possible options could be: none, max_abs, max_rms, loudness 
where the latter may employ some perceptual measure. To make this work, we need to analyze the 
sample on load and store the results somewhere. Maybe AudioFileStreamPreloaded could be an 
appropriate place. In a (potential future) direct-from-disk streaming mode, there perhaps needs to
be some metadata file next to the .wav - see below under "Misc Ideas".

-invert (or negate):  
Switches polarity of signal. Might be integrated into the Amplifier unit. But actually, that may be
redundant: negation should already be possible by a certain setting of width (may -100?), so maybe
we don't really need it. It could be convenient though. However, let's try to keep the number of
new opcodes low...

-filN_q_freqtrack (or resonanceN_freqtrack):  
Amount of tracking of the filter's quality factor Q (or resonance) of the filter's cutoff 
frequency. The default value should map to constant-Q behavior (maybe it should be 100%). The 
formula can be something based on what we have in RAPT where the filter's decay time constant 
scales with frequency.

-spectral_slope:  
A filter that applies an (approximate) slope to the spectrum of the signal. For example, 
spectral_slope=-3.01 would be a pinkening filter, i.e. a filter that turns white noise into pink
noise. A value of -6.02 would turn it into brown noise, +6.02 blue noise etc. Could also be called
color - but spectral_slope is probably less ambiguous. It could also benefit from key- and 
vel-tracking.

-param_mode (maybe find a better name):  
Switch the behavior of what happens when global, group and region all define values for the same 
opcode. The default mode specified by sfz is "override": group settings override global settings 
and region settings override global and/or group settings. We could also have "accumulate" and 
"levels_are_busses" modes. In the latter case, the outputs of the regions within a group are mixed
together (to form a sub-bus) and then *additional* effect processors are applied to that (sub)mix.
Then, the outputs of the groups are mixed together and to that master mix, the global effect 
processors are applied. So if there's a global setting for cutoff and also a group- and region 
setting for it, there will actually be 3 filters when a single layer is playing back a region.
When a chord of 3 notes (of the same region) are played, there will be 5 filters: 3 for the 3 
regions, an fourth one applied to their submix into a group ("sub-bus") and a fifth one to the 
global master mix ("master-bus"). We already have the behavior implemented for the override and 
levels_are_busses modes (but it's not yet opcode controlled). Maybe implement also the "accumulate" 
mode. In this mode, for example, there is no such thing as a group filter but the cutoff specified 
in the group gets (somehow) accumulated into the cutoff specified by the region. Maybe simple 
addition is indeed appropriate. But such an "accumulate" would place the burden on us to specify 
the accumulation behavior for each new opcode to be defined and it may not always be obvious, how 
this should be done. We'll see...for the time being, only "override" and "levels_are_busses" should 
be implemented.
 
-param_range:  
Decides, whether or not parameters should be restricted to the range specified by sfz. Possible 
values could be "sfz1" or "clipped", "free". Perhaps there could also be some mode that does clip 
but at different values than the sfz spec says.
 
-param_quantize:  
Some of the parameters that sfz wants to be integer can easily also admit float values. Examples 
for such parameters are tune, transpose, loop_start, loop_end.


Misc Ideas
==========

Meta Data for Audio Files:
--------------------------

Certain opcodes could take advantage of some metadata that is stored either directly in the .wav 
file itself or in a sidecar .meta file that sits next to the actual audio file. Examples for such
opcodes are: normalize, loop_start, loop_end, pitch_keycenter. Maybe the format could be xml based
but it could also be a simple textual representation similar to sfz. like:

max_abs=0 max_rms=-3.01 mean_pitch=69 max_pitch=75.3 min_pitch=67.8 loop_start=23.54 
loop_end=156.65 cycles_in_loop=1 category=single_note (others: single_drum, noise, polyphonic, 
speech, singing, melody/sequence/monophonic, drumloop, mixdown, effect, chord, environment)

Maybe write a little script that batch-generates such metadata files for a folder of samples for 
those parameters that are easily analyzed (such as max-abs, max_rms, maybe pitch, loudness)




