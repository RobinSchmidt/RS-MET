In this file, I collect some ideas for possible new features of the sampler engine.

Ideas for new opcodes
=====================

- reverse: Toggle reverse playback in the RegionPlayer. The user may to also want to use the offset 
  opcode to determine where to start. If not defined, it starts at the last sample by default. 
  Maybe loop_start and loop_end should also reverse roles in reverse mode.

- normalize: Apply an additional scale factor (directly at the RegionPlayer) that normalizes the 
  sample on the fly according to different criteria. Possible options could be: none, max_abs, 
  max_rms, loudness where the latter may employ some perceptual measure. To make this work, we need 
  to analyze the sample on load and store the results somewhere. Maybe AudioFileStreamPreloaded 
  could be an appropriate place. In a (potential future) direct-from-disk streaming mode, there 
  perhaps needs to be some metadata file next to the .wav - see below under "Misc Ideas".

- retrigger_mode: Controls what happens when the same note is retriggered that is currently still 
  playing (presumably in release phase). Available modes could be: new_layer: just triggers a new 
  layer and the old one also keeps playing. retrigger_envs: re-uses the currently playing layer but
  re-triggers its envelopes. retrigger_all: retriggers also LFOs and restarts the sample (maybe 
  there should be a very short crossfade between the playing and restarted sample). new_and_choke: 
  triggers a new layer and shuts off the old one immediately with quick fade-out over a few 
  milliseconds maybe. May also be called new_and_suppress, new_and_shutdown, new_and_brake, 
  new_and_stall, new_and_stop. ...Check, if there is already something similar in sfz. Maybe the 
  quick fade out could have an additional time parameter - like new_and_stop_10ms for a 10 ms 
  fade-out. new_layer could also be called overlap or self_overlap. The envelope retriggering could
  itself operate in different modes: start from current, start from start, additive. Maybe these 
  options could be appended to the retrigger_envs option..or maybe retigger=envs_from_current, 
  envs_from_start, envs_additive, overlap, restart, quick_fade_10ms

- invert (or negate): Switches polarity of signal. Might be integrated into the Amplifier unit. But
  actually, that may be redundant: negation should already be possible by a certain setting of 
  width (maybe -100?), so maybe we don't really need it. It could be convenient though. However, 
  let's try to keep the number of new opcodes low...

- filN_q_freqtrack (or resonanceN_freqtrack): Amount of tracking of the filter's quality factor Q
  (or resonance) of the filter's cutoff frequency. The default value should map to constant-Q
  behavior (maybe it should be 100%). The formula can be something based on what we have in RAPT 
  where the filter's decay time constant scales with frequency. It may perhaps be parametrized as
  ringing time tracking, i.e. filN_ringing_freqtrack or as decay time filN_decay_freqtrack.

- spectral_slope: A filter that applies an (approximate) slope to the spectrum of the signal. For
  example, spectral_slope=-3.01 would be a pinkening filter, i.e. a filter that turns white noise 
  into pink noise. A value of -6.02 would turn it into brown noise, +6.02 blue noise etc. Could also
  be called color - but spectral_slope is probably less ambiguous. It could also benefit from key-
  and vel-tracking.
  
- functionN: Let the user define an arbitrary function (as string) to be applied to the signal. For 
  example: ``functionN="tanh(a*x + b) / a - tanh(b)" functionN_a=1.7 functionN_b=-0.2``. It could 
  use  the same expression evaluator engine as in FuncShaper. Maybe it could also support a syntax 
  like: ``"m = x1 + x2; s = x1 - x2; m = a*m; s = b*s; y1 = m+s; y2 = m-s"`` where x1,x2 are the 
  two input channels and y1, y2 are the output channels (I deliberately do not use xL, xR for 
  "left" and "right" because maybe later, we want to support other multichannel configurations, so 
  it seems better to be generic). Maybe that functionality could also be integrated into the 
  waveshaper. But that would make it quite heavyweight. Maybe the engine should automatically 
  switch between loading a lightweight implementation when a predefined shape is used and only load 
  the heavyweight, expression evaluator based, variant when it's actually used in the sfz. We may 
  actually want to apply a similar strategy to load different filter implementations depending on 
  the selected filter type. Many filters will often be simple 1-poles but some will be 
  sophisticated virtual analog models (or maybe EngineersFilter) and it will waste resources to 
  load the big filter when a smaller one would suffice so we may have devise such a strategy 
  anyway. Also, some filter types actually have totally disjoint features (think, VA and elltiptic,
  for example), so it doesn't make sense to lump them all into one big uberfilter class when many
  of such objects are needed. The opcode may also be called formulaN. Maybe that's better.

- tableN: Works like functionN but instead of specifying the function as a string, it loads a 
  table. The table could be stored in a .wav file or maybe in a .txt or .dat file like the ones 
  generated by GnuPlotCPP, i.e. use a GnuPlot compatible textual syntax. Maybe we could also have
  tableN_interpolation=linear (other options: next_neighbor, cubic_spline, cubic_hermite, ...)

- The LFO could have additional "sequencing" features to control the speed. For example 
  lfo1_speeds=1,2,3,5 lfo1_repeats=1,2,4 would order the LFO to play one cycle a 1x speed, then 2 
  cycles at 2x speed, then 4 cycles at 3x speed, thene 1 cycle at 5x speed, then 2 cycles at 1x 
  speed, etc. Could be useful for dubstep wobbles.

- param_mode (maybe find a better name): Switch the behavior of what happens when global, group and 
  region all define values for the same opcode. The default mode specified by sfz is "override": 
  group settings override global settings and region settings override global and/or group 
  settings. We could also have "accumulate" and "levels_are_busses" modes. In the latter mode, the 
  3 hierarchy levels of the sfz specification (region, group, global) are mapped to busses: regions 
  form the single channels, groups form sub-busses mixing togther the regions/channels within the 
  group and the global instrument mixes together the groups/sub-busses into a master-bus. When 
  effect opcodes are specified for a group and/or globally, they don't merely provide fallback 
  values for the regions (as they would in normal sfz operation) but instead specify settings for
  *additional* effect processors to be applied to the (sub- or master) bus. So if there's a global 
  setting for cutoff and also a group- and region setting for it, there will actually be 3 filters 
  when a single layer is playing back a region. When a chord of 3 notes (of the same region) is 
  played, there will be 5 filters: 3 for the 3 regions, a fourth one applied to their submix into a
  group ("sub-bus") and a fifth one to the global master mix ("master-bus"). We already have the
  behavior implemented for the override and levels_are_busses modes (but it's not yet opcode 
  controlled). Maybe implement also the "accumulate" mode. In this mode, for example, there is no
  such thing as a group filter but the cutoff specified in the group gets (somehow) accumulated 
  into the cutoff specified by the region. Maybe simple addition is indeed appropriate. But such
  an "accumulate" would place the burden on us to specify the accumulation behavior for each new 
  opcode to be defined and it may not always be obvious, how this should be done. We'll see...for 
  the time being, only "override" and "levels_are_busses" are implemented. This feature totally 
  alters the signal flow and is quite a big deal. When a waveshaper is set up in a group or 
  globally, this will now also allow intermodulation distortion between the regions, such that
  playing of powerchords becomes a possibility.
 
- param_range: Decides, whether or not parameters should be restricted to the range specified by 
  sfz. Possible values could be "sfz1" or "clipped", "free". Perhaps there could also be some mode 
  that does clip but at different (i.e. extended) values than what the sfz spec says. For example, 
  the range for the width opcode may be extended to +-200% (sfz specifies the range to be +-100% 
  which allows only stereo narrowing but not widening - which is a bit sad).
 
- param_quantize: Some of the parameters that sfz specifies to be integer can easily also admit 
  float values. Examples for such parameters are tune, transpose, loop_start, loop_end. By 
  specifying param_quantize=none, we could lift these unnecessary restrictions. 

In the current implementation, we actually do not quantize the integer parameters anyway nor do we 
clip parameters to their specified ranges, so currently we already operate in that unrestricted 
mode. It would actually take some extra programming to enforce these restrictions and the "gain" 
would be *less* flexibility. Maybe we should or maybe we shouldn't do it. If users want to restrict 
themselves to using only integers for certain parameters, they can already do it. No additional 
code needed. The same goes for the ranges. However, in some cases we may indeed want to put some
limits on the ranges to ensure bibo stability (think of feedback) or sane resource requirements 
(think of delay) or safety for the user's equipment and/or ears (think of gain). But maybe such 
nannying of sfz authors is inappropriate - after all, if the author "programs" a "buggy" sfz, it 
may be argued that it's "their bug". I think, instead of safeguarding the user from mistakes by 
imposing (arbitrary) limits, I'd rather go with the "trust the programmer" approach. 
...but I'm not sure yet...we'll see...

Maybe the defined extensions should be prefixed with rs_ to avoid potential future name-clashes? 
SFZ v2 has a "type" opcode for "vendor-specific effect name" which can be, for example, 
"com.mda.Overdrive". Maybe we should use such a syntax for our special effects? There is also the
"vendor_specific" opcode defined by ARIA - but I have no idea how that is supposed to be used.


Sound Synthesis
===============

It may be nice to turn this into some sort of hybrid sampler/synthesizer thing. We could use the 
existing architecture as is and just define the synthesis units as new "effect" processors that 
actually generate some sound instead of (or in addition to) processing input signals. Maybe they 
could respond to a "gate" signal at their input which should be just a constant value of either 
0 or 1. They should also do something sensible in case of other values or arbitrary input signals, 
though. 

### What is possible already?

A couple of synthesizer-ish things can already be done with the standard sampler-ish features 
without actually adding any dedicated sound generator "effects". Among these are:

- A simple oscillator can be set up by just loading a single cycle waveform file and setting the 
  loop appropriately, so a simple oscillator with arbitrary waveform is already
  covered. Together with the filter(s), this should allow us to produce a wide range of typical 
  subtractive synthesizer sounds. Of course, we can also use multiple such oscillators, detune them
  to produce beating and supersaw sounds. A highpass filter may be beneficial to get the typical
  curved analog-like sawtooth (and square) shapes. A midi controller could be routed to the tune 
  opcodes of the oscillators to control the supersaw's "Detune"

- The filter (in bandpass mode) with keytracking of 100% and a white noise input sample can be 
  used to synthesize "whistle" sounds. A second bandpass with possibly a different resonance 
  setting can be used to further shape the attack. 

- The equalizers can be used to apply formant-like spectral peaks and anti-formant-like troughs to 
  any source signal. Ideally, the source signal should not already have a formant structure of its
  own. It could be an instrument sample that was pre-processed to flatten away the formants or it 
  could be rendered by some algorithm - maybe use something like Karplus-Strong to render a "raw" 
  pluck sound and then apply formants in the sampler to shape it further.

- A sample containing an impulse-like signal (possibly just a single value of 1 followed by
  all zeros) can be fed into the filter(s) with high resonance settings. The resonating filters 
  will produce decaying sinusoids. When keytracking is at 100%, these resonance "blips" are 
  musically playable. Putting a second such filter in series can be used for a smoother attack. 
  The lowpass mode may also be nice for that. This can also be used to beef up bassdrum sounds 
  with a low-frequency attack/decay sinusoid. In the higher ranges, it produces a sort of 
  "Popcorn" sound (https://www.youtube.com/watch?v=NjxNnqTcHhg) (...could use some reverb though 
  and maybe also waveshaping to give it some overtones)

### What else could be needed?

- Another interesting unit could be a comb filter with feedback and possibly a lowpass in the
  feedback patch. With an impulse like input, we could synthesize Karplus-Strong like "plucked 
  string" sounds. When we introduce a flanger unit, we should have an eye on enabling such 
  things in addition to the usual "flanging". It would primarily be a flanger but could be abused
  as Karplus-Strong synthesizer.
  
- We may need to allow the layer to have a ring-out phase where the RegionPlayer keeps playing 
  for some specified amount of time even after the sample has finished. Or it could be mocked by
  using a sample that has a tail containing some silence. By using a loop, that silence portion
  could be very short because the looping would lengthen it as needed.

- LFOs could be allowed to produce audio-rate signals and have keytracking of their frequencies. 
  When routed to pitch, we could mock frequency modulation synthesis. When the input sample is a 
  square wave, maybe we can also mock pulse-width modulation. We would need a linear mode of 
  operation of the pitch LFO - I'm not sure, what sfz wants. If sfz prescribes logarithmic pitch
  modulation for the pitch LFO, we may introduce another opcode to switch that or we could define
  separate opcodes lfoN_to_pitch and lfoN_to_freq to distinguish between these modes. Maybe the 
  modulators should produce stereo signals, too. Think of an LFO mapped to cutoff with a phase
  offset between left and right channel. That will easily give nice stereo movement. However,
  not for all parameters does it make sense to have different left/right values. Think, for 
  example, of a pan or width parameter - there is no meaningful concept of applying different
  pan values to left and right channels. I'm not yet sure, how to handle that. Maybe let's have two
  generic outputs, i.e. channel1/channel2 instead of left/right and how they are used is opcode
  defined - some (like cutoff) may interpret them as left/right, others (like pan) may just ignore 
  the second channel. 

- A reverb algo with an impulsive input sample could be used to create noise-burst like source 
  signals which could be shaped further by key-tracked comb and/or bandpass filters. Maybe a 
  series of narrow peak equalizers can be used instead of bandpasses.
    
 - Maybe we could have modal filter bank available as effect algorithm. With impulse and/or noisy
   inputs, we could do some modal synthesis. Maybe each modal filter could generate 8 or 16 modes.
   It would need some sort of dry/wet or pass-through parameter to run it in parallel with the 
   input. That could actually make sense for the regular filter as well. The opcodes could be 
   modalN_freqM, modalN_phaseM, modalN_gainM, modalN_decayM, modalN_attackM, modalN_through. But 
   no - even with the pass-through parameter we wouldn't get a parallel connection of the
   individual modal filters - so maybe a single one should support more than 8 modes...maybe up to
   64?


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
speech, singing, melody/sequence/monophonic or music_monophonic, music_polyphonic, 
music_mixdown/full_mix, drum_loop, instrument_loop, effect, chord, environmental, atmo/ambience, 
foley)

Maybe write a little script that batch-generates such metadata files for a folder of samples for 
those parameters that are easily analyzed (such as max-abs, max_rms, maybe pitch, loudness). Other
opcodes: tempo=140 (in bpm), freq_range=treble (also: bass, low_mid, mid, high_mid)...maybe that 
could also be a more quantitative spectral_centroid measure. Also interesting: spectral_slope - 
could be useful to the up the spectral_slope filter. Crest factors: mean_crest, max_crest, 
min_crest. Statistical measures: mean, variance, skew, kurtosis ..maybe also as short-time 
min/mean/max values. Maybe also analyze the stereo_width, dynamic_range...well, now we are getting 
into mastering territory...that's a bit beyond the scope of a sampler, but we are actually just 
talking about general metadata now. So...yeah - maybe we should define an sfz-like metadata format
for audio files. That could be useful for the sampler but also in other contexts. Maybe .amd for
"audio meta data" - not sure, if a certain CPU vendor would be happy with this choice, though...
how about .adi (audio data info), .afi (audio file info), .adp (audio data properties)




