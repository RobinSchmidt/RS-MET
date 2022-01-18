The RS-MET Sampler Engine
===========================

...is still in early stages and very much under construction...

The RS-MET sampler engine is built around sfz which is a popular text-based file format to specify 
sample-based instruments. In sfz, the instrument is specified in terms of so called opcodes where 
each opcode controls a certain playback parameter. For example, an opcode specification like 
cutoff=1000 controls the cutoff frequency of a filter and sets it to 1000 Hz. The original sfz 
specification (henceforth referred to as sfz v1) specifies more than 100 of such opcodes. On top of 
that, various extensions exist, for example sfz v2, Cakewalk's Rapture, the ARIA engine, the Linux
Sampler Project, Sfizz and many more. The RS-MET engine also defines some extensions of its 
own. It thereby tries to be as compatible as possible to existing extensions ...tbc...



Opcode Support
-----------------


Listed below are the opcodes that are currently supported by the RS-MET sampler engine, categorized
by the specification/extension that first introduced them.

### SFZ v1


#### Player

sample, lokey, hikey, lovel, hivel, delay, offset, loop_mode (partial), transpose, tune, 
pitch_keycenter, pitch_keytrack

#### Filters

fil_type (partial), cutoff, resonance, fil_keytrack, fil_keycenter, fil_veltrack, eqN_gain,
eqN_freq, eqN_bw


#### Amplifiers

volume, pan, width (preliminary), position, amp_keytrack, amp_keycenter, amp_veltrack


### RS-MET

#### General Extensions

In sfz v1, there were 3 equalizers whose settings were controlled by the opcodes eqN_gain, 
eqN_freq, eqN_bw where N would be replaced by 1,2,3. In the RS-MET engine, N can be arbitrary so 
you can have an arbitrary number of equalizer bands. The same goes for the filter opcodes. In sfz 
v2, there were two filters and in v1 only one which is why in sfz v1 and v2 the settings of the 
first filter are just set by opcodes cutoff, resonance etc. (without any index) and in sfz v2, 
those of the second filter were set by cutoff2, etc. The RS-MET engine allows N to be arbitrary 
here as well where in the case of N=1, the 1 is optional to support the indexless sfz v1 (and v2) syntax for filter 1. The other indexed opcodes for which an index of N=1 is optional are:

filN_type, cutoffN, resonanceN, filN_keytrack, filN_keycenter, filN_veltrack, volumeN, panN, 
widthN, positionN, ampN_keytrack, ampN_keycenter, ampN_veltrack

You may notice that we have also volumeN, panN, etc. That means, the amplifier is also realized as
an effect and you can have as many amplifier units as you want. All effects are applied in series
and the order of the effects in the chain is determined by the first opcode that applies to a given
effect. For example, if you write into your sfz-file "volume=-6 cutoff=500 pan=50" then the 
amplifier (to which volume and pan apply) will be placed before the filter (to which cutoff 
applies). For effects of the same kind with an index, their order will be dictated by the index, so
if you write "eq2_gain=3 eq1_gain=6" then eq1 will be before eq2 regardless of order of appearance. The rule is: if there's an index, that index determines the position among the sibling effects. For
effects of different kinds, the rule of the first opcode applies. If you write "eq2_gain=6" without specifying any settings for eq1, then your effect chain will nevertheless contain two equalizers, 
but the first one will have neutral default settings so it won't do anything to your signal.


#### Parameter Quantization

In sfz, some opcodes were specified to take integer values but which could by their nature just as
easily admit floating point values. So, the RS-MET engine lifts this restriction where it makes 
sense. The affected opcodes are: ...tbc...


#### Parameter Ranges

The sfz spec prescribes a range, i.e. a minimum and and maximum value for each parameter. We are a 
bit more liberal with respect to these range limits ...tbc...


#### New Opcodes

The RS-MET engine also introduces some entirely new opcodes for new effects and settings. These 
are ...tbc...



Misc
----

### Acknowledgments

I'm especially grateful to Niall McCallum (www.modeaudio.com) for funding the initial development 
of this engine. It will be used in one of their upcoming products...


### SFZ Resources

https://sfzformat.com/                  A website dedicated to the sfz format  
https://sfzformat.com/legacy/           Original description of the sfz fomat and sfz v1 opcodes  
https://sfzformat.com/opcodes/          Comprehensive list of opcodes including extensions  
https://sfzformat.com/software/tools/   List of sfz authoring tools  
https://github.com/sfz                  sfz related stuff  
https://github.com/sfztools             Open source sfz tools  
https://github.com/sfztools/sfizz       An open source sfz engine  
https://www.linuxsampler.org/           Another one  
https://sfzinstruments.github.io/       Free sfz instruments  

