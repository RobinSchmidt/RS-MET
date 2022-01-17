The RS-MET sampler engine is built around sfz which is a popular text-based file format to specify 
sample-based instruments. In sfz, the instrument is specified in terms of so called opcodes where 
each opcode controls a certain playback parameter. For example, an opcode specification like 
cutoff=1000 controls the cutoff frequency of a filter and sets it to 1000 Hz. The original sfz 
specification (henceforth referred to as sfz v1) specifies nearly 200 of such opcodes and various
extensions to this basic set exist, for example in sfz v2, the ARIA engine, Cakewalk's Rapture, the linux sampler project, sfizz and many more. The RS-MET engine defines also some extensions of its 
own. It thereby tries to be as compatible as possible to existing extensions ...tbc...



Opcode support
##############

SFZ v1
======

Player
------

sample, lokey, hikey, lovel, hivel, delay, offset, loop_mode (partial), transpose, tune, 
pitch_keycenter, pitch_keytrack

Filters
-------

filN_type (partial), cutoffN, resonanceN, filN_keytrack, filN_keycenter, filN_veltrack, eqN_gain,
eqN_freq, eqN_bw


Amplifiers
----------


Notes
-----

For some opcodes that include an index such as eqN_gain, cutoffN, etc. the N is optional when N=1 
to cater for the fact that sfz1 only allowed one of those devices. These are:

filN_type, cutoffN, resonanceN, filN_keytrack, filN_keycenter, filN_veltrack, volumeN, panN, 
widthN, positionN


