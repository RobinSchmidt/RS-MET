#pragma once






/* 
todo: 
-implement  an attack/decay envelope based on a difference of exponentials
-use the formulas from the modal filter - there's code to compute the time-constants and weights 
 from attack/decay settings
-later extend it to allow for a two-stage decay (maybe in a subclass)
-make it polyphonic - allow to build a very basic subtractive synth from such envelopes, a simple 
 osc class (maybe TriSaw osc?) and the rsLadder filter
 -the filter cutoff, osc-frequency and overall amplitude should be (polyphonically) modulated by 
  this simple envelope
 -maybe also have some simple LFO class
*/