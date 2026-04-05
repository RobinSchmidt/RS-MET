#ifndef RAPT_WAVEFORMS_H_INCLUDED
#define RAPT_WAVEFORMS_H_INCLUDED


/** A class for helping to produce various waveforms that are useful for sound synthesis. 

It includes various functions that take a normalized phasor value in the range 0..1 and converts
that phasor into an actual sample of the respective waveform. I have been deliberatery vague 
about the question whether we are talking about the half-open interval [0,1) or the closed interval
[0,1] here because both variants may make sense in different contexts. In the continuous time 
domain, there isn't really any material difference between both interpretations because the time 
interval between when the phasor is "just before 1" and when it reaches "actually 1 itself" is 
precisely zero because "just before" may get infinitesimally close to "precisely there". But in the
discrete time world, there is actually a tiny finite time interval, namely: one sample, between 
these two instants which will introduce a difference between both interpretations. Some waveforms, 
such as the sawtooth, prefer the phasor to be given in the closed interval while others, such as 
the sine, require it to be given in the half open interval. ...TBC...

It is planned to later add functions to render various single cycle waveforms into buffers suitable
for use in the context of table-lookup synthesis and to also add function to manipulate such 
waveforms. */

template<class T>
class rsWaveForms
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Conversions from phasors in [0,1] to waveforms

  /** Converts a phasor value into an upward sawtooth wave starting at -1 and ramping up to +1. */
  static inline T saw(T p) { return T(-1) + T(2) * p; }

  /** Converts a phasor value p into a pulse wave with adjustable pulse width given by pw. */
  static inline T pulse(T p, T pw = T(0.5)) { return p < pw ? T(-1) : T(+1); }


  //-----------------------------------------------------------------------------------------------
  // \name Conversions from phasors in [0,1) to waveforms

  static inline T sine(T p) { return rsSin(T(2*PI) * p); }

};


#endif