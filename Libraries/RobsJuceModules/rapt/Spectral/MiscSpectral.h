#ifndef RAPT_MISCSPECTRAL_H_INCLUDED
#define RAPT_MISCSPECTRAL_H_INCLUDED

/** Stuff that has to do with spectral analysis/processing/synthesis and doesn't fit anywhere 
else. */


/** Given a reference to a sinusoidal partial, this function assumes that the partial is of 
exponentially decaying type and synthesizes a perfectly expontially decaying signal that can be 
used to replace the tail of the original partial. The function will estimate the desired amplitude
and decay-rate from two points in the amplitude envelope in the partial: (1) the global maximum of
the envelope and (2) a user-defined time-instant given in terms of the frame-index within the 
partial. It is assumed that client code will do a crossfade between original partial and the 
synthesized decaying sinusoid centered around the time-instant that corresponds to the 
spliceIndex, i.e. the time-instant that is recorded in the time-array in the partial at the given 
index. The sinusoid will be synthesized to match the instantaneous amplitude and phase at that
time-instant. */
template<class T>
std::vector<T> rsExpDecayTail(
  const RAPT::rsSinusoidalPartial<T>& partial, int spliceIndex, T sampleRate);
// maybe this should go into some class dealing with modal synthesis

#endif