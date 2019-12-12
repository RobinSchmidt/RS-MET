#ifndef RAPT_MISCSPECTRAL_H_INCLUDED
#define RAPT_MISCSPECTRAL_H_INCLUDED

/** Stuff that has to do with spectral analysis/processing/synthesis and doesn't fit anywhere 
else. */

template<class T>
std::vector<T> rsExpDecayTail(
  const RAPT::rsSinusoidalPartial<T>& partial, int spliceIndex, T sampleRate);
// maybe this should go into some class dealing with modal synthesis

#endif