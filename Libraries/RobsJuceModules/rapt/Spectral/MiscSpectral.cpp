

template<class T>
std::vector<T> rsExpDecayTail(const rsSinusoidalPartial<T>& partial, int spliceIndex, T sampleRate)
{
  std::vector<T> t = partial.getTimeArray();
  std::vector<T> a = partial.getAmplitudeArray();
  int numFrames = (int) partial.getNumDataPoints();
  int maxIndex  = RAPT::rsArrayTools::maxIndex(&a[0], numFrames);

  rsAssert(maxIndex < spliceIndex);
  // splicing point is supposed to be somewhere in the decaying section after the global maximum

  return rsExpDecayTail(numFrames, &t[0], &a[0], maxIndex, spliceIndex, sampleRate, 
    partial.getFreq(spliceIndex), partial.getPhase(spliceIndex), spliceIndex);
  // maybe instead of using the instantaneous frequency at the splice index, we should use the 
  // average frequency ...hmm....well...that seems suitable for modal modeling of the whole partial
  // but maybe not so much for Elan's tail-splicing use case - maybe we should have both versions
}
// Maybe instead of using a second point (spliceIndex) to fit the exponential decay, use the total 
// energy (together with the maximum amplitude). That would free client code from having to specify
// this second point (for which it is not clear, how we could automate that decision)
// The total area under the function f(t) = A * exp(-a*t) is given by the definite integral from 0 
// to infinity of that function and comes out as A/a - maybe use that formula - or the area under 
// that function function squared (for the total energy). Here, "A" is the amplitude of the maximum 
// and "a" can be computed by using a = A / area where the area can be computed using numeric 
// integration (from the time-instant where the maximum occurs to the end) -> experiment with both 
// formulas (using area and energy) and see, which one gives the better fit
// -this can then be used to find optimal settings for the peak-shadowing algorithm

// maybe this idea can be refined to use the energy/area of an actualy attack/decay envelope - the 
// attack time is known - it's the time-instant of the maximum - with that, we may be able to get 
// more accurate values
