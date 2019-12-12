

template<class T>
std::vector<T> rsExpDecayTail(const rsSinusoidalPartial<T>& partial, int spliceIndex, T sampleRate)
{
  std::vector<T> t = partial.getTimeArray();
  std::vector<T> a = partial.getAmplitudeArray();
  int numFrames = (int) partial.getNumDataPoints();
  int maxIndex  = RAPT::rsArray::maxIndex(&a[0], numFrames);

  rsAssert(maxIndex < spliceIndex);
  // splicing point is supposed to be somewhere in the decaying section after the global maximum

  return rsExpDecayTail(numFrames, &t[0], &a[0], maxIndex, spliceIndex, sampleRate, 
    partial.getFreq(spliceIndex), partial.getPhase(spliceIndex), spliceIndex);
  // maybe instead of using the instantaneous frequency at the splice index, we should use the 
  // average frequency ...hmm....well...that seems suitable for modal modeling of the whole partial
  // but maybe not so much for Elan's tail-splicing use case - maybe we should have both versions
}