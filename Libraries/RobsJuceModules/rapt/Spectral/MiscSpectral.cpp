

template<class T>
std::vector<T> rsExpDecayTail(const rsSinusoidalPartial<T>& partial, int spliceIndex, T sampleRate)
{
  std::vector<T> t = partial.getTimeArray();
  std::vector<T> a = partial.getAmplitudeArray();
  int numFrames = (int) partial.getNumDataPoints();
  int maxIndex  = RAPT::rsArray::maxIndex(&a[0], numFrames);

  // todo: pass the maxIndex into the function call below:

  return rsExpDecayTail(numFrames, &t[0], &a[0], 
    spliceIndex, sampleRate, partial.getFreq(spliceIndex), partial.getPhase(spliceIndex));
  // maybe instead of using the instantaneous frequency at the splice index, we should use the 
  // average frequency ...hmm....well...that seems suitable for modal modeling of the whole partial
  // but maybe not so much for Elan's tail-splicing use case - maybe we should have both versions
}