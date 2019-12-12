

template<class T>
std::vector<T> rsExpDecayTail(const rsSinusoidalPartial<T>& partial, int spliceIndex, T sampleRate)
{
  return rsExpDecayTail(partial.getTimeArray(), partial.getAmplitudeArray(), spliceIndex, sampleRate, 
    partial.getFreq(spliceIndex), partial.getPhase(spliceIndex));
  // maybe instead of using the instantaneous frequency at the splice index, we should use the 
  // average frequency ...hmm....well...that seems suitable for modal modeling of the whole partial
  // but maybe not so much for Elan's tail-splicing use case - maybe we should have both versions
}