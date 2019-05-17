template<class T> 
rsEnvelopeFollower2<T>::rsEnvelopeFollower2()
{
  preFilter.setApproximationMethod(rosic::rsPrototypeDesignerD::BUTTERWORTH);
  //preFilter.setApproximationMethod(rosic::rsPrototypeDesignerD::BESSEL);
  preFilter.setOrder(6);

  slewLimiter.setAttackTime(1.0);
  slewLimiter.setReleaseTime(100.0);
  slewLimiter.setHoldTime(0.0);

  minMaxSmoother.setSlewRateLimiting(1.0);
  minMaxSmoother.setMinMaxMix(0.5);

  postFilter.setApproximationMethod(rosic::rsPrototypeDesignerD::BESSEL);
  postFilter.setOrder(6);

  setSampleRate(44100);
}

template<class T> 
void rsEnvelopeFollower2<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;

  preFilter.setSampleRate(sampleRate);
  preFilter.setFrequency(sampleRate/T(6));

  slewLimiter.setSampleRate(sampleRate);

  postFilter.setSampleRate(sampleRate);

  updateSmoothingFilters();
}

template<class T> 
void rsEnvelopeFollower2<T>::setInputFrequency(T newFreq)
{
  inputFreq = newFreq;
  updateSmoothingFilters();
}

template<class T> 
void rsEnvelopeFollower2<T>::reset()
{
  preFilter.reset();
  slewLimiter.reset();
  minMaxSmoother.reset();
  postFilter.reset();
}

template<class T> 
void rsEnvelopeFollower2<T>::updateSmoothingFilters()
{
  minMaxSmoother.setLength((int) ceil(sampleRate/inputFreq)); // one cycle
  postFilter.setFrequency(inputFreq/T(2));
}