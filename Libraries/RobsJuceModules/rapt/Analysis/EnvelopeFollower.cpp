template<class T> 
rsEnvelopeFollower2<T>::rsEnvelopeFollower2()
  : minMaxSmoother(10000) // constructor needs a max-length - todo: provide setMaxLength function
{
  typedef rsPrototypeDesigner<T> PD;

  //preFilter.setApproximationMethod(PD::BUTTERWORTH);
  preFilter.setApproximationMethod(PD::BESSEL);
  preFilter.setOrder(6);

  slewLimiter.setAttackTime(1.0);
  slewLimiter.setReleaseTime(100.0);
  slewLimiter.setHoldTime(0.0);

  minMaxSmoother.setSlewRateLimiting(1.0);
  minMaxSmoother.setMinMaxMix(0.5);

  postFilter.setApproximationMethod(PD::BESSEL);
  postFilter.setOrder(6);

  setSampleRate(44100);
}

template<class T> 
void rsEnvelopeFollower2<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;

  preFilter.setSampleRate(sampleRate);
  preFilter.setFrequency(sampleRate/T(6));
  // maybe instead of using a fixed fraction of the sample-rate, we should use a fixed absolute
  // frequency? ...but maybe not - the test uses a sample-rate of 2kHz, so a fixed freq would end 
  // up above fs/2

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

  // experimental - maybe attack and release of the slewrate limiter should be set proportionally
  // to the input period...
  T inputPeriod = T(1) / inputFreq;
  T att = T(1000 *  0.1) * inputPeriod;   // attack in milliseconds  (1/10 of a cycle)
  T rel = T(1000 * 10.0) * inputPeriod;   // release in milliseconds (10 cycles)
  slewLimiter.setAttackTime(att);
  slewLimiter.setReleaseTime(rel);

  // todo: make the hard-coded 0.1 and 10.0 factors accessible to client code, also use the hold
  // parameter of slewLimiter and make it accessible from client code
}