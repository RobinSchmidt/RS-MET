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


/*

Ideas:

-Detect intersample peaks:
 -Remember x[n-2] and x[n-1]. At time n, we have also x[n] available. I assume here that x[n] is 
  already positive, i.e. the absolute value of some input signal
 -If x[n-1] >= x[n-2] and x[n-1] >= x[n], we have detected a peak. It's located somewhere near 
  x[n-1].
 -The exact location of the peak could be anywhere between n-2 and n. Locate it by fitting a 
  parabola. xPeak is the value of the peak, frac is the fractional part of its position.
 -To compute output and update the smoothing filter states, we need to do 2 fractional filter steps 
  (using non-uniform filtering) and one normal step like (I think):
 -If the peak is between n-2 and n-1:
     smoother.setTimeConstant(attack);              //
     out1 = smoother.getSample(xPeak, fracPos);     // 1st fractional step x[n-2] is currently in the filter's state
     smoother.setTimeConstant(release);
     out2 = smoother.getSample(x[n-1], 1-fracPos);  // 2nd fractional step
     out3 = smoother.getSample(x[n]);               // this is a full-sample step
 -Else (peak is between n-1 and n): 
     smoother.setTimeConstant(release);             // may not be needed?
     out1 = smoother.getSample(x[n-1])              // full sample step
     smoother.setTimeConstant(attack);
     out2 = smoother.getSample(xPeak, fracPos);     // 1st fractional step
     smoother.setTimeConstant(release);
     out2 = smoother.getSample(x[n], 1-fracPos);    // 2nd fractional step
 -The scheme introduces a delay of 2 samples, i.e our output at time n applies to the input at n-2
  I think. Or is it n-1? 

-Try to produce release shapes other than the 1-pole exponential decay. For example, provide a 
 linear decay, maybe by using rsSlewRateLimiterLinear instead of rsSlewRateLimiter. But maybe 
 implement a class rsSlewRateLimiterFlexible that has a flexible shape. It could have members:
   user params: attackSamples, releaseSamples, shape
   algo params: coeffAttack, coeffRelease
 where the coeffs may mean different things depending on shape. Behaviors:
   exponential:  y = y * c
   linear:       y = max(0, y - c)
   gauss:        y = c * y * min(y, 1)
 where "gauss" behaves like exp(-x^2), I think. We multiply the output by the coeff c < 1 but also
 by the output y itself whenever y is less than 1 (this is needed to prevent it from blowing up).
 All of them can be combined with additional smoothing - maybe up to 4th order or something. This 
 might be especially nice for the linear shape. The cutoff of the smoothing filter should be high 
 enough to not become the factor that determines the release time.


*/