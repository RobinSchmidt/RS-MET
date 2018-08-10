template<class TSig, class TPar>
rsFakeResonanceFilter<TSig, TPar>::rsFakeResonanceFilter()
{
  fs     = 44100;
  cutoff = 1000;

  // resonance parameters:
  shift       = 0;
  gr          = 1.0;
  pr          = PI/2;
  dr          = 0.01;
  decayByFreq = 0.0;      // constant absolute bandwidth behavior
  ar          = 0.001;
  delay       = 1.0;      // use "optimal" delay 

  lpf1.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);
  lpf2.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);
  lpf3.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);
  lpf4.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT);

  hpf.setMode(rsOnePoleFilter<TSig, TPar>::HIGHPASS_MZT);

  //dl.setInterpolationMethod(rsInterpolator<TSig>::LINEAR); 
  //dl.setInterpolationMethod(rsInterpolator<TSig>::ALLPASS); 
  dl.setInterpolationMethod(rsInterpolator<TSig>::WARPED_ALLPASS); 
    // it seems to be a good idea to use the more sophisticated (possibly warped) allpass 
    // interpolation instead of linear - otherwise, the resonance amplitude will be reduced at 
    // high resonance frequencies whenever the delay is such that the fractional part of the 
    // delay (in samples) is 0.5 - the linear interpolator has its strongest lowpass character at 
    // this value for the fractional offset

  setSampleRate(fs);
  setLowpassCutoff(cutoff);
  setHighpassCutoff(20000);


  //reset();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  fs = newSampleRate;
  lpf1.setSampleRate(fs);
  lpf2.setSampleRate(fs);
  lpf3.setSampleRate(fs);
  lpf4.setSampleRate(fs);
  hpf.setSampleRate(fs);
  dl.setSampleRate(fs);
  setupResonance();
  setupResonanceDelay();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setLowpassCutoff(TPar newCutoff)
{
  cutoff = newCutoff;
  lpf1.setCutoff(cutoff);
  lpf2.setCutoff(cutoff);
  lpf3.setCutoff(cutoff);
  lpf4.setCutoff(cutoff);
  setupResonance();
  setupResonanceDelay();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setHighpassCutoff(TPar newCutoff)
{
  hpf.setCutoff(newCutoff);
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setResonanceShift(TPar newShift)
{
  shift = newShift;
  setupResonance();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setResonanceAmplitude(TPar newAmplitude)
{
  gr = newAmplitude;
  setupResonance();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setResonancePhase(TPar newPhase)
{
  pr = newPhase;
  setupResonance();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setResonanceDecay(TPar newDecay)
{
  dr = newDecay;
  setupResonance();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setDecayByFrequency(TPar newDecayByFreq)
{
  decayByFreq = newDecayByFreq;
  setupResonance();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setResonanceAttack(TPar newAttack)
{
  ar = newAttack;
  setupResonance();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setResonanceDelay(TPar newDelay)
{
  delay = newDelay;
  setupResonanceDelay();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::reset()
{
  lpf1.reset();
  lpf2.reset();
  lpf3.reset();
  lpf4.reset();
  hpf.reset();
  rf.reset();
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setupResonance()
{
  TPar fr  = cutoff * rsPitchOffsetToFreqFactor(shift);
  TPar scl = pow(2.0, -decayByFreq*rsLog2(0.001*fr));
    // scaler for the attack/decay times
    // i think, this may be optimized so as to use only exp/log instead of pow/log2 (the latter
    // pair of which is more expensive to compute) by using some rules for exponentials and 
    // logarithms

  rf.setModalParameters(fr, gr, scl*ar*dr, scl*dr, pr, fs, 1.0);
}

template<class TSig, class TPar>
void rsFakeResonanceFilter<TSig, TPar>::setupResonanceDelay()
{
  dl.setDelayTime(delay/cutoff);
}
