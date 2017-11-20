using namespace RSLib;

rsFakeResonanceFilter::rsFakeResonanceFilter()
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

  lpf1.setMode(rsOnePoleFilter::LOWPASS);
  lpf2.setMode(rsOnePoleFilter::LOWPASS);
  lpf3.setMode(rsOnePoleFilter::LOWPASS);
  lpf4.setMode(rsOnePoleFilter::LOWPASS);

  hpf.setMode(rsOnePoleFilter::HIGHPASS);

  //dl.setInterpolationMethod(rsInterpolator::LINEAR); 
  //dl.setInterpolationMethod(rsInterpolator::ALLPASS); 
  dl.setInterpolationMethod(rsInterpolator::WARPED_ALLPASS); 
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

void rsFakeResonanceFilter::setSampleRate(double newSampleRate)
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

void rsFakeResonanceFilter::setLowpassCutoff(double newCutoff)
{
  cutoff = newCutoff;
  lpf1.setCutoff(cutoff);
  lpf2.setCutoff(cutoff);
  lpf3.setCutoff(cutoff);
  lpf4.setCutoff(cutoff);
  setupResonance();
  setupResonanceDelay();
}

void rsFakeResonanceFilter::setHighpassCutoff(double newCutoff)
{
  hpf.setCutoff(newCutoff);
}

void rsFakeResonanceFilter::setResonanceShift(double newShift)
{
  shift = newShift;
  setupResonance();
}

void rsFakeResonanceFilter::setResonanceAmplitude(double newAmplitude)
{
  gr = newAmplitude;
  setupResonance();
}

void rsFakeResonanceFilter::setResonancePhase(double newPhase)
{
  pr = newPhase;
  setupResonance();
}

void rsFakeResonanceFilter::setResonanceDecay(double newDecay)
{
  dr = newDecay;
  setupResonance();
}

void rsFakeResonanceFilter::setDecayByFrequency(double newDecayByFreq)
{
  decayByFreq = newDecayByFreq;
  setupResonance();
}

void rsFakeResonanceFilter::setResonanceAttack(double newAttack)
{
  ar = newAttack;
  setupResonance();
}

void rsFakeResonanceFilter::setResonanceDelay(double newDelay)
{
  delay = newDelay;
  setupResonanceDelay();
}

void rsFakeResonanceFilter::reset()
{
  lpf1.reset();
  lpf2.reset();
  lpf3.reset();
  lpf4.reset();
  hpf.reset();
  rf.reset();
}

void rsFakeResonanceFilter::setupResonance()
{
  double fr  = cutoff * rsPitchOffsetToFreqFactor(shift);

  double scl = pow(2.0, -decayByFreq*rsLog2(0.001*fr));
    // scaler for the attack/decay times
    // i think, this may be optimized so as to use only exp/log instead of pow/log2 (the latter
    // pair of which is more expensive to compute) by using some rules for exponentials and 
    // logarithms

  rf.setModalParameters(fr, gr, scl*ar*dr, scl*dr, pr, fs, 1.0);
}

void rsFakeResonanceFilter::setupResonanceDelay()
{
  dl.setDelayTime(delay/cutoff);
}
