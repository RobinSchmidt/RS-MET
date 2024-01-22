// class FreqShifterHalfbandFilter:

FreqShifterHalfbandFilter::FreqShifterHalfbandFilter()
: stage1(8), stage2(8), stage3(8)
{
  static const int order = 24;

  // Design the filter as biquad cascade:
  rsEngineersFilterMono halfbandFilterBiquad;
  halfbandFilterBiquad.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC);
  halfbandFilterBiquad.setSampleRate(44100.0);
  halfbandFilterBiquad.setFrequency(11025);        // place passband edge exactly at sampleRate/4
  halfbandFilterBiquad.setRipple(0.1);
  halfbandFilterBiquad.setStopbandRejection(95.0);
  halfbandFilterBiquad.setPrototypeOrder(order);   // stopband edge will be <= 11040

  // Convert the biquad cascade filter into 3 8th order sections:
  double b1[order/3+1];
  double a1[order/3+1];
  double b2[order/3+1];
  double a2[order/3+1];
  double b3[order/3+1];
  double a3[order/3+1];
  rsFilterCoefficientConverterD::biquadCascadeToDirectForm(order/6, 
    halfbandFilterBiquad.getAddressB0(),                                                     
    halfbandFilterBiquad.getAddressB1(), halfbandFilterBiquad.getAddressB2(), 
    halfbandFilterBiquad.getAddressA1(), halfbandFilterBiquad.getAddressA2(),
    b1, a1);
  rsFilterCoefficientConverterD::biquadCascadeToDirectForm(order/6, 
    halfbandFilterBiquad.getAddressB0()+order/6, 
    halfbandFilterBiquad.getAddressB1()+order/6, halfbandFilterBiquad.getAddressB2()+order/6, 
    halfbandFilterBiquad.getAddressA1()+order/6, halfbandFilterBiquad.getAddressA2()+order/6,
    b2, a2);
  rsFilterCoefficientConverterD::biquadCascadeToDirectForm(order/6, 
    halfbandFilterBiquad.getAddressB0()+2*order/6, 
    halfbandFilterBiquad.getAddressB1()+2*order/6, halfbandFilterBiquad.getAddressB2()+2*order/6, 
    halfbandFilterBiquad.getAddressA1()+2*order/6, halfbandFilterBiquad.getAddressA2()+2*order/6,
    b3, a3);
  stage1.setCoefficients(a1, b1, order/3);
  stage2.setCoefficients(a2, b2, order/3);
  stage3.setCoefficients(a3, b3, order/3);
}

//=================================================================================================
// class FrequencyShifter:

// construction/destruction:

FrequencyShifter::FrequencyShifter()
{
  setSampleRate(44100.0);
  shiftInHz      = 0.0;
  feedbackFactor = 0.0;
  yOld           = 0.0;

  // Design parameters for the two elliptic halfband filters:
  static const int    protoFilterOrder = 24;
  static const double sampleRate       = 44100;  // in Hz
  static const double cutoff           = 11025;  // in Hz
  static const double ripple           = 0.1;    // in dB
  static const double rejection        = 95.0;   // in dB
  // These will result in a stopband frequency <= 11040.

  halfbandFilter1.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC);
  halfbandFilter1.setSampleRate(sampleRate);
  halfbandFilter1.setFrequency(cutoff);
  halfbandFilter1.setRipple(ripple);
  halfbandFilter1.setStopbandRejection(rejection);
  halfbandFilter1.setPrototypeOrder(protoFilterOrder);

  halfbandFilter2.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC);
  halfbandFilter2.setSampleRate(sampleRate);
  halfbandFilter2.setFrequency(cutoff);
  halfbandFilter2.setRipple(ripple);
  halfbandFilter2.setStopbandRejection(rejection);
  halfbandFilter2.setPrototypeOrder(protoFilterOrder);
  // Setting up the order last is actually a good idea because all the other calculation-triggering
  // setter calls will only calculate coeffs for a low order filter. The default order is 2. Of 
  // course, it would be even better, if we would not trigger intermediate calculations at all. 
  // Maybe provide a setup() method to set all the parameters at once and therefore also only 
  // triggers the calculation once. Also, the 2nd filter is identical to the 1st so it would be 
  // nice, if we could just copy the settings and coeffs. Maybe make a method 
  // halfbandFilter2.copySettingsFrom(halfbandFilter1). When we do it like this, we could actually
  // get rid of the variables again. We should replace the whole filter setup stuff with 2 lines 
  // like:
  //
  // halfbandFilter1.setup(LOWPASS, 24, ELLIPTIC, 44100, 11025, 0.1, 95.0);
  // halfbandFilter2.copySettingsFrom(halfbandFilter1);
  //
  // Actually, the halfbandFilter2 doesn't even need to be of class rsEngineersFilter. It can be of
  // type rsBiquadCascade which already has a copySettingsFrom() method. That would reduce the 
  // memory footprint of the freq-shifter, too. But before implementing all these optimizations, 
  // create a unit test. Maybe shift a 1 kHz sinewave by 200 Hz and check the result - maybe by 
  // comparing it to to a 1.2 kHz sinewave. But we don't really know, what the phase of that should
  // be - but we could find it out experimentally. Or we could just look at the magnitude spectrum.
  // It might also be nice to be able avoid the design procedure in the standard constructor when 
  // we at some point want to create arrays of freq-shifters. It would be wasteful if each shifter 
  // in the array designs the filter from scratch. Maybe we should just design the filter somewhere
  // at compile time and hardcode the coeffs.

  sinOsc1.setStartPhase(0.0);
  sinOsc2.setStartPhase(0.0);
  cosOsc1.setStartPhase(0.5*PI);
  cosOsc2.setStartPhase(0.5*PI);

  setupOscillators();
}

FrequencyShifter::~FrequencyShifter()
{

}

//-------------------------------------------------------------------------------------------------
// setup:
    
void FrequencyShifter::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    mutex.lock();
    sampleRate = newSampleRate;
    sinOsc1.setSampleRate(sampleRate);
    sinOsc2.setSampleRate(sampleRate);
    cosOsc1.setSampleRate(sampleRate);
    cosOsc2.setSampleRate(sampleRate);
    mutex.unlock();
  }
}

//-------------------------------------------------------------------------------------------------
// others:
    
void FrequencyShifter::reset()
{
  yOld = 0.0;
  halfbandFilter1.reset();
  halfbandFilter2.reset();
  nyquistBlocker.reset();
}

void FrequencyShifter::setupOscillators()
{    
  mutex.lock();

  // use always the variant which doesn't show the mirrored frequencies (see article):
  if( shiftInHz >= 0.0 )
  {
    sinOsc1.setFrequency(0.25*sampleRate);
    cosOsc1.setFrequency(0.25*sampleRate);
    sinOsc2.setFrequency(0.25*sampleRate+shiftInHz);
    cosOsc2.setFrequency(0.25*sampleRate+shiftInHz);
  }
  else
  {
    sinOsc1.setFrequency(0.25*sampleRate-shiftInHz);
    cosOsc1.setFrequency(0.25*sampleRate-shiftInHz);
    sinOsc2.setFrequency(0.25*sampleRate);
    cosOsc2.setFrequency(0.25*sampleRate);
  }

  // maybe this can be dragged out into the constructor:
  sinOsc1.trigger();
  sinOsc2.trigger();
  cosOsc1.trigger();
  cosOsc2.trigger();
  
  mutex.unlock();
}

//=================================================================================================
// class FrequencyShifterStereo:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FrequencyShifterStereo::FrequencyShifterStereo()
{
  dry    = 0.0;
  wet    = 1.0;
  shift  = 0.0;
  offset = 0.0;
}

FrequencyShifterStereo::~FrequencyShifterStereo()
{

}