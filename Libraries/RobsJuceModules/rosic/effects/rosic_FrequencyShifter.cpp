// class FreqShifterHalfbandFilter:

FreqShifterHalfbandFilter::FreqShifterHalfbandFilter()
: stage1(8), stage2(8), stage3(8)
{
  static const int order = 24;

  // design the filter as biquad cascade:
  rsEngineersFilterMono halfbandFilterBiquad;
  halfbandFilterBiquad.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC);
  halfbandFilterBiquad.setSampleRate(44100.0);
  halfbandFilterBiquad.setFrequency(11025);        // place passband edge exactly at sampleRate/4
  halfbandFilterBiquad.setRipple(0.1);
  halfbandFilterBiquad.setStopbandRejection(95.0);
  halfbandFilterBiquad.setPrototypeOrder(order);   // stopband edge will be <= 11040

  // convert the biquad cascade filter into 3 8th order sections:
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
  // There is some debugging code inserted here. See comments below.


  setSampleRate(44100.0);
  shiftInHz      = 0.0;
  feedbackFactor = 0.0;
  yOld           = 0.0;


  // Design parameters for the elliptic halfband filters:
  static const int protoFilterOrder = 24;
  // 24 was originally desired but causes access violations. 23 still causes access violations but
  // the address is different. With 22, we seem to be on the safe side.
  // Increasing rsEngineersFilter::maxNumBiquads to e.g. 30 does *not* seem to have an effect on 
  // how high we can go here. That makes it even more strange.



  halfbandFilter1.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC);
  halfbandFilter1.setSampleRate(44100.0);
  halfbandFilter1.setFrequency(11025);
  halfbandFilter1.setRipple(0.1);
  halfbandFilter1.setStopbandRejection(95.0);
  halfbandFilter1.setPrototypeOrder(protoFilterOrder);  // This call sets us up for desaster later



  halfbandFilter2.initBiquadCoeffs();  // For debug - triggers the same access violation as below.
  // It seems like the call  halfbandFilter1.setPrototypeOrder(24);  immediately before the 
  // call halfbandFilter2.initBiquadCoeffs();  messes up the address of a1 in halfbandFilter1. When
  // we then try to access a1 inside halfbandFilter2.initBiquadCoeffs(), we get the access 
  // violation. Apparently, we overwrite memory that isn't ours in rsEngineersFilterMono when we 
  // try to create a filter that is close to the maximum order. Trying to increase the default 
  // maxNumStages in rsBiquadCascade by setting rsEngineersFilter::maxNumBiquads to something 
  // higher (e.g. 32 instead of 25), doesn't actually seem to help in any way. On  the other hand,
  // using  protoFilterOrder = 22;  here seems to fix it (or rather work around it). 23 will still 
  // cause problems, but the messed address of a1 will be different.



  halfbandFilter2.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC); // Access violation!!!
  halfbandFilter2.setSampleRate(44100.0);
  halfbandFilter2.setFrequency(11025);
  halfbandFilter2.setRipple(0.1);
  halfbandFilter2.setStopbandRejection(95.0);
  halfbandFilter2.setPrototypeOrder(protoFilterOrder);
  // will result in a stopband frequency <= 11040


  sinOsc1.setStartPhase(0.0);
  sinOsc2.setStartPhase(0.0);
  cosOsc1.setStartPhase(0.5*PI);
  cosOsc2.setStartPhase(0.5*PI);

  setupOscillators();

  // ToDo:
  // -Figure out and dix the bug in rsEngineersFilterMono that we run into here. When that's done,
  //  clean up the code here, i.e. remove the stuff that we have inserted for debugging.
  //  There is now a function engineersFilterUnitTest() in the TestsRosicAndRapt project which
  //  successfully triggers this behavior in a simpler context.
  // -Define constants also for the other halfband filter design parameters to avoid using 
  //  magic numbers in the setters.
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