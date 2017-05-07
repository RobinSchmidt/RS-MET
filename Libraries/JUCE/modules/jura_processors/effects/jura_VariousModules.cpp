#include "rosof_VariousModulesAndEditors.h"
using namespace rosof;


//=================================================================================================
// Distortion Effects:

//-------------------------------------------------------------------------------------------------
// BitCrusher:

BitCrusherAudioModule::BitCrusherAudioModule(CriticalSection *newPlugInLock, rosic::BitCrusher *newBitCrusherToWrap) 
: AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);
  jassert( newBitCrusherToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedBitCrusher = newBitCrusherToWrap;
  moduleName = juce::String(T("BitCrusher"));

  // maybe these 2 calls can be absorbed into 1 initialize() call or something:
  //setActiveDirectory(getApplicationDirectory() + File::separatorString + juce::String(T("BitCrusherPresets")) );
  setActiveDirectory(getApplicationDirectory() + File::separatorString + moduleName + juce::String(T("Presets")) );
  createStaticParameters();
}

void BitCrusherAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Decimation", 1.0, 128.0, 1.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedBitCrusher, &BitCrusher::setDecimationFactor);

  p = new AutomatableParameter(plugInLock, "Quantization", 0.001, 1.0, 0.0, 0.00001, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedBitCrusher, &BitCrusher::setQuantizationInterval);

  p = new AutomatableParameter(plugInLock, "Amount", -200.0, 200.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedBitCrusher, &BitCrusher::setAmount); 

  for(int i=0; i < (int) observedParameters.size(); i++ )
    observedParameters[i]->resetToDefaultValue(true, true);
}

BitCrusherModuleEditor::BitCrusherModuleEditor(CriticalSection *newPlugInLock, BitCrusherAudioModule* newBitCrusherAudioModule) 
: AudioModuleEditor(newPlugInLock, newBitCrusherAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newBitCrusherAudioModule != NULL ); // you must pass a valid module here

  addWidget( decimationSlider = new RSlider (T("DecimationSlider")) );
  decimationSlider->assignParameter( moduleToEdit->getParameterByName("Decimation") );
  decimationSlider->setDescription(juce::String(T("Decimation factor for the sample-rate")));
  decimationSlider->setDescriptionField(infoField);
  decimationSlider->setStringConversionFunction(&rojue::valueToString0);

  addWidget( quantizationSlider = new RSlider (T("QuantizationSlider")) );
  quantizationSlider->assignParameter( moduleToEdit->getParameterByName("Quantization") );
  quantizationSlider->setDescription(juce::String(T("Quantization interval for the amplitude")));
  quantizationSlider->setDescriptionField(infoField);
  quantizationSlider->setStringConversionFunction(&rojue::valueToString4);

  addWidget( amountSlider = new RSlider (T("AmountSlider")) );
  amountSlider->assignParameter( moduleToEdit->getParameterByName("Amount") );
  amountSlider->setDescription(juce::String(T("Amount of the effect in percent")));
  amountSlider->setDescriptionField(infoField);
  amountSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  updateWidgetsAccordingToState();
}

void BitCrusherModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  decimationSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  quantizationSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  amountSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// ModulatedAllpass:

ModulatedAllpassAudioModule::ModulatedAllpassAudioModule(CriticalSection *newPlugInLock, 
                                                         rosic::ModulatedAllpass *newModulatedAllpassToWrap)
                                                          : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);
  jassert( newModulatedAllpassToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedModulatedAllpass = newModulatedAllpassToWrap;
  moduleName  = juce::String(T("ModulatedAllpass"));

  setActiveDirectory(getApplicationDirectory() + File::separatorString + juce::String(T("ModulatedAllpassPresets")) );
  createStaticParameters();
}

void ModulatedAllpassAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Factor", -10.0, 10.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedModulatedAllpass, &ModulatedAllpass::setFactor);

  p = new AutomatableParameter(plugInLock, "Offset", -1.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedModulatedAllpass, &ModulatedAllpass::setOffset);

  for(int i=0; i < (int) observedParameters.size(); i++ )
    observedParameters[i]->resetToDefaultValue(true, true);
}

ModulatedAllpassModuleEditor::ModulatedAllpassModuleEditor(CriticalSection *newPlugInLock, 
                                                           ModulatedAllpassAudioModule* newModulatedAllpassAudioModule) 
: AudioModuleEditor(newPlugInLock, newModulatedAllpassAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newModulatedAllpassAudioModule != NULL ); // you must pass a valid module here

  addWidget( factorSlider = new RSlider (T("FactorSlider")) );
  factorSlider->assignParameter( moduleToEdit->getParameterByName("Factor") );
  factorSlider->setDescription(juce::String(T("Factor for the modulating signal")));
  factorSlider->setDescriptionField(infoField);
  factorSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( offsetSlider = new RSlider (T("OffsetSlider")) );
  offsetSlider->assignParameter( moduleToEdit->getParameterByName("Offset") );
  offsetSlider->setDescription(juce::String(T("Offset for the modulating signal")));
  offsetSlider->setDescriptionField(infoField);
  offsetSlider->setStringConversionFunction(&rojue::valueToString2);

  updateWidgetsAccordingToState();
}

void ModulatedAllpassModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  factorSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  offsetSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// SlewRateLimiter:

SlewRateLimiterAudioModule::SlewRateLimiterAudioModule(CriticalSection *newPlugInLock, 
                                                       rosic::SlewRateLimiterStereo *newSlewRateLimiterToWrap)
                                                        : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newSlewRateLimiterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedSlewRateLimiter = newSlewRateLimiterToWrap;
  moduleName  = juce::String(T("SlewRateLimiter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("SlewRateLimiterPresets")) );
  createStaticParameters();
}

void SlewRateLimiterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Attack",  0.0, 10.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSlewRateLimiter, &SlewRateLimiterStereo::setAttackTime);

  p = new AutomatableParameter(plugInLock, "Release", 0.0, 10.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSlewRateLimiter, &SlewRateLimiterStereo::setReleaseTime);

  for(int i=0; i < (int) observedParameters.size(); i++ )
    observedParameters[i]->resetToDefaultValue(true, true);
}

SlewRateLimiterModuleEditor::SlewRateLimiterModuleEditor(CriticalSection *newPlugInLock, 
                                                         SlewRateLimiterAudioModule* newSlewRateLimiterAudioModule) 
: AudioModuleEditor(newPlugInLock, newSlewRateLimiterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newSlewRateLimiterAudioModule != NULL ); // you must pass a valid module here

  addWidget( attackSlider = new RSlider (T("AttackSlider")) );
  attackSlider->assignParameter( moduleToEdit->getParameterByName("Attack") );
  attackSlider->setDescription(juce::String(T("Slew rate for upward jumps")));
  attackSlider->setDescriptionField(infoField);
  attackSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( releaseSlider = new RSlider (T("ReleaseSlider")) );
  releaseSlider->assignParameter( moduleToEdit->getParameterByName("Release") );
  releaseSlider->setDescription(juce::String(T("Slew rate for dwonward jumps")));
  releaseSlider->setDescriptionField(infoField);
  releaseSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  updateWidgetsAccordingToState();
}

void SlewRateLimiterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  attackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  releaseSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// Harmonics:

HarmonicsAudioModule::HarmonicsAudioModule(CriticalSection *newPlugInLock, rosic::Harmonics *newHarmonicsToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newHarmonicsToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedHarmonics = newHarmonicsToWrap;
  moduleName  = juce::String(T("Harmonics"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("HarmonicsPresets")) );
  createStaticParameters();
}

void HarmonicsAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);

  // \todo get rid of this function and implement the new callback mechanism - maybe we need a class IndexedParameter or something

  if( wrappedHarmonics == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  //case  0: wrappedHarmonics->setDrive(                       value);   break;
  //case  1: wrappedHarmonics->setDryWetRatio(                 value);   break;
  //case  2: wrappedHarmonics->setInputHighpassCutoff(         value);   break;
  //case  3: wrappedHarmonics->setInputLowpassCutoff(          value);   break;
  //case  4: wrappedHarmonics->setOutputHighpassCutoff(        value);   break;
  //case  5: wrappedHarmonics->setOutputLowpassCutoff(         value);   break;
  case  6: wrappedHarmonics->setHarmonicGainFactor( 2,  0.01*value);   break;
  case  7: wrappedHarmonics->setHarmonicGainFactor( 3,  0.01*value);   break;
  case  8: wrappedHarmonics->setHarmonicGainFactor( 4,  0.01*value);   break;
  case  9: wrappedHarmonics->setHarmonicGainFactor( 5,  0.01*value);   break;
  case 10: wrappedHarmonics->setHarmonicGainFactor( 6,  0.01*value);   break;
  case 11: wrappedHarmonics->setHarmonicGainFactor( 7,  0.01*value);   break;
  case 12: wrappedHarmonics->setHarmonicGainFactor( 8,  0.01*value);   break;
  case 13: wrappedHarmonics->setHarmonicGainFactor( 9,  0.01*value);   break;
  case 14: wrappedHarmonics->setHarmonicGainFactor(10,  0.01*value);   break;
  case 15: wrappedHarmonics->setHarmonicGainFactor(11,  0.01*value);   break;
  case 16: wrappedHarmonics->setHarmonicGainFactor(12,  0.01*value);   break;
  } 
}

void HarmonicsAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  // ...here we either need a callback with two parameters, or we must stick to the old callback mechansim - for now, we'll do the latter

  p = new AutomatableParameter(plugInLock, juce::String(T("Drive")), -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedHarmonics, &Harmonics::setDrive);

  p = new AutomatableParameter(plugInLock, juce::String(T("DryWetRatio")), 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedHarmonics, &Harmonics::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, juce::String(T("InputHighpass")), 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedHarmonics, &Harmonics::setInputHighpassCutoff);

  p = new AutomatableParameter(plugInLock, juce::String(T("InputLowpass")), 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedHarmonics, &Harmonics::setInputLowpassCutoff);

  p = new AutomatableParameter(plugInLock, juce::String(T("OutputHighpass")), 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedHarmonics, &Harmonics::setOutputHighpassCutoff);

  p = new AutomatableParameter(plugInLock, juce::String(T("OutputLowpass")), 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedHarmonics, &Harmonics::setOutputLowpassCutoff);

  // these parameters are still used with the old callback mechanism (because they invoke two-parametric functions) - we may do something
  // better here later:
  p = new AutomatableParameter(plugInLock, juce::String(T("H02")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H03")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H04")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H05")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H06")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H07")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H08")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H09")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H10")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H11")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p = new AutomatableParameter(plugInLock, juce::String(T("H12")), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  for(int i=0; i < (int) observedParameters.size(); i++ )
  {   
    parameterChanged(observedParameters[i]); // because some parameters use the old callback mechanism
    observedParameters[i]->resetToDefaultValue(true, true);
  }
}

HarmonicsModuleEditor::HarmonicsModuleEditor(CriticalSection *newPlugInLock, HarmonicsAudioModule* newHarmonicsAudioModule) 
: AudioModuleEditor(newPlugInLock, newHarmonicsAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newHarmonicsAudioModule != NULL ); // you must pass a valid module here

  addWidget( globalLabel = new RTextField( juce::String(T("Global"))) );
  globalLabel->setJustification(Justification::centred);
  globalLabel->setDescription(T("Global parameters"));
  globalLabel->setDescriptionField(infoField);

  addWidget( harmonicsLabel = new RTextField( juce::String(T("Harmonics"))) );
  harmonicsLabel->setJustification(Justification::centred);
  harmonicsLabel->setDescription(T("Gains for the individual harmonics"));
  harmonicsLabel->setDescriptionField(infoField);

  addWidget( inFilterLabel = new RTextField( juce::String(T("Input Filter:"))) );
  inFilterLabel->setJustification(Justification::centredLeft);
  inFilterLabel->setDescription(T("Parameters of the filter at the input stage"));
  inFilterLabel->setDescriptionField(infoField);

  addWidget( outFilterLabel = new RTextField( juce::String(T("Output Filter:"))) );
  outFilterLabel->setJustification(Justification::centredLeft);
  outFilterLabel->setDescription(T("Parameters of the filter at the output stage"));
  outFilterLabel->setDescriptionField(infoField);

  addWidget( driveSlider = new RSlider (T("InLevelSlider")) );
  driveSlider->assignParameter( moduleToEdit->getParameterByName("Drive") );
  driveSlider->setDescription(juce::String(T("Gain for the input signal (pre waveshaper)")));
  driveSlider->setDescriptionField(infoField);
  driveSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);


  addWidget( inHighpassSlider = new RSlider (T("InHighpassSlider")) );
  inHighpassSlider->assignParameter( moduleToEdit->getParameterByName("InputHighpass") );
  inHighpassSlider->setSliderName(juce::String(T("Highpass")));
  inHighpassSlider->setDescription(juce::String(T("Cutoff frequency of the input highpass filter")));
  inHighpassSlider->setDescriptionField(infoField);
  inHighpassSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( inLowpassSlider = new RSlider (T("InLowpassSlider")) );
  inLowpassSlider->assignParameter( moduleToEdit->getParameterByName("InputLowpass") );
  inLowpassSlider->setSliderName(juce::String(T("Lowpass")));
  inLowpassSlider->setDescription(juce::String(T("Cutoff frequency of the input lowpass filter")));
  inLowpassSlider->setDescriptionField(infoField);
  inLowpassSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( outHighpassSlider = new RSlider (T("OutHighpassSlider")) );
  outHighpassSlider->assignParameter( moduleToEdit->getParameterByName("OutputHighpass") );
  outHighpassSlider->setSliderName(juce::String(T("Highpass")));
  outHighpassSlider->setDescription(juce::String(T("Cutoff frequency of the output highpass filter")));
  outHighpassSlider->setDescriptionField(infoField);
  outHighpassSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( outLowpassSlider = new RSlider (T("OutLowpassSlider")) );
  outLowpassSlider->assignParameter( moduleToEdit->getParameterByName("OutputLowpass") );
  outLowpassSlider->setSliderName(juce::String(T("Lowpass")));
  outLowpassSlider->setDescription(juce::String(T("Cutoff frequency of the output lowpass filter")));
  outLowpassSlider->setDescriptionField(infoField);
  outLowpassSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( h02Slider = new RSlider (T("H02Slider")) );
  h02Slider->assignParameter( moduleToEdit->getParameterByName("H02") );
  h02Slider->setDescription(juce::String(T("Gain (in percent) for the 2nd harmonic")));
  h02Slider->setDescriptionField(infoField);
  h02Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h03Slider = new RSlider (T("H03Slider")) );
  h03Slider->assignParameter( moduleToEdit->getParameterByName("H03") );
  h03Slider->setDescription(juce::String(T("Gain (in percent) for the 3rd harmonic")));
  h03Slider->setDescriptionField(infoField);
  h03Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h04Slider = new RSlider (T("H04Slider")) );
  h04Slider->assignParameter( moduleToEdit->getParameterByName("H04") );
  h04Slider->setDescription(juce::String(T("Gain (in percent) for the 4th harmonic")));
  h04Slider->setDescriptionField(infoField);
  h04Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h05Slider = new RSlider (T("H05Slider")) );
  h05Slider->assignParameter( moduleToEdit->getParameterByName("H05") );
  h05Slider->setDescription(juce::String(T("Gain (in percent) for the 5th harmonic")));
  h05Slider->setDescriptionField(infoField);
  h05Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h06Slider = new RSlider (T("H06Slider")) );
  h06Slider->assignParameter( moduleToEdit->getParameterByName("H06") );
  h06Slider->setDescription(juce::String(T("Gain (in percent) for the 6th harmonic")));
  h06Slider->setDescriptionField(infoField);
  h06Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h07Slider = new RSlider (T("H07Slider")) );
  h07Slider->assignParameter( moduleToEdit->getParameterByName("H07") );
  h07Slider->setDescription(juce::String(T("Gain (in percent) for the 7th harmonic")));
  h07Slider->setDescriptionField(infoField);
  h07Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h08Slider = new RSlider (T("H08Slider")) );
  h08Slider->assignParameter( moduleToEdit->getParameterByName("H08") );
  h08Slider->setDescription(juce::String(T("Gain (in percent) for the 8th harmonic")));
  h08Slider->setDescriptionField(infoField);
  h08Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h09Slider = new RSlider (T("H09Slider")) );
  h09Slider->assignParameter( moduleToEdit->getParameterByName("H09") );
  h09Slider->setDescription(juce::String(T("Gain (in percent) for the 9th harmonic")));
  h09Slider->setDescriptionField(infoField);
  h09Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h10Slider = new RSlider (T("H10Slider")) );
  h10Slider->assignParameter( moduleToEdit->getParameterByName("H10") );
  h10Slider->setDescription(juce::String(T("Gain (in percent) for the 10th harmonic")));
  h10Slider->setDescriptionField(infoField);
  h10Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h11Slider = new RSlider (T("H11Slider")) );
  h11Slider->assignParameter( moduleToEdit->getParameterByName("H11") );
  h11Slider->setDescription(juce::String(T("Gain (in percent) for the 11th harmonic")));
  h11Slider->setDescriptionField(infoField);
  h11Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( h12Slider = new RSlider (T("H12Slider")) );
  h12Slider->assignParameter( moduleToEdit->getParameterByName("H12") );
  h12Slider->setDescription(juce::String(T("Gain (in percent) for the 12th harmonic")));
  h12Slider->setDescriptionField(infoField);
  h12Slider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  updateWidgetsAccordingToState();
}

void HarmonicsModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth()/2;
  int h = getHeight()-y;  // usable height


  globalRect.setBounds(x, y, w, h);
  harmonicsRect.setBounds(x+w, y, w, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(globalRect);
  guiLayoutRectangles.add(harmonicsRect);

  x = globalRect.getX();
  y = globalRect.getY();
  w = globalRect.getWidth();
  h = globalRect.getHeight();
  globalLabel->setBounds(x+4, y+2, w-8, 16);
  y = globalLabel->getBottom(); 
  driveSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  inFilterLabel->setBounds(x+4, y+4, w-8, 16);
  y += 16;  
  inHighpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  inLowpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  outFilterLabel->setBounds(x+4, y+4, w-8, 16);
  y += 16;  
  outHighpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  outLowpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 

  x = harmonicsRect.getX();
  y = harmonicsRect.getY();
  w = harmonicsRect.getWidth();
  h = harmonicsRect.getHeight();
  harmonicsLabel->setBounds(x+4, y+2, w-8, 16);
  y = harmonicsLabel->getBottom()-4; 
  h02Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h03Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h04Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h05Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h06Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h07Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h08Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h09Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h10Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h11Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  h12Slider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  

}

//-------------------------------------------------------------------------------------------------
// WaveShaper:

WaveShaperAudioModule::WaveShaperAudioModule(CriticalSection *newPlugInLock, rosic::WaveShaper *newWaveShaperToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newWaveShaperToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedWaveShaper = newWaveShaperToWrap;
  moduleName  = juce::String(T("WaveShaper"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("WaveShaperPresets")) );
  createStaticParameters();
}

void WaveShaperAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, juce::String(T("CurveType")), 0.0, 3.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Linear")));
  p->addStringValue(juce::String(T("Tanh")));
  p->addStringValue(juce::String(T("Hardclip")));
  p->addStringValue(juce::String(T("Quintic")));
  //p->addStringValue(juce::String(T("Cubic")));
  //p->addStringValue(juce::String(T("Heptic")));
  //p->addStringValue(juce::String(T("Sin")));
  //p->setValue(1.0, false);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setTransferFunction);

  p = new AutomatableParameter(plugInLock, "Drive", -24.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setDrive);

  p = new AutomatableParameter(plugInLock, "DC", -5.0, 5.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setDcOffset);

  p = new AutomatableParameter(plugInLock, "Amount", -200.0, 200.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setAmount);

  p = new AutomatableParameter(plugInLock, "OutputLevel", -48.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setOutputLevel);

  p = new AutomatableParameter(plugInLock, "Oversampling", 1.0, 16.0, 1.0, 4.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setOversampling);

  p = new AutomatableParameter(plugInLock, "InterceptY", -0.5, 0.5, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setPenticInterceptY);

  p = new AutomatableParameter(plugInLock, "Slope", -1.0, 3.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWaveShaper, &WaveShaper::setPenticSlopeAtZero);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

WaveShaperModuleEditor::WaveShaperModuleEditor(CriticalSection *newPlugInLock, WaveShaperAudioModule* newWaveShaperAudioModule) 
: AudioModuleEditor(newPlugInLock, newWaveShaperAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newWaveShaperAudioModule != NULL ); // you must pass a valid module here
  waveShaperModuleToEdit = newWaveShaperAudioModule;

  addWidget( curveComboBox = new RComboBox(juce::String(T("CurveComboBox"))) );
  curveComboBox->assignParameter( moduleToEdit->getParameterByName(T("CurveType")) );
  curveComboBox->setDescription(T("Choose the curve-shape of the transfer function"));
  curveComboBox->setDescriptionField(infoField);
  curveComboBox->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( driveSlider = new RSlider (T("DriveSlider")) );
  driveSlider->assignParameter( moduleToEdit->getParameterByName("Drive") );
  driveSlider->setDescription(juce::String(T("Gain of the input signal before the waveshaper")));
  driveSlider->setDescriptionField(infoField);
  driveSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);
  driveSlider->addListener(this); // to update the plot

  addWidget( dcSlider = new RSlider (T("DCSlider")) );
  dcSlider->assignParameter( moduleToEdit->getParameterByName("DC") );
  dcSlider->setDescription(juce::String(T("DC offset for the input signal (after drive, before waveshaper)")));
  dcSlider->setDescriptionField(infoField);
  dcSlider->setStringConversionFunction(&rojue::valueToString2);
  dcSlider->addListener(this); // to update the plot

  addWidget( amountSlider = new RSlider (T("AmountSlider")) );
  amountSlider->assignParameter( moduleToEdit->getParameterByName("Amount") );
  amountSlider->setDescription(juce::String(T("Amount of the effect in percent")));
  amountSlider->setDescriptionField(infoField);
  amountSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);
  amountSlider->addListener(this); // to update the plot

  addWidget( outputLevelSlider = new RSlider (T("OutputLevelSlider")) );
  outputLevelSlider->assignParameter( moduleToEdit->getParameterByName("OutputLevel") );
  outputLevelSlider->setSliderName(juce::String(T("Level")));
  outputLevelSlider->setDescription(juce::String(T("Overall level of the output signal")));
  outputLevelSlider->setDescriptionField(infoField);
  outputLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( oversamplingSlider = new RSlider (T("OversamplingSlider")) );
  oversamplingSlider->assignParameter( moduleToEdit->getParameterByName("Oversampling") );
  oversamplingSlider->setDescription(juce::String(T("Oversampling factor (to avoid aliasing)")));
  oversamplingSlider->setDescriptionField(infoField);
  oversamplingSlider->setStringConversionFunction(&rojue::valueToString0);

  addWidget( slopeSlider = new RSlider (T("SlopeSlider")) );
  slopeSlider->assignParameter( moduleToEdit->getParameterByName("Slope") );
  slopeSlider->setDescription(juce::String(T("Slope of the transfer function at point x=0")));
  slopeSlider->setDescriptionField(infoField);
  slopeSlider->setStringConversionFunction(&rojue::valueToString2);
  slopeSlider->addListener(this); // to update the plot

  addWidget( interceptSlider = new RSlider (T("InterceptSlider")) );
  interceptSlider->assignParameter( moduleToEdit->getParameterByName("InterceptY") );
  interceptSlider->setSliderName(juce::String(T("y-Intercept")));
  interceptSlider->setDescription(juce::String(T("Interception of the y-axis (function value at x=0)")));
  interceptSlider->setDescriptionField(infoField);
  interceptSlider->setStringConversionFunction(&rojue::valueToString2);
  interceptSlider->addListener(this); // to update the plot

  numValues = 0;
  xValues   = NULL;
  yValues   = NULL;
  plot = new CurveFamilyPlotOld(juce::String(T("Plot")));
  plot->setDescription(juce::String(T("Shows the input-output transfer function")));
  plot->setAxisLabels(juce::String(T("")), juce::String(T("")));
  plot->setCurrentRange(-1.3, 1.3, -1.3, 1.3);
  plot->setVerticalCoarseGrid(1.0, true);
  plot->setHorizontalCoarseGrid(1.0, true);
  addPlot(plot);

  updateWidgetsAccordingToState();
}

WaveShaperModuleEditor::~WaveShaperModuleEditor()
{
  ScopedLock scopedLock(*plugInLock);
  if( xValues != NULL ) { delete[] xValues; xValues = NULL; }
  if( yValues != NULL ) { delete[] yValues; yValues = NULL; }
}

void WaveShaperModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = 2*getWidth()/5;
  int h = getHeight();
  y = getPresetSectionBottom();
  curveComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 24; 

  driveSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  dcSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  amountSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  outputLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y += 24;  

  slopeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  interceptSlider->setBounds(x+4, y+4, w-8, 16);
  //y += 28;  

  y = getHeight() - 16 - 8;
  oversamplingSlider->setBounds(x+4, y, w-8, 16);

  x = w;
  w = getWidth()-x;
  y = getPresetSectionBottom();
  plot->setBounds(x+4, y+4, w-8, w-8);


  numValues = plot->getWidth();
  if( xValues != NULL ) { delete[] xValues; xValues = NULL; }
  if( yValues != NULL ) { delete[] yValues; yValues = NULL; }
  xValues = new double[numValues];
  yValues = new double[numValues];
  updatePlot();
}

void WaveShaperModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updateWidgetEnablement();
  updatePlot();
}

void WaveShaperModuleEditor::rSliderValueChanged(rojue::RSlider *rSliderThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updatePlot();
}

void WaveShaperModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::updateWidgetsAccordingToState();
  updatePlot();
}

void WaveShaperModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);
  if( waveShaperModuleToEdit == NULL )
    return;
  if( waveShaperModuleToEdit->wrappedWaveShaper == NULL )
    return;
  rosic::WaveShaper* core = waveShaperModuleToEdit->wrappedWaveShaper;

  interceptSlider->setEnabled(core->doesTransferFunctionSupportInterceptY());
  slopeSlider->setEnabled(core->doesTransferFunctionSupportSlope());
}

void WaveShaperModuleEditor::updatePlot()
{
  ScopedLock scopedLock(*plugInLock);
  if( waveShaperModuleToEdit == NULL )
    return;
  if( waveShaperModuleToEdit->wrappedWaveShaper == NULL )
    return;
  rosic::WaveShaper* core = waveShaperModuleToEdit->wrappedWaveShaper;

  double xMin  = plot->getCurrentRangeMinX();
  double xMax  = plot->getCurrentRangeMaxX();
  double xStep = (xMax-xMin) / (numValues-1);
  for(int i=0; i<numValues; i++)
  {
    xValues[i] = xMin + i*xStep;
    yValues[i] = waveShaperModuleToEdit->wrappedWaveShaper->transferFunctionAt(xValues[i]);
  }
  plot->setCurveValues(numValues, xValues, yValues);
}

//-------------------------------------------------------------------------------------------------
// CompShaper:

CompShaperAudioModule::CompShaperAudioModule(CriticalSection *newPlugInLock, rosic::CompShaper *newCompShaperToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);
  jassert( newCompShaperToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedCompShaper = newCompShaperToWrap;
  moduleName  = juce::String(T("CompShaper"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("CompShaperPresets")) );
  createStaticParameters();
}

void CompShaperAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Threshold", -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCompShaper, &CompShaper::setThreshold);

  p = new AutomatableParameter(plugInLock, "Ratio", 1.0, 40.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCompShaper, &CompShaper::setRatio);

  p = new AutomatableParameter(plugInLock, "KneeWidth", 0.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCompShaper, &CompShaper::setKneeWidth);

  p = new AutomatableParameter(plugInLock, "Clip", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedCompShaper, &CompShaper::setToClipperMode);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

CompShaperModuleEditor::CompShaperModuleEditor(CriticalSection *newPlugInLock, CompShaperAudioModule* newCompShaperAudioModule) 
: AudioModuleEditor(newPlugInLock, newCompShaperAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newCompShaperAudioModule != NULL ); // you must pass a valid module here

  addWidget( curveLabel = new RTextField( juce::String(T("Transfer Curve"))) );
  curveLabel->setJustification(Justification::centred);
  curveLabel->setDescription(T("Parameters for the static transfer curve"));
  curveLabel->setDescriptionField(infoField);

  addWidget( timeLabel = new RTextField( juce::String(T("Level Detection"))) );
  timeLabel->setJustification(Justification::centred);
  timeLabel->setDescription(T("Parameters for the level detector"));
  timeLabel->setDescriptionField(infoField);

  addWidget( othersLabel = new RTextField( juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);

  addWidget( thresholdSlider = new RSlider (T("ThresholdSlider")) );
  thresholdSlider->assignParameter( moduleToEdit->getParameterByName("Threshold") );
  thresholdSlider->setDescription(juce::String(T("Threshold above which the signal will be attenuated")));
  thresholdSlider->setDescriptionField(infoField);
  thresholdSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( ratioSlider = new RSlider (T("RatioSlider")) );
  ratioSlider->assignParameter( moduleToEdit->getParameterByName("Ratio") );
  ratioSlider->setDescription(juce::String(T("Ratio .....find a good short explanation")));
  ratioSlider->setDescriptionField(infoField);
  ratioSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( kneeSlider = new RSlider (T("KneeSlider")) );
  kneeSlider->assignParameter( moduleToEdit->getParameterByName("KneeWidth") );
  kneeSlider->setDescription(juce::String(T("Transition width between the two slopes")));
  kneeSlider->setDescriptionField(infoField);
  kneeSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( clipButton = new RButton(juce::String(T("Clip"))) );
  clipButton->assignParameter( moduleToEdit->getParameterByName(T("Clip")) );
  clipButton->setDescription(juce::String(T("Clip signal at the threshold (infinite ratio)")));
  clipButton->setDescriptionField(infoField);
  clipButton->setClickingTogglesState(true);

  addWidget( driveSlider = new RSlider (T("InLevelSlider")) );
  driveSlider->assignParameter( moduleToEdit->getParameterByName("Drive") );
  driveSlider->setDescription(juce::String(T("Gain for the input signal (pre waveshaper)")));
  driveSlider->setDescriptionField(infoField);
  driveSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( outLevelSlider = new RSlider (T("OutLevelSlider")) );
  outLevelSlider->assignParameter( moduleToEdit->getParameterByName("OutLevel") );
  outLevelSlider->setDescription(juce::String(T("Gain for the output signal (post waveshaper)")));
  outLevelSlider->setDescriptionField(infoField);
  outLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  updateWidgetsAccordingToState();
}

void CompShaperModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight()-y;  // usable height

  int h3 = h/3;
  curveParametersRect.setBounds(x, y, w, h3);
  y += h3;
  timeParametersRect.setBounds(x, y, w, h3);
  y += h3;
  otherParametersRect.setBounds(x, y, w, h3);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(curveParametersRect);
  guiLayoutRectangles.add(timeParametersRect); 
  guiLayoutRectangles.add(otherParametersRect);

  x = curveParametersRect.getX();
  y = curveParametersRect.getY();
  w = curveParametersRect.getWidth();
  h = curveParametersRect.getHeight();
  curveLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  thresholdSlider->setBounds(x+4, y+4, 2*w/3-8, 16);
  y += 20;  
  ratioSlider->setBounds(x+4, y+4, 2*w/3-8, 16);
  y -= 20;  
  x = thresholdSlider->getRight();
  w = w-x;
  kneeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  clipButton->setBounds(x+4+16, y+4, w-8-32, 16);

  /*
  x = timeParametersRect.getX();
  y = timeParametersRect.getY();
  w = timeParametersRect.getWidth();
  h = timeParametersRect.getHeight();
  timeLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  attackSlider->setBounds(x+4, y+4, w/2-8, 16);
  releaseSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20; 
  x  = w/2 - w/4;
  lookAheadSlider->setBounds(x+4, y+4, w/2-8, 16);
  */

  x = otherParametersRect.getX();
  y = otherParametersRect.getY();
  w = otherParametersRect.getWidth();
  h = otherParametersRect.getHeight();
  othersLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  driveSlider->setBounds(x+4, y+4, w/2-8, 16);
  outLevelSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20; 
  //amountSlider->setBounds(x+4, y+4, w/2-8, 16);
  //autoGainButton->setBounds(x+w/2+4+32, y+4, w/2-8-64, 16);
}

//=================================================================================================
// Dynamics:

//-------------------------------------------------------------------------------------------------
// Compressor:

CompressorAudioModule::CompressorAudioModule(CriticalSection *newPlugInLock, rosic::SoftKneeCompressor *newCompressorToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newCompressorToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedCompressor = newCompressorToWrap;
  moduleName  = juce::String(T("Compressor"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("CompressorPresets")) );
  createStaticParameters();
}

void CompressorAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, T("Attack"), 1.0, 1000.0, 0.01, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setAttackTime);

  p = new AutomatableParameter(plugInLock, T("Release"), 10.0, 1000.0, 0.1, 100.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setReleaseTime);

  p = new AutomatableParameter(plugInLock, T("LookAhead"), 0.0, 50.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setLookAheadTime);

  p = new AutomatableParameter(plugInLock, T("InLevel"), -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setInputGain);

  p = new AutomatableParameter(plugInLock, T("OutLevel"), -24.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setOutputGain);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, T("Threshold"), -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setThreshold);

  p = new AutomatableParameter(plugInLock, T("Ratio"), 1.0, 40.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setRatio);

  p = new AutomatableParameter(plugInLock, T("KneeWidth"), 0.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setKneeWidth);

  p = new AutomatableParameter(plugInLock, T("Limit"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setLimiterMode);

  p = new AutomatableParameter(plugInLock, T("AutoGain"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback<SoftKneeCompressor>(wrappedCompressor, &SoftKneeCompressor::setAutoGain);

  //p = new AutomatableParameter(plugInLock, T("AntiAlias"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  //addObservedParameter(p);
  //p->setValueChangeCallback(wrappedCompressor, &SoftKneeCompressor::setAntiAlias);
  // only for test - usually on

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

//-------------------------------------------------------------------------------------------------
// CompressorModuleEditor:

CompressorModuleEditor::CompressorModuleEditor(CriticalSection *newPlugInLock, CompressorAudioModule* newCompressorAudioModule) 
: AudioModuleEditor(newPlugInLock, newCompressorAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newCompressorAudioModule != NULL ); // you must pass a valid module here

  addWidget( curveLabel = new RTextField( juce::String(T("Transfer Curve"))) );
  curveLabel->setJustification(Justification::centred);
  curveLabel->setDescription(T("Parameters for the static transfer curve"));
  curveLabel->setDescriptionField(infoField);

  addWidget( timeLabel = new RTextField( juce::String(T("Time Response"))) );
  timeLabel->setJustification(Justification::centred);
  timeLabel->setDescription(T("Parameters for the level detector"));
  timeLabel->setDescriptionField(infoField);

  addWidget( othersLabel = new RTextField( juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);

  addWidget( attackSlider = new RSlider (T("AttackSlider")) );
  attackSlider->assignParameter( moduleToEdit->getParameterByName("Attack") );
  attackSlider->setDescription(juce::String(T("Attack time for the envelope detector")));
  attackSlider->setDescriptionField(infoField);
  attackSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( releaseSlider = new RSlider (T("ReleaseSlider")) );
  releaseSlider->assignParameter( moduleToEdit->getParameterByName("Release") );
  releaseSlider->setDescription(juce::String(T("Release time for the envelope detector")));
  releaseSlider->setDescriptionField(infoField);
  releaseSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( lookAheadSlider = new RSlider (T("LookAheadSlider")) );
  lookAheadSlider->assignParameter( moduleToEdit->getParameterByName("LookAhead") );
  lookAheadSlider->setDescription(juce::String(T("LookAhead time for the envelope detector")));
  lookAheadSlider->setDescriptionField(infoField);
  lookAheadSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( inLevelSlider = new RSlider (T("InLevelSlider")) );
  inLevelSlider->assignParameter( moduleToEdit->getParameterByName("InLevel") );
  inLevelSlider->setDescription(juce::String(T("Gain for the input signal (pre compression)")));
  inLevelSlider->setDescriptionField(infoField);
  inLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( outLevelSlider = new RSlider (T("OutLevelSlider")) );
  outLevelSlider->assignParameter( moduleToEdit->getParameterByName("OutLevel") );
  outLevelSlider->setDescription(juce::String(T("Gain for the output signal (post compression)")));
  outLevelSlider->setDescriptionField(infoField);
  outLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Mix ratio between original and compressed signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( thresholdSlider = new RSlider (T("ThresholdSlider")) );
  thresholdSlider->assignParameter( moduleToEdit->getParameterByName("Threshold") );
  thresholdSlider->setDescription(juce::String(T("Threshold above which the signal will be attenuated")));
  thresholdSlider->setDescriptionField(infoField);
  thresholdSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( ratioSlider = new RSlider (T("RatioSlider")) );
  ratioSlider->assignParameter( moduleToEdit->getParameterByName("Ratio") );
  ratioSlider->setDescription(juce::String(T("Ratio .....find a good short explanation")));
  ratioSlider->setDescriptionField(infoField);
  ratioSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( kneeSlider = new RSlider (T("KneeSlider")) );
  kneeSlider->assignParameter( moduleToEdit->getParameterByName("KneeWidth") );
  kneeSlider->setSliderName(juce::String(T("Knee")));
  kneeSlider->setDescription(juce::String(T("Transition width between the two slopes")));
  kneeSlider->setDescriptionField(infoField);
  kneeSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( autoGainButton = new RButton(juce::String(T("AutoGain"))) );
  autoGainButton->assignParameter( moduleToEdit->getParameterByName(T("AutoGain")) );
  autoGainButton->setDescription(juce::String(T("Automatic gain compensation")));
  autoGainButton->setDescriptionField(infoField);
  autoGainButton->setClickingTogglesState(true);

  addWidget( limitButton = new RButton(juce::String(T("Limit"))) );
  limitButton->assignParameter( moduleToEdit->getParameterByName(T("Limit")) );
  limitButton->setDescription(juce::String(T("Limit signal at the threshold (infinite ratio)")));
  limitButton->setDescriptionField(infoField);
  limitButton->setClickingTogglesState(true);

  addWidget( antiAliasButton = new RButton(juce::String(T("AntiAlias"))) );
  antiAliasButton->assignParameter( moduleToEdit->getParameterByName(T("AntiAlias")) );
  antiAliasButton->setDescription(juce::String(T("AntiAliasing")));
  antiAliasButton->setDescriptionField(infoField);
  antiAliasButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void CompressorModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight()-y;  // usable height

  int h3 = h/3;
  curveParametersRect.setBounds(x, y, w/2, 2*h3);
  timeParametersRect.setBounds(x+w/2, y, w/2, 2*h3);
  y += 2*h3;
  otherParametersRect.setBounds(x, y, w, h3);

  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(curveParametersRect);
  guiLayoutRectangles.add(timeParametersRect); 
  guiLayoutRectangles.add(otherParametersRect);

  x = curveParametersRect.getX();
  y = curveParametersRect.getY();
  w = curveParametersRect.getWidth();
  h = curveParametersRect.getHeight();
  curveLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  thresholdSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  ratioSlider->setBounds(    x+4, y+4, w-8, 16);
  y += 28;  
  kneeSlider->setBounds(     x+4, y+4, w-8, 16);
  y += 20;  
  limitButton->setBounds(x+4+32, y+4, w-8-64, 16);

  x = timeParametersRect.getX();
  y = timeParametersRect.getY();
  w = timeParametersRect.getWidth();
  h = timeParametersRect.getHeight();
  timeLabel->setBounds(      x+4, y+2, w-8, 16);
  y += 20; 
  attackSlider->setBounds(   x+4, y+4, w-8, 16);
  y += 20; 
  releaseSlider->setBounds(  x+4, y+4, w-8, 16);
  y += 28; 
  lookAheadSlider->setBounds(x+4, y+4, w-8, 16);

  x = otherParametersRect.getX();
  y = otherParametersRect.getY();
  w = otherParametersRect.getWidth();
  h = otherParametersRect.getHeight();
  othersLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  inLevelSlider->setBounds(x+4, y+4, w/2-8, 16);
  outLevelSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20; 
  dryWetSlider->setBounds(x+4, y+4, w/2-8, 16);
  autoGainButton->setBounds(x+w/2+4+32, y+4, w/2-8-64, 16);
  //antiAliasButton->setBounds(0, 0, 40, 16);
}

//-------------------------------------------------------------------------------------------------
// Expander:

ExpanderAudioModule::ExpanderAudioModule(CriticalSection *newPlugInLock, rosic::SoftKneeExpander *newExpanderToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newExpanderToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedExpander = newExpanderToWrap;
  moduleName  = juce::String(T("Expander"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("ExpanderPresets")) );
  createStaticParameters();
}

void ExpanderAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, T("Attack"), 1.0, 1000.0, 0.01, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setAttackTime);

  p = new AutomatableParameter(plugInLock, T("Release"), 10.0, 1000.0, 0.1, 100.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setReleaseTime);

  p = new AutomatableParameter(plugInLock, T("LookAhead"), 0.0, 50.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setLookAheadTime);

  p = new AutomatableParameter(plugInLock, T("InLevel"), -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setInputGain);

  p = new AutomatableParameter(plugInLock, T("OutLevel"), -24.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setOutputGain);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, T("Threshold"), -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setThreshold);

  p = new AutomatableParameter(plugInLock, T("Ratio"), 1.0, 40.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setRatio);

  p = new AutomatableParameter(plugInLock, T("KneeWidth"), 0.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<SoftKneeExpander>(wrappedExpander, &SoftKneeExpander::setKneeWidth);

  //p = new Parameter("Gate", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  //addObservedParameter(p);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

ExpanderModuleEditor::ExpanderModuleEditor(CriticalSection *newPlugInLock, ExpanderAudioModule* newExpanderAudioModule) 
: AudioModuleEditor(newPlugInLock, newExpanderAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newExpanderAudioModule != NULL ); // you must pass a valid module here

  addWidget( curveLabel = new RTextField( juce::String(T("Transfer Curve"))) );
  curveLabel->setJustification(Justification::centred);
  curveLabel->setDescription(T("Parameters for the static transfer curve"));
  curveLabel->setDescriptionField(infoField);

  addWidget( timeLabel = new RTextField( juce::String(T("Time Response"))) );
  timeLabel->setJustification(Justification::centred);
  timeLabel->setDescription(T("Parameters for the level detector"));
  timeLabel->setDescriptionField(infoField);

  addWidget( othersLabel = new RTextField( juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);

  addWidget( attackSlider = new RSlider (T("AttackSlider")) );
  attackSlider->assignParameter( moduleToEdit->getParameterByName("Attack") );
  attackSlider->setDescription(juce::String(T("Attack time for the envelope detector")));
  attackSlider->setDescriptionField(infoField);
  attackSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( releaseSlider = new RSlider (T("ReleaseSlider")) );
  releaseSlider->assignParameter( moduleToEdit->getParameterByName("Release") );
  releaseSlider->setDescription(juce::String(T("Release time for the envelope detector")));
  releaseSlider->setDescriptionField(infoField);
  releaseSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( lookAheadSlider = new RSlider (T("LookAheadSlider")) );
  lookAheadSlider->assignParameter( moduleToEdit->getParameterByName("LookAhead") );
  lookAheadSlider->setDescription(juce::String(T("LookAhead time for the envelope detector")));
  lookAheadSlider->setDescriptionField(infoField);
  lookAheadSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( inLevelSlider = new RSlider (T("InLevelSlider")) );
  inLevelSlider->assignParameter( moduleToEdit->getParameterByName("InLevel") );
  inLevelSlider->setDescription(juce::String(T("Gain for the input signal (pre compression)")));
  inLevelSlider->setDescriptionField(infoField);
  inLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( outLevelSlider = new RSlider (T("OutLevelSlider")) );
  outLevelSlider->assignParameter( moduleToEdit->getParameterByName("OutLevel") );
  outLevelSlider->setDescription(juce::String(T("Gain for the output signal (post compression)")));
  outLevelSlider->setDescriptionField(infoField);
  outLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Mix ratio between original and compressed signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( thresholdSlider = new RSlider (T("ThresholdSlider")) );
  thresholdSlider->assignParameter( moduleToEdit->getParameterByName("Threshold") );
  thresholdSlider->setDescription(juce::String(T("Threshold above which the signal will be attenuated")));
  thresholdSlider->setDescriptionField(infoField);
  thresholdSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( ratioSlider = new RSlider (T("RatioSlider")) );
  ratioSlider->assignParameter( moduleToEdit->getParameterByName("Ratio") );
  ratioSlider->setDescription(juce::String(T("Ratio .....find a good short explanation")));
  ratioSlider->setDescriptionField(infoField);
  ratioSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( kneeSlider = new RSlider (T("KneeSlider")) );
  kneeSlider->assignParameter( moduleToEdit->getParameterByName("KneeWidth") );
  kneeSlider->setSliderName(juce::String(T("Knee")));
  kneeSlider->setDescription(juce::String(T("Transition width between the two slopes")));
  kneeSlider->setDescriptionField(infoField);
  kneeSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( gateButton = new RButton(juce::String(T("Gate"))), true, false); // not used
  gateButton->assignParameter( moduleToEdit->getParameterByName(T("Gate")) );
  gateButton->setDescription(juce::String(T("Gate signal at the threshold (infinite ratio)")));
  gateButton->setDescriptionField(infoField);
  gateButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void ExpanderModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight()-y;  // usable height

  int h3 = h/3;
  curveParametersRect.setBounds(x, y, w/2, 2*h3);
  timeParametersRect.setBounds(x+w/2, y, w/2, 2*h3);
  y += 2*h3;
  otherParametersRect.setBounds(x, y, w, h3);

  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(curveParametersRect);
  guiLayoutRectangles.add(timeParametersRect); 
  guiLayoutRectangles.add(otherParametersRect);

  x = curveParametersRect.getX();
  y = curveParametersRect.getY();
  w = curveParametersRect.getWidth();
  h = curveParametersRect.getHeight();
  curveLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  thresholdSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  ratioSlider->setBounds(    x+4, y+4, w-8, 16);
  y += 28;  
  kneeSlider->setBounds(     x+4, y+4, w-8, 16);
  y += 20;  
  gateButton->setBounds(x+4+32, y+4, w-8-64, 16);

  x = timeParametersRect.getX();
  y = timeParametersRect.getY();
  w = timeParametersRect.getWidth();
  h = timeParametersRect.getHeight();
  timeLabel->setBounds(      x+4, y+2, w-8, 16);
  y += 20; 
  attackSlider->setBounds(   x+4, y+4, w-8, 16);
  y += 20; 
  releaseSlider->setBounds(  x+4, y+4, w-8, 16);
  y += 28; 
  lookAheadSlider->setBounds(x+4, y+4, w-8, 16);

  x = otherParametersRect.getX();
  y = otherParametersRect.getY();
  w = otherParametersRect.getWidth();
  h = otherParametersRect.getHeight();
  othersLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  inLevelSlider->setBounds(x+4, y+4, w/2-8, 16);
  outLevelSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20; 
  dryWetSlider->setBounds(x+4, y+4, w/2-8, 16);
}

//-------------------------------------------------------------------------------------------------
// Limiter:

LimiterAudioModule::LimiterAudioModule(CriticalSection *newPlugInLock, rosic::Limiter *newLimiterToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newLimiterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedLimiter = newLimiterToWrap;
  moduleName  = juce::String(T("Limiter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("LimiterPresets")) );
  createStaticParameters();
}

void LimiterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, T("Attack"), 0.0, 10.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setAttackTime);

  p = new AutomatableParameter(plugInLock, T("Release"), 10.0, 1000.0, 0.1, 100.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setReleaseTime);

  p = new AutomatableParameter(plugInLock, T("LookAhead"), 0.0, 50.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setLookAheadTime);

  p = new AutomatableParameter(plugInLock, T("InLevel"), -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setInputGain);

  p = new AutomatableParameter(plugInLock, T("OutLevel"), -24.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setOutputGain);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, T("Limit"), -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Limiter>(wrappedLimiter, &Limiter::setLimit);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

LimiterModuleEditor::LimiterModuleEditor(CriticalSection *newPlugInLock, LimiterAudioModule* newLimiterAudioModule) 
: AudioModuleEditor(newPlugInLock, newLimiterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newLimiterAudioModule != NULL ); // you must pass a valid module here

  addWidget( curveLabel = new RTextField( juce::String(T("Transfer Curve"))) );
  curveLabel->setJustification(Justification::centred);
  curveLabel->setDescription(T("Parameters for the static transfer curve"));
  curveLabel->setDescriptionField(infoField);

  addWidget( timeLabel = new RTextField( juce::String(T("Time Response"))) );
  timeLabel->setJustification(Justification::centred);
  timeLabel->setDescription(T("Parameters for the level detector"));
  timeLabel->setDescriptionField(infoField);

  addWidget( othersLabel = new RTextField( juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);

  addWidget( attackSlider = new RSlider(T("AttackSlider")) );
  attackSlider->assignParameter( moduleToEdit->getParameterByName(T("Attack")) );
  attackSlider->setDescription(juce::String(T("Attack time for the envelope detector")));
  attackSlider->setDescriptionField(infoField);
  attackSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( releaseSlider = new RSlider(T("ReleaseSlider")) );
  releaseSlider->assignParameter( moduleToEdit->getParameterByName(T("Release")) );
  releaseSlider->setDescription(juce::String(T("Release time for the envelope detector")));
  releaseSlider->setDescriptionField(infoField);
  releaseSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( lookAheadSlider = new RSlider(T("LookAheadSlider")) );
  lookAheadSlider->assignParameter( moduleToEdit->getParameterByName(T("LookAhead")) );
  lookAheadSlider->setDescription(juce::String(T("LookAhead time for the envelope detector")));
  lookAheadSlider->setDescriptionField(infoField);
  lookAheadSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( inLevelSlider = new RSlider(T("InLevelSlider")) );
  inLevelSlider->assignParameter( moduleToEdit->getParameterByName(T("InLevel")) );
  inLevelSlider->setDescription(juce::String(T("Gain for the input signal (pre compression)")));
  inLevelSlider->setDescriptionField(infoField);
  inLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( outLevelSlider = new RSlider(T("OutLevelSlider")) );
  outLevelSlider->assignParameter( moduleToEdit->getParameterByName(T("OutLevel")) );
  outLevelSlider->setDescription(juce::String(T("Gain for the output signal (post compression)")));
  outLevelSlider->setDescriptionField(infoField);
  outLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider(T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName(T("DryWetRatio")) );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Mix ratio between original and compressed signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( limitSlider = new RSlider(T("LimitSlider")) );
  limitSlider->assignParameter( moduleToEdit->getParameterByName(T("Limit")) );
  limitSlider->setDescription(juce::String(T("Limit above which the signal will be attenuated")));
  limitSlider->setDescriptionField(infoField);
  limitSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  updateWidgetsAccordingToState();
}

void LimiterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight()-y;  // usable height

  int h3 = h/3;
  curveParametersRect.setBounds(x, y, w/2, 2*h3);
  timeParametersRect.setBounds(x+w/2, y, w/2, 2*h3);
  y += 2*h3;
  otherParametersRect.setBounds(x, y, w, h3);

  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(curveParametersRect);
  guiLayoutRectangles.add(timeParametersRect); 
  guiLayoutRectangles.add(otherParametersRect);

  x = curveParametersRect.getX();
  y = curveParametersRect.getY();
  w = curveParametersRect.getWidth();
  h = curveParametersRect.getHeight();
  curveLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  limitSlider->setBounds(x+4, y+4, w-8, 16);

  x = timeParametersRect.getX();
  y = timeParametersRect.getY();
  w = timeParametersRect.getWidth();
  h = timeParametersRect.getHeight();
  timeLabel->setBounds(      x+4, y+2, w-8, 16);
  y += 20; 
  attackSlider->setBounds(   x+4, y+4, w-8, 16);
  y += 20; 
  releaseSlider->setBounds(  x+4, y+4, w-8, 16);
  y += 28; 
  lookAheadSlider->setBounds(x+4, y+4, w-8, 16);

  x = otherParametersRect.getX();
  y = otherParametersRect.getY();
  w = otherParametersRect.getWidth();
  h = otherParametersRect.getHeight();
  othersLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  inLevelSlider->setBounds(x+4, y+4, w/2-8, 16);
  outLevelSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20; 
  dryWetSlider->setBounds(x+4, y+4, w/2-8, 16);
}

//-------------------------------------------------------------------------------------------------
// NoiseGate:

NoiseGateAudioModule::NoiseGateAudioModule(CriticalSection *newPlugInLock, rosic::NoiseGate *newNoiseGateToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newNoiseGateToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedNoiseGate = newNoiseGateToWrap;
  moduleName  = juce::String(T("NoiseGate"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("NoiseGatePresets")) );
  createStaticParameters();
}

void NoiseGateAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, T("Attack"), 0.1, 10.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setAttackTime);

  p = new AutomatableParameter(plugInLock, T("Hold"), 0.1, 100.0, 0.01, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setHoldTime);

  p = new AutomatableParameter(plugInLock, T("Release"), 1.0, 1000.0, 1.0, 10.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setReleaseTime);

  p = new AutomatableParameter(plugInLock, T("LookAhead"), 0.0, 50.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setLookAheadTime);

  p = new AutomatableParameter(plugInLock, T("InLevel"), -24.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setInputGain);

  p = new AutomatableParameter(plugInLock, T("OutLevel"), -24.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setOutputGain);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, T("Threshold"), -48.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setThreshold);

  p = new AutomatableParameter(plugInLock, T("Hysteresis"), 0.0, 12.0, 0.1, 3.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<NoiseGate>(wrappedNoiseGate, &NoiseGate::setHysteresis);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

NoiseGateModuleEditor::NoiseGateModuleEditor(CriticalSection *newPlugInLock, NoiseGateAudioModule* newNoiseGateAudioModule) 
: AudioModuleEditor(newPlugInLock, newNoiseGateAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newNoiseGateAudioModule != NULL ); // you must pass a valid module here

  addWidget( curveLabel = new RTextField( juce::String(T("Transfer Curve"))) );
  curveLabel->setJustification(Justification::centred);
  curveLabel->setDescription(T("Parameters for the static transfer curve"));
  curveLabel->setDescriptionField(infoField);

  addWidget( timeLabel = new RTextField( juce::String(T("Time Response"))) );
  timeLabel->setJustification(Justification::centred);
  timeLabel->setDescription(T("Parameters for the level detector"));
  timeLabel->setDescriptionField(infoField);

  addWidget( othersLabel = new RTextField( juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);

  addWidget( attackSlider = new RSlider (T("AttackSlider")) );
  attackSlider->assignParameter( moduleToEdit->getParameterByName("Attack") );
  attackSlider->setDescription(juce::String(T("Time for the gate to open")));
  attackSlider->setDescriptionField(infoField);
  attackSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( holdSlider = new RSlider (T("HoldSlider")) );
  holdSlider->assignParameter( moduleToEdit->getParameterByName("Hold") );
  holdSlider->setDescription(juce::String(T("Time for the gate to stay open after close threshold was undercut")));
  holdSlider->setDescriptionField(infoField);
  holdSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( releaseSlider = new RSlider (T("ReleaseSlider")) );
  releaseSlider->assignParameter( moduleToEdit->getParameterByName("Release") );
  releaseSlider->setDescription(juce::String(T("Release time for the envelope detector")));
  releaseSlider->setDescriptionField(infoField);
  releaseSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( lookAheadSlider = new RSlider (T("LookAheadSlider")) );
  lookAheadSlider->assignParameter( moduleToEdit->getParameterByName("LookAhead") );
  lookAheadSlider->setDescription(juce::String(T("LookAhead time for the envelope detector")));
  lookAheadSlider->setDescriptionField(infoField);
  lookAheadSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( inLevelSlider = new RSlider (T("InLevelSlider")) );
  inLevelSlider->assignParameter( moduleToEdit->getParameterByName("InLevel") );
  inLevelSlider->setDescription(juce::String(T("Gain for the input signal (pre compression)")));
  inLevelSlider->setDescriptionField(infoField);
  inLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( outLevelSlider = new RSlider (T("OutLevelSlider")) );
  outLevelSlider->assignParameter( moduleToEdit->getParameterByName("OutLevel") );
  outLevelSlider->setDescription(juce::String(T("Gain for the output signal (post compression)")));
  outLevelSlider->setDescriptionField(infoField);
  outLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Mix ratio between original and compressed signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( thresholdSlider = new RSlider (T("ThresholdSlider")) );
  thresholdSlider->assignParameter( moduleToEdit->getParameterByName("Threshold") );
  thresholdSlider->setDescription(juce::String(T("Threshold for the gate to open")));
  thresholdSlider->setDescriptionField(infoField);
  thresholdSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( hysteresisSlider = new RSlider (T("HysteresisSlider")) );
  hysteresisSlider->assignParameter( moduleToEdit->getParameterByName("Hysteresis") );
  hysteresisSlider->setDescription(juce::String(T("Difference between opening and closing threshold")));
  hysteresisSlider->setDescriptionField(infoField);
  hysteresisSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  updateWidgetsAccordingToState();
}

void NoiseGateModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight()-y;  // usable height

  int h3 = h/3;
  curveParametersRect.setBounds(x, y, w/2, 2*h3);
  timeParametersRect.setBounds(x+w/2, y, w/2, 2*h3);
  y += 2*h3;
  otherParametersRect.setBounds(x, y, w, h3);

  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(curveParametersRect);
  guiLayoutRectangles.add(timeParametersRect); 
  guiLayoutRectangles.add(otherParametersRect);

  x = curveParametersRect.getX();
  y = curveParametersRect.getY();
  w = curveParametersRect.getWidth();
  h = curveParametersRect.getHeight();
  curveLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  thresholdSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  hysteresisSlider->setBounds(    x+4, y+4, w-8, 16);

  x = timeParametersRect.getX();
  y = timeParametersRect.getY();
  w = timeParametersRect.getWidth();
  h = timeParametersRect.getHeight();
  timeLabel->setBounds(      x+4, y+2, w-8, 16);
  y += 20; 
  attackSlider->setBounds(   x+4, y+4, w-8, 16);
  y += 20; 
  holdSlider->setBounds(   x+4, y+4, w-8, 16);
  y += 20; 
  releaseSlider->setBounds(  x+4, y+4, w-8, 16);
  y += 28; 
  lookAheadSlider->setBounds(x+4, y+4, w-8, 16);

  x = otherParametersRect.getX();
  y = otherParametersRect.getY();
  w = otherParametersRect.getWidth();
  h = otherParametersRect.getHeight();
  othersLabel->setBounds(x+4, y+2, w-8, 16);
  y += 16; 
  inLevelSlider->setBounds(x+4, y+4, w/2-8, 16);
  outLevelSlider->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20; 
  dryWetSlider->setBounds(x+4, y+4, w/2-8, 16);
}

//=================================================================================================
// Filters:

//-------------------------------------------------------------------------------------------------
// CombBank:

CombBankAudioModule::CombBankAudioModule(CriticalSection *newPlugInLock, rosic::CombBank *newCombBankToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newCombBankToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedCombBank = newCombBankToWrap;
  moduleName  = juce::String(T("CombBank"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("CombBankPresets")) );
  createStaticParameters();
}

void CombBankAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, "Level", -48.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setLevel);

  p = new AutomatableParameter(plugInLock, "Frequency", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setReferenceFrequency);

  p = new AutomatableParameter(plugInLock, "Detune", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setDetune);

  p = new AutomatableParameter(plugInLock, "Pan1", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setPan1);

  p = new AutomatableParameter(plugInLock, "Pan2", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setPan2);

  p = new AutomatableParameter(plugInLock, "OddOnly", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setOddOnlyMode);

  p = new AutomatableParameter(plugInLock, "DecayTime", 0.1, 10.0, 0.1, 3.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setDecayTime);

  p = new AutomatableParameter(plugInLock, "HighDecayScale", 0.1, 10.0, 0.01, 0.3, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setHighDecayScale);

  p = new AutomatableParameter(plugInLock, "LowDecayScale", 0.1, 10.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setLowDecayScale);

  p = new AutomatableParameter(plugInLock, "HighCrossoverFrequency", 20.0, 20000.0, 0.0, 4000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setHighCrossoverFreq);

  p = new AutomatableParameter(plugInLock, "LowCrossoverFrequency", 20.0, 20000.0, 0.0, 250.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombBank, &CombBank::setLowCrossoverFreq);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

CombBankModuleEditor::CombBankModuleEditor(CriticalSection *newPlugInLock, CombBankAudioModule* newCombBankAudioModule) 
: AudioModuleEditor(newPlugInLock, newCombBankAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newCombBankAudioModule != NULL ); // you must pass a valid module here

  addWidget( toneLabel = new RTextField( juce::String(T("Tone"))) );
  toneLabel->setJustification(Justification::centred);
  toneLabel->setDescription(T("Parameters for the tone/timbre of the comb"));
  toneLabel->setDescriptionField(infoField);

  addWidget( decayLabel = new RTextField( juce::String(T("Decay"))) );
  decayLabel->setJustification(Justification::centred);
  decayLabel->setDescription(T("Parameters controlling the frequency dependent decay"));
  decayLabel->setDescriptionField(infoField);

  /*
  addWidget( othersLabel = new RTextField(juce::String(T("OthersLabel")), juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);
  */

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Fundamental frequency of the resonator")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( levelSlider = new RSlider (T("LevelSlider")) );
  levelSlider->assignParameter( moduleToEdit->getParameterByName("Level") );
  levelSlider->setDescription(juce::String(T("Overall output level")));
  levelSlider->setDescriptionField(infoField);
  levelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( detuneSlider = new RSlider (T("DetuneSlider")) );
  detuneSlider->assignParameter( moduleToEdit->getParameterByName("Detune") );
  detuneSlider->setDescription(juce::String(T("Detuning between the two resonators")));
  detuneSlider->setDescriptionField(infoField);
  detuneSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit2);

  addWidget( pan1Slider = new RSlider (T("Pan1Slider")) );
  pan1Slider->assignParameter( moduleToEdit->getParameterByName("Pan1") );
  pan1Slider->setDescription(juce::String(T("Panorama position of resonator 1")));
  pan1Slider->setDescriptionField(infoField);
  pan1Slider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( pan2Slider = new RSlider (T("Pan2Slider")) );
  pan2Slider->assignParameter( moduleToEdit->getParameterByName("Pan2") );
  pan2Slider->setDescription(juce::String(T("Panorama position of resonator 2")));
  pan2Slider->setDescriptionField(infoField);
  pan2Slider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( decayTimeSlider = new RSlider (T("DecayTimeSlider")) );
  decayTimeSlider->assignParameter( moduleToEdit->getParameterByName("DecayTime") );
  decayTimeSlider->setDescription(juce::String(T("Time for the tail to decay to -60 dB")));
  decayTimeSlider->setDescriptionField(infoField);
  decayTimeSlider->setStringConversionFunction(&rojue::secondsToStringWithUnit3);

  addWidget( highDecayScaleSlider = new RSlider (T("HighDecayScaleSlider")) );
  highDecayScaleSlider->assignParameter( moduleToEdit->getParameterByName("HighDecayScale") );
  highDecayScaleSlider->setDescription(juce::String(T("Scaler for the decay-time at high frequencies")));
  highDecayScaleSlider->setDescriptionField(infoField);
  highDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( lowDecayScaleSlider = new RSlider (T("LowDecayScaleSlider")) );
  lowDecayScaleSlider->assignParameter( moduleToEdit->getParameterByName("LowDecayScale") );
  lowDecayScaleSlider->setDescription(juce::String(T("Scaler for the decay-time at low frequencies")));
  lowDecayScaleSlider->setDescriptionField(infoField);
  lowDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( highFreqSlider = new RSlider (T("HighFreqSlider")) );
  highFreqSlider->assignParameter( moduleToEdit->getParameterByName("HighCrossoverFrequency") );
  highFreqSlider->setDescription(juce::String(T("Crossover frequency between mid and high frequencies")));
  highFreqSlider->setSliderName(juce::String(T("HighFreq")));
  highFreqSlider->setDescriptionField(infoField);
  highFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( lowFreqSlider = new RSlider (T("LowFreqSlider")) );
  lowFreqSlider->assignParameter( moduleToEdit->getParameterByName("LowCrossoverFrequency") );
  lowFreqSlider->setSliderName(juce::String(T("LowFreq")));
  lowFreqSlider->setDescription(juce::String(T("Crossover frequency between low and mid frequencies")));
  lowFreqSlider->setDescriptionField(infoField);
  lowFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( oddOnlyButton = new RButton(juce::String(T("OddOnly"))) );
  oddOnlyButton->assignParameter( moduleToEdit->getParameterByName(T("OddOnly")) );
  oddOnlyButton->setDescription(juce::String(T("Lets the comb create only odd harmonics.")));
  oddOnlyButton->setDescriptionField(infoField);
  oddOnlyButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void CombBankModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();

  dryWetSlider->setBounds(x+4,     y+4, w/2-8, 16);
  levelSlider->setBounds( x+w/2+4, y+4, w/2-8, 16);

  y = dryWetSlider->getBottom()+4;
  h = getHeight()-y;
  toneParametersRect.setBounds( x,     y, w/2, h);
  decayParametersRect.setBounds(x+w/2, y, w/2, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(toneParametersRect);
  guiLayoutRectangles.add(decayParametersRect); 

  x = toneParametersRect.getX();
  y = toneParametersRect.getY();
  w = toneParametersRect.getWidth();
  h = toneParametersRect.getHeight();
  toneLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  detuneSlider->setBounds(    x+4, y+4, w-8, 16);
  y += 20;  
  pan1Slider->setBounds(     x+4, y+4, w-8, 16);
  y += 20;  
  pan2Slider->setBounds(     x+4, y+4, w-8, 16);
  y += 20;  
  oddOnlyButton->setBounds(x+4+32, y+4, w-8-64, 16);

  x = decayParametersRect.getX();
  y = decayParametersRect.getY();
  w = decayParametersRect.getWidth();
  h = decayParametersRect.getHeight();
  decayLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  decayTimeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 28;  
  highDecayScaleSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowDecayScaleSlider->setBounds(x+4, y+4, w-8, 16);
  y += 28;  
  highFreqSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowFreqSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// CombResonator:

CombResonatorAudioModule::CombResonatorAudioModule(CriticalSection *newPlugInLock, rosic::CombResonatorStereo *newCombResonatorToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newCombResonatorToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedCombResonator = newCombResonatorToWrap;
  moduleName  = juce::String(T("CombResonator"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("CombResonatorPresets")) );
  createStaticParameters();
}

void CombResonatorAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, "Level", -48.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setLevel);

  p = new AutomatableParameter(plugInLock, "Frequency", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setFrequency);

  p = new AutomatableParameter(plugInLock, "Detune", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setDetune);

  p = new AutomatableParameter(plugInLock, "Pan1", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setPan1);

  p = new AutomatableParameter(plugInLock, "Pan2", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setPan2);

  p = new AutomatableParameter(plugInLock, "OddOnly", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setOddOnlyMode);

  p = new AutomatableParameter(plugInLock, "DecayTime", 0.1, 10.0, 0.1, 3.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setDecayTime);

  p = new AutomatableParameter(plugInLock, "HighDecayScale", 0.1, 10.0, 0.01, 0.3, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setHighDecayScale);

  p = new AutomatableParameter(plugInLock, "LowDecayScale", 0.1, 10.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setLowDecayScale);

  p = new AutomatableParameter(plugInLock, "HighCrossoverFrequency", 20.0, 20000.0, 0.0, 4000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setHighCrossoverFreq);

  p = new AutomatableParameter(plugInLock, "LowCrossoverFrequency", 20.0, 20000.0, 0.0, 250.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedCombResonator, &CombResonatorStereo::setLowCrossoverFreq);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

CombResonatorModuleEditor::CombResonatorModuleEditor(CriticalSection *newPlugInLock, CombResonatorAudioModule* newCombResonatorAudioModule) 
: AudioModuleEditor(newPlugInLock, newCombResonatorAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newCombResonatorAudioModule != NULL ); // you must pass a valid module here

  addWidget( toneLabel = new RTextField( juce::String(T("Tone"))) );
  toneLabel->setJustification(Justification::centred);
  toneLabel->setDescription(T("Parameters for the tone/timbre of the comb"));
  toneLabel->setDescriptionField(infoField);

  addWidget( decayLabel = new RTextField( juce::String(T("Decay"))) );
  decayLabel->setJustification(Justification::centred);
  decayLabel->setDescription(T("Parameters controlling the frequency dependent decay"));
  decayLabel->setDescriptionField(infoField);

  /*
  addWidget( othersLabel = new RTextField(juce::String(T("OthersLabel")), juce::String(T("Gain And Mix"))) );
  othersLabel->setJustification(Justification::centred);
  othersLabel->setDescription(T("Parameters for input-/output-gain and dry/wet mix"));
  othersLabel->setDescriptionField(infoField);
  */

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Fundamental frequency of the resonator")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( levelSlider = new RSlider (T("LevelSlider")) );
  levelSlider->assignParameter( moduleToEdit->getParameterByName("Level") );
  levelSlider->setDescription(juce::String(T("Overall output level")));
  levelSlider->setDescriptionField(infoField);
  levelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( detuneSlider = new RSlider (T("DetuneSlider")) );
  detuneSlider->assignParameter( moduleToEdit->getParameterByName("Detune") );
  detuneSlider->setDescription(juce::String(T("Detuning between the two resonators")));
  detuneSlider->setDescriptionField(infoField);
  detuneSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit2);

  addWidget( pan1Slider = new RSlider (T("Pan1Slider")) );
  pan1Slider->assignParameter( moduleToEdit->getParameterByName("Pan1") );
  pan1Slider->setDescription(juce::String(T("Panorama position of resonator 1")));
  pan1Slider->setDescriptionField(infoField);
  pan1Slider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( pan2Slider = new RSlider (T("Pan2Slider")) );
  pan2Slider->assignParameter( moduleToEdit->getParameterByName("Pan2") );
  pan2Slider->setDescription(juce::String(T("Panorama position of resonator 2")));
  pan2Slider->setDescriptionField(infoField);
  pan2Slider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( decayTimeSlider = new RSlider (T("DecayTimeSlider")) );
  decayTimeSlider->assignParameter( moduleToEdit->getParameterByName("DecayTime") );
  decayTimeSlider->setDescription(juce::String(T("Time for the tail to decay to -60 dB")));
  decayTimeSlider->setDescriptionField(infoField);
  decayTimeSlider->setStringConversionFunction(&rojue::secondsToStringWithUnit3);

  addWidget( highDecayScaleSlider = new RSlider (T("HighDecayScaleSlider")) );
  highDecayScaleSlider->assignParameter( moduleToEdit->getParameterByName("HighDecayScale") );
  highDecayScaleSlider->setDescription(juce::String(T("Scaler for the decay-time at high frequencies")));
  highDecayScaleSlider->setDescriptionField(infoField);
  highDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( lowDecayScaleSlider = new RSlider (T("LowDecayScaleSlider")) );
  lowDecayScaleSlider->assignParameter( moduleToEdit->getParameterByName("LowDecayScale") );
  lowDecayScaleSlider->setDescription(juce::String(T("Scaler for the decay-time at low frequencies")));
  lowDecayScaleSlider->setDescriptionField(infoField);
  lowDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( highFreqSlider = new RSlider (T("HighFreqSlider")) );
  highFreqSlider->assignParameter( moduleToEdit->getParameterByName("HighCrossoverFrequency") );
  highFreqSlider->setDescription(juce::String(T("Crossover frequency between mid and high frequencies")));
  highFreqSlider->setSliderName(juce::String(T("HighFreq")));
  highFreqSlider->setDescriptionField(infoField);
  highFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( lowFreqSlider = new RSlider(T("LowFreqSlider")) );
  lowFreqSlider->assignParameter( moduleToEdit->getParameterByName("LowCrossoverFrequency") );
  lowFreqSlider->setSliderName(juce::String(T("LowFreq")));
  lowFreqSlider->setDescription(juce::String(T("Crossover frequency between low and mid frequencies")));
  lowFreqSlider->setDescriptionField(infoField);
  lowFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( oddOnlyButton = new RButton(juce::String(T("OddOnly"))) );
  oddOnlyButton->assignParameter( moduleToEdit->getParameterByName(T("OddOnly")) );
  oddOnlyButton->setDescription(juce::String(T("Lets the comb create only odd harmonics.")));
  oddOnlyButton->setDescriptionField(infoField);
  oddOnlyButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void CombResonatorModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();

  dryWetSlider->setBounds(x+4,     y+4, w/2-8, 16);
  levelSlider->setBounds( x+w/2+4, y+4, w/2-8, 16);

  y = dryWetSlider->getBottom()+4;
  h = getHeight()-y;
  toneParametersRect.setBounds( x,     y, w/2, h);
  decayParametersRect.setBounds(x+w/2, y, w/2, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(toneParametersRect);
  guiLayoutRectangles.add(decayParametersRect); 

  x = toneParametersRect.getX();
  y = toneParametersRect.getY();
  w = toneParametersRect.getWidth();
  h = toneParametersRect.getHeight();
  toneLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  detuneSlider->setBounds(    x+4, y+4, w-8, 16);
  y += 20;  
  pan1Slider->setBounds(     x+4, y+4, w-8, 16);
  y += 20;  
  pan2Slider->setBounds(     x+4, y+4, w-8, 16);
  y += 20;  
  oddOnlyButton->setBounds(x+4+32, y+4, w-8-64, 16);

  x = decayParametersRect.getX();
  y = decayParametersRect.getY();
  w = decayParametersRect.getWidth();
  h = decayParametersRect.getHeight();
  decayLabel->setBounds(x+4, y+2, w-8, 16);
  y += 20; 
  decayTimeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 28;  
  highDecayScaleSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowDecayScaleSlider->setBounds(x+4, y+4, w-8, 16);
  y += 28;  
  highFreqSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowFreqSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// DualTwoPoleFilter:

DualTwoPoleFilterAudioModule::DualTwoPoleFilterAudioModule(CriticalSection *newPlugInLock, 
                                                           rosic::DualTwoPoleFilter *newDualTwoPoleFilterToWrap)
                                                            : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newDualTwoPoleFilterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedDualTwoPoleFilter = newDualTwoPoleFilterToWrap;
  moduleName  = juce::String(T("DualTwoPoleFilter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("DualTwoPoleFilterPresets")) );
  createStaticParameters();
}

void DualTwoPoleFilterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new ParameterTwoPoleFilterMode(plugInLock, juce::String(T("Mode1")));
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setMode1);

  p = new AutomatableParameter(plugInLock, "Frequency1", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setFrequency1);

  p = new AutomatableParameter(plugInLock, "Gain1", -48, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setGain1);

  p = new AutomatableParameter(plugInLock, "Bandwidth1", 0.2, 6.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setBandwidth1);

  p = new ParameterTwoPoleFilterMode(plugInLock, juce::String(T("Mode2")));
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setMode2);

  p = new AutomatableParameter(plugInLock, "Frequency2", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setFrequency2);

  p = new AutomatableParameter(plugInLock, "Gain2", -48, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setGain2);

  p = new AutomatableParameter(plugInLock, "Bandwidth2", 0.2, 6.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setBandwidth2);

  p = new AutomatableParameter(plugInLock, "SerialParallelBlend", 0.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setSerialParallelBlend);

  p = new AutomatableParameter(plugInLock, "FrequencyScale", 0.125, 8.0, 0.001, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setFrequencyScale);

  p = new AutomatableParameter(plugInLock, "GainScale", 0.125, 8.0, 0.001, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setGainScale);

  p = new AutomatableParameter(plugInLock, "BandwidthScale", 0.125, 8.0, 0.001, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedDualTwoPoleFilter, &DualTwoPoleFilter::setBandwidthScale);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

DualTwoPoleFilterModuleEditor::DualTwoPoleFilterModuleEditor(CriticalSection *newPlugInLock, 
                                                             DualTwoPoleFilterAudioModule* newDualTwoPoleFilterAudioModule) 
: AudioModuleEditor(newPlugInLock, newDualTwoPoleFilterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newDualTwoPoleFilterAudioModule != NULL ); // you must pass a valid module here
  twoPoleFilterModuleToEdit = newDualTwoPoleFilterAudioModule;

  addWidget( globalLabel = new RTextField( juce::String(T("Global"))) );
  globalLabel->setJustification(Justification::centredLeft);
  globalLabel->setDescription(T("Global parameters, affecting both filters"));
  globalLabel->setDescriptionField(infoField);

  addWidget( filter1Label = new RTextField( juce::String(T("Filter 1"))) );
  filter1Label->setJustification(Justification::centred);
  filter1Label->setDescription(T("Parameters for the first filter"));
  filter1Label->setDescriptionField(infoField);

  addWidget( filter2Label = new RTextField( juce::String(T("Filter 2"))) );
  filter2Label->setJustification(Justification::centred);
  filter2Label->setDescription(T("Parameters for the second filter"));
  filter2Label->setDescriptionField(infoField);

  addWidget( modeComboBox1 = new RComboBox(juce::String(T("ModeComboBox1"))) );
  modeComboBox1->assignParameter( moduleToEdit->getParameterByName(T("Mode1")) );
  modeComboBox1->setDescription(T("Choose the mode of the first filter"));
  modeComboBox1->setDescriptionField(infoField);
  modeComboBox1->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider1 = new RSlider (T("FrequencySlider1")) );
  frequencySlider1->assignParameter( moduleToEdit->getParameterByName("Frequency1") );
  frequencySlider1->setDescription(juce::String(T("Characteristic frequency of the first filter")));
  frequencySlider1->setDescriptionField(infoField);
  frequencySlider1->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( gainSlider1 = new RSlider (T("GainSlider1")) );
  gainSlider1->assignParameter( moduleToEdit->getParameterByName("Gain1") );
  gainSlider1->setDescription(juce::String(T("Gain of the first filter")));
  gainSlider1->setDescriptionField(infoField);
  gainSlider1->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( bandwidthSlider1 = new RSlider (T("BandwidthSlider1")) );
  bandwidthSlider1->assignParameter( moduleToEdit->getParameterByName("Bandwidth1") );
  bandwidthSlider1->setDescription(juce::String(T("Bandwidth of the first filter")));
  bandwidthSlider1->setDescriptionField(infoField);
  bandwidthSlider1->setStringConversionFunction(&rojue::octavesToStringWithUnit2);

  addWidget( modeComboBox2 = new RComboBox(juce::String(T("ModeComboBox2"))) );
  modeComboBox2->assignParameter( moduleToEdit->getParameterByName(T("Mode2")) );
  modeComboBox2->setDescription(T("Choose the mode of the second filter"));
  modeComboBox2->setDescriptionField(infoField);
  modeComboBox2->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider2 = new RSlider (T("FrequencySlider2")) );
  frequencySlider2->assignParameter( moduleToEdit->getParameterByName("Frequency2") );
  frequencySlider2->setDescription(juce::String(T("Characteristic frequency of the second filter")));
  frequencySlider2->setDescriptionField(infoField);
  frequencySlider2->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( gainSlider2 = new RSlider (T("GainSlider2")) );
  gainSlider2->assignParameter( moduleToEdit->getParameterByName("Gain2") );
  gainSlider2->setDescription(juce::String(T("Gain of the second filter")));
  gainSlider2->setDescriptionField(infoField);
  gainSlider2->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( bandwidthSlider2 = new RSlider (T("BandwidthSlider2")) );
  bandwidthSlider2->assignParameter( moduleToEdit->getParameterByName("Bandwidth2") );
  bandwidthSlider2->setDescription(juce::String(T("Bandwidth of the second filter")));
  bandwidthSlider2->setDescriptionField(infoField);
  bandwidthSlider2->setStringConversionFunction(&rojue::octavesToStringWithUnit2);

  addWidget( serialParallelBlendSlider = new RSlider (T("SerialParallelBlendSlider")) );
  serialParallelBlendSlider->assignParameter( moduleToEdit->getParameterByName("SerialParallelBlend") );
  serialParallelBlendSlider->setSliderName(juce::String(T("S/P")));
  serialParallelBlendSlider->setDescription(juce::String(T("Varies smoothly between serial and parallel connection")));
  serialParallelBlendSlider->setDescriptionField(infoField);
  serialParallelBlendSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( frequencyScaleSlider = new RSlider (T("FrequencyScaleSlider")) );
  frequencyScaleSlider->assignParameter( moduleToEdit->getParameterByName("FrequencyScale") );
  frequencyScaleSlider->setDescription(juce::String(T("Scales the frequency of both filters")));
  frequencyScaleSlider->setDescriptionField(infoField);
  frequencyScaleSlider->setStringConversionFunction(&rojue::valueToString3);

  addWidget( gainScaleSlider = new RSlider (T("GainScaleSlider")) );
  gainScaleSlider->assignParameter( moduleToEdit->getParameterByName("GainScale") );
  gainScaleSlider->setSliderName(juce::String(T("GnScl")));
  gainScaleSlider->setDescription(juce::String(T("Scales the gain of both filters")));
  gainScaleSlider->setDescriptionField(infoField);
  gainScaleSlider->setStringConversionFunction(&rojue::valueToString3);

  addWidget( bandwidthScaleSlider = new RSlider (T("BandwidthScaleSlider")) );
  bandwidthScaleSlider->assignParameter( moduleToEdit->getParameterByName("BandwidthScale") );
  bandwidthScaleSlider->setSliderName(juce::String(T("BwScl")));
  bandwidthScaleSlider->setDescription(juce::String(T("Scales the bandwidth of both filters")));
  bandwidthScaleSlider->setDescriptionField(infoField);
  bandwidthScaleSlider->setStringConversionFunction(&rojue::valueToString3);

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void DualTwoPoleFilterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = 48;

  globalRect.setBounds(x, y, w, h);

  y  = globalRect.getBottom();
  //h  = 106;
  h = getHeight()-y;
  w /= 2;
  filter1Rect.setBounds(x, y, w, h);
  filter2Rect.setBounds(x+w, y, w, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(globalRect);
  guiLayoutRectangles.add(filter1Rect);
  guiLayoutRectangles.add(filter2Rect);

  x = globalRect.getX();
  y = globalRect.getY();
  w = globalRect.getWidth();
  h = globalRect.getHeight();
  globalLabel->setBounds(x+4, y+4, 48, 16);
  frequencyScaleSlider->setBounds(globalLabel->getRight(), y+4, w-globalLabel->getRight()-4, 16);
  w /= 3;
  y += 20;
  serialParallelBlendSlider->setBounds(x+4,     y+4, w-8, 16);
  gainScaleSlider->setBounds(          x+w+4,   y+4, w-8, 16);
  bandwidthScaleSlider->setBounds(     x+2*w+4, y+4, w-8, 16);

  x = filter1Rect.getX();
  y = filter1Rect.getY();
  w = filter1Rect.getWidth();
  h = filter1Rect.getHeight();
  filter1Label->setBounds(x+4, y+2, w-8, 16);
  y = filter1Label->getBottom();
  modeComboBox1->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  frequencySlider1->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  gainSlider1->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  bandwidthSlider1->setBounds(x+4, y+4, w-8, 16);

  x = filter2Rect.getX();
  y = filter2Rect.getY();
  w = filter2Rect.getWidth();
  h = filter2Rect.getHeight();
  filter2Label->setBounds(x+4, y+2, w-8, 16);
  y = filter2Label->getBottom();
  modeComboBox2->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  frequencySlider2->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  gainSlider2->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  bandwidthSlider2->setBounds(x+4, y+4, w-8, 16);
}

void DualTwoPoleFilterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updateWidgetEnablement();
}

void DualTwoPoleFilterModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);

  if( twoPoleFilterModuleToEdit == NULL )
    return;
  if( twoPoleFilterModuleToEdit->wrappedDualTwoPoleFilter == NULL )
    return;
  rosic::DualTwoPoleFilter* core = twoPoleFilterModuleToEdit->wrappedDualTwoPoleFilter;

  frequencySlider1->setEnabled(core->doesModeSupportFrequency1());
  gainSlider1->setEnabled(     core->doesModeSupportGain1());
  bandwidthSlider1->setEnabled(core->doesModeSupportBandwidth1());

  frequencySlider2->setEnabled(core->doesModeSupportFrequency2());
  gainSlider2->setEnabled(     core->doesModeSupportGain2());
  bandwidthSlider2->setEnabled(core->doesModeSupportBandwidth2());
}

//-------------------------------------------------------------------------------------------------
// FourPoleFilter:

FourPoleFilterAudioModule::FourPoleFilterAudioModule(CriticalSection *newPlugInLock, rosic::FourPoleFilter *newFourPoleFilterToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newFourPoleFilterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedFourPoleFilter = newFourPoleFilterToWrap;
  moduleName  = juce::String(T("FourPoleFilter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("FourPoleFilterPresets")) );
  createStaticParameters();
}

void FourPoleFilterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new ParameterFourPoleFilterMode(plugInLock);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedFourPoleFilter, &FourPoleFilter::setMode);

  p = new AutomatableParameter(plugInLock, "Frequency", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFourPoleFilter, &FourPoleFilter::setFrequency);

  p = new AutomatableParameter(plugInLock, "Gain", -48, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFourPoleFilter, &FourPoleFilter::setGain);

  // not yet used:
  //p = new AutomatableParameter(plugInLock, "Bandwidth", 0.2, 6.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  //addObservedParameter(p); 
  //p->setValueChangeCallback(wrappedFourPoleFilter, &FourPoleFilter::setBandwidth);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

FourPoleFilterModuleEditor::FourPoleFilterModuleEditor(CriticalSection *newPlugInLock, 
                                                       FourPoleFilterAudioModule* newFourPoleFilterAudioModule) 
: AudioModuleEditor(newPlugInLock, newFourPoleFilterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newFourPoleFilterAudioModule != NULL ); // you must pass a valid module here
  fourPoleFilterModuleToEdit = newFourPoleFilterAudioModule;

  addWidget( modeComboBox = new FourPoleFilterModeComboBox(juce::String(T("ModeComboBox"))) );
  modeComboBox->assignParameter( moduleToEdit->getParameterByName(T("Mode")) );
  modeComboBox->setDescription(T("Choose the mode of the filter"));
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Characteristic frequency of the filter")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( gainSlider = new RSlider (T("GainSlider")) );
  gainSlider->assignParameter( moduleToEdit->getParameterByName("Gain") );
  gainSlider->setDescription(juce::String(T("Gain of the filter")));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( bandwidthSlider = new RSlider (T("BandwidthSlider")) );
  bandwidthSlider->assignParameter( moduleToEdit->getParameterByName("Bandwidth") );
  bandwidthSlider->setDescription(juce::String(T("Bandwidth of the filter")));
  bandwidthSlider->setDescriptionField(infoField);
  bandwidthSlider->setStringConversionFunction(&rojue::octavesToStringWithUnit2);

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void FourPoleFilterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  modeComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  gainSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  bandwidthSlider->setBounds(x+4, y+4, w-8, 16);
}

void FourPoleFilterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updateWidgetEnablement();
}

void FourPoleFilterModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);

  if( fourPoleFilterModuleToEdit == NULL )
    return;
  if( fourPoleFilterModuleToEdit->wrappedFourPoleFilter == NULL )
    return;
  rosic::FourPoleFilter* core = fourPoleFilterModuleToEdit->wrappedFourPoleFilter;

  //frequencySlider->setEnabled(core->doesModeSupportFrequency());
  //gainSlider->setEnabled(     core->doesModeSupportGain());
  //bandwidthSlider->setEnabled(core->doesModeSupportBandwidth());
}

//-------------------------------------------------------------------------------------------------
// LadderFilter:

LadderFilterAudioModule::LadderFilterAudioModule(CriticalSection *newPlugInLock, rosic::LadderFilter *newLadderFilterToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newLadderFilterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedLadderFilter = newLadderFilterToWrap;
  moduleName  = juce::String(T("LadderFilter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("LadderFilterPresets")) );
  createStaticParameters();
}

void LadderFilterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, juce::String(T("Mode")), 0.0, 1.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Lowpass")));
  p->addStringValue(juce::String(T("Highpass")));
  p->addStringValue(juce::String(T("Bandpass 6+12")));
  p->addStringValue(juce::String(T("Bandpass 6+18")));
  p->addStringValue(juce::String(T("Bandpass 12+12")));
  p->addStringValue(juce::String(T("Bandpass 18+6")));
  p->addStringValue(juce::String(T("Bandpass 12+6")));
  p->addStringValue(juce::String(T("Morph Low/Band/High")));
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setMode);

  p = new AutomatableParameter(plugInLock, juce::String(T("Frequency")), 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setCutoff);

  p = new AutomatableParameter(plugInLock, juce::String(T("Resonance")), 0.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setResonanceInPercent);

  p = new AutomatableParameter(plugInLock, juce::String(T("MakeUp")), 0.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setMakeUp);

  p = new AutomatableParameter(plugInLock, juce::String(T("Drive")), -24.0, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setDrive);

  p = new AutomatableParameter(plugInLock, juce::String(T("Order")), 0.0, 4.0, 1.0, 4.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setOutputStage);

  p = new AutomatableParameter(plugInLock, juce::String(T("Morph")), 0.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedLadderFilter, &LadderFilter::setMorph);

  /*
  p = new Parameter("Gain", -48, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new Parameter("Bandwidth", 0.2, 6.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  */

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

LadderFilterModuleEditor::LadderFilterModuleEditor(CriticalSection *newPlugInLock, LadderFilterAudioModule* newLadderFilterAudioModule) 
: AudioModuleEditor(newPlugInLock, newLadderFilterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newLadderFilterAudioModule != NULL ); // you must pass a valid module here
  ladderFilterModuleToEdit = newLadderFilterAudioModule;

  addWidget( modeComboBox = new RComboBox(juce::String(T("ModeComboBox"))) );
  modeComboBox->assignParameter( moduleToEdit->getParameterByName(T("Mode")) );
  modeComboBox->setDescription(T("Choose the mode of the filter"));
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Characteristic frequency of the filter")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( resonanceSlider = new RSlider (T("ResonanceSlider")) );
  resonanceSlider->assignParameter( moduleToEdit->getParameterByName("Resonance") );
  resonanceSlider->setDescription(juce::String(T("Resonance of the filter")));
  resonanceSlider->setDescriptionField(infoField);
  resonanceSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( makeUpSlider = new RSlider (T("MakeUpSlider")) );
  makeUpSlider->assignParameter( moduleToEdit->getParameterByName("MakeUp") );
  makeUpSlider->setDescription(juce::String(T("Make-up gain to compensate low frequency loss at high resonance")));
  makeUpSlider->setDescriptionField(infoField);
  makeUpSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( driveSlider = new RSlider (T("DriveSlider")) );
  driveSlider->assignParameter( moduleToEdit->getParameterByName("Drive") );
  driveSlider->setDescription(juce::String(T("Drive the filter into distortion")));
  driveSlider->setDescriptionField(infoField);
  driveSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( orderSlider = new RSlider (T("OrderSlider")) );
  orderSlider->assignParameter( moduleToEdit->getParameterByName("Order") );
  orderSlider->setDescription(juce::String(T("Order of the filter")));
  orderSlider->setDescriptionField(infoField);
  orderSlider->setStringConversionFunction(&rojue::valueToString0);


  addWidget( morphSlider = new RSlider (T("MorphSlider")) );
  morphSlider->assignParameter( moduleToEdit->getParameterByName("Morph") );
  morphSlider->setDescription(juce::String(T("Morph between highpass through bandpass to lowpass")));
  morphSlider->setDescriptionField(infoField);
  morphSlider->setStringConversionFunction(&rojue::valueToString2);

  /*
  addWidget( bandwidthSlider = new RSlider (T("BandwidthSlider")) );
  bandwidthSlider->assignParameter( moduleToEdit->getParameterByName("Bandwidth") );
  bandwidthSlider->setDescription(juce::String(T("Bandwidth of the filter")));
  bandwidthSlider->setDescriptionField(infoField);
  bandwidthSlider->setStringConversionFunction(&rojue::octavesToStringWithUnit2);
  */

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void LadderFilterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  modeComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  resonanceSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  makeUpSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  driveSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  orderSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  morphSlider->setBounds(x+4, y+4, w-8, 16);
}

void LadderFilterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updateWidgetEnablement();
}

void LadderFilterModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);
  if( ladderFilterModuleToEdit == NULL )
    return;
  if( ladderFilterModuleToEdit->wrappedLadderFilter == NULL )
    return;
  rosic::LadderFilter* core = ladderFilterModuleToEdit->wrappedLadderFilter;


  //frequencySlider->setEnabled(core->doesModeSupportFrequency());
  //gainSlider->setEnabled(     core->doesModeSupportGain());
  //bandwidthSlider->setEnabled(core->doesModeSupportBandwidth());
}

//-------------------------------------------------------------------------------------------------
// SlopeFilter:

SlopeFilterAudioModule::SlopeFilterAudioModule(CriticalSection *newPlugInLock, rosic::SlopeFilter *newSlopeFilterToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);
  jassert( newSlopeFilterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedSlopeFilter = newSlopeFilterToWrap;
  moduleName  = juce::String(T("SlopeFilter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString + juce::String(T("SlopeFilterPresets")) );
  createStaticParameters();
}

void SlopeFilterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Slope", -12.0, 12.0, 0.01, 0.0, Parameter::LINEAR);
  p->setValueChangeCallback(wrappedSlopeFilter, &SlopeFilter::setSlope);
  addObservedParameter(p); 

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

SlopeFilterModuleEditor::SlopeFilterModuleEditor(CriticalSection *newPlugInLock, SlopeFilterAudioModule* newSlopeFilterAudioModule) 
: AudioModuleEditor(newPlugInLock, newSlopeFilterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newSlopeFilterAudioModule != NULL ); // you must pass a valid module here
  slopeFilterModuleToEdit = newSlopeFilterAudioModule;

  addWidget( slopeSlider = new RSlider (T("SlopeSlider")) );
  slopeSlider->assignParameter( moduleToEdit->getParameterByName("Slope") );
  slopeSlider->setDescription(juce::String(T("Slope of the filter")));
  slopeSlider->setDescriptionField(infoField);
  slopeSlider->setStringConversionFunction(&rojue::decibelsPerOctaveToString2);

  updateWidgetsAccordingToState();
}

void SlopeFilterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom(); 
  slopeSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// TwoPoleFilter:

TwoPoleFilterAudioModule::TwoPoleFilterAudioModule(CriticalSection *newPlugInLock, rosic::TwoPoleFilter *newTwoPoleFilterToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newTwoPoleFilterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedTwoPoleFilter = newTwoPoleFilterToWrap;
  moduleName  = juce::String(T("TwoPoleFilter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("TwoPoleFilterPresets")) );
  createStaticParameters();
}

void TwoPoleFilterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new ParameterTwoPoleFilterMode(plugInLock);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedTwoPoleFilter, &TwoPoleFilter::setMode);

  p = new AutomatableParameter(plugInLock, "Frequency", 2.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedTwoPoleFilter, &TwoPoleFilter::setFrequency);

  p = new AutomatableParameter(plugInLock, "Gain", -48, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedTwoPoleFilter, &TwoPoleFilter::setGain);

  p = new AutomatableParameter(plugInLock, "Bandwidth", 0.2, 6.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedTwoPoleFilter, &TwoPoleFilter::setBandwidth);

  p = new AutomatableParameter(plugInLock, "Radius", -4.0, 4.0, 0.00001, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedTwoPoleFilter, &TwoPoleFilter::setRadius);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

TwoPoleFilterModuleEditor::TwoPoleFilterModuleEditor(CriticalSection *newPlugInLock, TwoPoleFilterAudioModule* newTwoPoleFilterAudioModule) 
: AudioModuleEditor(newPlugInLock, newTwoPoleFilterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newTwoPoleFilterAudioModule != NULL ); // you must pass a valid module here
  twoPoleFilterModuleToEdit = newTwoPoleFilterAudioModule;

  addWidget( modeComboBox = new RComboBox(juce::String(T("ModeComboBox"))) );
  modeComboBox->assignParameter( moduleToEdit->getParameterByName(T("Mode")) );
  modeComboBox->setDescription(T("Choose the mode of the filter"));
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Characteristic frequency of the filter")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( gainSlider = new RSlider (T("GainSlider")) );
  gainSlider->assignParameter( moduleToEdit->getParameterByName("Gain") );
  gainSlider->setDescription(juce::String(T("Gain of the filter")));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( bandwidthSlider = new RSlider (T("BandwidthSlider")) );
  bandwidthSlider->assignParameter( moduleToEdit->getParameterByName("Bandwidth") );
  bandwidthSlider->setDescription(juce::String(T("Bandwidth of the filter")));
  bandwidthSlider->setDescriptionField(infoField);
  bandwidthSlider->setStringConversionFunction(&rojue::octavesToStringWithUnit2);

  addWidget( radiusSlider = new RSlider (T("RadiusSlider")) );
  radiusSlider->assignParameter( moduleToEdit->getParameterByName("Radius") );
  radiusSlider->setDescription(juce::String(T("Radius of the pole or zero")));
  radiusSlider->setDescriptionField(infoField);
  radiusSlider->setStringConversionFunction(&rojue::valueToString5);
  radiusSlider->setVisible(false);

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void TwoPoleFilterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  modeComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  gainSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  bandwidthSlider->setBounds(x+4, y+4, w-8, 16);
  radiusSlider->setBounds(x+4, y+4, w-8, 16);
}

void TwoPoleFilterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  updateWidgetEnablement();
}

void TwoPoleFilterModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);

  if( twoPoleFilterModuleToEdit == NULL )
    return;
  if( twoPoleFilterModuleToEdit->wrappedTwoPoleFilter == NULL )
    return;
  rosic::TwoPoleFilter* core = twoPoleFilterModuleToEdit->wrappedTwoPoleFilter;

  frequencySlider->setEnabled(core->doesModeSupportFrequency());
  gainSlider->setEnabled(     core->doesModeSupportGain());
  bandwidthSlider->setEnabled(core->doesModeSupportBandwidth());
  if( core->doesModeSupportRadius() )
    radiusSlider->setVisible(true);
  else
    radiusSlider->setVisible(false);
}

//=================================================================================================
// Delay Effects:

//-------------------------------------------------------------------------------------------------
// PingPongEcho:

PingPongEchoAudioModule::PingPongEchoAudioModule(CriticalSection *newPlugInLock, rosic::PingPongEcho *newPingPongEchoToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newPingPongEchoToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedPingPongEcho = newPingPongEchoToWrap;
  moduleName  = juce::String(T("PingPongEcho"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("PingPongEchoPresets")) );
  createStaticParameters();
}

void PingPongEchoAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "DelayTime", 0.125, 1.0, 0.0125, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setDelayTime);

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, "Feedback", -100.0, 100.0, 0.1, 50.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setFeedbackInPercent);

  p = new AutomatableParameter(plugInLock, "Pan", -1.0, 1.0, 0.01, -0.4, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setPan);

  p = new AutomatableParameter(plugInLock, "HighDamp", 20.0, 20000.0, 0.0, 4000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setHighDamp);

  p = new AutomatableParameter(plugInLock, "LowDamp", 20.0, 20000.0, 0.0, 250.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setLowDamp);

  p = new AutomatableParameter(plugInLock, "PingPong", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setPingPongMode);

  p = new AutomatableParameter(plugInLock, "TrueStereo", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setTrueStereoMode);

  p = new AutomatableParameter(plugInLock, "TempoSync", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setSyncMode);

  p = new AutomatableParameter(plugInLock, "Activated", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedPingPongEcho, &PingPongEcho::setActive);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

PingPongEchoModuleEditor::PingPongEchoModuleEditor(CriticalSection *newPlugInLock, PingPongEchoAudioModule* newPingPongEchoAudioModule) 
: AudioModuleEditor(newPlugInLock, newPingPongEchoAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newPingPongEchoAudioModule != NULL ); // you must pass a valid module here

  addWidget( delayTimeSlider = new RSlider (T("DelayTimeSlider")) );
  delayTimeSlider->assignParameter( moduleToEdit->getParameterByName("DelayTime") );
  delayTimeSlider->setDescription(juce::String(T("Delay time in beats")));
  delayTimeSlider->setDescriptionField(infoField);
  delayTimeSlider->setStringConversionFunction(&rojue::beatsToStringWithUnit4);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( moduleToEdit->getParameterByName("Feedback") );
  feedbackSlider->setDescription(juce::String(T("Amount of feedback in percent")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( panSlider = new RSlider (T("PanSlider")) );
  panSlider->assignParameter( moduleToEdit->getParameterByName("Pan") );
  panSlider->setDescription(juce::String(T("Panorama position of the first echo")));
  panSlider->setDescriptionField(infoField);
  panSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( highDampSlider = new RSlider (T("HighDampSlider")) );
  highDampSlider->assignParameter( moduleToEdit->getParameterByName("HighDamp") );
  highDampSlider->setDescription(juce::String(T("Cutoff frequency for the high damping (lowpass) filter")));
  highDampSlider->setDescriptionField(infoField);
  highDampSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( lowDampSlider = new RSlider (T("LowDampSlider")) );
  lowDampSlider->assignParameter( moduleToEdit->getParameterByName("LowDamp") );
  lowDampSlider->setDescription(juce::String(T("Cutoff frequency for the low damping (highpass) filter")));
  lowDampSlider->setDescriptionField(infoField);
  lowDampSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( pingPongButton = new RButton(juce::String(T("PingPong"))) );
  pingPongButton->assignParameter( moduleToEdit->getParameterByName(T("PingPong")) );
  pingPongButton->setDescription(juce::String(T("Toggle ping-pong mode (alternating pan-positions) on/off")));
  pingPongButton->setDescriptionField(infoField);
  pingPongButton->setClickingTogglesState(true);

  addWidget( trueStereoButton = new RButton(juce::String(T("TrueStereo"))) );
  trueStereoButton->assignParameter( moduleToEdit->getParameterByName(T("TrueStereo")) );
  trueStereoButton->setDescription(juce::String(T("Toggle true-stereo mode on/off")));
  trueStereoButton->setDescriptionField(infoField);
  trueStereoButton->setClickingTogglesState(true);

  addWidget( tempoSyncButton = new RButton(juce::String(T("TempoSync"))) );
  tempoSyncButton->assignParameter( moduleToEdit->getParameterByName(T("TempoSync")) );
  tempoSyncButton->setButtonText(juce::String(T("Sync")));
  tempoSyncButton->setDescription(juce::String(T("Toggle tempo synchronization on/off")));
  tempoSyncButton->setDescriptionField(infoField);
  tempoSyncButton->setClickingTogglesState(true);
  tempoSyncButton->addRButtonListener(this);

  updateWidgetsAccordingToState();
}

void PingPongEchoModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  ScopedLock scopedLock(*plugInLock);
  if( buttonThatWasClicked == tempoSyncButton )
  {
    if( tempoSyncButton->getToggleState() == true )
      delayTimeSlider->setStringConversionFunction(&rojue::beatsToStringWithUnit4);
    else
      delayTimeSlider->setStringConversionFunction(&rojue::secondsToStringWithUnit4);
  }
}

void PingPongEchoModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth()/2;
  int h = getHeight();
  y = getPresetSectionBottom();
  delayTimeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowDampSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  highDampSlider->setBounds(x+4, y+4, w-8, 16);

  y = getPresetSectionBottom();
  dryWetSlider->setBounds(x+w+4, y+4, w-8, 16);
  y += 20;  
  panSlider->setBounds(x+w+4, y+4, w-8, 16);
  y += 20;  
  pingPongButton->setBounds(x+w+4, y+4, 64, 16);
  x = pingPongButton->getRight();
  trueStereoButton->setBounds(x+4, y+4, 72, 16);
}

//-------------------------------------------------------------------------------------------------
// Reverb:

ReverbAudioModule::ReverbAudioModule(CriticalSection *newPlugInLock, rosic::Reverb *newReverbToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newReverbToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedReverb = newReverbToWrap;
  moduleName  = juce::String(T("Reverb"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("ReverbPresets")) );
  createStaticParameters();
}

void ReverbAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, "FirstEcho", 10.0, 200.0, 0.1, 50.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setReferenceDelayTime);

  p = new AutomatableParameter(plugInLock, "PreDelay", 0.0, 250.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setPreDelay);

  p = new AutomatableParameter(plugInLock, "DecayTime", 0.1, 10.0, 0.01, 3.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setMidReverbTime);

  p = new AutomatableParameter(plugInLock, "HighDecayScale", 0.1, 10.0, 0.01, 0.3, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setHighReverbTimeScale);

  p = new AutomatableParameter(plugInLock, "LowDecayScale", 0.1, 10.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setLowReverbTimeScale);

  p = new AutomatableParameter(plugInLock, "HighCrossoverFrequency", 20.0, 20000.0, 0.0, 4000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setHighCrossoverFreq);

  p = new AutomatableParameter(plugInLock, "LowCrossoverFrequency", 20.0, 20000.0, 0.0, 250.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedReverb, &Reverb::setLowCrossoverFreq);

  p = new AutomatableParameter(plugInLock, "Pinking", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedReverb, &Reverb::setWetPinkingSwitch);

  p = new AutomatableParameter(plugInLock, "StereoSwap", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedReverb, &Reverb::setStereoSwapSwitch);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

ReverbModuleEditor::ReverbModuleEditor(CriticalSection *newPlugInLock, ReverbAudioModule* newReverbAudioModule) 
: AudioModuleEditor(newPlugInLock, newReverbAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newReverbAudioModule != NULL ); // you must pass a valid module here

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( firstEchoSlider = new RSlider (T("FirstEchoSlider")) );
  firstEchoSlider->assignParameter( moduleToEdit->getParameterByName("FirstEcho") );
  firstEchoSlider->setDescription(juce::String(T("Arrival time of the first reflection in seconds")));
  firstEchoSlider->setDescriptionField(infoField);
  firstEchoSlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( preDelaySlider = new RSlider (T("PreDelaySlider")) );
  preDelaySlider->assignParameter( moduleToEdit->getParameterByName("PreDelay") );
  preDelaySlider->setDescription(juce::String(T("Initial delay for the reverberation.")));
  preDelaySlider->setDescriptionField(infoField);
  preDelaySlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( decayTimeSlider = new RSlider (T("DecayTimeSlider")) );
  decayTimeSlider->assignParameter( moduleToEdit->getParameterByName("DecayTime") );
  decayTimeSlider->setDescription(juce::String(T("Time for the tail to decay to -60 dB")));
  decayTimeSlider->setDescriptionField(infoField);
  decayTimeSlider->setStringConversionFunction(&rojue::secondsToStringWithUnit2);

  addWidget( highDecayScaleSlider = new RSlider (T("HighDecayScaleSlider")) );
  highDecayScaleSlider->assignParameter( moduleToEdit->getParameterByName("HighDecayScale") );
  highDecayScaleSlider->setDescription(juce::String(T("Scaler for the decay-time at high frequencies")));
  highDecayScaleSlider->setDescriptionField(infoField);
  highDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( lowDecayScaleSlider = new RSlider (T("LowDecayScaleSlider")) );
  lowDecayScaleSlider->assignParameter( moduleToEdit->getParameterByName("LowDecayScale") );
  lowDecayScaleSlider->setDescription(juce::String(T("Scaler for the decay-time at low frequencies")));
  lowDecayScaleSlider->setDescriptionField(infoField);
  lowDecayScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( highFreqSlider = new RSlider (T("HighFreqSlider")) );
  highFreqSlider->assignParameter( moduleToEdit->getParameterByName("HighCrossoverFrequency") );
  highFreqSlider->setDescription(juce::String(T("Crossover frequency between mid and high frequencies")));
  highFreqSlider->setSliderName(juce::String(T("HighFreq")));
  highFreqSlider->setDescriptionField(infoField);
  highFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( lowFreqSlider = new RSlider (T("LowFreqSlider")) );
  lowFreqSlider->assignParameter( moduleToEdit->getParameterByName("LowCrossoverFrequency") );
  lowFreqSlider->setSliderName(juce::String(T("LowFreq")));
  lowFreqSlider->setDescription(juce::String(T("Crossover frequency between low and mid frequencies")));
  lowFreqSlider->setDescriptionField(infoField);
  lowFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( pinkButton = new RButton(juce::String(T("Pink"))) );
  pinkButton->assignParameter( moduleToEdit->getParameterByName(T("Pinking")) );
  pinkButton->setDescription(juce::String(T("Toggle pinking filter (3 dB/oct lowpass) for wet signal on/off")));
  pinkButton->setDescriptionField(infoField);
  pinkButton->setClickingTogglesState(true);

  addWidget( stereoSwapButton = new RButton(juce::String(T("StereoSwap"))) );
  stereoSwapButton->assignParameter( moduleToEdit->getParameterByName(T("StereoSwap")) );
  stereoSwapButton->setDescription(juce::String(T("Swap stereo channels (for wet signal) left for right")));
  stereoSwapButton->setButtonText(juce::String(T("Swap")));
  stereoSwapButton->setDescriptionField(infoField);
  stereoSwapButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void ReverbModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth()/2;
  int h = getHeight();
  y = getPresetSectionBottom();

  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  firstEchoSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  preDelaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  decayTimeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  highDecayScaleSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowDecayScaleSlider->setBounds(x+4, y+4, w-8, 16);

  y  = getPresetSectionBottom();
  x += w;
  pinkButton->setBounds(      x+4,     y+4, w/2-8, 16);
  stereoSwapButton->setBounds(x+w/2+4, y+4, w/2-8, 16);
  y += 20;  

  y += 20;  
  y += 20;  
  highFreqSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowFreqSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// SimpleDelay:

SimpleDelayAudioModule::SimpleDelayAudioModule(CriticalSection *newPlugInLock, rosic::FractionalDelayLineStereo *newSimpleDelayToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newSimpleDelayToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedSimpleDelay = newSimpleDelayToWrap;
  moduleName  = juce::String(T("SimpleDelay"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("SimpleDelayPresets")) );
  createStaticParameters();
}

void SimpleDelayAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "DelayTime", 0.0, 200.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSimpleDelay, &FractionalDelayLineStereo::setDelayTimeInMilliseconds);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

SimpleDelayModuleEditor::SimpleDelayModuleEditor(CriticalSection *newPlugInLock, SimpleDelayAudioModule* newSimpleDelayAudioModule) 
: AudioModuleEditor(newPlugInLock, newSimpleDelayAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newSimpleDelayAudioModule != NULL ); // you must pass a valid module here

  addWidget( delaySlider = new RSlider (T("DelaySlider")) );
  delaySlider->assignParameter( moduleToEdit->getParameterByName("DelayTime") );
  delaySlider->setDescription(juce::String(T("Delay in milliseconds")));
  delaySlider->setDescriptionField(infoField);
  delaySlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  updateWidgetsAccordingToState();
}

void SimpleDelayModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  delaySlider->setBounds(x+4, y+4, w-8, 16);
}

//=================================================================================================
// Modulation Effects:

//-------------------------------------------------------------------------------------------------
// ModulationEffect baseclass:

ModulationEffectAudioModule::ModulationEffectAudioModule(CriticalSection *newPlugInLock, 
                                                         rosic::ModulationEffect *newModulationEffectToWrap)
                                                          : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newModulationEffectToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedModulationEffect = newModulationEffectToWrap;
  moduleName  = juce::String(T("ModulationEffect"));
  setActiveDirectory(getApplicationDirectory() 
    + File::separatorString + juce::String(T("ModulationEffectPresets")) );
  lfoModule = new LowFrequencyOscillatorAudioModule(newPlugInLock, &wrappedModulationEffect->lfo);
  lfoModule->setModuleName(juce::String(T("LowFrequencyOscillator")));
  addChildAudioModule(lfoModule);
}

ModulationEffectModuleEditor::ModulationEffectModuleEditor(CriticalSection *newPlugInLock, 
                                                           ModulationEffectAudioModule* newModulationEffectAudioModule) 
: AudioModuleEditor(newPlugInLock, newModulationEffectAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newModulationEffectAudioModule != NULL ); // you must pass a valid module here
  modulationEffectModuleToEdit = newModulationEffectAudioModule;

  lfoEditor = new LowFrequencyOscillatorEditor(plugInLock, modulationEffectModuleToEdit->lfoModule);
  lfoEditor->setHeadlineText(T("LFO"));
  addChildEditor(lfoEditor );

  addWidget( lfoLabel = new RTextField( juce::String(T("LFO"))) );
  lfoLabel->setJustification(Justification::centred);
  lfoLabel->setDescription(T("Parameters for the low frequency oscillator"));
  lfoLabel->setDescriptionField(infoField);

  addWidget( effectLabel = new RTextField( juce::String(T("Effect"))) );
  effectLabel->setJustification(Justification::centred);
  effectLabel->setDescription(T("Parameters for the actual effect"));
  effectLabel->setDescriptionField(infoField);

  updateWidgetsAccordingToState();
}

/*
Rectangle ModulationEffectModuleEditor::getLfoEditButtonBounds() const
{
  Rectangle r = lfoEditor->getEditButtonBounds();
  r.setPosition(r.getX()+lfoEditor->getX(), r.getY()+lfoEditor->getY());
  return r;
}
*/

void ModulationEffectModuleEditor::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::updateWidgetsAccordingToState();
  lfoEditor->updateWidgetsAccordingToState();
}

void ModulationEffectModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth()/2;
  int h = getHeight()-y;

  lfoRect.setBounds(x, y, w, h);
  effectRect.setBounds(x+w, y, w, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(lfoRect);
  guiLayoutRectangles.add(effectRect);

  lfoEditor->setBounds(lfoRect);

  x = effectRect.getX();
  y = effectRect.getY();
  w = effectRect.getWidth();
  h = effectRect.getHeight();
  effectLabel->setBounds(x+4, y+2, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// Flanger:

FlangerAudioModule::FlangerAudioModule(CriticalSection *newPlugInLock, rosic::Flanger *newFlangerToWrap)
: ModulationEffectAudioModule(newPlugInLock, newFlangerToWrap)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newFlangerToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedFlanger = newFlangerToWrap;
  moduleName  = juce::String(T("Flanger"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("FlangerPresets")) );
  createStaticParameters();
}

void FlangerAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;
  p = dynamic_cast<AutomatableParameter*> (lfoModule->getParameterByName(T("CycleLength")));
  p->setRange(0.25, 16.0);
  p->setDefaultValue(8.0, true);
  p->setScaling(Parameter::EXPONENTIAL);

  p = new AutomatableParameter(plugInLock, "Depth", 0.0, 48.0, 0.1, 12.0, Parameter::LINEAR);  // #04
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFlanger, &Flanger::setDepth);

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFlanger, &Flanger::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, "Frequency", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFlanger, &Flanger::setFrequency);

  p = new AutomatableParameter(plugInLock, "Feedback", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFlanger, &Flanger::setFeedbackInPercent);

  p = new AutomatableParameter(plugInLock, "Invert", 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedFlanger, &Flanger::setNegativePolarity);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

FlangerModuleEditor::FlangerModuleEditor(CriticalSection *newPlugInLock, FlangerAudioModule* newFlangerAudioModule) 
: ModulationEffectModuleEditor(newPlugInLock, newFlangerAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newFlangerAudioModule != NULL ); // you must pass a valid module here
  flangerModuleToEdit = newFlangerAudioModule;

  addWidget( depthSlider = new RSlider (T("DepthSlider")) );
  depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  depthSlider->setDescription(juce::String(T("Modulation depth in semitones")));
  depthSlider->setDescriptionField(infoField);
  depthSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet (filtered) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Frequency of first notch or peak in the comb-filter")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( moduleToEdit->getParameterByName("Feedback") );
  feedbackSlider->setDescription(juce::String(T("Feedback around the delayline")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( invertButton = new RButton(juce::String(T("Invert"))) );
  invertButton->assignParameter( moduleToEdit->getParameterByName(T("Invert")) );
  invertButton->setDescription(juce::String(T("Invert polarity of wet signal")));
  invertButton->setDescriptionField(infoField);
  invertButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void FlangerModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  ModulationEffectModuleEditor::resized();

  int x = effectRect.getX();
  int y = effectLabel->getBottom()-2;
  int w = effectRect.getWidth();
  depthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  dryWetSlider->setBounds(x+4,    y+4,  w-8, 16);
  y += 20;
  invertButton->setBounds(x+w-64, y+4,   60, 16);
  y += 28;
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// Phaser:

PhaserAudioModule::PhaserAudioModule(CriticalSection *newPlugInLock, rosic::Phaser *newPhaserToWrap)
: ModulationEffectAudioModule(newPlugInLock, newPhaserToWrap)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newPhaserToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedPhaser = newPhaserToWrap;
  moduleName  = juce::String(T("Phaser"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("PhaserPresets")) );
  createStaticParameters();
}

void PhaserAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = dynamic_cast<AutomatableParameter*> (lfoModule->getParameterByName(T("CycleLength")));
  p->setRange(0.25, 16.0);
  p->setDefaultValue(8.0, true);
  p->setScaling(Parameter::EXPONENTIAL);

  p = new AutomatableParameter(plugInLock, "Depth", 0.0, 48.0, 0.1, 12.0, Parameter::LINEAR);  // #04
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setDepth);

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, juce::String(T("FilterMode")), 0.0, 1.0, 1.0, 0.0, Parameter::STRING);
  //p->addStringValue(juce::String(T("Bypass")));
  p->addStringValue(juce::String(T("Allpass 1st order")));
  p->addStringValue(juce::String(T("Allpass 2nd order")));
  //p->setValue(0.0, false);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setFilterMode);

  p = new AutomatableParameter(plugInLock, "Frequency", 2.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setFrequency);

  p = new AutomatableParameter(plugInLock, "Q", 0.1, 10.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setQ);

  p = new AutomatableParameter(plugInLock, "Feedback", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setFeedbackInPercent);

  p = new AutomatableParameter(plugInLock, "NumStages", 1.0, 24.0, 1.0, 4.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaser, &Phaser::setNumStages);

  //p = new AutomatableParameter(plugInLock, "Activated", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  //addObservedParameter(p);
  //p->setValueChangeCallback(wrappedPhaser, &Phaser::setActive);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

PhaserModuleEditor::PhaserModuleEditor(CriticalSection *newPlugInLock, PhaserAudioModule* newPhaserAudioModule) 
: ModulationEffectModuleEditor(newPlugInLock, newPhaserAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newPhaserAudioModule != NULL ); // you must pass a valid module here
  phaserModuleToEdit = newPhaserAudioModule;

  addWidget( filterLabel = new RTextField( juce::String(T("Filter:"))) );
  filterLabel->setDescription(T("Parameters for the filter"));
  filterLabel->setDescriptionField(infoField);

  addWidget( depthSlider = new RSlider (T("DepthSlider")) );
  depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  depthSlider->setDescription(juce::String(T("Modulation depth in semitones")));
  depthSlider->setDescriptionField(infoField);
  depthSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet(filtered) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( modeComboBox = new RComboBox(juce::String(T("ModeComboBox"))) );
  modeComboBox->assignParameter( moduleToEdit->getParameterByName(T("FilterMode")) );
  modeComboBox->setDescription(T("Choose the mode of the filter"));
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Characteristic frequency of the filter")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( qSlider = new RSlider (T("QSlider")) );
  qSlider->assignParameter( moduleToEdit->getParameterByName("Q") );
  qSlider->setDescription(juce::String(T("Q of the filter")));
  qSlider->setDescriptionField(infoField);
  qSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( moduleToEdit->getParameterByName("Feedback") );
  feedbackSlider->setDescription(juce::String(T("Feedback around the allpass chain")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( stagesSlider = new RSlider (T("StagesSlider")) );
  stagesSlider->assignParameter( moduleToEdit->getParameterByName("NumStages") );
  stagesSlider->setName(juce::String(T("Stages")));
  stagesSlider->setDescription(juce::String(T("Number of allpass stages in the chain")));
  stagesSlider->setDescriptionField(infoField);
  stagesSlider->setStringConversionFunction(&rojue::valueToString0);

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void PhaserModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  ModulationEffectModuleEditor::resized();

  int x = effectRect.getX();
  int y = effectLabel->getBottom()-2;
  int w = effectRect.getWidth();
  depthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 24;

  filterLabel->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  modeComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  stagesSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14; 
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;  
  qSlider->setBounds(x+4, y+4, w-8, 16);
}

void PhaserModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  if( rComboBoxThatHasChanged == modeComboBox )
    updateWidgetEnablement();
}

void PhaserModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);
  if( phaserModuleToEdit == NULL )
    return;
  if( phaserModuleToEdit->wrappedPhaser == NULL )
    return;
  rosic::Phaser* core = phaserModuleToEdit->wrappedPhaser;
  qSlider->setEnabled(core->doesModeSupportQ());
}

//-------------------------------------------------------------------------------------------------
// Tremolo:

TremoloAudioModule::TremoloAudioModule(CriticalSection *newPlugInLock, rosic::Tremolo *newTremoloToWrap)
: ModulationEffectAudioModule(newPlugInLock, newTremoloToWrap)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newTremoloToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedTremolo = newTremoloToWrap;
  moduleName  = juce::String(T("Tremolo"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("TremoloPresets")) );
  createStaticParameters();
}

void TremoloAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;
  p = dynamic_cast<AutomatableParameter*> (lfoModule->getParameterByName(T("CycleLength")));
  p->setRange(0.125, 4.0);
  p->setUpperAutomationLimit(4.0);

  p = new AutomatableParameter(plugInLock, "Depth", 0.0, 100.0, 1.0, 50.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Tremolo>(wrappedTremolo, &Tremolo::setDepthInPercent);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

TremoloModuleEditor::TremoloModuleEditor(CriticalSection *newPlugInLock, TremoloAudioModule* newTremoloAudioModule) 
: ModulationEffectModuleEditor(newPlugInLock, newTremoloAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newTremoloAudioModule != NULL ); // you must pass a valid module here
  tremoloModuleToEdit = newTremoloAudioModule;

  addWidget( depthSlider = new RSlider (T("DepthSlider")) );
  depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  depthSlider->setDescription(juce::String(T("Modulation depth in percent")));
  depthSlider->setDescriptionField(infoField);
  depthSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  updateWidgetsAccordingToState();
}

void TremoloModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  ModulationEffectModuleEditor::resized();

  int x = effectRect.getX();
  int y = effectLabel->getBottom();
  int w = effectRect.getWidth();
  depthSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// Vibrato:

VibratoAudioModule::VibratoAudioModule(CriticalSection *newPlugInLock, rosic::Vibrato *newVibratoToWrap)
: ModulationEffectAudioModule(newPlugInLock, newVibratoToWrap)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newVibratoToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedVibrato = newVibratoToWrap;
  moduleName  = juce::String(T("Vibrato"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("VibratoPresets")) );
  createStaticParameters();
}

void VibratoAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Depth", 0.0, 2.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Vibrato>(wrappedVibrato, &Vibrato::setDepthInPercent);

  p = new AutomatableParameter(plugInLock, "DryWet", 0.0, 1.0, 0.01, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback<Vibrato>(wrappedVibrato, &Vibrato::setDryWetRatio);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

VibratoModuleEditor::VibratoModuleEditor(CriticalSection *newPlugInLock, VibratoAudioModule* newVibratoAudioModule) 
: ModulationEffectModuleEditor(newPlugInLock, newVibratoAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newVibratoAudioModule != NULL ); // you must pass a valid module here
  vibratoModuleToEdit = newVibratoAudioModule;

  addWidget( depthSlider = new RSlider (T("DepthSlider")) );
  depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  depthSlider->setDescription(juce::String(T("Modulation depth in semitones")));
  depthSlider->setDescriptionField(infoField);
  depthSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit2);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWet") );
  dryWetSlider->setDescription(juce::String(T("Mix ratio between original and vibrato'ed signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  updateWidgetsAccordingToState();
}

void VibratoModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  ModulationEffectModuleEditor::resized();

  int x = effectRect.getX();
  int y = effectLabel->getBottom();
  int w = effectRect.getWidth();
  depthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// WahWah:

WahWahAudioModule::WahWahAudioModule(CriticalSection *newPlugInLock, rosic::WahWah *newWahWahToWrap)
: ModulationEffectAudioModule(newPlugInLock, newWahWahToWrap)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newWahWahToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedWahWah = newWahWahToWrap;
  moduleName  = juce::String(T("WahWah"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("WahWahPresets")) );
  createStaticParameters();
}

void WahWahAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = dynamic_cast<AutomatableParameter*> (lfoModule->getParameterByName(T("CycleLength")));
  p->setRange(0.25, 8.0);
  p->setDefaultValue(0.5, true);
  p->setScaling(Parameter::EXPONENTIAL);

  p = new AutomatableParameter(plugInLock, "Depth", 0.0, 48.0, 0.1, 12.0, Parameter::LINEAR);  // #04
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWahWah, &WahWah::setDepth);

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWahWah, &WahWah::setDryWetRatio);

  p = new ParameterTwoPoleFilterMode(plugInLock);
  p->setDefaultValue(rosic::TwoPoleFilter::BANDPASS, true);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedWahWah, &WahWah::setFilterMode);

  p = new AutomatableParameter(plugInLock, "Frequency", 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWahWah, &WahWah::setFrequency);

  p = new AutomatableParameter(plugInLock, "Gain", -48, 48.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWahWah, &WahWah::setGain);

  p = new AutomatableParameter(plugInLock, "Bandwidth", 0.25, 2.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedWahWah, &WahWah::setBandwidth);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

WahWahModuleEditor::WahWahModuleEditor(CriticalSection *newPlugInLock, WahWahAudioModule* newWahWahAudioModule) 
: ModulationEffectModuleEditor(newPlugInLock, newWahWahAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newWahWahAudioModule != NULL ); // you must pass a valid module here
  wahWahModuleToEdit = newWahWahAudioModule;

  addWidget( filterLabel = new RTextField( juce::String(T("Filter:"))) );
  filterLabel->setDescription(T("Parameters for the filter"));
  filterLabel->setDescriptionField(infoField);

  addWidget( depthSlider = new RSlider (T("DepthSlider")) );
  depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  depthSlider->setDescription(juce::String(T("Modulation depth in semitones")));
  depthSlider->setDescriptionField(infoField);
  depthSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry and wet(filtered) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( modeComboBox = new RComboBox(juce::String(T("ModeComboBox"))) );
  modeComboBox->assignParameter( moduleToEdit->getParameterByName(T("Mode")) );
  modeComboBox->setDescription(T("Choose the mode of the filter"));
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update enablement of the sliders

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Characteristic frequency of the filter")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( gainSlider = new RSlider (T("GainSlider")) );
  gainSlider->assignParameter( moduleToEdit->getParameterByName("Gain") );
  gainSlider->setDescription(juce::String(T("Gain of the filter")));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( bandwidthSlider = new RSlider (T("BandwidthSlider")) );
  bandwidthSlider->assignParameter( moduleToEdit->getParameterByName("Bandwidth") );
  bandwidthSlider->setDescription(juce::String(T("Bandwidth of the filter")));
  bandwidthSlider->setDescriptionField(infoField);
  bandwidthSlider->setStringConversionFunction(&rojue::octavesToStringWithUnit2);

  updateWidgetsAccordingToState();
  updateWidgetEnablement();
}

void WahWahModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  ModulationEffectModuleEditor::resized();

  int x = effectRect.getX();
  int y = effectLabel->getBottom();
  int w = effectRect.getWidth();
  depthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 28;
  filterLabel->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  modeComboBox->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  gainSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  bandwidthSlider->setBounds(x+4, y+4, w-8, 16);
}

void WahWahModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);
  if( rComboBoxThatHasChanged == modeComboBox )
    updateWidgetEnablement();
}

void WahWahModuleEditor::updateWidgetEnablement()
{
  ScopedLock scopedLock(*plugInLock);

  if( wahWahModuleToEdit == NULL )
    return;
  if( wahWahModuleToEdit->wrappedWahWah == NULL )
    return;
  rosic::WahWah* core = wahWahModuleToEdit->wrappedWahWah;

  frequencySlider->setEnabled(core->doesModeSupportFrequency());
  gainSlider->setEnabled(     core->doesModeSupportGain());
  bandwidthSlider->setEnabled(core->doesModeSupportBandwidth());
}


//=================================================================================================
// Spectral Effects:

//-------------------------------------------------------------------------------------------------
// FormantShifter:

FormantShifterAudioModule::FormantShifterAudioModule(CriticalSection *newPlugInLock, rosic::FormantShifterStereo *newFormantShifterToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newFormantShifterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedFormantShifter = newFormantShifterToWrap;
  moduleName  = juce::String(T("FormantShifter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("FormantShifterPresets")) );
  createStaticParameters();
}

void FormantShifterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  //p = new AutomatableParameter(plugInLock, T("BlockSize"), 32.0, 8192.0, 0.0, 2048.0, Parameter::EXPONENTIAL);
  //addObservedParameter(p); 
  //p->setValueChangeCallback(wrappedFormantShifter, &FormantShifterStereo::setBlockSize);

  p = new AutomatableParameter(plugInLock, T("FormantScale"), 0.25, 4.0, 0.01, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFormantShifter, &FormantShifterStereo::setFormantScale);

  p = new AutomatableParameter(plugInLock, T("FormantOffset"), -200.0, 200.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFormantShifter, &FormantShifterStereo::setFormantOffset);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFormantShifter, &FormantShifterStereo::setDryWetRatio);

  //p = new AutomatableParameter(plugInLock, T("Mono"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  //addObservedParameter(p);
  //p->setValueChangeCallback(wrappedFormantShifter, &FormantShifterStereo::setMono);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

FormantShifterModuleEditor::FormantShifterModuleEditor(CriticalSection *newPlugInLock, 
                                                       FormantShifterAudioModule* newFormantShifterAudioModule) 
: AudioModuleEditor(newPlugInLock, newFormantShifterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newFormantShifterAudioModule != NULL ); // you must pass a valid module here

  addWidget( formantScaleSlider = new RSlider (T("FormantScaleSlider")) );
  formantScaleSlider->assignParameter( moduleToEdit->getParameterByName(T("FormantScale")) );
  formantScaleSlider->setDescription(juce::String(T("Formant scale factor")));
  formantScaleSlider->setDescriptionField(infoField);
  formantScaleSlider->setStringConversionFunction(&rojue::valueToString2);

  addWidget( formantOffsetSlider = new RSlider (T("FormantOffsetSlider")) );
  formantOffsetSlider->assignParameter( moduleToEdit->getParameterByName(T("FormantOffset")) );
  formantOffsetSlider->setDescription(juce::String(T("Formant offset in Hz")));
  formantOffsetSlider->setDescriptionField(infoField);
  formantOffsetSlider->setStringConversionFunction(&rojue::hertzToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName(T("DryWetRatio")) );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry (original) and wet (ringmodulated) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  /*
  addWidget( monoButton = new RButton(juce::String(T("Mono"))) );
  monoButton->assignParameter( moduleToEdit->getParameterByName(T("Mono")) );
  monoButton->setDescription(juce::String(T("Toggle mono-mode (saves CPU power)")));
  monoButton->setDescriptionField(infoField);
  monoButton->setClickingTogglesState(true);
  */

  updateWidgetsAccordingToState();
}

void FormantShifterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  formantScaleSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  formantOffsetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  //antiAliasButton->setBounds(x+4, y+4, 60, 16);
  //y += 20;  
}

//=================================================================================================
// Other Effects:

//-------------------------------------------------------------------------------------------------
// Chorus:

ChorusAudioModule::ChorusAudioModule(CriticalSection *newPlugInLock, rosic::Chorus *newChorusToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newChorusToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedChorus = newChorusToWrap;
  moduleName  = juce::String(T("Chorus"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("ChorusPresets")) );
  createStaticParameters();
}

void ChorusAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*plugInLock);

  // \todo get rid of this function and implement the new callback mechanism - maybe we need a class IndexedParameter or something

  if( wrappedChorus == NULL )
    return;
  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case   0: wrappedChorus->setAverageDelayTime(               value);  break;
  case   1: wrappedChorus->setCycleLength(                    value);  break;
  case   2: wrappedChorus->setDepth(                          value);  break;
  case   3: wrappedChorus->setStereoPhaseOffsetInDegrees(     value);  break;
  case   4: wrappedChorus->setGlobalFeedback(            0.01*value);  break;
  //case   5: wrappedChorus->setCrossMixFactor(                0.01*value);  break;
  //case   6: wrappedChorus->setFeedbackFactor2(               0.01*value);  break;
  case   7: wrappedChorus->setDryWetRatio(                    value       );  break;

  case   8: wrappedChorus->setVoiceDelayScale(0,         0.01*value);  break;
  case   9: wrappedChorus->setVoiceDepthScale(0,         0.01*value);  break;
  case  10: wrappedChorus->setVoiceAmpScale(  0,         0.01*value);  break;
  case  11: wrappedChorus->activateVoice(     0,      value >= 0.5f);  break;

  case  12: wrappedChorus->setVoiceDelayScale(1,         0.01*value);  break;
  case  13: wrappedChorus->setVoiceDepthScale(1,         0.01*value);  break;
  case  14: wrappedChorus->setVoiceAmpScale(  1,         0.01*value);  break;
  case  15: wrappedChorus->activateVoice(     1,      value >= 0.5f);  break;

  case  16: wrappedChorus->setVoiceDelayScale(2,         0.01*value);  break;
  case  17: wrappedChorus->setVoiceDepthScale(2,         0.01*value);  break;
  case  18: wrappedChorus->setVoiceAmpScale(  2,         0.01*value);  break;
  case  19: wrappedChorus->activateVoice(     2,      value >= 0.5f);  break;

  case  20: wrappedChorus->setVoiceDelayScale(3,         0.01*value);  break;
  case  21: wrappedChorus->setVoiceDepthScale(3,         0.01*value);  break;
  case  22: wrappedChorus->setVoiceAmpScale(  3,         0.01*value);  break;
  case  23: wrappedChorus->activateVoice(     3,      value >= 0.5f);  break;
  } 
}

void ChorusAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  // ...here we either need a callback with two parameters, or we must stick to the old callback mechansim - for now, we'll do the latter

  p = new AutomatableParameter(plugInLock, "Delay", 10.0, 50.0, 1.0, 20.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "CycleLength", 0.5, 8.0, 0.25, 4.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "Depth", 0.0, 2.0, 0.01, 0.25, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "StereoPhase", 0.0, 180.0, 1.0, 90.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "GlobalFeedback", -99.0, 99.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "CrossMix", -100.0, 100.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "FeedbackPostCrossMix", -99.0, 99.0, 1.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 

  p = new AutomatableParameter(plugInLock, "DelayScaleVoice1", 1.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "DepthScaleVoice1", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "AmpScaleVoice1", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "OnOffVoice1", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "DelayScaleVoice2", 1.0, 100.0, 1.0, 83.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "DepthScaleVoice2", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "AmpScaleVoice2", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "OnOffVoice2", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "DelayScaleVoice3", 1.0, 100.0, 1.0, 69.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "DepthScaleVoice3", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "AmpScaleVoice3", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "OnOffVoice3", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  p = new AutomatableParameter(plugInLock, "DelayScaleVoice4", 1.0, 100.0, 1.0, 58.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "DepthScaleVoice4", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "AmpScaleVoice4", -100.0, 100.0, 1.0, 100.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p = new AutomatableParameter(plugInLock, "OnOffVoice4", 0.0, 1.0, 1.0, 1.0, Parameter::BOOLEAN);
  addObservedParameter(p);

  for(int i=0; i < (int) observedParameters.size(); i++ )
  {   
    parameterChanged(observedParameters[i]); // because some parameters use the old callback mechanism
    observedParameters[i]->resetToDefaultValue(true, true);
  }
}

ChorusModuleEditor::ChorusModuleEditor(CriticalSection *newPlugInLock, ChorusAudioModule* newChorusAudioModule) 
: AudioModuleEditor(newPlugInLock, newChorusAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newChorusAudioModule != NULL ); // you must pass a valid module here

  addWidget( globalLabel = new RTextField( juce::String(T("Global"))) );
  globalLabel->setJustification(Justification::centred);
  globalLabel->setDescription(T("Global parameters for the chorus effect"));
  globalLabel->setDescriptionField(infoField);

  addWidget( delaySlider = new RSlider (T("DelaySlider")) );
  delaySlider->assignParameter( moduleToEdit->getParameterByName("Delay") );
  delaySlider->setDescription(juce::String(T("Average delay in the delaylines")));
  delaySlider->setDescriptionField(infoField);
  delaySlider->setStringConversionFunction(&rojue::millisecondsToStringWithUnit2);

  addWidget( cycleLengthSlider = new RSlider (T("CycleLengthSlider")) );
  cycleLengthSlider->assignParameter( moduleToEdit->getParameterByName("CycleLength") );
  cycleLengthSlider->setSliderName(juce::String(T("Cycle")));
  cycleLengthSlider->setDescription(juce::String(T("Length of one cycle (in beats)")));
  cycleLengthSlider->setDescriptionField(infoField);
  cycleLengthSlider->setStringConversionFunction(&rojue::beatsToStringWithUnit4);

  addWidget( depthSlider = new RSlider (T("DepthSlider")) );
  depthSlider->assignParameter( moduleToEdit->getParameterByName("Depth") );
  depthSlider->setDescription(juce::String(T("Depth of the modulation")));
  depthSlider->setDescriptionField(infoField);
  depthSlider->setStringConversionFunction(&rojue::semitonesToStringWithUnit2);

  addWidget( globalFeedbackSlider = new RSlider (T("GlobalFeedbackSlider")) );
  globalFeedbackSlider->assignParameter( moduleToEdit->getParameterByName("GlobalFeedback") );
  globalFeedbackSlider->setSliderName(juce::String(T("Feedback")));
  globalFeedbackSlider->setDescription(juce::String(T("Feedback around the chorus effect")));
  globalFeedbackSlider->setDescriptionField(infoField);
  globalFeedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( crossMixSlider = new RSlider (T("CrossMixSlider")) );
  crossMixSlider->assignParameter( moduleToEdit->getParameterByName("CrossMix") );
  crossMixSlider->setDescription(juce::String(T("Amount by which left wet signal goes to right channel and vice versa")));
  crossMixSlider->setDescriptionField(infoField);
  crossMixSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( feedback2Slider = new RSlider (T("Feedback2Slider")) );
  feedback2Slider->assignParameter( moduleToEdit->getParameterByName("FeedbackPostCrossMix") );
  feedback2Slider->setDescription(juce::String(T("Feedback around the chorus effect after channel cross-mix")));
  feedback2Slider->setDescriptionField(infoField);
  feedback2Slider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( stereoPhaseSlider = new RSlider (T("StereoPhaseSlider")) );
  stereoPhaseSlider->assignParameter( moduleToEdit->getParameterByName("StereoPhase") );
  stereoPhaseSlider->setDescription(juce::String(T("Phase offset between LFOs for left and right channel")));
  stereoPhaseSlider->setDescriptionField(infoField);
  stereoPhaseSlider->setStringConversionFunction(&rojue::degreesToStringWithUnit0);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setDescription(juce::String(T("Ratio between dry (original) and wet (chorused) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);


  addWidget( voice1Button = new RButton(juce::String(T("Voice 1"))) );
  voice1Button->assignParameter( moduleToEdit->getParameterByName(T("OnOffVoice1")) );
  voice1Button->setDescription(juce::String(T("Switch voice 1 on/off")));
  voice1Button->setDescriptionField(infoField);
  voice1Button->setClickingTogglesState(true);

  addWidget( voice1DelaySlider = new RSlider (T("Voice1DelaySlider")) );
  voice1DelaySlider->assignParameter( moduleToEdit->getParameterByName("DelayScaleVoice1") );
  voice1DelaySlider->setSliderName(juce::String(T("Dly")));
  voice1DelaySlider->setDescription(juce::String(T("Scales the delay-time for voice 1")));
  voice1DelaySlider->setDescriptionField(infoField);
  voice1DelaySlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice1DepthSlider = new RSlider (T("Voice1DepthSlider")) );
  voice1DepthSlider->assignParameter( moduleToEdit->getParameterByName("DepthScaleVoice1") );
  voice1DepthSlider->setSliderName(juce::String(T("Dpt")));
  voice1DepthSlider->setDescription(juce::String(T("Scales the modulation depth for voice 1")));
  voice1DepthSlider->setDescriptionField(infoField);
  voice1DepthSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice1AmpSlider = new RSlider (T("Voice1AmpSlider")) );
  voice1AmpSlider->assignParameter( moduleToEdit->getParameterByName("AmpScaleVoice1") );
  voice1AmpSlider->setSliderName(juce::String(T("Amp")));
  voice1AmpSlider->setDescription(juce::String(T("Scales the amplitude for voice 1")));
  voice1AmpSlider->setDescriptionField(infoField);
  voice1AmpSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);


  addWidget( voice2Button = new RButton(juce::String(T("Voice 2"))) );
  voice2Button->assignParameter( moduleToEdit->getParameterByName(T("OnOffVoice2")) );
  voice2Button->setDescription(juce::String(T("Switch voice 2 on/off")));
  voice2Button->setDescriptionField(infoField);
  voice2Button->setClickingTogglesState(true);

  addWidget( voice2DelaySlider = new RSlider (T("Voice2DelaySlider")) );
  voice2DelaySlider->assignParameter( moduleToEdit->getParameterByName("DelayScaleVoice2") );
  voice2DelaySlider->setSliderName(juce::String(T("Dly")));
  voice2DelaySlider->setDescription(juce::String(T("Scales the delay-time for voice 2")));
  voice2DelaySlider->setDescriptionField(infoField);
  voice2DelaySlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice2DepthSlider = new RSlider (T("Voice2DepthSlider")) );
  voice2DepthSlider->assignParameter( moduleToEdit->getParameterByName("DepthScaleVoice2") );
  voice2DepthSlider->setSliderName(juce::String(T("Dpt")));
  voice2DepthSlider->setDescription(juce::String(T("Scales the modulation depth for voice 2")));
  voice2DepthSlider->setDescriptionField(infoField);
  voice2DepthSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice2AmpSlider = new RSlider (T("Voice2AmpSlider")) );
  voice2AmpSlider->assignParameter( moduleToEdit->getParameterByName("AmpScaleVoice2") );
  voice2AmpSlider->setSliderName(juce::String(T("Amp")));
  voice2AmpSlider->setDescription(juce::String(T("Scales the amplitude for voice 2")));
  voice2AmpSlider->setDescriptionField(infoField);
  voice2AmpSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);


  addWidget( voice3Button = new RButton(juce::String(T("Voice 3"))) );
  voice3Button->assignParameter( moduleToEdit->getParameterByName(T("OnOffVoice3")) );
  voice3Button->setDescription(juce::String(T("Switch voice 3 on/off")));
  voice3Button->setDescriptionField(infoField);
  voice3Button->setClickingTogglesState(true);

  addWidget( voice3DelaySlider = new RSlider (T("Voice3DelaySlider")) );
  voice3DelaySlider->assignParameter( moduleToEdit->getParameterByName("DelayScaleVoice3") );
  voice3DelaySlider->setSliderName(juce::String(T("Dly")));
  voice3DelaySlider->setDescription(juce::String(T("Scales the delay-time for voice 3")));
  voice3DelaySlider->setDescriptionField(infoField);
  voice3DelaySlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice3DepthSlider = new RSlider (T("Voice3DepthSlider")) );
  voice3DepthSlider->assignParameter( moduleToEdit->getParameterByName("DepthScaleVoice3") );
  voice3DepthSlider->setSliderName(juce::String(T("Dpt")));
  voice3DepthSlider->setDescription(juce::String(T("Scales the modulation depth for voice 3")));
  voice3DepthSlider->setDescriptionField(infoField);
  voice3DepthSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice3AmpSlider = new RSlider (T("Voice3AmpSlider")) );
  voice3AmpSlider->assignParameter( moduleToEdit->getParameterByName("AmpScaleVoice3") );
  voice3AmpSlider->setSliderName(juce::String(T("Amp")));
  voice3AmpSlider->setDescription(juce::String(T("Scales the amplitude for voice 3")));
  voice3AmpSlider->setDescriptionField(infoField);
  voice3AmpSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);


  addWidget( voice4Button = new RButton(juce::String(T("Voice 4"))) );
  voice4Button->assignParameter( moduleToEdit->getParameterByName(T("OnOffVoice4")) );
  voice4Button->setDescription(juce::String(T("Switch voice 4 on/off")));
  voice4Button->setDescriptionField(infoField);
  voice4Button->setClickingTogglesState(true);

  addWidget( voice4DelaySlider = new RSlider (T("Voice4DelaySlider")) );
  voice4DelaySlider->assignParameter( moduleToEdit->getParameterByName("DelayScaleVoice4") );
  voice4DelaySlider->setSliderName(juce::String(T("Dly")));
  voice4DelaySlider->setDescription(juce::String(T("Scales the delay-time for voice 4")));
  voice4DelaySlider->setDescriptionField(infoField);
  voice4DelaySlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice4DepthSlider = new RSlider (T("Voice4DepthSlider")) );
  voice4DepthSlider->assignParameter( moduleToEdit->getParameterByName("DepthScaleVoice4") );
  voice4DepthSlider->setSliderName(juce::String(T("Dpt")));
  voice4DepthSlider->setDescription(juce::String(T("Scales the modulation depth for voice 4")));
  voice4DepthSlider->setDescriptionField(infoField);
  voice4DepthSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  addWidget( voice4AmpSlider = new RSlider (T("Voice4AmpSlider")) );
  voice4AmpSlider->assignParameter( moduleToEdit->getParameterByName("AmpScaleVoice4") );
  voice4AmpSlider->setSliderName(juce::String(T("Amp")));
  voice4AmpSlider->setDescription(juce::String(T("Scales the amplitude for voice 4")));
  voice4AmpSlider->setDescriptionField(infoField);
  voice4AmpSlider->setStringConversionFunction(&rojue::percentToStringWithUnit0);

  updateWidgetsAccordingToState();
}

void ChorusModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth()/2;
  int h = getHeight()-y;

  globalRect.setBounds(x, y, w, h);
  x += w;
  w /= 2;
  h /= 2;
  voice1Rect.setBounds(x,   y,   w, h);
  voice2Rect.setBounds(x+w, y,   w, h);
  voice3Rect.setBounds(x,   y+h, w, h);
  voice4Rect.setBounds(x+w, y+h, w, h);
  guiLayoutRectangles.clear();
  guiLayoutRectangles.add(globalRect);
  guiLayoutRectangles.add(voice1Rect);
  guiLayoutRectangles.add(voice2Rect);
  guiLayoutRectangles.add(voice3Rect);
  guiLayoutRectangles.add(voice4Rect);

  x = globalRect.getX();
  y = globalRect.getY(); 
  w = globalRect.getWidth();
  h = globalRect.getHeight();
  globalLabel->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 28;
  delaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  cycleLengthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  depthSlider->setBounds(x+4, y+4, w-8, 16);
  //y += 20;  
  //globalFeedbackSlider->setBounds(x+4, y+4, w-8, 16);
  //y += 20;  
  //crossMixSlider->setBounds(x+4, y+4, w-8, 16);
  //y += 20;  
  //feedback2Slider->setBounds(x+4, y+4, w-8, 16);
  //y += 20;  
  //stereoPhaseSlider->setBounds(x+4, y+4, w-8, 16);
  //y += 20;  


  x = voice1Rect.getX();
  y = voice1Rect.getY(); 
  w = voice1Rect.getWidth();
  h = voice1Rect.getHeight();
  voice1Button->setBounds(x+12, y+4, w-24, 16);
  y += 24;
  voice1DelaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice1DepthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice1AmpSlider->setBounds(x+4, y+4, w-8, 16);

  x = voice2Rect.getX();
  y = voice2Rect.getY(); 
  w = voice2Rect.getWidth();
  h = voice2Rect.getHeight();
  voice2Button->setBounds(x+12, y+4, w-24, 16);
  y += 24;
  voice2DelaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice2DepthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice2AmpSlider->setBounds(x+4, y+4, w-8, 16);

  x = voice3Rect.getX();
  y = voice3Rect.getY(); 
  w = voice3Rect.getWidth();
  h = voice3Rect.getHeight();
  voice3Button->setBounds(x+12, y+4, w-24, 16);
  y += 24;
  voice3DelaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice3DepthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice3AmpSlider->setBounds(x+4, y+4, w-8, 16);

  x = voice4Rect.getX();
  y = voice4Rect.getY(); 
  w = voice4Rect.getWidth();
  h = voice4Rect.getHeight();
  voice4Button->setBounds(x+12, y+4, w-24, 16);
  y += 24;
  voice4DelaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice4DepthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  voice4AmpSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// FrequencyShifter:

FrequencyShifterAudioModule::FrequencyShifterAudioModule(CriticalSection *newPlugInLock, 
                                                         rosic::FrequencyShifterStereo *newFrequencyShifterToWrap)
                                                          : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newFrequencyShifterToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedFrequencyShifter = newFrequencyShifterToWrap;
  moduleName  = juce::String(T("FrequencyShifter"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("FrequencyShifterPresets")) );
  createStaticParameters();
}

void FrequencyShifterAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "FrequencyShift", -200.0, 200.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFrequencyShifter, &FrequencyShifterStereo::setFrequencyShift);

  p = new AutomatableParameter(plugInLock, "Feedback", -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFrequencyShifter, &FrequencyShifterStereo::setFeedbackInPercent);

  p = new AutomatableParameter(plugInLock, "StereoOffset", -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedFrequencyShifter, &FrequencyShifterStereo::setStereoOffset);

  //p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  //addObservedParameter(p); 
  //p = new AutomatableParameter(plugInLock, "MidSideRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  //addObservedParameter(p); 

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

FrequencyShifterModuleEditor::FrequencyShifterModuleEditor(CriticalSection *newPlugInLock, 
                                                           FrequencyShifterAudioModule* newFrequencyShifterAudioModule) 
: AudioModuleEditor(newPlugInLock, newFrequencyShifterAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newFrequencyShifterAudioModule != NULL ); // you must pass a valid module here

  addWidget( shiftSlider = new RSlider (T("ShiftSlider")) );
  shiftSlider->assignParameter( moduleToEdit->getParameterByName("FrequencyShift") );
  shiftSlider->setDescription(juce::String(T("Frequency shift in Hz")));
  shiftSlider->setDescriptionField(infoField);
  shiftSlider->setStringConversionFunction(&rojue::hertzToStringWithUnit1);

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( moduleToEdit->getParameterByName("Feedback") );
  feedbackSlider->setDescription(juce::String(T("Feedback around the shifter")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( stereoOffsetSlider = new RSlider (T("StereoOffsetSlider")) );
  stereoOffsetSlider->assignParameter( moduleToEdit->getParameterByName("StereoOffset") );
  stereoOffsetSlider->setDescription(juce::String(T("Stereo offset of the shift between left and right in Hz")));
  stereoOffsetSlider->setDescriptionField(infoField);
  stereoOffsetSlider->setStringConversionFunction(&rojue::hertzToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry (original) and wet (frequency shifted) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  /*
  addWidget( midSideSlider = new RSlider (T("MidSideSlider")) );
  midSideSlider->assignParameter( moduleToEdit->getParameterByName("MidSideRatio") );
  midSideSlider->setSliderName(juce::String(T("Mid/Side")));
  midSideSlider->setDescription(juce::String(T("Ratio between mid and side wet signal")));
  midSideSlider->setDescriptionField(infoField);
  midSideSlider->setStringConversionFunction(&rojue::ratioToString0);
  */

  updateWidgetsAccordingToState();
}

void FrequencyShifterModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  shiftSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  stereoOffsetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
}

//-------------------------------------------------------------------------------------------------
// PhaseStereoizer:

PhaseStereoizerAudioModule::PhaseStereoizerAudioModule(CriticalSection *newPlugInLock, rosic::PhaseStereoizer *newPhaseStereoizerToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newPhaseStereoizerToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedPhaseStereoizer = newPhaseStereoizerToWrap;
  moduleName  = juce::String(T("PhaseStereoizer"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("PhaseStereoizerPresets")) );
  createStaticParameters();
}

void PhaseStereoizerAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "StereoPhaseOffset", 0.0, 180.0, 1.0, 90.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaseStereoizer, &PhaseStereoizer::setPhaseOffset);

  p = new AutomatableParameter(plugInLock, "DryWetRatio", 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaseStereoizer, &PhaseStereoizer::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, "MidSideRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaseStereoizer, &PhaseStereoizer::setMidSideRatio);

  //p = new AutomatableParameter(plugInLock, "Gain", -6.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  //addObservedParameter(p); 
  //p->setValueChangeCallback(wrappedPhaseStereoizer, &PhaseStereoizer::setGain);

  p = new AutomatableParameter(plugInLock, "Lowpass", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaseStereoizer, &PhaseStereoizer::setLowpassCutoff);

  p = new AutomatableParameter(plugInLock, "Highpass", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedPhaseStereoizer, &PhaseStereoizer::setHighpassCutoff);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

PhaseStereoizerModuleEditor::PhaseStereoizerModuleEditor(CriticalSection *newPlugInLock, 
                                                         PhaseStereoizerAudioModule* newPhaseStereoizerAudioModule) 
: AudioModuleEditor(newPlugInLock, newPhaseStereoizerAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newPhaseStereoizerAudioModule != NULL ); // you must pass a valid module here

  addWidget( phaseOffsetSlider = new RSlider (T("PhaseOffsetSlider")) );
  phaseOffsetSlider->assignParameter( moduleToEdit->getParameterByName("StereoPhaseOffset") );
  phaseOffsetSlider->setDescription(juce::String(T("Phase offset between left and right (wet) signal in degrees")));
  phaseOffsetSlider->setDescriptionField(infoField);
  phaseOffsetSlider->setStringConversionFunction(&rojue::degreesToStringWithUnit0);

  addWidget( dryWetRatioSlider = new RSlider (T("DryWetRatioSlider")) );
  dryWetRatioSlider->assignParameter( moduleToEdit->getParameterByName("DryWetRatio") );
  dryWetRatioSlider->setDescription(juce::String(T("Ratio between dry (original) and wet (phase-shifted) signal")));
  dryWetRatioSlider->setDescriptionField(infoField);
  dryWetRatioSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( sideLowpassSlider = new RSlider (T("SideLowpassSlider")) );
  sideLowpassSlider->assignParameter( moduleToEdit->getParameterByName("Lowpass") );
  sideLowpassSlider->setSliderName(juce::String(T("Lowpass")));
  sideLowpassSlider->setDescription(juce::String(T("Cutoff frequency of the lowpass filter for the wet side signal")));
  sideLowpassSlider->setDescriptionField(infoField);
  sideLowpassSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( sideHighpassSlider = new RSlider (T("SideHighpassSlider")) );
  sideHighpassSlider->assignParameter( moduleToEdit->getParameterByName("Highpass") );
  sideHighpassSlider->setSliderName(juce::String(T("Highpass")));
  sideHighpassSlider->setDescription(juce::String(T("Cutoff frequency of the highpass filter for the wet side signal")));
  sideHighpassSlider->setDescriptionField(infoField);
  sideHighpassSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( midSideRatioSlider = new RSlider (T("MidSideRatioSlider")) );
  midSideRatioSlider->assignParameter( moduleToEdit->getParameterByName("MidSideRatio") );
  midSideRatioSlider->setDescription(juce::String(T("Ratio between mid and side signal (stereo-width)")));
  midSideRatioSlider->setDescriptionField(infoField);
  midSideRatioSlider->setStringConversionFunction(&rojue::ratioToString0);

  /*
  addWidget( gainSlider = new RSlider (T("GainSlider")) );
  gainSlider->assignParameter( moduleToEdit->getParameterByName("Gain") );
  gainSlider->setDescription(juce::String(T("Global gain for compensation of gain changes")));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);
  */

  updateWidgetsAccordingToState();
}

void PhaseStereoizerModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth()/2;
  int h = getHeight();
  y = getPresetSectionBottom();

  phaseOffsetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
  midSideRatioSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
  dryWetRatioSlider->setBounds(x+4, y+4, w-8, 16);


  x += w;
  y  = getPresetSectionBottom();
  sideHighpassSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  sideLowpassSlider->setBounds(x+4, y+4, w-8, 16);
  //gainSlider->setBounds(x+4, y+4, w-8, 16);
}

//-------------------------------------------------------------------------------------------------
// RingModulator:

RingModulatorAudioModule::RingModulatorAudioModule(CriticalSection *newPlugInLock, rosic::RingModulatorStereo *newRingModulatorToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newRingModulatorToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedRingModulator = newRingModulatorToWrap;
  moduleName  = juce::String(T("RingModulator"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("RingModulatorPresets")) );
  createStaticParameters();
}

void RingModulatorAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, T("Frequency"), 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedRingModulator, &RingModulatorStereo::setModulatorFrequency);

  p = new AutomatableParameter(plugInLock, T("Feedback"), -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedRingModulator, &RingModulatorStereo::setFeedbackInPercent);

  p = new AutomatableParameter(plugInLock, T("StereoOffset"), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedRingModulator, &RingModulatorStereo::setStereoOffset);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedRingModulator, &RingModulatorStereo::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, T("AntiAlias"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedRingModulator, &RingModulatorStereo::setAntiAliasing);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

RingModulatorModuleEditor::RingModulatorModuleEditor(CriticalSection *newPlugInLock, RingModulatorAudioModule* newRingModulatorAudioModule) 
: AudioModuleEditor(newPlugInLock, newRingModulatorAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newRingModulatorAudioModule != NULL ); // you must pass a valid module here

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName(T("Frequency")) );
  frequencySlider->setDescription(juce::String(T("Frequency in Hz")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( moduleToEdit->getParameterByName(T("Feedback")) );
  feedbackSlider->setDescription(juce::String(T("Feedback around the ringmodulator")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( stereoOffsetSlider = new RSlider (T("StereoOffsetSlider")) );
  stereoOffsetSlider->assignParameter( moduleToEdit->getParameterByName(T("StereoOffset")) );
  stereoOffsetSlider->setDescription(juce::String(T("Stereo offset of the frequency between left and right in Hz")));
  stereoOffsetSlider->setDescriptionField(infoField);
  stereoOffsetSlider->setStringConversionFunction(&rojue::hertzToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName(T("DryWetRatio")) );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry (original) and wet (ringmodulated) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( antiAliasButton = new RButton(juce::String(T("AntiAlias"))) );
  antiAliasButton->assignParameter( moduleToEdit->getParameterByName(T("AntiAlias")) );
  antiAliasButton->setDescription(juce::String(T("Toggle anti-aliasing by oversampling on/off")));
  antiAliasButton->setDescriptionField(infoField);
  antiAliasButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void RingModulatorModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  stereoOffsetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  antiAliasButton->setBounds(x+4, y+4, 60, 16);
  y += 20;  
}

//-------------------------------------------------------------------------------------------------
// SingleSidebandModulator:

SingleSidebandModulatorAudioModule::SingleSidebandModulatorAudioModule(CriticalSection *newPlugInLock, 
  rosic::SingleSidebandModulatorStereo *newSingleSidebandModulatorToWrap)
   : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newSingleSidebandModulatorToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedSingleSidebandModulator = newSingleSidebandModulatorToWrap;
  moduleName  = juce::String(T("SingleSidebandModulator"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("SingleSidebandModulatorPresets")) );
  createStaticParameters();
}

void SingleSidebandModulatorAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, T("Frequency"), 20.0, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setModulatorFrequency);

  p = new AutomatableParameter(plugInLock, T("UpperSidebandLevel"), -60.0, 0.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setUpperSidebandLevel);

  p = new AutomatableParameter(plugInLock, T("LowerSidebandLevel"), -60.0, 0.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setLowerSidebandLevel);

  p = new AutomatableParameter(plugInLock, T("Feedback"), -99.0, 99.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setFeedbackInPercent);

  p = new AutomatableParameter(plugInLock, T("StereoOffset"), -100.0, 100.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setStereoOffset);

  p = new AutomatableParameter(plugInLock, T("DryWetRatio"), 0.0, 1.0, 0.01, 1.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setDryWetRatio);

  p = new AutomatableParameter(plugInLock, T("AntiAlias"), 0.0, 1.0, 1.0, 0.0, Parameter::BOOLEAN);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedSingleSidebandModulator, &SingleSidebandModulatorStereo::setAntiAliasing);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

SingleSidebandModulatorModuleEditor::SingleSidebandModulatorModuleEditor(CriticalSection *newPlugInLock, 
                                                                         SingleSidebandModulatorAudioModule* newSingleSidebandModulatorAudioModule) 
: AudioModuleEditor(newPlugInLock, newSingleSidebandModulatorAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newSingleSidebandModulatorAudioModule != NULL ); // you must pass a valid module here

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName(T("Frequency")) );
  frequencySlider->setDescription(juce::String(T("Frequency in Hz")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( upperSidebandLevelSlider = new RSlider (T("UpperSidebandLevelSlider")) );
  upperSidebandLevelSlider->assignParameter( moduleToEdit->getParameterByName(T("UpperSidebandLevel")) );
  upperSidebandLevelSlider->setDescription(juce::String(T("Upper sideband level in dB")));
  upperSidebandLevelSlider->setDescriptionField(infoField);
  upperSidebandLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( lowerSidebandLevelSlider = new RSlider (T("LowerSidebandLevelSlider")) );
  lowerSidebandLevelSlider->assignParameter( moduleToEdit->getParameterByName(T("LowerSidebandLevel")) );
  lowerSidebandLevelSlider->setDescription(juce::String(T("Lower sideband level in dB")));
  lowerSidebandLevelSlider->setDescriptionField(infoField);
  lowerSidebandLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( feedbackSlider = new RSlider (T("FeedbackSlider")) );
  feedbackSlider->assignParameter( moduleToEdit->getParameterByName(T("Feedback")) );
  feedbackSlider->setDescription(juce::String(T("Feedback around the SSB-modulator")));
  feedbackSlider->setDescriptionField(infoField);
  feedbackSlider->setStringConversionFunction(&rojue::percentToStringWithUnit1);

  addWidget( stereoOffsetSlider = new RSlider (T("StereoOffsetSlider")) );
  stereoOffsetSlider->assignParameter( moduleToEdit->getParameterByName(T("StereoOffset")) );
  stereoOffsetSlider->setDescription(juce::String(T("Stereo offset of the frequency between left and right in Hz")));
  stereoOffsetSlider->setDescriptionField(infoField);
  stereoOffsetSlider->setStringConversionFunction(&rojue::hertzToStringWithUnit1);

  addWidget( dryWetSlider = new RSlider (T("DryWetSlider")) );
  dryWetSlider->assignParameter( moduleToEdit->getParameterByName(T("DryWetRatio")) );
  dryWetSlider->setSliderName(juce::String(T("Dry/Wet")));
  dryWetSlider->setDescription(juce::String(T("Ratio between dry (original) and wet (modulated) signal")));
  dryWetSlider->setDescriptionField(infoField);
  dryWetSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( antiAliasButton = new RButton(juce::String(T("AntiAlias"))) );
  antiAliasButton->assignParameter( moduleToEdit->getParameterByName(T("AntiAlias")) );
  antiAliasButton->setDescription(juce::String(T("Toggle anti-aliasing by oversampling on/off")));
  antiAliasButton->setDescriptionField(infoField);
  antiAliasButton->setClickingTogglesState(true);

  updateWidgetsAccordingToState();
}

void SingleSidebandModulatorModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  upperSidebandLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  lowerSidebandLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  feedbackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  stereoOffsetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  dryWetSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  antiAliasButton->setBounds(x+4, y+4, 60, 16);
  y += 20;  
}

//-------------------------------------------------------------------------------------------------
// StereoPan:

StereoPanAudioModule::StereoPanAudioModule(CriticalSection *newPlugInLock, rosic::StereoPan *newStereoPanToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newStereoPanToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedStereoPan = newStereoPanToWrap;
  moduleName  = juce::String(T("StereoPan"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("StereoPanPresets")) );
  createStaticParameters();
}

void StereoPanAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, juce::String(T("PanLaw")), 0.0, 8.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String(T("Linear")));
  p->addStringValue(juce::String(T("Trigonometric")));
  p->addStringValue(juce::String(T("Square Root")));
  p->addStringValue(juce::String(T("Linear Clipped")));
  p->addStringValue(juce::String(T("Linear Renormalized")));
  p->addStringValue(juce::String(T("Linear Clipped, Renormalized")));
  p->addStringValue(juce::String(T("X-Mix Linear")));
  p->addStringValue(juce::String(T("X-Mix Trigonometric")));
  p->addStringValue(juce::String(T("X-Mix Square Root")));
  //p->setValue(1.0, false);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedStereoPan, &StereoPan::setPanLaw);

  p = new AutomatableParameter(plugInLock, "Pan", -1.0, 1.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedStereoPan, &StereoPan::setPanoramaPosition);

  p = new AutomatableParameter(plugInLock, "Gain", -12.0, 12.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedStereoPan, &StereoPan::setGain);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

StereoPanModuleEditor::StereoPanModuleEditor(CriticalSection *newPlugInLock, StereoPanAudioModule* newStereoPanAudioModule) 
: AudioModuleEditor(newPlugInLock, newStereoPanAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newStereoPanAudioModule != NULL ); // you must pass a valid module here
  stereoPanModuleToEdit = newStereoPanAudioModule;

  addWidget( panLawLabel = new RTextField( juce::String(T("Pan Law:"))) );
  panLawLabel->setDescription(T("Select the pan-law"));
  panLawLabel->setDescriptionField(infoField);

  addWidget( panLawComboBox = new RComboBox(juce::String(T("PanLawComboBox"))) );
  panLawComboBox->assignParameter( moduleToEdit->getParameterByName(T("PanLaw")) );
  panLawComboBox->setDescription(panLawLabel->getDescription());
  panLawComboBox->setDescriptionField(infoField);
  //panLawComboBox->addListener(this); // to update the plot

  addWidget( panSlider = new RSlider (T("PanSlider")) );
  panSlider->assignParameter( moduleToEdit->getParameterByName("Pan") );
  panSlider->setDescription(juce::String(T("Panorama position")));
  panSlider->setDescriptionField(infoField);
  panSlider->setStringConversionFunction(&rojue::valueToString2);
  //panSlider->addListener(this); // to update the plot

  addWidget( gainSlider = new RSlider (T("GainSlider")) );
  gainSlider->assignParameter( moduleToEdit->getParameterByName("Gain") );
  gainSlider->setDescription(juce::String(T("Gain")));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit2);

  numValues = 0;
  xValues   = NULL;
  yValues   = NULL;
  plot = new StereoPanPlotEditor(juce::String(T("Plot")));
  plot->setStereoPanToEdit(stereoPanModuleToEdit->wrappedStereoPan);
  plot->assignParameterPan(    moduleToEdit->getParameterByName("Pan"));
  plot->assignParameterPanLaw(moduleToEdit->getParameterByName("PanLaw"));
  addPlot(plot);

  updateWidgetsAccordingToState();
}

StereoPanModuleEditor::~StereoPanModuleEditor()
{
  ScopedLock scopedLock(*plugInLock);
  if( xValues != NULL ) { delete[] xValues; xValues = NULL; }
  if( yValues != NULL ) { delete[] yValues; yValues = NULL; }
}

void StereoPanModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight()-y;

  plot->setBounds(x+4, y+4, w-8, h-48);

  y = plot->getBottom();
  w = w/2;
  panLawLabel->setBounds(x+4, y+4, 56, 16); 
  panLawComboBox->setBounds(panLawLabel->getRight()+4, y+4, 180, 16);
  y += 20;  
  panSlider->setBounds( x+4,   y+4, w-8, 16);
  gainSlider->setBounds(x+w+4, y+4, w-8, 16);
  // maybe add an routing-diagram....
}

//-------------------------------------------------------------------------------------------------
// StereoWidth:

StereoWidthAudioModule::StereoWidthAudioModule(CriticalSection *newPlugInLock, rosic::StereoWidth *newStereoWidthToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newStereoWidthToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedStereoWidth = newStereoWidthToWrap;
  moduleName  = juce::String(T("StereoWidth"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("StereoWidthPresets")) );
  createStaticParameters();
}

void StereoWidthAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "MidSideRatio", 0.0, 1.0, 0.01, 0.5, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedStereoWidth, &StereoWidth::setMidSideRatio);

  p = new AutomatableParameter(plugInLock, "Gain", -6.0, 24.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedStereoWidth, &StereoWidth::setGlobalGain);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

StereoWidthModuleEditor::StereoWidthModuleEditor(CriticalSection *newPlugInLock, StereoWidthAudioModule* newStereoWidthAudioModule) 
: AudioModuleEditor(newPlugInLock, newStereoWidthAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newStereoWidthAudioModule != NULL ); // you must pass a valid module here

  addWidget( midSideRatioSlider = new RSlider (T("MidSideRatioSlider")) );
  midSideRatioSlider->assignParameter( moduleToEdit->getParameterByName("MidSideRatio") );
  midSideRatioSlider->setDescription(juce::String(T("Ratio between mid and side signal (stereo-width)")));
  midSideRatioSlider->setDescriptionField(infoField);
  midSideRatioSlider->setStringConversionFunction(&rojue::ratioToString0);

  addWidget( gainSlider = new RSlider (T("GainSlider")) );
  gainSlider->assignParameter( moduleToEdit->getParameterByName("Gain") );
  gainSlider->setDescription(juce::String(T("Global gain for compensation of gain changes")));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  updateWidgetsAccordingToState();
}

void StereoWidthModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  midSideRatioSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
  gainSlider->setBounds(x+4, y+4, w-8, 16);
}

//=================================================================================================
// Signal Generators:

//-------------------------------------------------------------------------------------------------
// SineOscillator:

SineOscillatorAudioModule::SineOscillatorAudioModule(CriticalSection *newPlugInLock, rosic::SineOscillator *newSineOscillatorToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newSineOscillatorToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedSineOscillator = newSineOscillatorToWrap;
  moduleName  = juce::String(T("SineOscillator"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("SineOscillatorPresets")) );
  createStaticParameters();
}

void SineOscillatorAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "Frequency", 0.2, 20000.0, 0.0, 1000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedSineOscillator, &SineOscillator::setFrequency);

  // \todo: add gain parameter

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

SineOscillatorModuleEditor::SineOscillatorModuleEditor(CriticalSection *newPlugInLock, 
                                                       SineOscillatorAudioModule* newSineOscillatorAudioModule) 
: AudioModuleEditor(newPlugInLock, newSineOscillatorAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newSineOscillatorAudioModule != NULL ); // you must pass a valid module here

  addWidget( frequencySlider = new RSlider (T("FrequencySlider")) );
  frequencySlider->assignParameter( moduleToEdit->getParameterByName("Frequency") );
  frequencySlider->setDescription(juce::String(T("Frequency of the sinusoid")));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  updateWidgetsAccordingToState();
}

void SineOscillatorModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();
  frequencySlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;  
}

//-------------------------------------------------------------------------------------------------
// Noisifier:

NoisifierAudioModule::NoisifierAudioModule(CriticalSection *newPlugInLock, rosic::Noisifier *newNoisifierToWrap)
 : AudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*plugInLock);

  jassert( newNoisifierToWrap != NULL ); // you must pass a valid rosic-object 
  wrappedNoisifier = newNoisifierToWrap;
  moduleName       = juce::String(T("Noisifier"));
  setActiveDirectory(getApplicationDirectory() + File::separatorString 
    + juce::String(T("NoisifierPresets")) );
  createStaticParameters();
}

void NoisifierAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*plugInLock);

  AutomatableParameter* p;

  p = new AutomatableParameter(plugInLock, "PassLevel", -96.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p); 
  p->setValueChangeCallback(wrappedNoisifier, &Noisifier::setPassThroughLevel);

  p = new AutomatableParameter(plugInLock, "NoiseLevel", -96.0, 6.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedNoisifier, &Noisifier::setNoiseLevel);

  p = new AutomatableParameter(plugInLock, "SpectralSlope", -12.0, 12.0, 0.1, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedNoisifier, &Noisifier::setNoiseSpectralSlope);

  p = new AutomatableParameter(plugInLock, "LowestFrequency", 20.0, 20000.0, 0.0, 20.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedNoisifier, &Noisifier::setLowestFrequency);

  p = new AutomatableParameter(plugInLock, "HighestFrequency", 20.0, 20000.0, 0.0, 20000.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedNoisifier, &Noisifier::setHighestFrequency);

  for(int i=0; i < (int) observedParameters.size(); i++)
    observedParameters[i]->resetToDefaultValue(true, true);
}

NoisifierModuleEditor::NoisifierModuleEditor(CriticalSection *newPlugInLock, NoisifierAudioModule* newNoisifierAudioModule) 
: AudioModuleEditor(newPlugInLock, newNoisifierAudioModule)
{
  ScopedLock scopedLock(*plugInLock);

  jassert(newNoisifierAudioModule != NULL ); // you must pass a valid module here

  addWidget( passLevelSlider = new RSlider (T("PassLevelSlider")) );
  passLevelSlider->assignParameter( moduleToEdit->getParameterByName("PassLevel") );
  passLevelSlider->setDescription(juce::String(T("Level with which the input signal is passed through")));
  passLevelSlider->setDescriptionField(infoField);
  passLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( noiseLevelSlider = new RSlider (T("NoiseLevelSlider")) );
  noiseLevelSlider->assignParameter( moduleToEdit->getParameterByName("NoiseLevel") );
  noiseLevelSlider->setDescription(juce::String(T("Level with which the generated noise is mixed in")));
  noiseLevelSlider->setDescriptionField(infoField);
  noiseLevelSlider->setStringConversionFunction(&rojue::decibelsToStringWithUnit1);

  addWidget( spectralSlopeSlider = new RSlider (T("SpectralSlopeSlider")) );
  spectralSlopeSlider->assignParameter( moduleToEdit->getParameterByName("SpectralSlope") );
  spectralSlopeSlider->setDescription(juce::String(T("Spectral slope of the generated noise")));
  spectralSlopeSlider->setDescriptionField(infoField);
  spectralSlopeSlider->setStringConversionFunction(&rojue::decibelsPerOctaveToString2);

  addWidget( lowestFreqSlider = new RSlider (T("LowestFreqSlider")) );
  lowestFreqSlider->assignParameter( moduleToEdit->getParameterByName("LowestFrequency") );
  lowestFreqSlider->setDescription(juce::String(T("Lowest frequency present in the noise")));
  lowestFreqSlider->setDescriptionField(infoField);
  lowestFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);

  addWidget( highestFreqSlider = new RSlider (T("HighestFreqSlider")) );
  highestFreqSlider->assignParameter( moduleToEdit->getParameterByName("HighestFrequency") );
  highestFreqSlider->setDescription(juce::String(T("Highest frequency present in the noise")));
  highestFreqSlider->setDescriptionField(infoField);
  highestFreqSlider->setStringConversionFunction(&rojue::hertzToStringWithUnitTotal5);


  updateWidgetsAccordingToState();
}

void NoisifierModuleEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);

  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  y = getPresetSectionBottom();

  passLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
  noiseLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
  spectralSlopeSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
  lowestFreqSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
  highestFreqSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20; 
}
