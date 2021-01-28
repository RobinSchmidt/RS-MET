
// construction/destruction:

AciDevilAudioModule::AciDevilAudioModule(CriticalSection *newPlugInLock, 
  rosic::AciDevil *aciDevilToWrap) : AudioModuleWithMidiIn(newPlugInLock) 
{
  jassert(aciDevilToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedAciDevil = aciDevilToWrap;
  init();
}

AciDevilAudioModule::AciDevilAudioModule(CriticalSection *newPlugInLock) 
  : AudioModuleWithMidiIn(newPlugInLock) 
{
  wrappedAciDevil = new rosic::AciDevil;
  wrappedAciDevilIsOwned = true;
  init();
}

void AciDevilAudioModule::init()
{
  setModuleTypeName("AcidDevil");
  createParameters();
  sequencerModule = new AcidSequencerAudioModule(lock, &wrappedAciDevil->sequencer);
  sequencerModule->setModuleName(juce::String("Sequencer"));
  addChildAudioModule(sequencerModule);
}

AciDevilAudioModule::~AciDevilAudioModule()
{
  if(wrappedAciDevilIsOwned)
    delete wrappedAciDevil;
}

AudioModuleEditor* AciDevilAudioModule::createEditor(int type)
{
  return new jura::AciDevilModuleEditor(lock, this); // get rid of passing the lock
}

// internal functions:

void AciDevilAudioModule::createParameters()
{
  //typedef MetaControlledParameter Param;
  typedef ModulatableParameter Param;
  Param* p;

  typedef rosic::AciDevil AD;
  AD* ad = wrappedAciDevil;

  p = new Param("MasterLevel", -60.0, 0.0, -12.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setMasterLevel);
  addObservedParameter(p);

  p = new Param("Accent", 0.0, 100.0, 50.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setAccent);
  addObservedParameter(p);

  p = new Param("SlideTime", 1.0, 500.0, 60.0, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setSlideTime);
  addObservedParameter(p);

  p = new Param("Waveform", 0.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<AD>(ad, &AD::setWaveform);
  addObservedParameter(p);

  p = new Param("PulseWidth", 1.0, 100.0, 45.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setPulseWidth);
  addObservedParameter(p);
  // 45 is the default value because that's roughly what i have measured in real 303 samples
  // ...but the DSP object does not to respond to it...why? is it not yet implemented?

  p = new Param("SubOscLevel", -60.0, 0.0, -60.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setSubOscLevel);
  addObservedParameter(p);

  p = new Param("SubOscWaveform", 0.0, 1.0, 1.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<AD>(ad, &AD::setSubOscWaveform);
  addObservedParameter(p);

  p = new Param("Cutoff", 200.0, 10000.0, 300.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<AD>(ad, &AD::setCutoff);
  addObservedParameter(p);

  p = new Param("Resonance", 0.0, 100.0, 50.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setResonance);
  addObservedParameter(p);

  p = new Param("FilterMode", 0.0, 14.0, 1.0, Parameter::STRING, 1.0);
  p->addStringValue("Flat");
  p->addStringValue("Lowpass 6");
  p->addStringValue("Lowpass 12");
  p->addStringValue("Lowpass 18");
  p->addStringValue("Lowpass 24");
  p->addStringValue("Highpass 6");
  p->addStringValue("Highpass 12");
  p->addStringValue("Highpass 18");
  p->addStringValue("Highpass 24");
  p->addStringValue("Bandpass 12+12");
  p->addStringValue("Bandpass 6+18");
  p->addStringValue("Bandpass 18+6");
  p->addStringValue("Bandpass 6+12");
  p->addStringValue("Bandpass 12+6");
  p->addStringValue("Bandpass 6+6");
  p->setValue(3.0, false, false);
  p->setValueChangeCallback<AD>(ad, &AD::setFilterMode);
  addObservedParameter(p);

  p = new Param("EnvMod", 0.0, 80.0, 12.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setEnvMod);
  addObservedParameter(p);

  p = new Param("NormalDecay", 30.0, 3000.0, 200.0, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setNormalDecay);
  addObservedParameter(p);

  p = new Param("AccentDecay", 30.0, 300.0, 60.0, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setAccentDecay);
  addObservedParameter(p);

  p = new Param("NormalAttack", 3.0, 50.0, 3.0, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setNormalAttack);
  addObservedParameter(p);

  p = new Param("AccentAttack", 3.0, 50.0, 10.0, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setAccentAttack);
  addObservedParameter(p);

  //p = new Param("UpwardFraction", 0.0, 100.0, 66.6, Parameter::LINEAR, 0.1);
  //p->setValueChangeCallback<AD>(ad, &AD::setUpwardFraction);
  //addObservedParameter(p);

  p = new Param("AmpDecay", 3.0, 3000.0, 1230.0, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setAmpDecay);
  addObservedParameter(p);

  p = new Param("AmpSustain", -60.0, 0.0, -60.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setAmpSustain);
  addObservedParameter(p);

  p = new Param("AmpRelease", 0.3, 50.0, 0.5, Parameter::EXPONENTIAL, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setAmpRelease);
  addObservedParameter(p);

  p = new Param("DistortionDrive", -24.0, 60.0, 0.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<AD>(ad, &AD::setClipperDrive);
  addObservedParameter(p);
}

//=================================================================================================

AciDevilModuleEditor::AciDevilModuleEditor(CriticalSection *newPlugInLock, 
  AciDevilAudioModule* newAciDevilAudioModule) 
  : AudioModuleEditor(newAciDevilAudioModule)
{
  setHeadlineStyle(MAIN_HEADLINE);

  // assign the pointer to the rosic::AciDevil object to be used as aduio engine:
  jassert(newAciDevilAudioModule != NULL ); // you must pass a valid module here
  aciDevilModuleToEdit = newAciDevilAudioModule;

  createWidgets();
  updateWidgetsAccordingToState();
  //setSize(772, 394);
  setSize(634, 394);
}

void AciDevilModuleEditor::createWidgets()
{
  typedef rsAutomatableSlider Sld;
  typedef rsAutomatableComboBox Box;
  //typedef rsAutomatableButton Btn;
  typedef RTextField Lbl;
  Sld* s;
  Box* c;
  //Btn* b;
  Lbl* l;

  addWidget( globalLabel = l = new Lbl("Global"));
  //l->setJustificationType(Justification::centred);
  l->setDescription("Global parameters");
  l->setDescriptionField(infoField);

  addWidget( masterLevelSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("MasterLevel") );
  s->setSliderName("Level");
  s->setDescription("Master level in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( accentSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Accent") );
  s->setDescription("Accent in percent");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( slideTimeSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("SlideTime") );
  s->setSliderName("Slide");
  s->setDescription("Slide time in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( oscLabel = new Lbl("Oscillator") );
  oscLabel->setJustification(Justification::centred);
  oscLabel->setDescription("Oscillator parameters");
  oscLabel->setDescriptionField(infoField);

  addWidget( waveformSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Waveform") );
  s->setSliderName("Saw/Pulse");
  s->setDescription("Mix between saw- and pulse-wave for main oscillator");
  s->setStringConversionFunction(ratioToString0);
  s->setDescriptionField(infoField);

  addWidget( pulseWidthSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("PulseWidth") );
  s->setSliderName("PulseWidth");
  s->setDescription("Width of high section of the rectangular pulse waveform");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( subOscLabel = new Lbl("SubOsc:") );
  subOscLabel->setJustification(Justification::centredLeft);
  subOscLabel->setDescription("Sub-oscillator settings");
  subOscLabel->setDescriptionField(infoField);

  addWidget( subOscLevelSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("SubOscLevel") );
  s->setSliderName("Level");
  s->setDescription("Sub-oscillator level in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( subOscWaveformSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("SubOscWaveform") );
  s->setSliderName("Saw/Pulse");
  s->setDescription("Mix between saw- and pulse-wave for suboscillator");
  s->setStringConversionFunction(ratioToString0);
  s->setDescriptionField(infoField);


  addWidget( filterLabel = new Lbl("Filter") );
  filterLabel->setJustification(Justification::centred);
  filterLabel->setDescription("Filter");
  filterLabel->setDescriptionField(infoField);

  addWidget( cutoffSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Cutoff") );
  s->setDescription("Filter cutoff frequency in Hz");
  s->setStringConversionFunction(hertzToStringWithUnitTotal5);
  s->setDescriptionField(infoField);

  addWidget( resonanceSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Resonance") );
  s->setDescription("Resonance in percent");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( filterModeLabel = new Lbl("Mode:") );
  filterModeLabel->setJustification(Justification::centredLeft);
  filterModeLabel->setDescription("Choose the filter mode");
  filterModeLabel->setDescription("Mode:");
  filterModeLabel->setDescriptionField(infoField);

  addWidget( filterModeBox = c = new Box );
  c->assignParameter( moduleToEdit->getParameterByName("FilterMode") );
  c->setDescription(filterModeLabel->getDescription());
  c->setDescriptionField(infoField);

  addWidget( envModSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("EnvMod") );
  s->setDescription("Amount of modulation of cutoff frequency in semitones");
  s->setStringConversionFunction(semitonesToStringWithUnit1);
  s->setDescriptionField(infoField);


  addWidget( filterEnvLabel = new Lbl("Filter Envelope") );
  filterEnvLabel->setJustification(Justification::centred);
  filterEnvLabel->setDescription("Filter envelope parameters");
  filterEnvLabel->setDescriptionField(infoField);

  addWidget( normalLabel = new Lbl("Normal:") );
  normalLabel->setJustification(Justification::centredLeft);
  normalLabel->setDescription("Time values for normal (un-accented) notes");
  normalLabel->setDescriptionField(infoField);

  addWidget( normalDecaySlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("NormalDecay") );
  s->setSliderName("Decay");
  s->setDescription("Decay time for normal (un-accented) notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( normalAttackSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("NormalAttack") );
  s->setSliderName("Attack");
  s->setDescription("Attack time for normal (un-accented) notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);


  addWidget( accentLabel = new Lbl("Accent:") );
  accentLabel->setJustification(Justification::centredLeft);
  accentLabel->setDescription("Time values for accented notes");
  accentLabel->setDescriptionField(infoField);

  addWidget( accentDecaySlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AccentDecay") );
  s->setSliderName("Decay");
  s->setDescription("Decay time for accented notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( accentAttackSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AccentAttack") );
  s->setSliderName("Attack");
  s->setDescription("Attack time for accented notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( ampLabel = new Lbl("Amp Envelope") );
  ampLabel->setJustification(Justification::centred);
  ampLabel->setDescription("Amplide envelope and distortion parameters");
  ampLabel->setDescriptionField(infoField);

  addWidget( ampDecaySlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AmpDecay") );
  s->setSliderName("Decay");
  s->setDescription("Decay time for amplitude envelope in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( ampSustainSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AmpSustain") );
  s->setSliderName("Sustain");
  s->setDescription("Sustain level for amplitude envelope in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( ampReleaseSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AmpRelease") );
  s->setSliderName("Release");
  s->setDescription("Release time for amplitude envelope in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);


  addWidget( distLabel = new Lbl("Distortion") );
  distLabel->setJustification(Justification::centred);
  distLabel->setDescription("Distortion Settings");
  distLabel->setDescriptionField(infoField);

  addWidget( distortionDriveSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("DistortionDrive") );
  s->setSliderName("Drive");
  s->setDescription("Drive for distortion unit in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  sequencerEditor = new AcidSequencerModuleEditor(lock, aciDevilModuleToEdit->sequencerModule);
  addChildEditor( sequencerEditor );
  sequencerEditor->setDescriptionField(infoField, true);
}

//-------------------------------------------------------------------------------------------------
// setup:





//-------------------------------------------------------------------------------------------------
// callbacks:

void AciDevilModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  sequencerEditor->updateWidgetsAccordingToState();
}

void AciDevilModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);
  // maybe write into the empty area that a more sophisticated distortion is to come
  // or maybe use it for a filter-plot
}

void AciDevilModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = getHeadlineBottom()+4;
  w = 220;
  h = 120;

  guiLayoutRectangles.clear();
  globalRectangle.setBounds(x, y, w, h);
  x = globalRectangle.getRight()-2;
  w = 140;
  oscRectangle.setBounds(x, y, w, h);
  x = oscRectangle.getRight()-2;
  w = 140;
  filterRectangle.setBounds(x, y, w, h);
  x = filterRectangle.getRight()-2;
  w = 140;
  filterEnvRectangle.setBounds(x, y, w, h);
  w = 140;
  y = filterEnvRectangle.getBottom()-2;
  h = 80;  // test
  ampRectangle.setBounds(x, y, w, h); 

  y = ampRectangle.getBottom()-2;
  h = getHeight() - y;
  distRectangle.setBounds(x, y, w, h);


  guiLayoutRectangles.add(globalRectangle);
  guiLayoutRectangles.add(oscRectangle);
  guiLayoutRectangles.add(filterRectangle);
  guiLayoutRectangles.add(filterEnvRectangle);
  guiLayoutRectangles.add(ampRectangle);
  guiLayoutRectangles.add(distRectangle);


  x = globalRectangle.getX();
  y = globalRectangle.getY();
  w = globalRectangle.getWidth();

  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setBounds(x+4, y+4, w-8, 32);
  y = stateWidgetSet->getBottom()+4;

  y = stateWidgetSet->getBottom() + 4 + 32; // leave space for tuning widgets

  masterLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y = masterLevelSlider->getBottom();
  accentSlider->setBounds(    x+4,     y+4, w/2-8, 16);
  slideTimeSlider->setBounds( x+w/2+4, y+4, w/2-8, 16);


  x = oscRectangle.getX();
  y = oscRectangle.getY();
  w = oscRectangle.getWidth();
  oscLabel->setBounds(x, y+2, w, 16);
  y = oscLabel->getBottom();
  waveformSlider->setBounds(x+4, y+4, w-8, 16);
  //y += 16;
  //pulseWidthSlider->setBounds(x+4, y+4, w-8, 16);
  y += 24;
  subOscLabel->setBounds(x+4, y+4, w-8, 16);
  y += 16;

  subOscLevelSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  subOscWaveformSlider->setBounds(x+4, y+4, w-8, 16);

  x = filterRectangle.getX();
  y = filterRectangle.getY();
  w = filterRectangle.getWidth();
  filterLabel->setBounds(x, y+2, w, 16);
  y = filterLabel->getBottom();
  cutoffSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  resonanceSlider->setBounds(x+4, y+4, w-8, 16);
  y += 20;
  filterModeLabel->setBounds(x+4,    y+4, 40,     16);
  filterModeBox->setBounds(  x+40+4, y+4, w-40-8, 16);
  y += 20;  // maybe use 4 or 8 pixels more - it should have a greate distance
  envModSlider->setBounds(x+4, y+4, w-8, 16);

  x = filterEnvRectangle.getX();
  y = filterEnvRectangle.getY();
  w = filterEnvRectangle.getWidth();
  filterEnvLabel->setBounds(x, y+2, w, 16);
  y = filterEnvLabel->getBottom();
  normalLabel->setBounds(x+4, y+4, w-8, 16);
  y += 16;
  normalDecaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  normalAttackSlider->setBounds(x+4, y+4, w-8, 16);
  y += 16;
  accentLabel->setBounds(x+4, y+4, w-8, 16);
  y += 16;
  accentDecaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  accentAttackSlider->setBounds(x+4, y+4, w-8, 16);

  x = ampRectangle.getX();
  y = ampRectangle.getY();
  w = ampRectangle.getWidth();
  ampLabel->setBounds(x, y+2, w, 16);
  y = ampLabel->getBottom();
  ampDecaySlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  ampSustainSlider->setBounds(x+4, y+4, w-8, 16);
  y += 14;
  ampReleaseSlider->setBounds(x+4, y+4, w-8, 16);

  y = distRectangle.getY();
  distLabel->setBounds(x, y+2, w, 16);
  y = distLabel->getBottom();
  distortionDriveSlider->setBounds(x+4, y+4, w-8, 16);


  // todo: set up dist label


  y = globalRectangle.getBottom()-2;
  w = filterRectangle.getRight();
  sequencerEditor->setBounds(0, y, w, 252);


  int noteSize    = sizeof(rosic::AcidNote);     // 12/5    with int/byte
  int patternSize = sizeof(rosic::AcidPattern);  // 208/96

  // rename to "Amplifier" to "Amplifier Envelope", pack the A/D/R sliders densely
  // make a distortion section below the amp env: parameters: drive, shape, DC etc.
  // shape could be several parameters
  // the filte mode and envmod should be a little lower, maybe by 8 pixels ..or maybe only
  // the envmod slider should be lowered

  // what can we do below the preset section? maybe a little scope? or some meters?
}

/*
Ideas:
-the shift functionality for the sequencer should be available separately for accent, glide, 
 octave and notes - currently everything is shifted togther

Todo: move the amplifier secstion under the filter env, maybe ad some more distortion options

*/