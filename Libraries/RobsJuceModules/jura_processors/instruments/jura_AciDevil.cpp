
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
  //p = new Param("Cutoff", 20.0, 10000.0, 300.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<AD>(ad, &AD::setCutoff);
  addObservedParameter(p);
  // The cutoff should have a lower minimum setting. Just changing the parameter's minimum from 200
  // to 20, for example, doesn't seem to work though. Somewhere below 200, the cutoff stops 
  // responding to further changes as if there's some sort of cutoff = max(cutoff, 180) or 
  // something going on. Check in the DSP code, if we limit the cutoff range there and try to do 
  // something about it. When we extend the range in an update, we really need to think about how 
  // to keep the plugin backward compatible with older automation data. Maybe we could have 
  // different behavioral modes that can be set up in the global section (like "old behavior" and
  // new behavior"). Maybe we could also create a general framework that allows the user to modify
  // parameter ranges. Maybe with a right click on a slider, we could give options like 
  // setMin/setMax and we store these values in the patch data as e.g. CutoffMin, CutoffMax, etc.
  // Maybe we should create a subclass of Parameter for that - maybe VariableRangeParameter. It 
  // should save the min/max values only when they are different from their defaults because 
  // otherwise we would blow up the patch data too much. It would imply that we can't change these
  // defaults in later updates but that seems to be a reasonable compromise. Implementing such a 
  // solution is considerable work but doing it would also solve any similar problems in the 
  // future. The desire to extend parameter ranges in later updates might be a common enough problem
  // to warrant for such a solution.


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
  //setHeadlineStyle(MAIN_HEADLINE);  // old

  setHeadlineStyle(SUB_HEADLINE); // Aletnative - allows the GUI to be smaller - maybe we should do it
  setPresetSectionPosition(BELOW_HEADLINE);

  // Assign the pointer to the rosic::AciDevil object to be used as aduio engine:
  jassert(newAciDevilAudioModule != NULL ); // you must pass a valid module here
  aciDevilModuleToEdit = newAciDevilAudioModule;


  createWidgets();
  updateWidgetsAccordingToState();
  //setSize(772, 394);
  //setSize(634, 394);
  setSize(634, 370);
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
  l->setItemDescription("Global parameters");
  l->setDescriptionField(infoField);

  addWidget( masterLevelSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("MasterLevel") );
  s->setSliderName("Level");
  s->setItemDescription("Master level in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( accentSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Accent") );
  s->setItemDescription("Accent in percent");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( slideTimeSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("SlideTime") );
  s->setSliderName("Slide");
  s->setItemDescription("Slide time in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( oscLabel = new Lbl("Oscillator") );
  oscLabel->setJustification(Justification::centred);
  oscLabel->setItemDescription("Oscillator parameters");
  oscLabel->setDescriptionField(infoField);

  addWidget( waveformSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Waveform") );
  s->setSliderName("Saw/Pulse");
  s->setItemDescription("Mix between saw- and pulse-wave for main oscillator");
  s->setStringConversionFunction(ratioToString0);
  s->setDescriptionField(infoField);

  addWidget( pulseWidthSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("PulseWidth") );
  s->setSliderName("PulseWidth");
  s->setItemDescription("Width of high section of the rectangular pulse waveform");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( subOscLabel = new Lbl("SubOsc:") );
  subOscLabel->setJustification(Justification::centredLeft);
  subOscLabel->setItemDescription("Sub-oscillator settings");
  subOscLabel->setDescriptionField(infoField);

  addWidget( subOscLevelSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("SubOscLevel") );
  s->setSliderName("Level");
  s->setItemDescription("Sub-oscillator level in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( subOscWaveformSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("SubOscWaveform") );
  s->setSliderName("Saw/Pulse");
  s->setItemDescription("Mix between saw- and pulse-wave for suboscillator");
  s->setStringConversionFunction(ratioToString0);
  s->setDescriptionField(infoField);


  addWidget( filterLabel = new Lbl("Filter") );
  filterLabel->setJustification(Justification::centred);
  filterLabel->setItemDescription("Filter");
  filterLabel->setDescriptionField(infoField);

  addWidget( cutoffSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Cutoff") );
  s->setItemDescription("Filter cutoff frequency in Hz");
  s->setStringConversionFunction(hertzToStringWithUnitTotal5);
  s->setDescriptionField(infoField);

  addWidget( resonanceSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("Resonance") );
  s->setItemDescription("Resonance in percent");
  s->setStringConversionFunction(percentToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( filterModeLabel = new Lbl("Mode:") );
  filterModeLabel->setJustification(Justification::centredLeft);
  filterModeLabel->setItemDescription("Choose the filter mode");
  filterModeLabel->setItemDescription("Mode:");
  filterModeLabel->setDescriptionField(infoField);

  addWidget( filterModeBox = c = new Box );
  c->assignParameter( moduleToEdit->getParameterByName("FilterMode") );
  c->setItemDescription(filterModeLabel->getItemDescription());
  c->setDescriptionField(infoField);

  addWidget( envModSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("EnvMod") );
  s->setItemDescription("Amount of modulation of cutoff frequency in semitones");
  s->setStringConversionFunction(semitonesToStringWithUnit1);
  s->setDescriptionField(infoField);


  addWidget( filterEnvLabel = new Lbl("Filter Envelope") );
  filterEnvLabel->setJustification(Justification::centred);
  filterEnvLabel->setItemDescription("Filter envelope parameters");
  filterEnvLabel->setDescriptionField(infoField);

  addWidget( normalLabel = new Lbl("Normal:") );
  normalLabel->setJustification(Justification::centredLeft);
  normalLabel->setItemDescription("Time values for normal (un-accented) notes");
  normalLabel->setDescriptionField(infoField);

  addWidget( normalDecaySlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("NormalDecay") );
  s->setSliderName("Decay");
  s->setItemDescription("Decay time for normal (un-accented) notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( normalAttackSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("NormalAttack") );
  s->setSliderName("Attack");
  s->setItemDescription("Attack time for normal (un-accented) notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);


  addWidget( accentLabel = new Lbl("Accent:") );
  accentLabel->setJustification(Justification::centredLeft);
  accentLabel->setItemDescription("Time values for accented notes");
  accentLabel->setDescriptionField(infoField);

  addWidget( accentDecaySlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AccentDecay") );
  s->setSliderName("Decay");
  s->setItemDescription("Decay time for accented notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( accentAttackSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AccentAttack") );
  s->setSliderName("Attack");
  s->setItemDescription("Attack time for accented notes in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( ampLabel = new Lbl("Amp Envelope") );
  ampLabel->setJustification(Justification::centred);
  ampLabel->setItemDescription("Amplide envelope and distortion parameters");
  ampLabel->setDescriptionField(infoField);

  addWidget( ampDecaySlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AmpDecay") );
  s->setSliderName("Decay");
  s->setItemDescription("Decay time for amplitude envelope in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);

  addWidget( ampSustainSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AmpSustain") );
  s->setSliderName("Sustain");
  s->setItemDescription("Sustain level for amplitude envelope in decibels");
  s->setStringConversionFunction(decibelsToStringWithUnit1);
  s->setDescriptionField(infoField);

  addWidget( ampReleaseSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("AmpRelease") );
  s->setSliderName("Release");
  s->setItemDescription("Release time for amplitude envelope in milliseconds");
  s->setStringConversionFunction(millisecondsToStringWithUnit2);
  s->setDescriptionField(infoField);


  addWidget( distLabel = new Lbl("Distortion") );
  distLabel->setJustification(Justification::centred);
  distLabel->setItemDescription("Distortion Settings");
  distLabel->setDescriptionField(infoField);

  addWidget( distortionDriveSlider = s = new Sld );
  s->assignParameter( aciDevilModuleToEdit->getParameterByName("DistortionDrive") );
  s->setSliderName("Drive");
  s->setItemDescription("Drive for distortion unit in decibels");
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

  //y = getHeadlineBottom()+4;  // results in y = 24


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
  //y = globalRectangle.getY();
  y = getHeadlineBottom()+4;
  w = globalRectangle.getWidth();
  stateWidgetSet->setLayout(StateLoadSaveWidgetSet::LABEL_AND_BUTTONS_ABOVE);
  stateWidgetSet->setBounds(x+4, y+4, w-8, 32);


  y = stateWidgetSet->getBottom()+4;

  //y = stateWidgetSet->getBottom() + 4 + 32; // leave space for tuning widgets

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

-Give the oscillators a start-phase parameter and also a polarity switch (both for main and subosc 
 separately). That's potentially important when mixing the 303 basslines with bassdrums and/or bass
 sounds. Maybe the subosc should have and adjustable detune. That would bring us actually into 2 
 osc subtractive synth territory.

-Maybe instead of just blending between Saw/Pulse, let the user load custom waveforms.

-Provide a pulse-width parameter. Maybe it should depend on osc frequency via some other parameter
 like PulseWidthByFreq.

-The cutoff should have a lower minimum setting. See comment in createParameters().

-The filter needs some more work: the resonance drops too much towards higher cutoffs for the
 mid/high resonance range - at the upper end, it seems better, but slightly below, the behavior is
 suboptimal.

-Maybe implement a simple undo/redo mechanism by keeping track of the applied transformations

-Maybe make a smaller headline not at the top but at the top-left - move preset section some 20 
 pixels down and say Slot1-AcidDevil - similar to FuncShaper's GUI and all the other, "smaller"
 plugins, see AudioModuleEditor::setHeadlineStyle

-The sequencer should have a pair of Undo/Redo buttons They should be labeled like 
 "Undo (3)" "Redo (5)" indicating that currently there are 3 possible undo steps available and
  5 possible redo steps.

- The sequencer should have an "Export" button that let's the user write a .mid file with the 
  current pattern. See:
  https://docs.juce.com/master/classMidiFile.html
  https://docs.juce.com/master/classMidiMessageSequence.html
  https://www.youtube.com/watch?v=P27ml4M3V7A
 
- Maybe the Export function can be integrated into the the "Save" button. If the user chooses to
  save to a .mid file, it will be exported. We don't need to modify the GUI for that. But then, we
  also need to implement midi import in Load.
  
- Accented and non-accented notes should use two well defined Velocity levels. 
  Maybe 64 and 127. Or maybe 50 and 100 to have some headroom. Or maybe 40 and 80. Maybe AcidDevil 
  should have a switch for midi interpeting accents binary (on off based on velocity threshold) or 
  continuously.

- Maybe make important sliders bigger (20 or 24 pixel high): Level (maybe), Saw/Pulse, Cutoff, 
  Drive, EnvMod, Resonance ...those that the user is likely to automate. Accent is actually also a
  good automation target

 -The distortion unit should get a mode, maybe pre/post filters and and some manipluators for the 
  transfer function and perhaps a little display for the function...or maybe not. But it should 
  have a DC parameter. And maye it should be a softclipper with afjustable hardness. Maybe for 
  that, the smooth-crossfade function in the research codebase in testSmoothCrossFade() could be
  used. Maybe experiment with that in the context of the sampler engine. The next milestone is to
  get some more serious DSP algos going anyway.

- It would be nice if the filter could be morphable. I'd really like to be able to morph through: 
  LP24 -> LP18 -> LP12 -> LP6 -> Flat -> HP6 -> HP12 -> HP18 -> HP24. Especially the range
  LP18 -> LP12 -> LP6 is interesting. I tend to like the sound of LP12 most. LP 18 is too dull, LP6
  too bright. 12 is the sweet spot - but it might be nice to adjust it more finely. Maybe something
  like LP11, LP10 etc. may sound even better in some cases. How can we do this? Maybe crossfade 
  between the coefficient sets? Experiment a bit with a LP6/LP12 crossfade. Or maybe we could try 
  post filtering with an adjustable slope/tilt filter?


*/