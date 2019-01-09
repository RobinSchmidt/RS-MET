ModalSynthAudioModule::ModalSynthAudioModule(CriticalSection *lockToUse)
  : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("ModalSynth"); // find a better name - maybe Modacous?
  createParameters();
}

void ModalSynthAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef rosic::rsModalSynth MS;
  typedef Parameter Param;
  typedef ParameterWithKeyVelScaling ParamKV;
  Param* p;
  //ParamKV* pkv;



  // mode frequency parameters:

  //maxNumModes = p = new Param("MaxNumModes", 10.0, 1024.0, 1024.0, Parameter::EXPONENTIAL, 1.0);
  maxNumModes = p = new Param("MaxNumModes", 10.0, 1024.0, 32.0, Parameter::EXPONENTIAL, 1.0);
  p->setValueChangeCallback<MS>(&core, &MS::setMaxNumPartials);
  addObservedParameter(p);
  // maybe use 32 as defualt value only for debug builds ...or use 1024 always and in debug mode
  // just set it initially to something lower

  ratioProfileTopLeft = p = new Param("RatiosTopLeft", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile1);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("RatiosTopRight", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile2);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("RatiosBottomLeft", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile3);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("RatiosBottomRight", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfile4);
  addObservedParameter(p);

  freqRatiosX = p = new Param("FreqRatiosX", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioMixX);
  addObservedParameter(p);

  freqRatiosY = p = new Param("FreqRatiosY", 0.0, 1.0, 0.5, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioMixY);
  addObservedParameter(p);


  // mode amplitude and envelope parameters:

  amp = p = new Param("Amplitude", 0.01, 10.0, 1.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<MS>(&core, &MS::setAmplitude);
  addObservedParameter(p);

  ampByRatio = p = new Param("AmplitudeByRatio", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setAmplitudeByRatio);
  addObservedParameter(p);

  ampByKey = p = new Param("AmplitudeByKey", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setAmplitudeByKey);
  addObservedParameter(p);

  ampByVel = p = new Param("AmplitudeByVel", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setAmplitudeByVel);
  addObservedParameter(p);





  //attack = pkv = new ParamKV("Attack", 0.1, 100.0, 5.0, Parameter::EXPONENTIAL, 0.0);
  attack = p = new Param("Attack", 0.0, 100.0, 5.0, Parameter::LINEAR, 0.0); // linear seems better for attack
  p->setValueChangeCallback<MS>(&core, &MS::setAttack);
  addObservedParameter(p);

  decay = p = new Param("Decay", 10.0, 10000.0, 5.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<MS>(&core, &MS::setDecay);
  addObservedParameter(p);

}

/*
AudioModuleEditor* ModalSynthAudioModule::createEditor(int type)
{
  return new ModalSynthEditor(this);
}
*/

void ModalSynthAudioModule::processBlock(double **buf, int numChannels, int numSamples)
{
  for(int n = 0; n < numSamples; n++)
    core.getSampleFrameStereo(&buf[0][n], &buf[1][n]);
}

void ModalSynthAudioModule::processStereoFrame(double *left, double *right)
{
  core.getSampleFrameStereo(left, right);
}

void ModalSynthAudioModule::setSampleRate(double newSampleRate)
{
  core.setSampleRate(newSampleRate);
}

void ModalSynthAudioModule::reset()
{
  core.reset();
}

void ModalSynthAudioModule::populateFreqRatioProfileParam(Parameter* p)
{
  // make sure that the order is the same as in the enumerations in rsModalSynth

  // 0-parametric settings:
  p->addStringValue("All Harmonics");
  p->addStringValue("Odd Harmonics");
  p->addStringValue("Pseudo Harmonic 12-TET");
  p->addStringValue("12-TET");
  p->addStringValue("Rod Free/Free");
  p->addStringValue("Rod Free/Clamped");
}

//=================================================================================================

ModalSynthEditor::ModalSynthEditor(jura::ModalSynthAudioModule *newModalSynthModuleToEdit)
  : AudioModuleEditor(newModalSynthModuleToEdit)
{
  ScopedLock scopedLock(*lock);
  modalModule = newModalSynthModuleToEdit;
  createWidgets();
  setSize(500, 460);
}

void ModalSynthEditor::resized()
{

}

void ModalSynthEditor::createWidgets()
{
  typedef RSlider Sld;
  typedef RComboBox Box;
  //typedef RButton Btn;

  Sld* sld;
  Box* box;

  addWidget( box = boxTopLeftRatios = new Box );
  box->assignParameter( modalModule->getParameterByName("RatiosTopLeft") );
  box->setDescription("Mode frequency ratios in top left corner");
  box->setDescriptionField(infoField);
  //box->registerComboBoxObserver(this);  // may be needed later

  addWidget( box = boxTopRightRatios = new Box );
  box->assignParameter( modalModule->getParameterByName("RatiosTopRight") );
  box->setDescription("Mode frequency ratios in top right corner");
  box->setDescriptionField(infoField);
  //box->registerComboBoxObserver(this);  // may be needed later

  addWidget( box = boxBottomLeftRatios = new Box );
  box->assignParameter( modalModule->getParameterByName("RatiosBottomLeft") );
  box->setDescription("Mode frequency ratios in bottom left corner");
  box->setDescriptionField(infoField);
  //box->registerComboBoxObserver(this);  // may be needed later

  addWidget( box = boxBottomRightRatios = new Box );
  box->assignParameter( modalModule->getParameterByName("RatiosBottomRight") );
  box->setDescription("Mode frequency ratios in bottom right corner");
  box->setDescriptionField(infoField);
  //box->registerComboBoxObserver(this);  // may be needed later


  addWidget(xyPadRatios = new rsVectorPad);
  xyPadRatios->assignParameterX(modalModule->getParameterByName("FreqRatiosX"));
  xyPadRatios->assignParameterY(modalModule->getParameterByName("FreqRatiosY"));
  xyPadRatios->setDescription("Vector morph of the frequency ratios in the 4 corners");
  xyPadRatios->setDescriptionField(infoField);
  // it would be nice, if the xy-pad would show a plot of the resulting freq-ratios
  // maybe show the 4 edge freq-ratios in plots next to the middle xy-pad-plot

  addWidget( sld = sldMaxNumModes = new Sld );
  sld->assignParameter( modalModule->getParameterByName("MaxNumModes") );
  sld->setSliderName("MaxNumModes");
  sld->setDescription("Maximum number of modes to be produced (to control CPU load)");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&valueToString0);

}