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
  maxNumModes = p = new Param("MaxNumModes", 1.0, 1024.0, 32.0, Parameter::EXPONENTIAL, 1.0);
  p->setValueChangeCallback<MS>(&core, &MS::setMaxNumPartials);
  addObservedParameter(p);
  // maybe use 32 as defualt value only for debug builds ...or use 1024 always and in debug mode
  // just set it initially to something lower

  ratioProfileTopLeft = p = new Param("RatiosTopLeft", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfileTopLeft);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("RatiosTopRight", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfileTopRight);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("RatiosBottomLeft", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfileBottomLeft);
  addObservedParameter(p);

  ratioProfileTopLeft = p = new Param("RatiosBottomRight", 0, 1, 0, Parameter::STRING, 1);
  populateFreqRatioProfileParam(p);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioProfileBottomRight);
  addObservedParameter(p);

  freqRatiosX = p = new Param("FreqRatiosX", -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setFreqRatioMixX);
  addObservedParameter(p);

  freqRatiosY = p = new Param("FreqRatiosY", -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
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


  attack = p = new Param("Attack", 0.0, 100.0, 10.0, Parameter::LINEAR, 0.0); // linear seems better for attack
  p->setValueChangeCallback<MS>(&core, &MS::setAttack);
  addObservedParameter(p);

  attackByRatio = p = new Param("AttackByRatio", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setAttackByRatio);
  addObservedParameter(p);

  attackByKey = p = new Param("AttackByKey", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setAttackByKey);
  addObservedParameter(p);

  attackByVel = p = new Param("AttackByVel", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setAttackByVel);
  addObservedParameter(p);


  decay = p = new Param("Decay", 10.0, 10000.0, 1000.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<MS>(&core, &MS::setDecay);
  addObservedParameter(p);

  decayByRatio = p = new Param("DecayByRatio", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setDecayByRatio);
  addObservedParameter(p);

  decayByKey = p = new Param("DecayByKey", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setDecayByKey);
  addObservedParameter(p);

  decayByVel = p = new Param("DecayByVel", -200, 100, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MS>(&core, &MS::setDecayByVel);
  addObservedParameter(p);






}

AudioModuleEditor* ModalSynthAudioModule::createEditor(int type)
{
  return new ModalSynthEditor(this);
}

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
  setSize(480, 400);
}

void ModalSynthEditor::resized()
{
  AudioModuleEditor::resized();
  int m  = 4; // margin
  int y  = getPresetSectionBottom() + m;
  int wh = 16;     // widget height
  int w  = getWidth() / 3;
  int x  = 0;
  int xyPadSize = w;

  // freq-ratio widgets:
  boxTopLeftRatios->setBounds(      0, y, w, wh);
  boxTopRightRatios->setBounds(   2*w, y, w, wh);
  xyPadRatios->setBounds(           w, y, w, w);
  y += w-wh;
  boxBottomLeftRatios->setBounds(   0, y, w, wh);
  boxBottomRightRatios->setBounds(2*w, y, w, wh);
  y = xyPadRatios->getY() + (w-wh)/2;
  sldRatiosX->setBounds(    0+m, y, w-2*m, wh);
  sldRatiosY->setBounds(  2*w+m, y, w-2*m, wh);
  x = xyPadRatios->getX();
  y = xyPadRatios->getBottom();
  sldMaxNumModes->setBounds(x, y, w, wh);

  // amplitude:
  y = sldMaxNumModes->getBottom() + m;
  x = m; 
  int sw = getWidth()/2 - 2*m;   // slider width
  int sw2 = (sw-12) / 3;   // width of the attached Rat/Key/Vel sliders
  int xk = x + sw/3 + 4;
  int xv = x + sw - sw2;
  int dy = wh-2;

  sldAmp       ->setBounds(x,  y, sw,  wh);  y += dy;
  sldAmpByRatio->setBounds(x,  y, sw2, wh);
  sldAmpByKey  ->setBounds(xk, y, sw2, wh);
  sldAmpByVel  ->setBounds(xv, y, sw2, wh);
  int dy2 = dy+2*m;
  y += dy2;
  sldAttack       ->setBounds(x,  y, sw,  wh);  y += dy;
  sldAttackByRatio->setBounds(x,  y, sw2, wh);
  sldAttackByKey  ->setBounds(xk, y, sw2, wh);
  sldAttackByVel  ->setBounds(xv, y, sw2, wh);
  y += dy2;
  sldDecay       ->setBounds(x,  y, sw,  wh);  y += dy;
  sldDecayByRatio->setBounds(x,  y, sw2, wh);
  sldDecayByKey  ->setBounds(xk, y, sw2, wh);
  sldDecayByVel  ->setBounds(xv, y, sw2, wh);








  
  /* Layout ideas for the amplitude/envelope widgets:

   Attack    Amplitude   AttackScale    
   R K V     R K V       DecayScale
   Decay     Blend       FreqDelta 
   R K V     R K V       PhaseDelta


   Attack      Amplitude   Decay
   R K V       R K V       R K V

   AttScl      Blend       DecScl
   FreqDelta   R K V       PhaseDelta


   Amplitude       FreqScale
   R K V           R K V
   Attack          Decay
   R K V           R K V
   AttScl          DecScl
            Blend


  PhaseRandomness

  PhaseSeed

  PhaseAlternate 

  mmhh...but actually, there are much more amplitude and decay- related widgets to come, so it
  doesn't make much sense to think too much about layout now

  above vector pad:

  Level   Detune   Tuning
  K   V   K    V   12-TET

  below vector pad:

  Attack     AmpSlope    AttScl
  R  K V     K      V    DecScl
  Decay       Blend      FreqDelta
  R K V       R K V      R     K V





  AmpCutoffHP     AmpCutoffSlopeHP
  K  V            K  V
  AmpCutoffLP     AmpCutoffSlopeLP
  K  V            K  V


  PhaseRandomSeed, PhaseRandomness, PhaseDelta, etc....


  */



  //int size = jmin(getWidth(), getHeight()-y);
  //xyPad->setBounds(0, y, xyPadSize, xyPadSize);

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

  addWidget( sld = sldRatiosX = new Sld );
  sld->assignParameter( modalModule->getParameterByName("FreqRatiosX") );
  sld->setSliderName("RatiosX");
  sld->setDescription("X-coordinate of the ratio vector morph");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&valueToStringTotal5);

  addWidget( sld = sldRatiosY = new Sld );
  sld->assignParameter( modalModule->getParameterByName("FreqRatiosY") );
  sld->setSliderName("RatiosY");
  sld->setDescription("Y-coordinate of the ratio vector morph");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&valueToStringTotal5);

  addWidget( sld = sldMaxNumModes = new Sld );
  sld->assignParameter( modalModule->getParameterByName("MaxNumModes") );
  sld->setSliderName("MaxNumModes");
  sld->setDescription("Maximum number of modes to be produced (to control CPU load)");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&valueToString0);


  // amplitude:

  addWidget( sld = sldAmp = new Sld );
  sld->assignParameter( modalModule->getParameterByName("Amplitude") );
  sld->setSliderName("Amplitude");
  sld->setDescription("Overall amplitude");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&valueToStringTotal5);

  addWidget( sld = sldAmpByRatio = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AmplitudeByRatio") );
  sld->setSliderName("R");
  sld->setDescription("Mode amplitude dependency on frequency ratio");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldAmpByKey = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AmplitudeByKey") );
  sld->setSliderName("K");
  sld->setDescription("Amplitude dependency on key");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldAmpByVel = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AmplitudeByVel") );
  sld->setSliderName("V");
  sld->setDescription("Amplitude dependency on velocity");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);



  addWidget( sld = sldAttack = new Sld );
  sld->assignParameter( modalModule->getParameterByName("Attack") );
  sld->setSliderName("Attack");
  sld->setDescription("Attack time");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&millisecondsToStringWithUnit2);

  addWidget( sld = sldAttackByRatio = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AttackByRatio") );
  sld->setSliderName("R");
  sld->setDescription("Attack dependency on frequency ratio");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldAttackByKey = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AttackByKey") );
  sld->setSliderName("K");
  sld->setDescription("Attack dependency on key");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldAttackByVel = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AttackByVel") );
  sld->setSliderName("V");
  sld->setDescription("Attack dependency on velocity");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);



  addWidget( sld = sldDecay = new Sld );
  sld->assignParameter( modalModule->getParameterByName("Decay") );
  sld->setSliderName("Decay");
  sld->setDescription("Decay time");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&millisecondsToStringWithUnit2);

  addWidget( sld = sldDecayByRatio = new Sld );
  sld->assignParameter( modalModule->getParameterByName("DecayByRatio") );
  sld->setSliderName("R");
  sld->setDescription("Decay dependency on frequency ratio");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldDecayByKey = new Sld );
  sld->assignParameter( modalModule->getParameterByName("DecayByKey") );
  sld->setSliderName("K");
  sld->setDescription("Decay dependency on key");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldDecayByVel = new Sld );
  sld->assignParameter( modalModule->getParameterByName("DecayByVel") );
  sld->setSliderName("V");
  sld->setDescription("Decay dependency on velocity");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);






  // envelope:
  /*
  addWidget( sld = sldAmp = new Sld );
  sld->assignParameter( modalModule->getParameterByName("Amplitude") );
  sld->setSliderName("Amplitude");
  sld->setDescription("Overall amplitude");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&valueToStringTotal5);

  addWidget( sld = sldAmpByRatio = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AmplitudeByRatio") );
  sld->setSliderName("R");
  sld->setDescription("Mode amplitude dependency on frequency ratio");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldAmpByKey = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AmplitudeByKey") );
  sld->setSliderName("K");
  sld->setDescription("Amplitude dependency on key");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);

  addWidget( sld = sldAmpByVel = new Sld );
  sld->assignParameter( modalModule->getParameterByName("AmplitudeByVel") );
  sld->setSliderName("V");
  sld->setDescription("Amplitude dependency on velocity");
  sld->setDescriptionField(infoField);
  sld->setStringConversionFunction(&percentToStringWithUnit2);
  */




}