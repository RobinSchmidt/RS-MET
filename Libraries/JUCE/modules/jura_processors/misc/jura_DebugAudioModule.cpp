DebugAudioModule::DebugAudioModule(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("DebugAudioModule");
  createParameters();
}

void DebugAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  typedef MetaControlledParameter Param;
  Param* p;

  leftParam = p = new Param("Left" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setLeftValue);
  addObservedParameter(p);

  rightParam = p = new Param("Right" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setRightValue);
  addObservedParameter(p);

  smoothParam = p = new Param("Smoothing" , 0.0, 100.0, 0.0, Parameter::LINEAR, 1.0);
  p->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setSmoothingTime);
  addObservedParameter(p);
}

AudioModuleEditor* DebugAudioModule::createEditor()
{
  //return AudioModuleWithMidiIn::createEditor();
  return new DebugModuleEditor(this);
}

void DebugAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int i = 0; i < numChannels; i++)
    for(int n = 0; n < numSamples; n++)
      inOutBuffer[i][n] = values[i];
  // todo: loop over child modules and let them add their value
  // or have one child module with scale and offset parameters
}

void DebugAudioModule::processStereoFrame(double *left, double *right)
{
  *left  = values[0];
  *right = values[1];;
}

void DebugAudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  AudioModule::setSampleRate(newSampleRate);
}

void DebugAudioModule::reset()
{
  ScopedLock scopedLock(*lock);
  AudioModule::reset();
}

void DebugAudioModule::setMidiController(int controllerNumber, float controllerValue)
{
  double v = controllerValue / 127.0; //  0..1
  v = 2*v-1;                          // -1..+1
  if(controllerNumber == 74)
    getParameterByName("Left")->setValue(v, true, true);
  else if(controllerNumber == 71)
    getParameterByName("Right")->setValue(v, true, true);
}

void DebugAudioModule::setLeftValue( double newValue) 
{ 
  values[0] = newValue;  
}

void DebugAudioModule::setRightValue(double newValue) 
{ 
  values[1] = newValue;  
}

void DebugAudioModule::setSmoothingTime(double newTime) 
{
  if(smoothingManager) 
    smoothingManager->setSmoothingTime(newTime);
}

//=================================================================================================

DebugModuleEditor::DebugModuleEditor(jura::DebugAudioModule *newDebugModuleToEdit) 
  : AudioModuleEditor(newDebugModuleToEdit)
{
  ScopedLock scopedLock(*lock);
  debugModule = newDebugModuleToEdit;
  createWidgets();
  setSize(500, 300);
}

void DebugModuleEditor::createWidgets()
{
  typedef AutomatableSlider Sld;
  Sld* s;
  Parameter* p;

  // create vector-pad:
  addWidget(xyPad = new rsVectorPad);

  addWidget( leftSlider = s = new Sld );
  s->assignParameter( p = debugModule->getParameterByName("Left") );
  //xyPad->assignParameterX(p);
  s->setSliderName("Left");
  s->setDescription("Left channel output value");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToStringTotal5);

  addWidget( rightSlider = s = new Sld );
  s->assignParameter( p = debugModule->getParameterByName("Right") );
  //xyPad->assignParameterY(p);
  s->setSliderName("Right");
  s->setDescription("Right channel output value");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToStringTotal5);

  addWidget( smoothSlider = s = new Sld );
  s->assignParameter( p = debugModule->getParameterByName("Smoothing") );
  s->setSliderName("Smoothing");
  s->setDescription("Parameter smoothing in milliseconds");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToStringTotal5);
}

void DebugModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int m  = 4; // margin
  int x  = m;
  int y  = getPresetSectionBottom() + m;
  int w  = getWidth() - 2*m;
  int h  = getHeight();
  int wh = 16;     // widget height
  int dy = wh-2;   // delta-y between widgets

  int size = jmin(getWidth(), getHeight()-y);
  xyPad->setBounds(0, y, size, size);


  x = xyPad->getRight() + m;
  w = getWidth() - x - m;
  y = getPresetSectionBottom() + m;
  leftSlider  ->setBounds(x, y, w, wh); y += dy;
  rightSlider ->setBounds(x, y, w, wh); y += dy;
  smoothSlider->setBounds(x, y, w, wh); y += dy;
}