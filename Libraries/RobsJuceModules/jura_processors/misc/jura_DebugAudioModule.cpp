DebugAudioModule::DebugAudioModule(CriticalSection *lockToUse) : AudioModuleWithMidiIn(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("DebugAudioModule");
  addChildAudioModule(eqModule = new EqualizerAudioModule(lock));
  createParameters();
}

void DebugAudioModule::createParameters()
{
  ScopedLock scopedLock(*lock);

  //typedef Parameter Param;
  //typedef rsSmoothableParameter Param;
  //typedef MetaControlledParameter Param;
  typedef ModulatableParameter Param; // has wrong value when being smoothed
  Param* p;

  leftParam = p = new Param("Left" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setLeftValue);
  addObservedParameter(p);

  rightParam = p = new Param("Right" , -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setRightValue);
  addObservedParameter(p);

  // Smoothing itself should not be smoothed:
  smoothParam = new Parameter("Smoothing" , 0.0, 1000.0, 0.0, Parameter::LINEAR, 1.0);
  smoothParam->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setSmoothingTime);
  addObservedParameter(smoothParam);


  testParam = p = new Param("Test", -0.99999, +0.99999, 0);  // workaround
  //testParam = p = new Param("Test", -1, +1, 0);
  p->setMapper(new rsParameterMapperTanh(-1, +1, 5));
  p->setValueChangeCallback<DebugAudioModule>(this, &DebugAudioModule::setTestParam);
  addObservedParameter(p);
}

AudioModuleEditor* DebugAudioModule::createEditor(int type)
{
  //return AudioModuleWithMidiIn::createEditor();
  return new DebugModuleEditor(this);
}

void DebugAudioModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  for(int i = 0; i < numChannels; i++)
    for(int n = 0; n < numSamples; n++)
      inOutBuffer[i][n] += values[i];
  if(eqModule)
    eqModule->processBlock(inOutBuffer, numChannels, numSamples);
}

void DebugAudioModule::processStereoFrame(double *left, double *right)
{
  *left  += values[0];
  *right += values[1];
  if(eqModule)
    eqModule->processStereoFrame(left, right);
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
  setSmoothingTime(leftParam,  newTime);
  setSmoothingTime(rightParam, newTime);
  //if(smoothingManager) 
  //  smoothingManager->setSmoothingTime(newTime);
}

void DebugAudioModule::setSmoothingTime(Parameter* p, double newTime)
{
  rsSmoothableParameter *sp = dynamic_cast<rsSmoothableParameter*>(p);
  if(sp)
    sp->setSmoothingTime(newTime);
}

void DebugAudioModule::setTestParam(double newValue)
{
  double value = newValue;
}

//=================================================================================================

DebugModuleEditor::DebugModuleEditor(jura::DebugAudioModule *newDebugModuleToEdit) 
  : AudioModuleEditor(newDebugModuleToEdit)
{
  ScopedLock scopedLock(*lock);
  debugModule = newDebugModuleToEdit;
  if(debugModule->eqModule != nullptr)
    addChildEditor(eqEditor = new EqualizerModuleEditor(lock, debugModule->eqModule));
  createWidgets();
  setSize(500, 460);
}

void DebugModuleEditor::createWidgets()
{
  typedef rsModulatableSlider Sld;
  Sld* s;
  Parameter* p;

  // create vector-pad and node-editor:
  addWidget(xyPad = new rsVectorPad);
  //addWidget(nodeEditor = new rsNodeEditor);

  addWidget( leftSlider = s = new Sld );
  s->assignParameter( p = debugModule->getParameterByName("Left") );
  xyPad->assignParameterX(p);
  s->setSliderName("Left");
  s->setDescription("Left channel output value");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToStringTotal5);

  addWidget( rightSlider = s = new Sld );
  s->assignParameter( p = debugModule->getParameterByName("Right") );
  xyPad->assignParameterY(p);
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

  addWidget( testSlider = s = new Sld );
  s->assignParameter( p = debugModule->getParameterByName("Test") );
  s->setSliderName("Test");
  s->setDescription("Test Parameter");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToStringTotal5);

  addWidget( popupButton = new RButton("Popup") );
  popupButton->addRButtonListener(this);
  popupButton->setDescription("Open/close popup menu");
  popupButton->setClickingTogglesState(false);
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

  int xyPadSize = 200;
  //int size = jmin(getWidth(), getHeight()-y);
  xyPad->setBounds(0, y, xyPadSize, xyPadSize);

  x = xyPad->getRight() + m;
  w = getWidth() - x - m;
  y = getPresetSectionBottom() + m;
  leftSlider  ->setBounds(x, y, w, wh); y += dy;
  rightSlider ->setBounds(x, y, w, wh); y += dy;
  testSlider  ->setBounds(x, y, w, wh); y += dy;
  //testSlider  ->setBounds(x, y, 30, wh); y += dy; // for testing the text entry field
  smoothSlider->setBounds(x, y, w, wh); y += dy;

  popupButton->setBounds(x, y, 64, wh); y += dy;

  if(nodeEditor)
  {
    // preliminary - it overlaps with eq-editor
    x = xyPad->getRight() + m;
    w = getWidth() - x - m;
    y = smoothSlider->getBottom() + m;
    h = w;
    nodeEditor->setBounds(x, y, w, h);
  }

  if(eqEditor)
  {
    y = xyPad->getBottom();
    h = getHeight() - y;
    eqEditor->setBounds(0, y, getWidth(), h);
  }
}

void DebugModuleEditor::rButtonClicked(RButton* button)
{
  // We use this to test, if it works even for a plugin GUI - because my oscillator context menu 
  // does not work - at least not in Tracktion - it's hidden behind the main GUI. 
  if(button == popupButton)
  {
    // See: https://docs.juce.com/master/classPopupMenu.html
    // https://docs.juce.com/master/classPopupMenu.html#details
    PopupMenu m;
    m.addItem(1, "Item 1");
    m.addItem(2, "Item 2");
    m.addItem(3, "Item 3");
    m.addItem(4, "Item 4");
    m.showMenuAsync (PopupMenu::Options(), [](int result)
      {
        if(result == 0)
        {
          // user dismissed the menu without picking anything
        }
        else if(result == 1)
        {
          // user picked item 1
        }
        else if(result == 2)
        {
          // user picked item 2
        }
      }
    );
  }


  int dummy = 0;
}