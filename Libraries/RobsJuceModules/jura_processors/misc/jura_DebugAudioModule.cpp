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

DebugModuleEditor::~DebugModuleEditor()
{
  delete popupRect;
  delete popupComponent;
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




  addWidget( popupButton1 = new RButton("Popup 1") );
  popupButton1->addRButtonListener(this);
  popupButton1->setDescription("Open a stack-allocated juce::PopupMenu");
  popupButton1->setClickingTogglesState(false);



  popupRect = new jura::RectangleComponent();  // is explicitly deleted in our destructor
  popupRect->setAlwaysOnTop(true);
  popupRect->setSize(400, 300);
  popupRect->setFillColour(Colours::black);

  addWidget( popupButton2 = new RButton("Popup 2") );
  popupButton2->addRButtonListener(this);
  popupButton2->setDescription("Open/close an owned jura::RectangleComponent");
  popupButton2->setClickingTogglesState(true);



  popupContent = new jura::RectangleComponent();      // will be owned by popupComponent
  popupContent->setSize(400, 300);
  popupContent->setFillColour(Colours::black);

  popupComponent = new jura::RPopUpComponent();       // is explicitly deleted in our destructor
  popupComponent->setContentComponent(popupContent);
  // We do *not* add popupComponent as child component via addChildComponent(popupComponent) 
  // because we want to add it to the Desktop later and doing so would remove it from the child 
  // components anyway. We need to explicitly delete it in our destructor because we cannot rely on
  // the deleteAllChildren call that our basclass constructor does. The popupContent, however,
  // does not need to be explicitly deleted because calling 
  // popupComponent->setContentComponent(popupContent); transfers ownership of the popupContent to
  // the popupComponent. Also, using addChildComponent would make the popup appear in the wrong 
  // place when it's opened for the first time.


  addWidget( popupButton3 = new RButton("Popup 3") );
  popupButton3->addRButtonListener(this);
  popupButton3->setDescription("Open/close an owned jura::RPopUpComponent");
  popupButton3->setClickingTogglesState(true);



  int dummy = 0;
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
  int bw = 64;  // button width

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

  popupButton1->setBounds(x, y, bw, wh);
  x += bw + 4;
  popupButton2->setBounds(x, y, bw, wh);
  x += bw + 4;
  popupButton3->setBounds(x, y, bw, wh);
  y += dy;



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
  if(button == popupButton1)
  {
    // See: https://docs.juce.com/master/classPopupMenu.html
    // https://docs.juce.com/master/classPopupMenu.html#details
    PopupMenu m;
    m.addItem(1, "Item 1");
    m.addItem(2, "Item 2");
    m.addItem(3, "Item 3");
    m.addItem(4, "Item 4");
    m.showMenuAsync(PopupMenu::Options(), [](int result)
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
        // ...
      }
    );
  }
  // OK - the popup menu shows up in Tracktion just fine. We now need to figure out, how that works
  // and adapt WaveOscEditorContextMenu accrodingly. But there's a difference. the context menu for
  // the osc whould not be modal. It should be persistent until the user closes it manually by 
  // either using the close button on the menu iteself or the "More" button on the osc editor which
  // toggles the context menu on and and off. Let's see what happens:
  //
  //   showMenuAsync  calls  showWithOptionalCallback  calls  createWindow
  //
  // These are all member functions of juce::PopupMenu. createWindow returns a 
  //
  //   new HelperClasses::MenuWindow(...)
  //
  // and I currently don't see where this is ever deleted. Naively looking at the code, this looks
  // like a memory leak to me. However, MenuWindow is a subclass of Component.
  //
  // Try it with RPopupMenu - I actually already know that this works - figure out how that handles
  // the addToDesktop/show/toFront etc. business. We have class  RPopUpComponent  which is a general
  // component that can pop up. Maybe try to show a rectangle component using that, i.e. set the
  // content component of a RPopUpComponent to a RecatngelComponent

  // OK - this rectangle component exposes the same behavior as the context menu in the osc editor:
  if(button == popupButton2)
  {
    if(popupButton2->getToggleState() == true)
    {
      //int x = getScreenX() + getWidth();
      //int y = getScreenY() - 102;
      int x = popupButton2->getScreenX();
      int y = popupButton2->getScreenY() + 16;
      popupRect->setTopLeftPosition(x, y);
      popupRect->addToDesktop(
        ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);


      //popupRect->setBroughtToFrontOnMouseClick(false);  // Test
      // Interesting: Even when calling  setBroughtToFrontOnMouseClick(false)   mouse clicks will 
      // still bring it to the front.


      popupRect->setVisible(true);
      popupRect->toFront(true);


      // Next test:
      juce::ComponentPeer* popupPeer = popupRect->getPeer();
      if(popupPeer)
      {
        popupPeer->setAlwaysOnTop(true);
        popupPeer->toFront(false);  // false: do not take keyboard focus
        // ...nope - this also doesn't help. 
        // There are also methods like ComponentPeer::setAlwaysOnTop etc. - Check them out, too.
      }



      //popupRect->grabKeyboardFocus();   // Test - nope - doesn't help

      //juce::MouseEvent me;
      //popupRect->mouseDown(me);
    }
    else
      popupRect->setVisible(false);
  }
  // Due to the fact that it is large enough, parts of it are visible - and when clicking on this 
  // part, it actually comes to the front. Could we perhaps somehow trigger the same action that is
  // triggered from mouse-clicks? Maybe some sort of grabFocus()? Check
  // setBroughtToFrontOnMouseClick(), setExplicitFocusOrder(), setFocusContainerType(), 
  // findFocusContainer()  - maybe we could just send a dummy mouseDown event to it? That would be 
  // an ugly hack, though. Generating mock mouse events seems not to be easy, though:
  // https://forum.juce.com/t/creating-a-mouseevent/32635/11
  // https://forum.juce.com/t/simulate-mouse-actions/29362/4
  // https://forum.juce.com/t/my-solution-for-key-mouse-remote-events-handling/37698

  // try setting  setBroughtToFrontOnMouseClick(false)  and check, if mouse click will then not 
  // bring it to front anymore. If not, check what actually happens inside the mouse-handler that 
  // brings it to front. OK - done - actually, calling setBroughtToFrontOnMouseClick(false); does 
  // not prevent it from coming to the front on mouse clicks


  if(button == popupButton3)
  {
    if(popupButton3->getToggleState() == true)
    {
      bool showModally = false;
      int x = popupButton3->getScreenX();
      int y = popupButton3->getScreenY() + 16;
      popupComponent->showAt(showModally, x, y, 400, 300);
      // ToDo: try to let it determine the desired size by itself by using the size of the 
      // contentComponent. Has it to do with the poisiton of the content component? Nope. Using
      // setBounds instead of setSize during creation doe not help.

      int dummy = 0;
    }
    else
    {

    }
  }
  // When clicking the button the 1st time, the window appears in a wrong position.
  // In Tracktion, this popup doesn not appear at all - not even in the background!
  // When using showModally = true; we get a memory leak on exit. But: when doing that, the 
  // rectangle actually does show up in Tracktion - but in the background.
  // Oh - the leak seems to happen even in non-modal mode. Could it be that adding a Component
  // to the desktop, removes it from the child components? I think so.
  // Oh - not calling addChildComponent also fixes the problem with the position of apperance.




  int dummy = 0;
}