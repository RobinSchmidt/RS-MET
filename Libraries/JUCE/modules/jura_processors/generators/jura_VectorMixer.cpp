
//-------------------------------------------------------------------------------------------------
// construction/destruction:

VectorMixerAudioModule::VectorMixerAudioModule(CriticalSection *newPlugInLock, 
  rosic::VectorMixer *newVectorMixerToWrap)
: AudioModule(newPlugInLock)
{
  jassert( newVectorMixerToWrap != NULL ); // you must pass a valid rosic-object to the constructor
  wrappedVectorMixer = newVectorMixerToWrap;
  setModuleTypeName("VectorMixer");
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// automation:

void VectorMixerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedVectorMixer == NULL )
    return;

  double value = parameterThatHasChanged->getValue();
  switch( getIndexOfParameter(parameterThatHasChanged) )
  {
  case  0: wrappedVectorMixer->setX(value);  break;
  case  1: wrappedVectorMixer->setY(value);  break;
  } // end of switch( parameterIndex )
}

/*
//-------------------------------------------------------------------------------------------------
// state saving and recall:

XmlElement* VectorMixerAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  // store the parameters of the underlying core object:
  if( wrappedVectorMixer != NULL )
    xmlState = vectorPadStateToXml(wrappedVectorMixer, xmlState);

  return xmlState;
}

void VectorMixerAudioModule::setStateFromXml(const XmlElement& xmlState,
                                             const juce::String& stateName, bool markAsClean)
{
  // restore the inherited controller mappings:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  // restore the parameters of the underlying core object:
  if( wrappedVectorMixer != NULL )
    vectorPadStateFromXml(wrappedVectorMixer, xmlState);
}
*/

//-------------------------------------------------------------------------------------------------
// internal functions:

void VectorMixerAudioModule::initializeAutomatableParameters()
{
  // create the automatable parameters and add them to the list - note that the order of the adds
  // is important because in parameterChanged(), the index (position in the array) will be used to
  // identify which particlua parameter has changed.

  // this pointer will be used to temporarily store the addresses of the created 
  // AutomatableParameter-objects:
  AutomatableParameter* p;

  // #000:
  p = new AutomatableParameter(lock, "X", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #001:
  p = new AutomatableParameter(lock, "Y", -1.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to setValue for each parameter in order to set up all the slave voices:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================


VectorMixerPad::VectorMixerPad(rosic::VectorMixer* newVectorMixerToEdit, const juce::String& name) 
  : CoordinateSystemOld(name)
{
  setDescription("Drag around the dot to adjust the mix between the 4 signals");

  // indicate that this ParameterObserver is a GUI element:
  ParameterObserver::setIsGuiElement(true);

  // assign the pointer to the VectorMixer to be edited:
  vectorMixerToEdit = newVectorMixerToEdit;

  // set up the plot range:
  setMaximumRange(-1.1, 1.1, -1.1, 1.1);
  setCurrentRange(-1.1, 1.1, -1.1, 1.1);
  setHorizontalCoarseGrid(1.0,  true);
  //setHorizontalFineGrid(  0.25, true);
  setVerticalCoarseGrid(  1.0,  true);
  //setVerticalFineGrid(    0.25, true);

  // use a normal mouse cursor:
  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  // choose some nice radius for the handle:
  dotRadius = 5.f;
  dotColour = Colours::blue;

  // initialize the two pointers to Parameter objects with NULL:
  xParameter  = NULL;
  yParameter  = NULL;

  // activate automation for this ParameterObserver:
  ParameterObserver::setLocalAutomationSwitch(true);
}

VectorMixerPad::~VectorMixerPad()
{
  // remove ourselves as listeners from the Parameter objects, such that they do 
  // not try to notify a nonexistent listener:
  ParameterObserver::setLocalAutomationSwitch(false);
  if( xParameter != NULL )
    xParameter->deRegisterParameterObserver(this);
  if( yParameter != NULL )
    yParameter->deRegisterParameterObserver(this);
}

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void VectorMixerPad::setVectorMixerToEdit(rosic::VectorMixer* newVectorMixerToEdit)
{
  vectorMixerToEdit = newVectorMixerToEdit;
}

void VectorMixerPad::assignParameterX(Parameter *parameterToAssign)
{
  xParameter = parameterToAssign;
  if( xParameter != NULL )
    xParameter->registerParameterObserver(this);
}

void VectorMixerPad::assignParameterY(Parameter *parameterToAssign)
{
  yParameter = parameterToAssign;
  if( yParameter != NULL )
    yParameter->registerParameterObserver(this);
}

void VectorMixerPad::unAssignParameterX()
{
  if( xParameter != NULL )
    xParameter->deRegisterParameterObserver(this);
  xParameter = NULL;
}

void VectorMixerPad::unAssignParameterY()
{
  if( yParameter != NULL )
    yParameter->deRegisterParameterObserver(this);
  yParameter = NULL;
}

void VectorMixerPad::parameterChanged(Parameter* parameterThatHasChanged)
{
  // send out a change-message: 
  sendChangeMessage();
  //updatePlot();
}

void VectorMixerPad::updateWidgetFromAssignedParameter(bool sendMessage)
{
  updateBackgroundImage();
  if( sendMessage == true )
    sendChangeMessage();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void VectorMixerPad::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  setLocalAutomationSwitch(false);

  // call the method which updates the widget:
  updateBackgroundImage();
  //updateWidgetFromAssignedParameter(false);


  // switch the wantsAutomationNotification flag on again:  
  setLocalAutomationSwitch(true);
}

void VectorMixerPad::mouseDown(const MouseEvent &e)
{
  if( vectorMixerToEdit == NULL )
    return;

  // preliminray: do not open the MIDI-learn menu on right-button - show the image export menu 
  // instead (inherited behaviour from CoordinateSytem):
  if( e.mods.isRightButtonDown() )
    CoordinateSystemOld::mouseDown(e);
  else
  {
    // get the position of the event in components coordinates
    double x = e.getMouseDownX();
    double y = e.getMouseDownY();

    setupVectorMixerAccordingToMousePosition(x, y);

    // send out a change-message:
    sendChangeMessage();
  }

  /*
  if( e.mods.isRightButtonDown() && xParameter != NULL && yParameter != NULL )
  {
  // prepare some strings for the popup menu:
  int freqCC = xParameter->getAssignedMidiController();
  juce::String freqjuce::String;
  if( freqCC > -1 )
  freqjuce::String = juce::String(T("(currently CC#")) + juce::String(freqCC) + juce::String(T(")"));
  else
  freqjuce::String = juce::String(T("(currently none)")); 
  juce::String minFreqjuce::String = hertzToStringWithUnitTotal5(xParameter->getLowerAutomationLimit());
  juce::String maxFreqjuce::String = hertzToStringWithUnitTotal5(xParameter->getUpperAutomationLimit());

  int resoCC = yParameter->getAssignedMidiController();
  juce::String resojuce::String;
  if( resoCC > -1 )
  resojuce::String = juce::String(T("(currently CC#")) + juce::String(resoCC) + juce::String(T(")"));
  else
  resojuce::String = juce::String(T("(currently none)")); 

  // ToDo: different cases - y can be reso, q or gain
  juce::String minResojuce::String = percentToStringWithUnit1(yParameter->getLowerAutomationLimit());
  juce::String maxResojuce::String = percentToStringWithUnit1(yParameter->getUpperAutomationLimit());

  // create a context menu to allow for MIDI learn and setting up min and max automation values:
  PopupMenu menu;
  menu.addItem(1, juce::String(T("MIDI learn frequency ") + freqjuce::String)  );
  menu.addItem(2, juce::String(T("MIDI learn resonance ") + resojuce::String)  );
  menu.addItem(3, juce::String(T("use current values as lower limits"))  );
  menu.addItem(4, juce::String(T("use current values as upper limits"))  );
  menu.addItem(5, juce::String(T("revert to defaults"))                  );

  const int result = menu.show();

  // retrieve current characteristic frequency and resonance:
  double freq = filterToEdit->getFrequencyNominal(); // frequency
  double reso = filterToEdit->getResonance();        // resonance

  if (result == 0)
  {
  // user dismissed the menu without picking anything - do nothing
  }
  else if (result == 1)
  {
  // user picked the frequency learn item:
  xParameter->switchIntoMidiLearnMode();
  }
  else if (result == 2)
  {
  // user picked the resonance learn item:
  yParameter->switchIntoMidiLearnMode();
  }
  else if (result == 3)
  {
  // user picked the lower-limit item:
  xParameter->setLowerAutomationLimit(freq);
  yParameter->setLowerAutomationLimit(reso);
  }
  else if (result == 4)
  {
  // user picked the upper-limit item:
  xParameter->setUpperAutomationLimit(freq);
  yParameter->setUpperAutomationLimit(reso);
  }
  else if (result == 5)
  {
  // user picked the revert to defaults item:
  xParameter->revertToDefaults();
  yParameter->revertToDefaults();
  }
  } // end of  if( e.mods.isRightButtonDown() )
  */
}

void VectorMixerPad::mouseDrag(const juce::MouseEvent &e)
{
  if( vectorMixerToEdit == NULL )
    return;

  /*
  if( e.mods.isRightButtonDown() && xParameter != NULL && yParameter != NULL )
  {
  // ignore mouse drags whne the right button is down and we have assigned automatable 
  // parameters because in that case, the right click is used for opening the MIDI-learn popup
  }
  else...
  */

  // get the position of the event in components coordinates
  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  double y = e.getMouseDownY() + e.getDistanceFromDragStartY();

  setupVectorMixerAccordingToMousePosition(x, y);

  // send out a change-message:
  sendChangeMessage();
}

void VectorMixerPad::setupVectorMixerAccordingToMousePosition(double mouseX, double mouseY)
{
  if( vectorMixerToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  double x = mouseX;
  double y = mouseY;

  // convert them into the vectorMixer's coordinates:
  transformFromComponentsCoordinates(x, y);

  // set up the VectorMixer and raise automation events to update other widgets that represent the
  // parameters:
  vectorMixerToEdit->setXY(x, y);
  if( xParameter != NULL )
    xParameter->setValue(x, true, true);
  if( yParameter != NULL )
    yParameter->setValue(y, true, true);
}

//-------------------------------------------------------------------------------------------------
// drawing:

void VectorMixerPad::resized()
{
  CoordinateSystemOld::resized();
  updateBackgroundImage();
}

void VectorMixerPad::drawCoordinateSystem(Graphics &g)
{
  CoordinateSystemOld::drawCoordinateSystem(g);

  if( vectorMixerToEdit == NULL )
    return;

  // draw the dot-handle and a crosshair:
  //double x = vectorMixerToEdit->getX();
  //double y = vectorMixerToEdit->getY();
  //transformToComponentsCoordinates(x, y);

  double x = coordinateMapper.mapX(vectorMixerToEdit->getX());
  double y = coordinateMapper.mapY(vectorMixerToEdit->getY());

  g.setColour(dotColour);
  g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius), 
    (float) (2*dotRadius), (float) (2*dotRadius) );
  g.setColour(dotColour.withAlpha(0.4f));

  // test this:
  float w = (float) getWidth();
  float h = (float) getHeight();
  g.drawLine(       0,(float)y,        w, (float)y, 1.f);  // horizontal
  g.drawLine((float)x,       0, (float)x,        h, 1.f);  // vertical
}

//=================================================================================================

// construction/destruction:

VectorMixerModuleEditor::VectorMixerModuleEditor(CriticalSection *newPlugInLock, 
  VectorMixerAudioModule* newVectorMixerAudioModule) 
  : AudioModuleEditor(newVectorMixerAudioModule) //, VectorMixerPad(NULL)
{
  // assign the pointer to the rosic::VectorMixer object to be used as audio engine:
  jassert( newVectorMixerAudioModule != NULL ); // you must pass a valid module here
  //VectorMixerPad::setVectorMixerToEdit(newVectorMixerAudioModule->wrappedVectorMixer);

  vectorPad = new VectorMixerPad(newVectorMixerAudioModule->wrappedVectorMixer);
  vectorPad->setDescriptionField(infoField);
  vectorPad->addChangeListener(this);
  addPlot( vectorPad );

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// overrides:

/*
void VectorMixerModuleEditor::updateWidgetsAccordingToState()
{
VectorMixerPad::updateWidgetsAccordingToState();
}
*/

/*
void VectorMixerModuleEditor::mouseDown(const MouseEvent &e)
{
VectorMixerPad::mouseDown(e);
}

void VectorMixerModuleEditor::mouseDrag(const MouseEvent &e)
{
VectorMixerPad::mouseDrag(e);
}

void VectorMixerModuleEditor::mouseMove(const MouseEvent &e)
{
VectorMixerPad::mouseMove(e);
}

void VectorMixerModuleEditor::mouseEnter(const MouseEvent &e)
{
VectorMixerPad::mouseEnter(e);
}

void VectorMixerModuleEditor::mouseExit(const MouseEvent &e)
{
VectorMixerPad::mouseExit(e);
}

void VectorMixerModuleEditor::mouseUp(const MouseEvent &e)
{
VectorMixerPad::mouseUp(e);
}

void VectorMixerModuleEditor::mouseDoubleClick(const MouseEvent &e)
{
VectorMixerPad::mouseDoubleClick(e);
}

void VectorMixerModuleEditor::mouseWheelMove(const MouseEvent &e, 
float wheelIncrementX, float wheelIncrementY)
{
VectorMixerPad::mouseWheelMove(e, wheelIncrementX, wheelIncrementY);
}

void VectorMixerModuleEditor::paint(Graphics &g)
{
VectorMixerPad::paint(g);
}
*/
void VectorMixerModuleEditor::resized()
{
  //vectorPad->setBounds(getBounds());
  vectorPad->setBounds(0, 0, getWidth(), getHeight());
}

/*
XmlElement* VectorMixerModuleEditor::getStateAsXml(const juce::String& stateName) const
{
return AudioModuleEditor::getStateAsXml(stateName);
}

bool VectorMixerModuleEditor::setStateFromXml(const XmlElement& xmlState, const juce::String& name)
{
return AudioModuleEditor::setStateFromXml(xmlState, name);
}
*/