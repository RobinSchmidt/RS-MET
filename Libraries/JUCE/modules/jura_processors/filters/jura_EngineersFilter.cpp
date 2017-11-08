//#include "rosof_EngineersFilterAudioModule.h"
//using namespace rosof;

// construction/destruction:

EngineersFilterAudioModule::EngineersFilterAudioModule(CriticalSection *newPlugInLock, 
  rosic::rsEngineersFilter *sciFilterToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(sciFilterToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedEngineersFilter = sciFilterToWrap;

  moduleName = juce::String("EngineersFilter");
  setActiveDirectory(getApplicationDirectory() + juce::String("/Presets/EngineersFilter") );
  initializeAutomatableParameters();
}

EngineersFilterAudioModule::EngineersFilterAudioModule(CriticalSection *newPlugInLock)
  : AudioModule(newPlugInLock)
{
  wrappedEngineersFilter = new rosic::rsEngineersFilter;
  wrappedEngineersFilterIsOwned = true;

  // todo: factor out this code (duplicated from the other constuctor) into an init() function that 
  // can be called from both constructors:
  moduleName = juce::String("EngineersFilter");
  setActiveDirectory(getApplicationDirectory() + juce::String("/Presets/EngineersFilter") );
  initializeAutomatableParameters();
}

EngineersFilterAudioModule::~EngineersFilterAudioModule()
{
  if(wrappedEngineersFilterIsOwned)
    delete wrappedEngineersFilter;
}

AudioModuleEditor* EngineersFilterAudioModule::createEditor()
{
  return new jura::EngineersFilterModuleEditor(lock, this); // get rid of passing the lock
}

// automation:

void EngineersFilterAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  if( wrappedEngineersFilter == NULL )
    return;

  // find out the index in the vector of the parameter that has been changed:
  int    index = getIndexOfParameter(parameterThatHasChanged);
  double value = parameterThatHasChanged->getValue();
  switch( index )
  {
  case 0: wrappedEngineersFilter->setMode(               roundToInt(value)     );   break;
  case 1: wrappedEngineersFilter->setApproximationMethod(roundToInt(value) + 1 );   break;
  case 2: wrappedEngineersFilter->setFrequency(          value);                    break;
  case 3: wrappedEngineersFilter->setPrototypeOrder(     roundToInt(value)     );   break;
  case 4: wrappedEngineersFilter->setBandwidth(          value);                    break;
  case 5: wrappedEngineersFilter->setGain(               value);                    break;
  case 6: wrappedEngineersFilter->setRipple(             value);                    break;
  case 7: wrappedEngineersFilter->setStopbandRejection(  value);                    break;
  /*

  case 3: wrappedEngineersFilter->setUpperFrequency(     value);                    break;

  */

  } // end of switch( parameterIndex )

  markStateAsDirty();
}

// internal functions:

void EngineersFilterAudioModule::initializeAutomatableParameters()
{
  juce::Array<double> defaultValues;
  AutomatableParameter* p;

  // #000:	
  p = new AutomatableParameter(lock, "Mode", 0.0, 7.0, 1.0, 1.0, Parameter::STRING);
  p->addStringValue(juce::String("Bypass"));
  p->addStringValue(juce::String("Lowpass"));
  p->addStringValue(juce::String("Highpass"));
  p->addStringValue(juce::String("Bandpass"));
  p->addStringValue(juce::String("Bandreject"));
  p->addStringValue(juce::String("Low Shelf"));
  p->addStringValue(juce::String("High Shelf"));
  p->addStringValue(juce::String("Peak/Dip"));
  addObservedParameter(p);

  // #001:	
  p = new AutomatableParameter(lock, "Method", 0.0, 5.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue(juce::String("Butterworth"));
  p->addStringValue(juce::String("Chebychev"));
  p->addStringValue(juce::String("Inverse Chebychev"));
  p->addStringValue(juce::String("Elliptic"));
  p->addStringValue(juce::String("Bessel"));
  p->addStringValue(juce::String("Papoulis"));
  //p->addStringValue(juce::String(T("Gauss")));
  //p->addStringValue(juce::String(T("Halpern")));
  //p->addStringValue(juce::String(T("Linkwitz/Riley")));
  //p->addStringValue(juce::String(T("Coincident Poles")));
  addObservedParameter(p);
  // maybe sort the differently: Bessel->Butter->Papoulis->Cheby->Ellip
  // or better: Coincident, Gauss, Bessel, Linkwitz, Butter, Papoulis, Halpern, Inv Cheby, Cheby
  // elliptic

  // #002:	
  p = new AutomatableParameter(lock, "Frequency" , 20.0, 20000.0, 0.0, 1000.0, 
    Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #003:	
  //p = new Parameter("Slope", 6.0, 120.0, 6.0, 24.0, Parameter::LINEAR);
  //addObservedParameter(p);
  p = new AutomatableParameter(lock, "Order", 1.0, 20.0, 1.0, 4.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #004:	
  p = new AutomatableParameter(lock, "Bandwidth", 0.1, 10.0, 0.01, 1.0, 
    Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #005:	
  p = new AutomatableParameter(lock, "Gain", -48, 24.0, 0.01, 0.0, Parameter::LINEAR);
  //p = new AutomatableParameter(plugInLock, "Gain", -108, 24.0, 0.01, 0.0, Parameter::LINEAR);
  addObservedParameter(p);

  // #006:	
  p = new AutomatableParameter(lock, "Ripple", 0.1, 12.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);

  // #007:	
  p = new AutomatableParameter(lock, "Rejection", 20, 120.0, 0.01, 60.0, Parameter::LINEAR);
  addObservedParameter(p);

  // make a call to parameterChanged for each parameter in order to set up the DSP-core to reflect 
  // the values the automatable parameters:
  for(int i=0; i < (int) parameters.size(); i++ )
    parameterChanged(parameters[i]);
}

//=================================================================================================

EngineersFilterPlotEditor::EngineersFilterPlotEditor(const juce::String& name) 
  : SpectrumDisplayOld(name)
{
  setDescription("Drag vertical lines at the corner frequencies to adjust them.");

  //ParameterObserver::isGuiElement = true;
  sciFilterToEdit = NULL;

  // set up the plot range:
  setMaximumRange(15.625, 32000.0, -120.0, 30.0);
  setCurrentRange(15.625, 32000.0, -120.0, 30.0);
  setHorizontalCoarseGrid(12.0, true);
  setHorizontalFineGrid(   3.0, false);
  setVerticalCoarseGridVisible( true);
  setVerticalFineGridVisible(   false);
  CoordinateSystemOld::setAxisValuesPositionX(CoordinateSystemOld::ABOVE_AXIS);
  CoordinateSystemOld::setAxisValuesPositionY(CoordinateSystemOld::RIGHT_TO_AXIS);
  showPositionAsDescription = true;

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  /*
  // initialize the two pointers to Parameter objects with NULL:
  lowFreqParameter   = NULL;
  lowSlopeParameter  = NULL;
  highFreqParameter  = NULL;
  highSlopeParameter = NULL;
  */

  // this stuff will be (re-) assigned in resized():
  numBins     = 0;
  frequencies = NULL;
  magnitudes  = NULL;

  currentlyDraggedHandle = NONE;

  // activate automation for this ParameterObserver:
  //ParameterObserver::localAutomationSwitch = true;
}

EngineersFilterPlotEditor::~EngineersFilterPlotEditor(void)
{
  /*
  // remove ourselves as listeners from the Parameter objects, such that they do not try to 
  // notify a nonexistent listener:
  ParameterObserver::localAutomationSwitch = false;
  if( lowFreqParameter != NULL )
  lowFreqParameter->deRegisterParameterObserver(this);
  if( lowSlopeParameter != NULL )
  lowSlopeParameter->deRegisterParameterObserver(this);
  if( highFreqParameter != NULL )
  highFreqParameter->deRegisterParameterObserver(this);
  if( highSlopeParameter != NULL )
  highSlopeParameter->deRegisterParameterObserver(this);
  */

  deleteAndZero(frequencies);
  deleteAndZero(magnitudes);
}

// parameter-settings:

void EngineersFilterPlotEditor::setEngineersFilterToEdit(
  rosic::rsEngineersFilter* newEngineersFilterToEdit)
{
  sciFilterToEdit = newEngineersFilterToEdit;
}

/*
void EngineersFilterPlotEditor::assignParameterLowFreq(Parameter *parameterToAssign)
{
lowFreqParameter = parameterToAssign;
if( lowFreqParameter != NULL )
lowFreqParameter->registerParameterObserver(this);
}

void EngineersFilterPlotEditor::assignParameterLowSlope(Parameter *parameterToAssign)
{
lowSlopeParameter = parameterToAssign;
if( lowSlopeParameter != NULL )
lowSlopeParameter->registerParameterObserver(this);
}

void EngineersFilterPlotEditor::assignParameterHighFreq(Parameter *parameterToAssign)
{
highFreqParameter = parameterToAssign;
if( highFreqParameter != NULL )
highFreqParameter->registerParameterObserver(this);
}

void EngineersFilterPlotEditor::assignParameterHighSlope(Parameter *parameterToAssign)
{
highSlopeParameter = parameterToAssign;
if( highSlopeParameter != NULL )
highSlopeParameter->registerParameterObserver(this);
}

void EngineersFilterPlotEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
sendChangeMessage();
}
*/

void EngineersFilterPlotEditor::updateWidgetFromAssignedParameter(bool sendMessage)
{
  updatePlot();
  if( sendMessage == true )
    sendChangeMessage();
}

// callbacks:

void EngineersFilterPlotEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  //localAutomationSwitch = false;

  // call the method which updates the widget:
  updatePlot();
  //updateWidgetFromAssignedParameter(false);

  // switch the wantsAutomationNotification flag on again:  
  //localAutomationSwitch = true;
}

void EngineersFilterPlotEditor::mouseMove(const MouseEvent &e)
{
  SpectrumDisplayOld::mouseMove(e);

  // ...we need some stuff from EasyQ here - maybe factor out a basclass 
  // FrequencyResponsePlotEditor

  /*
  if( getDragHandleAt(e.x, e.y) != NONE )
  setMouseCursor(MouseCursor::LeftRightResizeCursor);
  else
  setMouseCursor(MouseCursor::NormalCursor);
  */
}

void EngineersFilterPlotEditor::mouseDown(const MouseEvent &e)
{
  if( sciFilterToEdit == NULL )
    return;

  currentlyDraggedHandle = getDragHandleAt(e.x, e.y);
}

void EngineersFilterPlotEditor::mouseDrag(const juce::MouseEvent &e)
{
  if( sciFilterToEdit == NULL )
    return;

  // get the position of the event in components coordinates:
  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  double y = e.getMouseDownY() + e.getDistanceFromDragStartY();

  setupFilterAccordingToMousePosition(x, y);
  sendChangeMessage();
}

void EngineersFilterPlotEditor::mouseUp(const juce::MouseEvent &e)
{
  currentlyDraggedHandle = NONE;
}

int EngineersFilterPlotEditor::getDragHandleAt(int x, int y)
{
  return NONE;

  /*
  if( sciFilterToEdit == NULL )
  return NONE;

  // get the x-positions of the two vertical lines (in pixels):
  double xLowMid  = sciFilterToEdit->getLowMidFreq();
  double yDummy   = 0.0;
  transformToComponentsCoordinates(xLowMid, yDummy);
  double xMidHigh = sciFilterToEdit->getMidHighFreq();
  yDummy = 0.0;
  transformToComponentsCoordinates(xMidHigh, yDummy);

  // ckeck, if the passed x,y coordinates are near one of these vertical lines:
  double margin = 4.0;
  if( fabs( (double)x - xLowMid ) < margin )
  return LOW_MID;
  else if( fabs( (double)x - xMidHigh ) < margin && !sciFilterToEdit->isInTwoWayMode() )
  return MID_HIGH;
  else
  return NONE;
  */
}

void EngineersFilterPlotEditor::setupFilterAccordingToMousePosition(double mouseX, double mouseY)
{
  double freq  = mouseX;
  double dummy = mouseY;
  transformFromComponentsCoordinates(freq, dummy);

  /*
  if( currentlyDraggedHandle == LOW_MID && lowFreqParameter != NULL )
  lowFreqParameter->setValue(freq, true, true);
  if( currentlyDraggedHandle == MID_HIGH && highFreqParameter != NULL  )
  highFreqParameter->setValue(freq, true, true);
  */
}

// drawing:

void EngineersFilterPlotEditor::resized()
{
  SpectrumDisplayOld::resized();

  // (re) allocate and fill the arrays for the magnitude plot
  numBins = getWidth();
  if( frequencies == NULL )
    delete[] frequencies;
  if( magnitudes == NULL )
    delete[] magnitudes;
  frequencies = new double[numBins];
  magnitudes  = new double[numBins];

  getDisplayedFrequencies(frequencies, numBins);
  if( sciFilterToEdit != NULL )
    sciFilterToEdit->getMagnitudeResponse(frequencies, magnitudes, numBins, true, false);
  else
  {
    int k;
    for(k=0; k<numBins; k++)
    {
      frequencies[k] = 15.625;
      magnitudes[k]  = 0.0;
    }
  }
  setSpectrum(numBins, frequencies, magnitudes);
}

void EngineersFilterPlotEditor::updatePlot()
{
  if( sciFilterToEdit == NULL )
    return;
  sciFilterToEdit->getMagnitudeResponse(frequencies, magnitudes, numBins, true, false);
  setSpectrum(numBins, frequencies, magnitudes);
}

void EngineersFilterPlotEditor::plotCurveFamily(Graphics &g, juce::Image* targetImage, 
  XmlElement *targetSVG)
{
  if( sciFilterToEdit == NULL )
    return;

  CurveFamilyPlotOld::plotCurveFamily(g, targetImage, targetSVG);

  /*
  if( colourScheme.plotColours.size() < 3 )
  return;

  // draw the draggable vertical lines
  //g.setColour(Colour(0xffa0a0a0)); // a somewhat bright grey
  g.setColour(Colour(0xff707070));
  double x = sciFilterToEdit->getLowMidFreq();
  double y = 0.0;
  if( targetImage == NULL )
  transformToComponentsCoordinates(x, y);
  else
  transformToImageCoordinates(x, y, targetImage);
  //g.setColour(rojue::getMixedColour(colourScheme.plotColours[0], colourScheme.plotColours[1]));
  g.drawLine((float) x, 0.f, (float) x, (float) getHeight(), 2.0);

  if( !sciFilterToEdit->isInTwoWayMode() )
  {
  x = sciFilterToEdit->getMidHighFreq();
  y = 0.0;
  if( targetImage == NULL )
  transformToComponentsCoordinates(x, y);
  else
  transformToImageCoordinates(x, y, targetImage);
  //g.setColour(rojue::getMixedColour(colourScheme.plotColours[1], colourScheme.plotColours[2]));
  g.drawLine((float) x, 0.f, (float) x, (float) getHeight(), 2.0);
  }
  */
}

//=================================================================================================

EngineersFilterModuleEditor::EngineersFilterModuleEditor(CriticalSection *newPlugInLock, 
  EngineersFilterAudioModule* newEngineersFilterAudioModule) 
  : AudioModuleEditor(newEngineersFilterAudioModule)
{
  setHeadlineText( juce::String("EngineersFilter") );

  isTopLevelEditor = true;

  jassert(newEngineersFilterAudioModule != NULL ); // you must pass a valid module here
  sciFilterModuleToEdit = newEngineersFilterAudioModule;

  addWidget( modeComboBox = new RNamedComboBox("modeComboBox", "Mode:") );
  modeComboBox->assignParameter( sciFilterModuleToEdit->getParameterByName("Mode") );
  modeComboBox->setDescription("Mode or type of the filter");
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update visibility of the sliders

  addWidget( methodComboBox = new RNamedComboBox("methodComboBox", "Method:") );
  methodComboBox->assignParameter( sciFilterModuleToEdit->getParameterByName("Method") );
  methodComboBox->setDescription("Approximation method for the filter design");
  methodComboBox->setDescriptionField(infoField);
  methodComboBox->registerComboBoxObserver(this); // to update visibility of the sliders

  modeComboBox->setNameLabelWidth(methodComboBox->getNameLabelWidth()); // to align the actual boxes

  addWidget( frequencySlider = new RSlider ("FrequencySlider") );
  frequencySlider->assignParameter( sciFilterModuleToEdit->getParameterByName("Frequency") );
  frequencySlider->setSliderName(juce::String("Frequency"));
  frequencySlider->setDescription(juce::String("Characteristic frequency"));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);
  frequencySlider->addListener(this);

  addWidget( orderSlider = new RSlider ("OrderSlider") );
  orderSlider->assignParameter( sciFilterModuleToEdit->getParameterByName("Order") );
  orderSlider->setSliderName(juce::String("Order"));
  //orderSlider->setDescription(juce::String(T("Order of the prototype filter (bandpass-/reject and peak filters have actual order of twice this value)")));
  orderSlider->setDescription(juce::String("Order of the prototype filter"));
  orderSlider->setDescriptionField(infoField);
  orderSlider->setStringConversionFunction(&valueToString);
  orderSlider->addListener(this);

  addWidget( bandwidthSlider = new RSlider ("BandwidthSlider") );
  bandwidthSlider->assignParameter( sciFilterModuleToEdit->getParameterByName("Bandwidth") );
  bandwidthSlider->setSliderName(juce::String("Bandwidth"));
  bandwidthSlider->setDescription(juce::String("Bandwidth for bandpass-/bandreject-/peak-filters (in octaves)"));
  bandwidthSlider->setDescriptionField(infoField);
  bandwidthSlider->setStringConversionFunction(&octavesToStringWithUnit2);
  bandwidthSlider->addListener(this);

  addWidget( gainSlider = new RSlider("GainSlider") );
  gainSlider->assignParameter( sciFilterModuleToEdit->getParameterByName("Gain") );
  gainSlider->setSliderName(juce::String("Gain"));
  gainSlider->setDescription(juce::String("Gain for shelving- and peak-filters"));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  gainSlider->addListener(this);

  addWidget( rippleSlider = new RSlider("RippleSlider") );
  rippleSlider->assignParameter( sciFilterModuleToEdit->getParameterByName("Ripple") );
  rippleSlider->setSliderName(juce::String("Ripple"));
  rippleSlider->setDescription(juce::String("Ripple insider the filter's passband"));
  rippleSlider->setDescriptionField(infoField);
  rippleSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  rippleSlider->addListener(this);

  addWidget( rejectionSlider = new RSlider("RejectionSlider") );
  rejectionSlider->assignParameter( sciFilterModuleToEdit->getParameterByName("Rejection") );
  rejectionSlider->setSliderName(juce::String("Rejection"));
  rejectionSlider->setDescription(juce::String("Rejection level for the filter's stopband"));
  rejectionSlider->setDescriptionField(infoField);
  rejectionSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  rejectionSlider->addListener(this);

  plotEditor = new EngineersFilterPlotEditor(juce::String("SpectrumEditor"));
  plotEditor->setDescription(juce::String("Frequency response plot"));
  plotEditor->setDescriptionField(infoField);
  plotEditor->addChangeListener(this);
  plotEditor->setEngineersFilterToEdit(newEngineersFilterAudioModule->wrappedEngineersFilter);
  addPlot( plotEditor );
  juce::Array<Colour> graphColours;
  //graphColours.add(Colour(0xffa00000));   
  //graphColours.add(Colour(0xff30b030));
  //graphColours.add(Colour(0xff6060ff));
  setGraphColours(graphColours);

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(640, 300);  
}

// callbacks:

void EngineersFilterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  updateWidgetVisibility();
  plotEditor->updatePlot();
}

void EngineersFilterModuleEditor::rSliderValueChanged(RSlider *rSliderThatHasChanged)
{
  plotEditor->updatePlot();
}

void EngineersFilterModuleEditor::updateWidgetsAccordingToState()
{
  AudioModuleEditor::updateWidgetsAccordingToState();
  updateWidgetVisibility();
  plotEditor->updatePlot();
}

void EngineersFilterModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  modeComboBox->setBounds(   x+4, y+4,  w/3-8, 16);
  methodComboBox->setBounds( x+4, y+24, w/3-8, 16);

  x = modeComboBox->getRight() + 4;

  w = getWidth() - x; // remaining width
  w = w/3;

  frequencySlider->setBounds(x+4, y+4,  w-8, 16);
  orderSlider->setBounds(    x+4, y+24, w-8, 16);
  x += w;
  bandwidthSlider->setBounds(x+4, y+4,  w-8, 16);
  gainSlider->setBounds(     x+4, y+24, w-8, 16);
  x += w;
  rippleSlider->setBounds(   x+4, y+4,  w-8, 16);
  rejectionSlider->setBounds(x+4, y+24, w-8, 16);

  y = methodComboBox->getBottom();
  h = infoField->getY()-y;

  if( isTopLevelEditor )
    plotEditor->setBounds(4, y+4, getWidth()-8, h-4);
  else
    plotEditor->setBounds(4, y+4, getWidth()-8, h-8);
}

void EngineersFilterModuleEditor::updateWidgetVisibility()
{
  if( sciFilterModuleToEdit == NULL )
    return;
  if( sciFilterModuleToEdit->wrappedEngineersFilter == NULL )
    return;

  rosic::rsEngineersFilter* sf = sciFilterModuleToEdit->wrappedEngineersFilter;
  bandwidthSlider->setEnabled(sf->hasCurrentModeBandwidthParameter());
  gainSlider->setEnabled(     sf->hasCurrentModeGainParameter());
  rippleSlider->setEnabled(   sf->hasCurrentModeRippleParameter());
  rejectionSlider->setEnabled(sf->hasCurrentModeRejectionParameter());
  if( sf->hasCurrentModeGainParameter() )
  {
    if( sf->getApproximationMethod() == PrototypeDesigner::CHEBYCHEV )
      rippleSlider->setDescription(juce::String("Ripple inside the boost/cut band in percent of dB-peak-gain"));
    if( sf->getApproximationMethod() == PrototypeDesigner::INVERSE_CHEBYCHEV )
      rippleSlider->setDescription(juce::String("Ripple outside the boost/cut band in percent of dB-peak-gain"));
    if( sf->getApproximationMethod() == PrototypeDesigner::ELLIPTIC )
      rippleSlider->setDescription(juce::String("Ripple in- and outside the boost/cut band in percent of dB-peak-gain"));
    rippleSlider->setStringConversionFunction(&percentToStringWithUnit2);
  }
  else
  {
    rippleSlider->setDescription(juce::String("Ripple inside the filter's passband in dB"));
    rippleSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  }
}
