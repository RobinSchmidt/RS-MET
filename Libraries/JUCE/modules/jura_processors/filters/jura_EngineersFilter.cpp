// construction/destruction:

EngineersFilterAudioModule::EngineersFilterAudioModule(CriticalSection *newPlugInLock, 
  rosic::rsEngineersFilterStereo *sciFilterToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(sciFilterToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedEngineersFilter = sciFilterToWrap;
  init();
}

EngineersFilterAudioModule::EngineersFilterAudioModule(CriticalSection *newPlugInLock)
  : AudioModule(newPlugInLock)
{
  wrappedEngineersFilter = new rosic::rsEngineersFilterStereo;
  wrappedEngineersFilterIsOwned = true;
  init();
}

void EngineersFilterAudioModule::init()
{
  setModuleTypeName("EngineersFilter");
  createParameters();
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

// internal functions:

void EngineersFilterAudioModule::createParameters()
{
  typedef MetaControlledParameter Param;
  Param* p;

  typedef rosic::rsEngineersFilterStereo EF;
  EF* ef = wrappedEngineersFilter;

  //juce::Array<double> defaultValues;

  p = new Param("Mode", 0.0, 7.0, 1.0, Parameter::STRING, 1.0);
  p->setValueChangeCallback<EF>(ef, &EF::setMode);
  p->addStringValue("Bypass");
  p->addStringValue("Lowpass");
  p->addStringValue("Highpass");
  p->addStringValue("Bandpass");
  p->addStringValue("Bandreject");
  p->addStringValue("Low Shelf");
  p->addStringValue("High Shelf");
  p->addStringValue("Peak/Dip");
  addObservedParameter(p);

  p = new Param("Method", 0.0, 5.0, 0.0, Parameter::STRING, 1.0);
  p->setValueChangeCallback<EF>(ef, &EF::setApproximationMethod);
  p->addStringValue("Butterworth");
  p->addStringValue("Chebychev");
  p->addStringValue("Inverse Chebychev");
  p->addStringValue("Elliptic");
  p->addStringValue("Bessel");
  p->addStringValue("Papoulis");
  p->addStringValue("Halpern");
  //p->addStringValue("Gauss");
  //p->addStringValue("Linkwitz/Riley");
  //p->addStringValue("Coincident Poles");
  addObservedParameter(p);
  // maybe sort the differently: Bessel->Butter->Papoulis->Cheby->Ellip
  // or better: Coincident, Gauss, Bessel, Linkwitz, Butter, Papoulis, Halpern, Inv Cheby, Cheby
  // elliptic

  p = new Param("Frequency" , 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<EF>(ef, &EF::setFrequency);
  addObservedParameter(p);

  p = new Param("Order", 1.0, 20.0, 4.0, Parameter::LINEAR, 1.0);
  p->setValueChangeCallback<EF>(ef, &EF::setPrototypeOrder);
  addObservedParameter(p);

  p = new Param("Bandwidth", 0.1, 10.0, 1.0, Parameter::EXPONENTIAL, 0.01);
  p->setValueChangeCallback<EF>(ef, &EF::setBandwidth);
  addObservedParameter(p);

  p = new Param("Gain", -48, 24.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<EF>(ef, &EF::setGain);
  addObservedParameter(p);
	
  p = new Param("Ripple", 0.1, 12.0, 1.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<EF>(ef, &EF::setRipple);
  addObservedParameter(p);

  p = new Param("Rejection", 20, 120.0, 60.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<EF>(ef, &EF::setStopbandRejection);
  addObservedParameter(p);
}

//=================================================================================================

EngineersFilterPlotEditor::EngineersFilterPlotEditor(const juce::String& name) 
  : rsSpectrumPlot(name)
{
  //setDescription("Drag vertical lines at the corner frequencies to adjust them.");
  setIsGuiElement(true);

  sciFilterToEdit = NULL;

  // set up the plot range:
  setMaximumRange(15.625, 32000.0, -120.0, 30.0);
  setCurrentRange(15.625, 32000.0, -120.0, 30.0);
  setHorizontalCoarseGrid(12.0, true);
  setHorizontalFineGrid(   3.0, false);
  setVerticalCoarseGridVisible( true);
  setVerticalFineGridVisible(   false);
  rsPlot::setAxisValuesPositionX(rsPlotSettings::ABOVE_AXIS);
  rsPlot::setAxisValuesPositionY(rsPlotSettings::RIGHT_TO_AXIS);
  showPositionAsDescription = true;

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

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
  modeParam     ->deRegisterParameterObserver(this);
  methodParam   ->deRegisterParameterObserver(this);
  orderParam    ->deRegisterParameterObserver(this);
  freqParam     ->deRegisterParameterObserver(this);
  bandwidthParam->deRegisterParameterObserver(this);
  gainParam     ->deRegisterParameterObserver(this);
  rippleParam   ->deRegisterParameterObserver(this);
  rejectionParam->deRegisterParameterObserver(this);
  deleteAndZero(frequencies);
  deleteAndZero(magnitudes);
}

// parameter-settings:

void EngineersFilterPlotEditor::setEngineersFilterToEdit(
  rosic::rsEngineersFilterStereo* newEngineersFilterToEdit)
{
  sciFilterToEdit = newEngineersFilterToEdit;
}

void EngineersFilterPlotEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  updatePlot();  // should call repaintOnMessageThread();
}

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
  rsSpectrumPlot::mouseMove(e);

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
  toPixelCoordinates(xLowMid, yDummy);
  double xMidHigh = sciFilterToEdit->getMidHighFreq();
  yDummy = 0.0;
  toPixelCoordinates(xMidHigh, yDummy);

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
  fromPixelCoordinates(freq, dummy);

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
  rsSpectrumPlot::resized();

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

  rsDataPlot::plotCurveFamily(g, targetImage, targetSVG);

  /*
  if( colourScheme.plotColours.size() < 3 )
  return;

  // draw the draggable vertical lines
  //g.setColour(Colour(0xffa0a0a0)); // a somewhat bright grey
  g.setColour(Colour(0xff707070));
  double x = sciFilterToEdit->getLowMidFreq();
  double y = 0.0;
  if( targetImage == NULL )
  toPixelCoordinates(x, y);
  else
  transformToImageCoordinates(x, y, targetImage);
  //g.setColour(rojue::getMixedColour(colourScheme.plotColours[0], colourScheme.plotColours[1]));
  g.drawLine((float) x, 0.f, (float) x, (float) getHeight(), 2.0);

  if( !sciFilterToEdit->isInTwoWayMode() )
  {
  x = sciFilterToEdit->getMidHighFreq();
  y = 0.0;
  if( targetImage == NULL )
  toPixelCoordinates(x, y);
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
  jassert(newEngineersFilterAudioModule != nullptr ); // you must pass a valid module here
  sciFilterModuleToEdit = newEngineersFilterAudioModule;
  createWidgets();
  updateWidgetsAccordingToState();
  //setSize(640, 300);


  //setSize(627, 291);
   // w x h should be (N*11) x (M*25+66) for homogenous visual appearance of the grid. The plot will
   // then be (N*11) x (M*25). We choose N=57, M=9 here.

  // the formula doesn't work anymore after refactoring the drawing code - figure out the new optimal
  // setting..
  // width 606, 617, 628 looks good - so it should should be N*11+1, for N=55, it's 606

  setSize(628, 292);  // yep - looks perfect - figure out general formula

}

EngineersFilterModuleEditor::~EngineersFilterModuleEditor()
{
  sciFilterModuleToEdit->getParameterByName("Mode")->deRegisterParameterObserver(this);
  sciFilterModuleToEdit->getParameterByName("Method")->deRegisterParameterObserver(this);
}

// callbacks:

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

  y = methodComboBox->getBottom()+4;
  h = infoField->getY()-y;

  plotEditor->setBounds(0, y, getWidth(), h);
  //plotEditor->setBounds(0, 0, getWidth(), getHeight()); // for figuring out ideal size
}

void EngineersFilterModuleEditor::parameterChanged(Parameter* param)
{
  updateWidgetVisibility();
}

void EngineersFilterModuleEditor::createWidgets()
{
  typedef AutomatableSlider Sld;
  //typedef AutomatableComboBox Box;
  Sld* s;
  //Box* c;
  Parameter* p;

  plotEditor = new EngineersFilterPlotEditor("SpectrumEditor");
  plotEditor->setDescription("Frequency response plot");
  plotEditor->setDescriptionField(infoField);
  plotEditor->addChangeListener(this);
  plotEditor->setEngineersFilterToEdit(sciFilterModuleToEdit->wrappedEngineersFilter);
  addPlot( plotEditor );
  juce::Array<Colour> graphColours;
  //graphColours.add(Colour(0xffa00000));   
  //graphColours.add(Colour(0xff30b030));
  //graphColours.add(Colour(0xff6060ff));
  setGraphColours(graphColours);

  addWidget( modeComboBox = new RNamedComboBox("modeComboBox", "Mode:") );
  modeComboBox->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Mode") );
  p->registerParameterObserver(this);
  plotEditor->assignParameterMode(p);
  modeComboBox->setDescription("Mode or type of the filter");
  modeComboBox->setDescriptionField(infoField);

  addWidget( methodComboBox = new RNamedComboBox("methodComboBox", "Method:") );
  methodComboBox->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Method") );
  p->registerParameterObserver(this);
  plotEditor->assignParameterMethod(p);
  methodComboBox->setDescription("Approximation method for the filter design");
  methodComboBox->setDescriptionField(infoField);

  modeComboBox->setNameLabelWidth(methodComboBox->getNameLabelWidth()); // to align the actual boxes

  addWidget( frequencySlider = s = new Sld );
  s->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Frequency") );
  plotEditor->assignParameterFrequency(p);
  s->setSliderName("Frequency");
  s->setDescription("Characteristic frequency");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( orderSlider = s = new Sld );
  s->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Order") );
  plotEditor->assignParameterOrder(p);
  s->setSliderName("Order");
  s->setDescription("Order of the prototype filter (peak/band filters have twice this order)");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString);

  addWidget( bandwidthSlider = s = new Sld );
  s->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Bandwidth") );
  plotEditor->assignParameterBandwidth(p);
  s->setSliderName("Bandwidth");
  s->setDescription("Bandwidth for bandpass-/bandreject-/peak-filters (in octaves)");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&octavesToStringWithUnit2);

  addWidget( gainSlider = s = new Sld );
  s->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Gain") );
  plotEditor->assignParameterGain(p);
  s->setSliderName("Gain");
  s->setDescription("Gain for shelving- and peak-filters");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( rippleSlider = s = new Sld );
  s->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Ripple") );
  plotEditor->assignParameterRipple(p);
  s->setSliderName("Ripple");
  s->setDescription("Ripple insider the filter's passband");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( rejectionSlider = s = new Sld );
  s->assignParameter( p = sciFilterModuleToEdit->getParameterByName("Rejection") );
  plotEditor->assignParameterRejection(p);
  s->setSliderName("Rejection");
  s->setDescription("Rejection level for the filter's stopband");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&decibelsToStringWithUnit2);
}

void EngineersFilterModuleEditor::updateWidgetVisibility()
{
  if( sciFilterModuleToEdit == NULL )
    return;
  if( sciFilterModuleToEdit->wrappedEngineersFilter == NULL )
    return;

  rosic::rsEngineersFilterStereo* sf = sciFilterModuleToEdit->wrappedEngineersFilter;
  bandwidthSlider->setEnabled(sf->hasCurrentModeBandwidthParameter());
  gainSlider->setEnabled(     sf->hasCurrentModeGainParameter());
  rippleSlider->setEnabled(   sf->hasCurrentModeRippleParameter());
  rejectionSlider->setEnabled(sf->hasCurrentModeRejectionParameter());
  if( sf->hasCurrentModeGainParameter() )
  {
    if( sf->getApproximationMethod() == rsPrototypeDesigner::CHEBYCHEV )
      rippleSlider->setDescription("Ripple inside the boost/cut band in percent of dB-peak-gain");
    if( sf->getApproximationMethod() == rsPrototypeDesigner::INVERSE_CHEBYCHEV )
      rippleSlider->setDescription("Ripple outside the boost/cut band in percent of dB-peak-gain");
    if( sf->getApproximationMethod() == rsPrototypeDesigner::ELLIPTIC )
      rippleSlider->setDescription("Ripple in- and outside the boost/cut band in percent of dB-peak-gain");
    rippleSlider->setStringConversionFunction(&percentToStringWithUnit2);
  }
  else
  {
    rippleSlider->setDescription("Ripple inside the filter's passband in dB");
    rippleSlider->setStringConversionFunction(&decibelsToStringWithUnit2);
  }
}
