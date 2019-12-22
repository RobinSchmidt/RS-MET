
//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiModeFilterAudioModule::MultiModeFilterAudioModule(CriticalSection *newPlugInLock, 
  rosic::MultiModeFilter *newMultiModeFilterToWrap) : AudioModule(newPlugInLock)
{
  jassert( newMultiModeFilterToWrap != NULL ); // you must pass a valid rosic-object
  wrappedMultiModeFilter = newMultiModeFilterToWrap;
  setModuleTypeName("MultiModeFilter");
  createParameters();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MultiModeFilterAudioModule::createParameters()
{
  typedef MetaControlledParameter Param;
  Param* p;

  typedef rosic::MultiModeFilter MMF;
  MMF* mmf = wrappedMultiModeFilter;

  std::vector<double> defaultValues;

  p = new Param("TwoStages", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0);
  p->setValueChangeCallback<MMF>(mmf, &MMF::useTwoStages);
  addObservedParameter(p);

  p = new Param("Mode", 0.0, 14.0, 0.0, Parameter::STRING, 1.0);
  //p = new Param("Mode", 0.0, 14.0, 1.0, Parameter::STRING, 1.0);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setMode);
  p->addStringValue("Bypass");
  p->addStringValue("Moogish Lowpass");
  p->addStringValue("Lowpass 6 dB/oct");
  p->addStringValue("Lowpass 12 dB/oct");
  p->addStringValue("Highpass 6 dB/oct");
  p->addStringValue("Highpass 12 dB/oct");
  p->addStringValue("Bandpass 2*6 dB/oct");
  p->addStringValue("Bandstop 2*6 dB/oct");
  p->addStringValue("Peak/Dip");
  p->addStringValue("Low Shelv 1st order");
  p->addStringValue("Low Shelv 2nd order");
  p->addStringValue("High Shelv 1st order");
  p->addStringValue("High Shelv 2nd order");
  p->addStringValue("Allpass 1st order");
  p->addStringValue("Allpass 2nd order");
  //p->addStringValue("Morph Low/Band/High");
  //p->addStringValue("Morph Low/Peak/High");
  addObservedParameter(p);

  p = new Param("Frequency", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setFrequencyNominal);
  defaultValues.push_back(20.0);
  defaultValues.push_back(200.0);
  defaultValues.push_back(2000.0);
  defaultValues.push_back(20000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
    // todo: define more meaningful default values here - for example tune the frequency to 
    // harmonics (assuming keytrack==100%)

  p = new Param("FrequencyByKey", -200.0, 200.0, 0.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setFrequencyByKey);
  defaultValues.clear(); 
  defaultValues.push_back(0.0);
  defaultValues.push_back(25.0);
  defaultValues.push_back(50.0);
  defaultValues.push_back(75.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("FrequencyByVel", -200.0, 200.0, 0.0,  Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setFrequencyByVel);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Resonance", 0.0, 100.0, 10.0, Parameter::LINEAR, 0.1);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setResonance);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Q", 0.5, 50.0, sqrt(0.5), Parameter::EXPONENTIAL, 0.001);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setQ);
  defaultValues.clear(); 
  defaultValues.push_back(0.5);
  defaultValues.push_back(sqrt(0.5));
  defaultValues.push_back(1.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("PreAllpass", 20.0, 20000.0, 20000.0, Parameter::EXPONENTIAL, 0.0);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setAllpassFreq);
  defaultValues.clear(); 
  defaultValues.push_back(20.0);
  defaultValues.push_back(200.0);
  defaultValues.push_back(2000.0);
  defaultValues.push_back(20000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("MakeUp", 0.0, 100.0, 0.0,  Parameter::LINEAR, 1.0);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setMakeUp);
  defaultValues.clear(); 
  defaultValues.push_back(0.0);
  defaultValues.push_back(25.0);
  defaultValues.push_back(50.0);
  defaultValues.push_back(75.0);
  defaultValues.push_back(100.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Drive", -24.0, 24.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setDrive);
  defaultValues.clear(); 
  defaultValues.push_back(-24.0);
  defaultValues.push_back(-18.0);
  defaultValues.push_back(-12.0);
  defaultValues.push_back(-6.0);
  defaultValues.push_back(-3.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(18.0);
  defaultValues.push_back(24.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Dc", -1.0, 1.0, 0.0, Parameter::LINEAR, 0.01);
  //p->setValueChangeCallback<MMF>(mmf, &MMF::setDc);
  addObservedParameter(p);

  p = new Param("Gain", -60.0, 30.0, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setGain);
  defaultValues.clear(); 
  defaultValues.push_back(-24.0);
  defaultValues.push_back(-18.0);
  defaultValues.push_back(-12.0);
  defaultValues.push_back(-6.0);
  defaultValues.push_back(-3.0);
  defaultValues.push_back(0.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(6.0);
  defaultValues.push_back(12.0);
  defaultValues.push_back(18.0);
  defaultValues.push_back(24.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Morph", -0.99, 0.99, 0.0, Parameter::LINEAR, 0.01);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setMorph);
  defaultValues.clear(); 
  defaultValues.push_back(-0.99);
  defaultValues.push_back(0.5);
  defaultValues.push_back(0.99);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);

  p = new Param("Order", 0.0, 4.0, 4.0, Parameter::LINEAR, 1.0);
  p->setValueChangeCallback<MMF>(mmf, &MMF::setOrder);
  addObservedParameter(p);
}

//=================================================================================================

MultiModeFreqResponseEditor::MultiModeFreqResponseEditor(const juce::String& name) 
  : rsSpectrumPlot(name)
{
  setDescription("Drag around the node to adjust the filter's frequency and resonance, Q or gain");

  ParameterObserver::setIsGuiElement(true);
  filterToEdit = NULL;

  // set up the plot range:
  setMaximumRange(15.625, 32000.0, -60.0, 30.0);
  setCurrentRange(15.625, 32000.0, -60.0, 30.0);
  setHorizontalCoarseGrid(12.0, false);
  setHorizontalFineGrid(   3.0, false);
  setVerticalCoarseGridVisible( false);
  setVerticalFineGridVisible(   false);

  plotColourScheme.setCurveColouringStrategy(PlotColourScheme::UNIFORM);

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  dotRadius = 5.f;

  freqParameter  = NULL;
  resoParameter  = NULL;
  qParameter     = NULL;
  gainParameter  = NULL;
  morphParameter = NULL;

  // this stuff will be (re-) allocated in resized():
  numBins     = 0;
  frequencies = NULL;
  magnitudes  = NULL;

  // activate automation for this ParameterObserver:
  ParameterObserver::setLocalAutomationSwitch(true);
}

MultiModeFreqResponseEditor::~MultiModeFreqResponseEditor(void)
{
  // remove ourselves as listeners from the Parameter objects, such that they do not try to notify 
  // a nonexistent listener:
  ParameterObserver::setLocalAutomationSwitch(false);
  if( freqParameter != NULL )
    freqParameter->deRegisterParameterObserver(this);
  if( resoParameter != NULL )
    resoParameter->deRegisterParameterObserver(this);
  if( qParameter != NULL )
    qParameter->deRegisterParameterObserver(this);
  if( gainParameter != NULL )
    gainParameter->deRegisterParameterObserver(this);
  if( morphParameter != NULL )
    morphParameter->deRegisterParameterObserver(this);

  deleteAndZero(frequencies);
  deleteAndZero(magnitudes);
}

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void MultiModeFreqResponseEditor::setFilterToEdit(rosic::MultiModeFilter* newFilterToEdit)
{
  filterToEdit = newFilterToEdit;
}

void MultiModeFreqResponseEditor::assignParameterFreq(Parameter *parameterToAssign)
{
  freqParameter = parameterToAssign;
  if( freqParameter != NULL )
    freqParameter->registerParameterObserver(this);
}

void MultiModeFreqResponseEditor::assignParameterReso(Parameter *parameterToAssign)
{
  resoParameter = parameterToAssign;
  if( resoParameter != NULL )
    resoParameter->registerParameterObserver(this);
}

void MultiModeFreqResponseEditor::assignParameterQ(Parameter *parameterToAssign)
{
  qParameter = parameterToAssign;
  if( qParameter != NULL )
    qParameter->registerParameterObserver(this);
}

void MultiModeFreqResponseEditor::assignParameterGain(Parameter *parameterToAssign)
{
  gainParameter = parameterToAssign;
  if( gainParameter != NULL )
    gainParameter->registerParameterObserver(this);
}

void MultiModeFreqResponseEditor::assignParameterMorph(Parameter *parameterToAssign)
{
  morphParameter = parameterToAssign;
  if( morphParameter != NULL )
    morphParameter->registerParameterObserver(this);
}

void MultiModeFreqResponseEditor::unAssignParameterFreq()
{
  if( freqParameter != NULL )
    freqParameter->deRegisterParameterObserver(this);
  freqParameter = NULL;
}

void MultiModeFreqResponseEditor::unAssignParameterReso()
{
  if( resoParameter != NULL )
    resoParameter->deRegisterParameterObserver(this);
  resoParameter = NULL;
}

void MultiModeFreqResponseEditor::unAssignParameterQ()
{
  if( qParameter != NULL )
    qParameter->deRegisterParameterObserver(this);
  qParameter = NULL;
}

void MultiModeFreqResponseEditor::unAssignParameterGain()
{
  if( gainParameter != NULL )
    gainParameter->deRegisterParameterObserver(this);
  gainParameter = NULL;
}

void MultiModeFreqResponseEditor::unAssignParameterMorph()
{
  if( morphParameter != NULL )
    morphParameter->deRegisterParameterObserver(this);
  morphParameter = NULL;
}

void MultiModeFreqResponseEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  sendChangeMessage();
  updatePlot();
}

void MultiModeFreqResponseEditor::parameterWillBeDeleted(Parameter* p)
{

  // clear reference to parameter that will be deleted

}

void MultiModeFreqResponseEditor::updateWidgetFromAssignedParameter(bool sendMessage)
{
  updatePlot();
  if( sendMessage == true )
    sendChangeMessage();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void MultiModeFreqResponseEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  setLocalAutomationSwitch(false);

  // call the method which updates the widget:
  updatePlot();
  //updateWidgetFromAssignedParameter(false);


  // switch the wantsAutomationNotification flag on again:  
  setLocalAutomationSwitch(true);
}

void MultiModeFreqResponseEditor::mouseDown(const MouseEvent &e)
{

  if( filterToEdit == NULL )
    return;

  // preliminray: do not open the MIDI-learn menu on right-button - show the image export menu 
  // instead (inherited behaviour from CoordinateSytem):
  if( e.mods.isRightButtonDown() )
    rsPlot::mouseDown(e);
  else
  {
    // get the position of the event in components coordinates
    double x = e.getMouseDownX();
    double y = e.getMouseDownY();

    setupFilterAccordingToMousePosition(x, y);

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

void MultiModeFreqResponseEditor::mouseDrag(const juce::MouseEvent &e)
{
  if( filterToEdit == NULL )
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

  setupFilterAccordingToMousePosition(x, y);

  // send out a change-message:
  sendChangeMessage();
}

void MultiModeFreqResponseEditor::setupFilterAccordingToMousePosition(double mouseX, 
  double mouseY)
{
  if( filterToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  double x = mouseX;
  double y = mouseY;

  // convert them into a frequency and resonance/q/gain value:
  double gain = y;
  fromPixelCoordinates(x, gain);
  gain = RAPT::rsClip(gain, -60.0, 30.0);
  double freq = x;
  double reso = yToReso(y);
  double q    = yToQ(y);

  // restrict ranges (ToDo: actually the filter should take care of the itself....):
  freq = RAPT::rsClip(freq, 20.0, 20000.0);
  reso = RAPT::rsClip(reso, 0.0, 100.0);
  q    = RAPT::rsClip(q, 0.5, 50.0);

  // set up the filter and raise automation events to update other widgets that represent the
  // parameters:
  filterToEdit->setFrequencyNominal(freq);
  if( freqParameter != NULL )
    freqParameter->setValue(freq, true, true);
  if( filterToEdit->getMode() == MultiModeFilterParameters::MOOGISH_LOWPASS )
  {
    filterToEdit->setResonance(reso);
    if( resoParameter != NULL )
      resoParameter->setValue(reso, true, true);
  }
  else if( filterToEdit->currentModeSupportsGain() )
  {
    filterToEdit->setGain(gain);
    if( gainParameter != NULL )
      gainParameter->setValue(gain, true, true);
  }
  else
  {
    filterToEdit->setQ(q);
    if( qParameter != NULL )
      qParameter->setValue(q, true, true);
  }
}

//-------------------------------------------------------------------------------------------------
// drawing:

void MultiModeFreqResponseEditor::resized()
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
  for(int i=0; i<numBins; i++)
    magnitudes[i]  = 0.0;
  updatePlot();
}

void MultiModeFreqResponseEditor::updatePlot()
{
  if( filterToEdit == NULL )
    return;

  // fill the magnitude array with the magnitudes:
  filterToEdit->setFrequencyInstantaneous(filterToEdit->getFrequencyNominal(), true);
  filterToEdit->getMagnitudeResponse(frequencies, magnitudes, numBins, true, false);
  RAPT::rsArrayTools::clip(magnitudes, numBins, -120.0, 120.0);

  // overwrite the magnitude value at the bin closest to the cutoff-frequency with the magnitude at 
  // the exact cutoff frequency:
  double freq  = filterToEdit->getFrequencyNominal();
  double level = RAPT::rsAmpToDbWithCheck(filterToEdit->getMagnitudeAt(freq), 0.000001);
  level        = RAPT::rsClip(level, -120.0, +120.0);

  for(int bin=0; bin < (numBins-1); bin++)
  {
    if( frequencies[bin] < freq && frequencies[bin+1] > freq )
    {
      if( fabs(frequencies[bin]-freq) <= fabs(frequencies[bin+1]-freq) ) // lower bin is closer
        magnitudes[bin]   = level;
      else                                                               // upper bin is closer
        magnitudes[bin+1] = level;
    }
  }

  setSpectrum(numBins, frequencies, magnitudes);
  //repaint();
}

void MultiModeFreqResponseEditor::plotCurveFamily(Graphics &g, juce::Image* targetImage, 
  XmlElement *targetSVG)
{
  if( filterToEdit == NULL )
    return;

  rsDataPlot::plotCurveFamily(g, targetImage, targetSVG);


  //Colour graphColour = colourScheme.curves; // preliminary
  //if( colourScheme.plotColours.size() > 0 )
  //  graphColour = colourScheme.plotColours[0];
  Colour graphColour = plotColourScheme.getCurveColour(0);  
  //g.setColour(graphColour); 

  //double freq = filterToEdit->getFrequencyNominal(); // frequency
  //double reso = filterToEdit->getResonance();        // resonance

  // retrieve characteristic frequency and gain in order to draw the handle (y = gain will be used 
  // only if the mode dictates that - otherwise it will serve as dummy and the actual y-value will
  // be calculated seperately:
  double x = filterToEdit->getFrequencyNominal();
  double y = filterToEdit->getGain();  

  // determine the coordinates of the handle into image component coordinates (for export) or 
  // components coordinates for plot:
  toPixelCoordinates(x, y);

  // y is now the gain in component's coordinates - if we do not have a peaking or shelving type,
  // we need to re-assign it to some value related to resonance or Q:
  if( filterToEdit->getMode() == MultiModeFilterParameters::MOOGISH_LOWPASS )
    y = (float) resoToY( filterToEdit->getResonance());
  else if( filterToEdit->currentModeSupportsGain() ) 
  {
    // keep y to be the transformed gain
  }
  else if( filterToEdit->currentModeSupportsQ() ) 
    y = (float) qToY( filterToEdit->getQ());
  else  
    y = (1.0/3.0) * getHeight();

  // draw the handle and a crosshair:
  g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius), 
    (float) (2*dotRadius), (float) (2*dotRadius) );
  g.setColour(graphColour.withAlpha(0.4f));
  float w = (float) getWidth();
  float h = (float) getHeight();
  g.drawLine(       0,(float)y,        w, (float)y, 1.f);  // horizontal
  g.drawLine((float)x,       0, (float)x,        h, 1.f);  // vertical
}

//-------------------------------------------------------------------------------------------------
// internal functions:

double MultiModeFreqResponseEditor::resoToY(double reso)
{
  return (1.0-0.01*reso) * getHeight();
  //return dotRadius + (1.0-0.01*reso) * (getHeight()-2.f*dotRadius);
}

double MultiModeFreqResponseEditor::yToReso(double y)
{
  return 100.0 * ( 1.0 - y/getHeight() );
  //return 100.0 * ( 1.0 + (dotRadius-y) / (getHeight()-2.0*dotRadius) );
}

double MultiModeFreqResponseEditor::qToY(double q)
{
  return -RAPT::rsExpToLin(q, 0.5, 50.0, -(double)getHeight(), 0.0);
  //return 100.0;
}

double MultiModeFreqResponseEditor::yToQ(double y)
{
  return RAPT::rsLinToExp(-y, -(double)getHeight(), 0.0, 0.5, 50.0);
}

//=================================================================================================


MultiModeFilterModuleEditor::MultiModeFilterModuleEditor(CriticalSection *newPlugInLock, 
  MultiModeFilterAudioModule* newMultiModeFilterAudioModule) 
  : AudioModuleEditor(newMultiModeFilterAudioModule)
{
  jassert(newMultiModeFilterAudioModule != NULL ); // you must pass a valid module here
  filterToEdit = newMultiModeFilterAudioModule->wrappedMultiModeFilter;
  setHeadlineText("Filter");
  stateWidgetSet->stateLabel->setVisible(false);
  webLink->setVisible(false);
  infoField->setVisible(false);
  createWidgets();
  isTopLevelEditor = false;
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void MultiModeFilterModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( filterToEdit == NULL )
    return;
  if( buttonThatWasClicked == twoStagesButton )
  {
    filterToEdit->useTwoStages(twoStagesButton->getToggleState());
    frequencyResponseDisplay->updatePlot();
  }
  moduleToEdit->markStateAsDirty();
  //sendChangeMessage();
}

void MultiModeFilterModuleEditor::rComboBoxChanged(RComboBox *rComboBoxThatHasChanged)
{
  if( filterToEdit == NULL )
    return;

  if( rComboBoxThatHasChanged == modeComboBox )
  {
    //filterToEdit->setMode(modeComboBox->getSelectedId());
    //filterToEdit->setMode( stringToFilterModeIndex(modeComboBox->getText()) );
    updateWidgetArrangement();
    updateWidgetsAccordingToState();  
  }

  moduleToEdit->markStateAsDirty();
  //setPresetDirty();
  //frequencyResponseDisplay->updatePlot(); // already called from  updateWidgetsAccordingToState()
  //sendChangeMessage();
}

void MultiModeFilterModuleEditor::rSliderValueChanged(RSlider *sliderThatHasChanged)
{
  // there are some sliders for parameters which are not observed by the plot itself - when one of 
  // these changes, we do the update here:
  if(  sliderThatHasChanged == orderSlider || sliderThatHasChanged == makeUpSlider )
    frequencyResponseDisplay->updatePlot(); 
}


void MultiModeFilterModuleEditor::resized()
{
  linkPosition = INVISIBLE;
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  y = getPresetSectionBottom();
  frequencyResponseDisplay->setBounds(x, y+4, w, h-92);
  y = frequencyResponseDisplay->getBottom();

  modeComboBox->setBounds(x+4,  y+4, w/2-4, 20);

  y = modeComboBox->getBottom();
  freqSlider->setBounds(modeComboBox->getX(), y+4, w/2-4 , 16);

  y = freqSlider->getBottom()-2;
  w = freqSlider->getWidth();
  freqByKeySlider->setBounds(freqSlider->getX(),           y, w/2-2, 16);
  freqByVelSlider->setBounds(freqSlider->getX() + w/2 + 2, y, w/2-2, 16);

  y = modeComboBox->getY();
  w = getWidth();
  x = w/2;
  resoSlider->setBounds(x+4, y, w/2-8, 16);
  qSlider->setBounds(resoSlider->getBounds());

  y = resoSlider->getBottom();

  preAllpassSlider->setBounds(x+4,     y+4, w/4-4, 16);
  orderSlider->setBounds(     x+w/4+4, y+4, w/4-8, 16);
  twoStagesButton->setBounds(orderSlider->getBounds());

  gainSlider->setBounds(x+4, y+4, w/2-8, 16);
  morphSlider->setBounds(x+4, y+4, w/2-8, 16);

  y += 20;

  driveSlider->setBounds( x+4,     y+4, w/4-4, 16);
  makeUpSlider->setBounds(x+w/4+4, y+4, w/4-8, 16);
}

void MultiModeFilterModuleEditor::updateWidgetsAccordingToState()
{
  if( filterToEdit == NULL )
    return;

  //modeComboBox->setText(filterModeIndexToString(filterToEdit->getMode()), true); //old

  //modeComboBox->selectItemFromText(filterModeIndexToString(filterToEdit->getMode()), false);

  // something to do here...?

  freqSlider->setValue(filterToEdit->getFrequencyNominal(),               false, false);
  freqByKeySlider->setValue(filterToEdit->getFrequencyByKey(),            false, false);
  freqByVelSlider->setValue(filterToEdit->getFrequencyByVel(),            false, false);
  resoSlider->setValue(filterToEdit->getResonance(),                      false, false);
  qSlider->setValue(filterToEdit->getQ(),                                 false, false);
  gainSlider->setValue(filterToEdit->getGain(),                           false, false);
  driveSlider->setValue(filterToEdit->getDrive(),                         false, false);
  morphSlider->setValue(filterToEdit->getMorph(),                         false, false);
  /*
  freq2ScaleSlider->setValue(
  filterToEdit->twoStageBiquad.getSecondFrequencyScaleFactor(),         false, false);
  freq2OffsetSlider->setValue(
  filterToEdit->twoStageBiquad.getSecondFrequencyOffset(),              false, false);
  q2ScaleSlider->setValue(
  filterToEdit->twoStageBiquad.getSecondQScaleFactor(),                 false, false);
  gain2ScaleSlider->setValue(
  filterToEdit->twoStageBiquad.getSecondGainScaleFactor(),              false, false);
  */
  orderSlider->setValue(filterToEdit->getOrder(),                         false, false);
  twoStagesButton->setToggleState(filterToEdit->usesTwoStages(),          false);
  preAllpassSlider->setValue(filterToEdit->getAllpassFreq(),              false, false);
  makeUpSlider->setValue(filterToEdit->getMakeUp(),                       false, false);

  updateWidgetArrangement();
  frequencyResponseDisplay->updatePlot();

  stateWidgetSet->updateStateNameField();
}

void MultiModeFilterModuleEditor::createWidgets()
{
  typedef rsAutomatableSlider Sld;  Sld* s;
  typedef rsAutomatableButton Btn;  Btn* b;
  //typedef RTextField        Txf;  Txf* t;

  addWidget( modeComboBox = new RNamedComboBox("ModeComboBox", "Type:") );
  modeComboBox->assignParameter(moduleToEdit->getParameterByName("Mode") );
  modeComboBox->registerComboBoxObserver(this);
  modeComboBox->setDescription("Type/Mode of the filter.");

  addWidget( twoStagesButton = b = new Btn("2 Stages") );
  b->assignParameter(moduleToEdit->getParameterByName("TwoStages"));
  b->addRButtonListener(this);
  b->setDescription("Switch second filter stage on/off");
  b->setClickingTogglesState(true);

  addWidget( freqSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Frequency"));
  s->setDescription("Characteristic frequency of the filter");
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( freqByKeySlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("FrequencyByKey"));
  s->setSliderName("Key");
  s->setDescription("Key tracking of the filter's frequency");
  s->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( freqByVelSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("FrequencyByVel"));
  s->setSliderName("Vel");
  s->setDescription("Velocity tracking of the filter's frequency");
  s->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( resoSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Resonance"));
  s->setDescription("Resonance amount of the filter");
  s->setStringConversionFunction(&percentToStringWithUnit1);

  addWidget( qSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Q"));
  s->setSliderName("Q");
  s->setDescription("Quality factor of the filter");
  s->setStringConversionFunction(&valueToString3);

  addWidget( driveSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Drive"));
  s->setDescription("Drives the filter into distortion");
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( orderSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Order"));
  s->setDescription("Selects the order of the filter - this affects the slope");
  s->setStringConversionFunction(&valueToString0);
  s->addListener(this);

  addWidget( gainSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Gain"));
  s->setSliderName("Gain");
  s->setDescription("Gain for peak and shelving filter types");
  s->setStringConversionFunction(&decibelsToStringWithUnit2);

  addWidget( morphSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("Morph") );
  s->setDescription("Morph between filter types");
  s->setStringConversionFunction(&valueToString2);
  //s->setRange(-0.99, 0.99, 0.01, -0.99);
  //s->setScaling(Parameter::LINEAR_BIPOLAR);
  //automatableSliders.addIfNotAlreadyThere(morphSlider);
  //morphSlider->addListener(this);

  //addWidget( transitionSlider = s = new Sld );
  //s->assignParameter(moduleToEdit->getParameterByName("Transition") );
  //s->setSliderName("Transition");
  //s->setDescription("Determines the transition when morphing between filter types");
  //s->setStringConversionFunction(&valueToString3);

  addWidget( preAllpassSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("PreAllpass") );
  s->setSliderName("Allpass");
  s->setDescription("First order allpass before actual filter to pre-shape waveform");
  s->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( makeUpSlider = s = new Sld );
  s->assignParameter(moduleToEdit->getParameterByName("MakeUp") );
  s->setDescription("Compensates low frequency losses at high resonance via gain");
  s->setStringConversionFunction(&percentToStringWithUnit0);
  s->addListener(this);

  frequencyResponseDisplay = new MultiModeFreqResponseEditor("SpectrumEditor");
  frequencyResponseDisplay->setFilterToEdit(filterToEdit);
  frequencyResponseDisplay->addChangeListener(this);
  frequencyResponseDisplay->assignParameterFreq( moduleToEdit->getParameterByName("Frequency"));
  frequencyResponseDisplay->assignParameterReso( moduleToEdit->getParameterByName("Resonance"));
  frequencyResponseDisplay->assignParameterQ(    moduleToEdit->getParameterByName("Q"));
  frequencyResponseDisplay->assignParameterGain( moduleToEdit->getParameterByName("Gain"));
  frequencyResponseDisplay->assignParameterMorph(moduleToEdit->getParameterByName("Morph") );
  addPlot( frequencyResponseDisplay );

  // customize the descriptions for the load/save buttons:
  stateWidgetSet->stateLoadButton->setDescription(   "Load filter settings from file");
  stateWidgetSet->stateSaveButton->setDescription(   "Save filter settings to file");
  stateWidgetSet->statePlusButton->setDescription(   "Skip to next filter settings file in current directory");
  stateWidgetSet->stateMinusButton->setDescription(  "Skip to previous filter settings file in current directory");
  stateWidgetSet->stateFileNameLabel->setDescription("Name of current preset for the filter section (if any)");
}

void MultiModeFilterModuleEditor::updateWidgetArrangement()
{
  if( filterToEdit == NULL )
    return;

  //int x, y, w, h;

  setWidgetsVisible(false);
  modeComboBox->setVisible(true);

  int mode = filterToEdit->getMode();


  // in bypass-mode, we dont need any widgets so we return:
  if( mode == rosic::MultiModeFilterParameters::BYPASS )
    return;

  // the frequency related sliders must be visible for all modes:
  freqSlider->setVisible(true);
  freqByKeySlider->setVisible(true);
  freqByVelSlider->setVisible(true);

  if( mode == rosic::MultiModeFilterParameters::MOOGISH_LOWPASS )
  {
    arrangeWidgetsForMoogishLowpassMode();
  }
  else  // twoStageBiquad
  {
    // all types of the TwoStageBiquad support the two-stages switch, so make the button visible:
    twoStagesButton->setVisible(true);

    if(    mode == rosic::MultiModeFilterParameters::LOWPASS_6 
      || mode == rosic::MultiModeFilterParameters::HIGHPASS_6 
      || mode == rosic::MultiModeFilterParameters::ALLPASS_1ST )
    {
      arrangeWidgetsForFirstOrderWithoutGain();
    }
    else if(    mode == rosic::MultiModeFilterParameters::LOWPASS_RBJ
      || mode == rosic::MultiModeFilterParameters::HIGHPASS_RBJ
      || mode == rosic::MultiModeFilterParameters::BANDPASS_RBJ
      || mode == rosic::MultiModeFilterParameters::BANDREJECT_RBJ
      || mode == rosic::MultiModeFilterParameters::ALLPASS_RBJ )
    {
      arrangeWidgetsForSecondOrderWithoutGain();
    }
    else if(    mode == rosic::MultiModeFilterParameters::LOW_SHELV_1ST
      || mode == rosic::MultiModeFilterParameters::HIGH_SHELV_1ST )
    {
      arrangeWidgetsForFirstOrderWithGain();
    }
    else if(    mode == rosic::MultiModeFilterParameters::LOW_SHELV_RBJ
      || mode == rosic::MultiModeFilterParameters::HIGH_SHELV_RBJ
      || mode == rosic::MultiModeFilterParameters::PEAK_OR_DIP_RBJ )
    {
      arrangeWidgetsForSecondOrderWithGain();
    }
    else if(    mode == rosic::MultiModeFilterParameters::MORPH_LP_PK_HP )
    {
      arrangeWidgetsForMorphableMode();
    }

    /*
    // the first order types don't support 'Q', for the others we must make the slider visible:
    if(    filterToEdit->getMode() != rosic::MultiModeFilterParameters::LOWPASS_6 
    && filterToEdit->getMode() != rosic::MultiModeFilterParameters::HIGHPASS_6 
    && filterToEdit->getMode() != rosic::MultiModeFilterParameters::LOW_SHELV_1ST 
    && filterToEdit->getMode() != rosic::MultiModeFilterParameters::HIGH_SHELV_1ST
    && filterToEdit->getMode() != rosic::MultiModeFilterParameters::ALLPASS_1ST )
    {
    qSlider->setVisible(true);
    }
    // the peak and shelving types support the gain parameter:
    if(    filterToEdit->getMode() == rosic::MultiModeFilterParameters::PEAK_OR_DIP_RBJ 
    || filterToEdit->getMode() == rosic::MultiModeFilterParameters::LOW_SHELV_1ST
    || filterToEdit->getMode() == rosic::MultiModeFilterParameters::LOW_SHELV_RBJ
    || filterToEdit->getMode() == rosic::MultiModeFilterParameters::HIGH_SHELV_1ST
    || filterToEdit->getMode() == rosic::MultiModeFilterParameters::HIGH_SHELV_RBJ   )
    {
    gainSlider->setVisible(true);
    }

    if(    filterToEdit->getMode() == rosic::MultiModeFilterParameters::MORPH_LP_BP_HP
    || filterToEdit->getMode() == rosic::MultiModeFilterParameters::MORPH_LP_PK_HP )
    {
    morphSlider->setVisible(true);
    //twoStagesButton->setVisible(true);
    //transitionSlider->setVisible(true);
    y = morphSlider->getBottom();
    x = morphSlider->getRight()-64;
    w = 64;
    h = 16;
    twoStagesButton->setBounds(x, y+4, w, h);
    int dummy = 0;
    }
    */

  }

}

/*
void MultiModeFilterEditor::arrangeCommonWidgets()
{

}
*/

void MultiModeFilterModuleEditor::arrangeWidgetsForMoogishLowpassMode()
{
  resoSlider->setVisible(true);
  driveSlider->setVisible(true);
  orderSlider->setVisible(true);
  preAllpassSlider->setVisible(true);
  makeUpSlider->setVisible(true);
}

void MultiModeFilterModuleEditor::arrangeWidgetsForFirstOrderWithoutGain()
{
  int x = modeComboBox->getRight()+4;
  int y = modeComboBox->getY();
  twoStagesButton->setVisible(true);
  twoStagesButton->setBounds(x+4, y, 80, 16);
}

void MultiModeFilterModuleEditor::arrangeWidgetsForSecondOrderWithoutGain()
{
  int x = modeComboBox->getRight()+4;
  int y = modeComboBox->getY();
  int w = freqSlider->getWidth();

  twoStagesButton->setVisible(true);
  twoStagesButton->setBounds(x+4, y, 80, 16);
  y += 20;
  w  = getWidth()-x;
  qSlider->setVisible(true);
  qSlider->setBounds(x+4, y, w-8, 16);
}

void MultiModeFilterModuleEditor::arrangeWidgetsForFirstOrderWithGain()
{
  arrangeWidgetsForFirstOrderWithoutGain();
  gainSlider->setVisible(true);
  gainSlider->setBounds(qSlider->getX(), qSlider->getY(), qSlider->getWidth(), 16);
}

void MultiModeFilterModuleEditor::arrangeWidgetsForSecondOrderWithGain()
{
  arrangeWidgetsForSecondOrderWithoutGain();
  gainSlider->setVisible(true);
  gainSlider->setBounds(qSlider->getX(), qSlider->getY()+20, qSlider->getWidth(), 16);
}

void MultiModeFilterModuleEditor::arrangeWidgetsForMorphableMode()
{
  arrangeWidgetsForSecondOrderWithoutGain();
  morphSlider->setVisible(true);
  morphSlider->setBounds(qSlider->getX(), qSlider->getY()+20, qSlider->getWidth(), 16);
}
