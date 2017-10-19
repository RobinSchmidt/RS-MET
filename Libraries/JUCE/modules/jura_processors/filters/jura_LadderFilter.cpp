
Ladder::Ladder(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse,
  ModulationManager* modManagerToUse) 
  : ModulatableAudioModule(lockToUse, metaManagerToUse, modManagerToUse)
{
  ScopedLock scopedLock(*lock);
  moduleName = "Ladder";
  setActiveDirectory(getApplicationDirectory() + "/LadderPresets");

  createStaticParameters();
}

void Ladder::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  std::vector<double> defaultValues;
  //typedef MetaControlledParameter Param;
  typedef ModulatableParameter Param;
  Param* p;

  p = new Param("Cutoff", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL);
  p->setDefaultModParameters(20, 20000, -2, 3, true);
  //p->setModulationRangeMin(20.0);
  //p->setModulationRangeMax(20000.0);
  //p->setDefaultModulationDepthMin(-2.0);
  //p->setDefaultModulationDepthMax(+3.0);
  //p->setDefaultModulationModeRelative(true);
  defaultValues.clear();
  defaultValues.push_back(125.0);
  defaultValues.push_back(250.0);
  defaultValues.push_back(500.0);
  defaultValues.push_back(1000.0);
  defaultValues.push_back(2000.0);
  defaultValues.push_back(4000.0);
  defaultValues.push_back(8000.0);
  defaultValues.push_back(16000.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p); 
    // maybe we should have a subclass ModulatableAudioModule of AudioModule and a function
    // addModulatableParameter - this function could register the ModulationTarget for the
    // passed ModulatableParameter
  p->setValueChangeCallback<Ladder>(this, &Ladder::setCutoff);

  p = new Param("Resonance", 0.0, 1.0, 0.2, Parameter::LINEAR, 0.01);
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setResonance);

  p = new Param("StereoSpread", -24.0, +24.0, 0.0, Parameter::LINEAR_BIPOLAR);
  //p->setMapper(new rsParameterMapperSinh(-24, +24, 2.0)); // for test - something is wrong - maybe the slider class can't deal with it
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setStereoSpread);

  p = new Param("Mode", 0.0, 14.0, 4.0, Parameter::STRING);
  p->addStringValue("Flat");
  p->addStringValue("Lowpass 6 dB/oct");
  p->addStringValue("Lowpass 12 dB/oct");
  p->addStringValue("Lowpass 18 dB/oct");
  p->addStringValue("Lowpass 24 dB/oct");
  p->addStringValue("Highpass 6 dB/oct");
  p->addStringValue("Highpass 12 dB/oct");
  p->addStringValue("Highpass 18 dB/oct");
  p->addStringValue("Highpass 24 dB/oct");
  p->addStringValue("Bandpass 6/6 dB/oct");
  p->addStringValue("Bandpass 6/12 dB/oct");
  p->addStringValue("Bandpass 6/18 dB/oct");
  p->addStringValue("Bandpass 12/6 dB/oct");
  p->addStringValue("Bandpass 12/12 dB/oct");
  p->addStringValue("Bandpass 18/6 dB/oct");
  addObservedParameter(p);
  p->setValueChangeCallback<Ladder>(this, &Ladder::setMode);

  //// make sure that the parameters are initially in sync with the audio engine:
  //for(int i = 0; i < (int)parameters.size(); i++)
  //  parameters[i]->resetToDefaultValue(true, true); // do we need this?
  // replace this loop by a call to:
  //resetParametersToDefaultValues();
}

//-------------------------------------------------------------------------------------------------
// Editor creation:

AudioModuleEditor* Ladder::createEditor()
{
  return new jura::LadderEditor(this);
}

//-------------------------------------------------------------------------------------------------
// the audio processing callbacks:

void Ladder::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
  {
    inOutBuffer[0][n] = ladderL.getSample(inOutBuffer[0][n]);
    inOutBuffer[1][n] = ladderR.getSample(inOutBuffer[1][n]);
  }
}

void Ladder::processStereoFrame(double *left, double *right)
{
  *left  = ladderL.getSample(*left);
  *right = ladderR.getSample(*right);
}

//-------------------------------------------------------------------------------------------------
// parameter setters (callback targets for the Parameter objects):

void Ladder::setSampleRate(double newSampleRate){
  ladderL.setSampleRate(newSampleRate); 
  ladderR.setSampleRate(newSampleRate);
}
void Ladder::setCutoff(double newCutoff){
  cutoff = newCutoff;
  ladderL.setCutoff(freqFactorL * cutoff); 
  ladderR.setCutoff(freqFactorR * cutoff);
}
void Ladder::setResonance(double newResonance){
  ladderL.setResonance(newResonance);
  ladderR.setResonance(newResonance);
}
void Ladder::setMode(int newMode){
  ladderL.setMode(newMode);
  ladderR.setMode(newMode);
}
void Ladder::reset(){
  ladderL.reset();
  ladderR.reset();
}

inline double pitchOffsetToFreqFactor(double pitchOffset){ // move to RAPT library, eventually
  return exp(0.057762265046662109118102676788181 * pitchOffset);
}
void Ladder::setStereoSpread(double newSpread){
  freqFactorL = pitchOffsetToFreqFactor(+0.5*newSpread);
  freqFactorR = pitchOffsetToFreqFactor(-0.5*newSpread);  // == 1 / freqFactorL -> optimize later
  setCutoff(cutoff);
  // By spreading the cutoff frequencies of left and right channel, we not only get a stereo 
  // effect, but also a formant-filter like sound. Maybe the effectiveness can be improved, if we
  // have another parameter that makes this spreading dependent on the cutoff frequency. This way, 
  // we could control distance of the formants as function of the nominal center frequency between 
  // them. I think that should allow for more realistic formant sweeps. The spread parameter is 
  // also an interesting target for an LFO.
}
void Ladder::setMidSideMode(bool shouldBeInMidSideMode){
  midSideMode = shouldBeInMidSideMode;
}

// inquiry:

double Ladder::getMagnitudeAt(double frequency)
{
  return ladderL.getMagnitudeResponseAt(frequency);
}

void Ladder::getMagnitudeResponse(const double *frequencies, double *magnitudes, int numBins,
  bool inDecibels)
{
  int i;
  for(i = 0; i < numBins; i++)
    magnitudes[i] = getMagnitudeAt(frequencies[i]);
  if(inDecibels)
    for(i = 0; i < numBins; i++)
      magnitudes[i] = amp2dBWithCheck(magnitudes[i], 1.e-6);
}

//=================================================================================================
// LadderSpectrumEditor:

// construction/destruction:

LadderSpectrumEditor::LadderSpectrumEditor(const juce::String& name) 
  : SpectrumDisplayOld(name)
{
  // nulling these pointers must occur before setting up the plot range:
  freqParameter = nullptr;
  resoParameter = nullptr;
  filterToEdit  = nullptr;

  // this stuff will be (re-) allocated in resized():
  numBins     = 0;
  frequencies = NULL;
  magnitudes  = NULL;


  setDescription("Drag around the node to adjust the filter's frequency and resonance");
  ParameterObserver::isGuiElement = true;

  // set up the plot range:
  setMaximumRange(15.625, 32000.0, -60.0, 60.0);
  setCurrentRange(15.625, 32000.0, -60.0, 60.0);
  setHorizontalCoarseGrid(12.0, false);
  setHorizontalFineGrid(   3.0, false);
  setVerticalCoarseGridVisible( false);
  setVerticalFineGridVisible(   false);

  plotColourScheme.setCurveColouringStrategy(PlotColourScheme::UNIFORM);

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  dotRadius = 5.f;

  // activate automation for this ParameterObserver:
  ParameterObserver::localAutomationSwitch = true;
}

LadderSpectrumEditor::~LadderSpectrumEditor(void)
{
  // remove ourselves as listeners from the Parameter objects, such that they do not try to notify 
  // a nonexistent listener:
  ParameterObserver::localAutomationSwitch = false;
  if( freqParameter != NULL )
    freqParameter->deRegisterParameterObserver(this);
  if( resoParameter != NULL )
    resoParameter->deRegisterParameterObserver(this);

  deleteAndZero(frequencies);
  deleteAndZero(magnitudes);
}

//-------------------------------------------------------------------------------------------------
// parameter-settings:

void LadderSpectrumEditor::setFilterToEdit(jura::Ladder* newFilterToEdit)
{
  filterToEdit = newFilterToEdit;
}

void LadderSpectrumEditor::assignParameterFreq(Parameter *parameterToAssign)
{
  freqParameter = parameterToAssign;
  if( freqParameter != NULL )
    freqParameter->registerParameterObserver(this);
}

void LadderSpectrumEditor::assignParameterReso(Parameter *parameterToAssign)
{
  resoParameter = parameterToAssign;
  if( resoParameter != NULL )
    resoParameter->registerParameterObserver(this);
}

void LadderSpectrumEditor::unAssignParameterFreq()
{
  if( freqParameter != NULL )
    freqParameter->deRegisterParameterObserver(this);
  freqParameter = NULL;
}

void LadderSpectrumEditor::unAssignParameterReso()
{
  if( resoParameter != NULL )
    resoParameter->deRegisterParameterObserver(this);
  resoParameter = NULL;
}

void LadderSpectrumEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  //// old code:
  //// send out a change-message: 
  //sendChangeMessage(); // ? why that?
  ////updatePlot();

  // new:
  updatePlot();
}

void LadderSpectrumEditor::parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted)
{
  // clear reference to parameter that will be deleted
}

void LadderSpectrumEditor::updateWidgetFromAssignedParameter(bool sendMessage)
{
  updatePlot();
  if( sendMessage == true )
    sendChangeMessage();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void LadderSpectrumEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  // temporarily switch the wantsAutomationNotification flag from the ParameterObserver base 
  // class off to avoid circular notifications and updates:
  localAutomationSwitch = false;

    // ?? what is this good for ?? we don't send out any messages here... fihure out in debugger

  // call the method which updates the widget:
  updatePlot();
  //updateWidgetFromAssignedParameter(false);

  // switch the wantsAutomationNotification flag on again:  
  localAutomationSwitch = true;
}

void LadderSpectrumEditor::mouseDown(const MouseEvent &e)
{
  if( filterToEdit == NULL )
    return;

  // preliminary: do not open the MIDI-learn menu on right-button - show the image export menu 
  // instead (inherited behaviour from CoordinateSytem):
  if( e.mods.isRightButtonDown() )
    CoordinateSystemOld::mouseDown(e);
  else
  {
    // get the position of the event in components coordinates
    double x = e.getMouseDownX();
    double y = e.getMouseDownY();

    setupFilterAccordingToMousePosition(x, y);

    // send out a change-message:
    sendChangeMessage();
  }
}

void LadderSpectrumEditor::mouseDrag(const juce::MouseEvent &e)
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

void LadderSpectrumEditor::setupFilterAccordingToMousePosition(double mouseX, double mouseY)
{
  if( filterToEdit == NULL )
    return;

  // get the position of the event in components coordinates
  double x = mouseX;
  double y = mouseY;

  // convert them into a frequency and resonance/q/gain value:
  double gain = y;
  transformFromComponentsCoordinates(x, gain);
  gain = clip(gain, -60.0, 30.0);
  double freq = x;
  double reso = yToReso(y);

  // restrict ranges (ToDo: actually the filter should take care of the itself....):
  //freq = clip(freq, 20.0, 20000.0);
  //reso = clip(reso,  0.0,     1.0);   // change the upper limit later...

  // set up the filter and raise automation events to update other widgets that represent the
  // parameters:
  if(freqParameter != NULL)
    freqParameter->setValue(freq, true, true);
  if(resoParameter != NULL)
    resoParameter->setValue(reso, true, true);
}

//-------------------------------------------------------------------------------------------------
// drawing:

void LadderSpectrumEditor::resized()
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
  for(int i = 0; i < numBins; i++)
    magnitudes[i]  = 0.0;
  updatePlot();
}

inline void clipBuffer(double *buffer, int length, double min, double max)
{
  for(int i = 0; i < length; i++)
    buffer[i] = clip(buffer[i], min, max);
}
void LadderSpectrumEditor::updatePlot()
{
  if( filterToEdit == nullptr || freqParameter == nullptr )
    return;

  // fill the magnitude array with the magnitudes:
  filterToEdit->getMagnitudeResponse(frequencies, magnitudes, numBins, true);
  clipBuffer(magnitudes, numBins, -120.0, 120.0);

  // We overwrite the magnitude value at the bin closest to the cutoff-frequency with the magnitude 
  // at the exact cutoff frequency. This avoids some graphic artifacts when sweeping cutoff at high
  // resonance (the peak jumps up and down, depending on how close the bin center is to the actual
  // cutoff/resonance frequency):
  double freq  = freqParameter->getValue();
  double level = amp2dBWithCheck(filterToEdit->getMagnitudeAt(freq), 0.000001);
  level        = clip(level, -120.0, +120.0);
  for(int bin = 0; bin < numBins-1; bin++)
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

void LadderSpectrumEditor::plotCurveFamily(Graphics &g, juce::Image* targetImage, 
  XmlElement *targetSVG)
{
  if(freqParameter == nullptr || resoParameter == nullptr) // pointers need to be assigned
    return;

  CurveFamilyPlotOld::plotCurveFamily(g, targetImage, targetSVG);

  //Colour graphColour = colourScheme.curves; // preliminary
  //if( colourScheme.plotColours.size() > 0 )
  //  graphColour = colourScheme.plotColours[0];
  Colour graphColour = plotColourScheme.getCurveColour(0);  
  //g.setColour(graphColour); 

  double x = freqParameter->getValue();    // frequency in Hz
  double y = 0;                            // dummy, for the moment

  // determine the coordinates of the handle into image component coordinates (for export) or 
  // components coordinates for plot:
  if( targetImage == NULL )
    transformToComponentsCoordinates(x, y);
  else
    transformToImageCoordinates(x, y, targetImage);

  y = resoToY(resoParameter->getValue());  // ..now we assign also y


  // draw the handle and a crosshair:
  g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius), 
    (float) (2*dotRadius), (float) (2*dotRadius) );
  g.setColour(graphColour.withAlpha(0.4f));
  float w = (float) getPlotWidth(targetImage);
  float h = (float) getPlotHeight(targetImage);
  g.drawLine(       0,(float)y,        w, (float)y, 1.f);  // horizontal
  g.drawLine((float)x,       0, (float)x,        h, 1.f);  // vertical
}

// internal functions:

double LadderSpectrumEditor::resoToY(double reso, juce::Image *targetImage)
{
  double maxRes = resoParameter->getMaxValue();
  double y = (1 - reso/maxRes) * getPlotHeight(targetImage);
  return y;
}

double LadderSpectrumEditor::yToReso(double y, juce::Image *targetImage)
{
  double maxRes = resoParameter->getMaxValue();
  double reso = maxRes * ( 1.0 - y/getPlotHeight(targetImage) );    
  return reso;
}
// maybe use not the full plot height but plotHeight - dotRadius or something - at very high and
// low resonances, we have to go to the border of the plot with the mouse and half of the dot goes
// outside the visible area, which is not nice

//=================================================================================================
// the GUI editor class for the Ladder:

LadderEditor::LadderEditor(jura::Ladder *newLadderToEdit) : AudioModuleEditor(newLadderToEdit)
{
  ScopedLock scopedLock(*lock);
  // maybe we should avoid this lock here and instead have a function that connects the widgets 
  // with the parameters where we acquire the lock - but maybe not

  // set the plugIn-headline:
  setHeadlineText("Ladder");

  // assign the pointer to the edited object:
  jassert(newLadderToEdit != nullptr ); // you must pass a valid module here
  ladderToEdit = newLadderToEdit;

  // create the widgets and assign the automatable parameters to them:

  frequencyResponseDisplay = new LadderSpectrumEditor("SpectrumEditor");
  frequencyResponseDisplay->setFilterToEdit(ladderToEdit);
  frequencyResponseDisplay->addChangeListener(this); // do we need this?
  frequencyResponseDisplay->assignParameterFreq( moduleToEdit->getParameterByName("Cutoff"));
  frequencyResponseDisplay->assignParameterReso( moduleToEdit->getParameterByName("Resonance"));
  addPlot( frequencyResponseDisplay );

  addWidget( cutoffSlider = new ModulatableSlider() );
  cutoffSlider->assignParameter( ladderToEdit->getParameterByName("Cutoff") );
  cutoffSlider->setSliderName("Cutoff");
  cutoffSlider->setDescription("Cutoff frequency in Hz");
  cutoffSlider->setDescriptionField(infoField);
  cutoffSlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( resonanceSlider = new ModulatableSlider() );
  resonanceSlider->assignParameter( ladderToEdit->getParameterByName("Resonance") );
  resonanceSlider->setSliderName("Resonance");
  resonanceSlider->setDescription("Amount of feedback");
  resonanceSlider->setDescriptionField(infoField);
  resonanceSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( spreadSlider = new ModulatableSlider() );
  spreadSlider->assignParameter( ladderToEdit->getParameterByName("StereoSpread") );
  spreadSlider->setSliderName("Spread");
  spreadSlider->setDescription("Detunes cutoff frequencies of channels");
  spreadSlider->setDescriptionField(infoField);
  spreadSlider->setStringConversionFunction(&valueToStringTotal5);

  addWidget( modeComboBox = new AutomatableComboBox() );
  modeComboBox->assignParameter( ladderToEdit->getParameterByName("Mode") );
  modeComboBox->setDescription("Select frequency response type");
  modeComboBox->setDescriptionField(infoField);
  modeComboBox->registerComboBoxObserver(this); // to update plot when mode is switched
   // todo: pass the mode parameter to the frequency response plot, too - when use does 
   // right-click, the popup opens and the mode can be selected. then we don't need to observe
   // the combobox here...
   // maybe, we can allow also different modes for the two filters and different resonances

  // change to midSideButton
  //addWidget( invertButton = new RButton(juce::String(T("Invert"))) );
  //invertButton->assignParameter( pitchShifterModuleToEdit->getParameterByName(T("Invert")) );
  //invertButton->setDescription(juce::String(T("Invert polarity of wet (shifted) signal")));
  //invertButton->setDescriptionField(infoField);
  //invertButton->setClickingTogglesState(true);

  // set up the widgets:
  updateWidgetsAccordingToState();

  setSize(420, 240);  
  // There's a bug somewhere - if we request a size of 400 x 300, the widgets become invisible.
  // I guess, 400 x 300 is the size, it already has by default, and when we set it again to this 
  // size, it will be recognized that nothing changed and some update function is not being called
  // ...or something. ...need to figure that out
  // ..i think, resized() is not being called and thus, the widget positions won't be set up.
  // Maybe, the default size should not be 400 x 300 but something like 0 x 0 or (1 x 1, if that's
  // not allowed) - something that is assured to be never requested as actual size of an editor.
}

void LadderEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int x  = 0;
  int y  = 0;
  int w  = getWidth();
  int w2 = w/2;
  int h  = getHeight();

  y = getPresetSectionBottom();

  frequencyResponseDisplay->setBounds(x, y+4, w, h-y-2*20-8); // 2*20 for 2 widget-rows below it
  y = frequencyResponseDisplay->getBottom();

  cutoffSlider->setBounds(      4, y+4, w2-4, 16);
  resonanceSlider->setBounds(w2+4, y+4, w2-8, 16);

  y = cutoffSlider->getBottom();  

  modeComboBox->setBounds(   4, y+4, w2-4, 16);
  spreadSlider->setBounds(w2+4, y+4, w2-8, 16);

  //y = modeComboBox->getBottom(); 

  // maybe here, we somehow have to also resize our wrapper, if any
}

void LadderEditor::rComboBoxChanged(RComboBox* comboBoxThatHasChanged)
{
  frequencyResponseDisplay->updatePlot();
}