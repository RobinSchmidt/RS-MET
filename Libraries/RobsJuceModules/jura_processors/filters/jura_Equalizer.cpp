
//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

EqualizerAudioModule::EqualizerAudioModule(CriticalSection *newPlugInLock,
  rosic::EqualizerStereo *equalizerStereoToWrap) : ModulatableAudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*lock);
  jassert(equalizerStereoToWrap != NULL); // you must pass a valid rosic-object to the constructor
  //  ---nah, this Module admits NULL pointers for use inside EchoLab ...why?
  wrappedEqualizerStereo = equalizerStereoToWrap;
  init();
}

EqualizerAudioModule::EqualizerAudioModule(CriticalSection *newPlugInLock)
  : ModulatableAudioModule(newPlugInLock)
{
  ScopedLock scopedLock(*lock);
  wrappedEqualizerStereo = new rosic::EqualizerStereo;
  wrappedEqualizerIsOwned = true;
  init();
}

void EqualizerAudioModule::init()
{
  setModuleTypeName("Equalizer");
  selectedChannel  = 0;
  selectedIndex    = -1;
  patchFormatIndex = 2;
  createStaticParameters();
}

EqualizerAudioModule::~EqualizerAudioModule()
{
  if(wrappedEqualizerIsOwned)
    delete wrappedEqualizerStereo;
}

AudioModuleEditor* EqualizerAudioModule::createEditor(int type)
{
  return new EqualizerModuleEditor(lock, this);
}

//-------------------------------------------------------------------------------------------------
// automation and state management:

void EqualizerAudioModule::createStaticParameters()
{
  ScopedLock scopedLock(*lock);

  typedef ModulatableParameter Param;
  Param* p;

  std::vector<double> defaultValues;

  addObservedParameter( p = new Param("Bypass", 0.0, 1.0, 0.0, Parameter::BOOLEAN, 1.0) );
  p->setValueChangeCallback(wrappedEqualizerStereo, &EqualizerStereo::setBypass);

  p = new Param("StereoMode", 0.0, 3.0, 0.0, Parameter::STRING, 1.0);
  p->addStringValue("Stereo Linked");
  p->addStringValue("Left/Right");
  p->addStringValue("Mid/Side");
  p->addStringValue("Mono");
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedEqualizerStereo, &EqualizerStereo::setStereoMode);

  p = new Param("GainRange", 0.0, 4.0, 3.0, Parameter::STRING, 1.0);
  p->addStringValue("3 dB");
  //p->addStringValue(juce::String(T("\xF1 3 dB")));  // xF1 == hex for 241 -> plusminus in extended ASCII
  p->addStringValue("6 dB");
  p->addStringValue("12 dB");
  p->addStringValue("24 dB");
  p->addStringValue("48 dB");
  addObservedParameter(p);
  //p->setValueChangeCallback(wrappedEqualizerStereo, &EqualizerStereo::setStereoMode);

  p = new Param("GlobalGain", -48.0, 48.0, 0.0, Parameter::LINEAR_BIPOLAR, 0.1);
  defaultValues.clear();
  defaultValues.push_back(-24.0);
  defaultValues.push_back(-18.0);
  defaultValues.push_back(-12.0);
  defaultValues.push_back( -6.0);
  defaultValues.push_back(  0.0);
  defaultValues.push_back( +6.0);
  defaultValues.push_back(+12.0);
  defaultValues.push_back(+18.0);
  defaultValues.push_back(+24.0);
  p->setDefaultValues(defaultValues);
  addObservedParameter(p);
  p->setValueChangeCallback(wrappedEqualizerStereo, &EqualizerStereo::setGlobalGain);
}

void EqualizerAudioModule::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedLock scopedLock(*lock);

  if( getIndexOfParameter(parameterThatHasChanged) == 1 ) // StereoMode
  {
    if( getNumChannelsToPlot() == 1 )
    {
      selectChannel(0);
      return; // to avoid sending a change-message twice
    }
  }

  markStateAsDirty();
}

XmlElement* EqualizerAudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  ScopedLock scopedLock(*lock);
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if(  wrappedEqualizerStereo->getStereoMode() == rosic::EqualizerStereo::LEFT_RIGHT
    || wrappedEqualizerStereo->getStereoMode() == rosic::EqualizerStereo::MID_SIDE    )
  {
    for(int c = 0; c < 2; c++)
      xmlState->addChildElement(getChannelStateAsXml(c));
  }
  else
  {
    for(int i = 0; i < filterModeParameters[0].size(); i++)
      xmlState->addChildElement(getBandStateAsXml(0, i));
  }

  markStateAsClean();
  return xmlState;
}

void EqualizerAudioModule::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean)
{
  XmlElement convertedState = convertXmlStateIfNecessary(xmlState);

  ScopedLock scopedLock(*lock);
  selectedIndex = -1;
  AudioModule::setStateFromXml(convertedState, stateName, markAsClean);

  removeAllBands();

  int m = 1;
  double f, g, b;

  int c = 0;
  forEachXmlChildElementWithTagName(convertedState, channelChild, "Channel")
  {
    //int i = 0;

    forEachXmlChildElementWithTagName(*channelChild, bandChild, "Band")
    {
      juce::String modeString = bandChild->getStringAttribute("Mode",      "Bypass");
      f = bandChild->getDoubleAttribute("Frequency", 1000.0);
      g = bandChild->getDoubleAttribute("Gain",      0.0);
      b = bandChild->getDoubleAttribute("Bandwidth", 2.0*asinh(1.0/sqrt(2.0))/log(2.0));

      if( modeString == juce::String("Peak/Dip") )
        m = TwoPoleFilter::PEAK;
      else if( modeString == juce::String("Low Shelving") )
        m = TwoPoleFilter::LOW_SHELF;
      else if( modeString == juce::String("High Shelving") )
        m = TwoPoleFilter::HIGH_SHELF;
      else if( modeString == juce::String("Lowpass 6 dB/oct") )
        m = TwoPoleFilter::LOWPASS6;
      else if( modeString == juce::String("Lowpass 12 dB/oct") )
        m = TwoPoleFilter::LOWPASS12;
      else if( modeString == juce::String("Highpass 6 dB/oct") )
        m = TwoPoleFilter::HIGHPASS6;
      else if( modeString == juce::String("Highpass 12 dB/oct") )
        m = TwoPoleFilter::HIGHPASS12;
      else if( modeString == juce::String("Notch 2*6 dB/oct") )
        m = TwoPoleFilter::BANDREJECT;
      // this string-to-int conversion should be somehow automated...

      addBand(c, m, f, g, b, false, true);
    }
    c++;
  }

  sendParameterSetChangeNotification(this);
  markStateAsClean();
}

XmlElement EqualizerAudioModule::convertXmlStateIfNecessary(const XmlElement& xmlState)
{
  ScopedLock scopedLock(*lock);

  //int xmlPatchFormatIndex = xmlState.getIntAttribute("PatchFormat", patchFormatIndex);

  if( xmlState.getChildByName(juce::String("Channel")) == NULL )
  {
    // wrap the band-data into a "Channel" child element:
    XmlElement convertedState(xmlState);
    convertedState.deleteAllChildElements();

    XmlElement* channelState = new XmlElement(juce::String("Channel"));
    forEachXmlChildElementWithTagName(xmlState, bandChild, "Band")
    {
      XmlElement* newBandChild = new XmlElement("Band");
      *newBandChild = *bandChild;
      channelState->addChildElement(newBandChild);
    }

    convertedState.addChildElement(channelState);
    return convertedState;
  }
  else
    return xmlState;
}

XmlElement* EqualizerAudioModule::getChannelStateAsXml(int c)
{
  XmlElement* channelState = new XmlElement(juce::String("Channel"));

  for(int i = 0; i < filterModeParameters[c].size(); i++)
    channelState->addChildElement(getBandStateAsXml(c, i));

  return channelState;
}

XmlElement* EqualizerAudioModule::getBandStateAsXml(int c, int b)
{
  XmlElement* bandState = new XmlElement(juce::String("Band"));

  bandState->setAttribute(filterModeParameters[c][b]->getName(), filterModeParameters[c][b]->getStringValue() );
  bandState->setAttribute(frequencyParameters [c][b]->getName(), frequencyParameters [c][b]->getValue()       );
  bandState->setAttribute(gainParameters      [c][b]->getName(), gainParameters      [c][b]->getValue()       );
  bandState->setAttribute(bandwidthParameters [c][b]->getName(), bandwidthParameters [c][b]->getValue()       );

  return bandState;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void EqualizerAudioModule::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  wrappedEqualizerStereo->setSampleRate((float)newSampleRate);
}

void EqualizerAudioModule::selectChannel(int channelToSelect)
{
  ScopedLock scopedLock(*lock);
  jassert(channelToSelect >= 0 && channelToSelect <= 1); // there are only 2 channels
  selectedChannel = channelToSelect;
  deSelectBand();
}

int EqualizerAudioModule::selectBand(int channel, int indexToSelect)
{
  ScopedLock scopedLock(*lock);
  selectedChannel = channel;
  if( indexToSelect < 0 || indexToSelect >= wrappedEqualizerStereo->getNumBands(channel) )
    selectedIndex = -1;
  else
    selectedIndex = indexToSelect;
  sendParameterSetChangeNotification(this);
  return selectedIndex;
}

void EqualizerAudioModule::deSelectBand()
{
  ScopedLock scopedLock(*lock);
  selectedIndex = -1;
  sendParameterSetChangeNotification(this);
}

void EqualizerAudioModule::addBand(int channel, int mode, double frequency, double gain, double bandwidth,
                                   bool selectNewBand, bool suppressNotification)
{
  ScopedLock scopedLock(*lock);
  int newIndex = wrappedEqualizerStereo->addBand(channel, mode, frequency, gain, bandwidth);

  std::vector<double> defaultValues;
  typedef Parameter Param;
  //typedef ModulatableParameter Param; // experimental, does not yet work right
  Param* p;

  p = new Param("Mode", 0.0, 8.0, 0.0, Parameter::STRING, 1.0);
  p->setMutexToUse(lock);
  p->addStringValue("Bypass");
  p->addStringValue("Peak/Dip");
  p->addStringValue("Low Shelving");
  p->addStringValue("High Shelving");
  p->addStringValue("Lowpass 6 dB/oct");
  p->addStringValue("Lowpass 12 dB/oct");
  p->addStringValue("Highpass 6 dB/oct");
  p->addStringValue("Highpass 12 dB/oct");
  p->addStringValue("Notch 2*6 dB/oct");
  p->setValue(mode, false, false);
  p->registerParameterObserver(this);
  setupManagers(p);
  filterModeParameters[channel].add(p);

  p = new Param("Frequency", 20.0, 20000.0, 1000.0, Parameter::EXPONENTIAL, 0.0);
  p->setMutexToUse(lock);
  defaultValues.clear();
  defaultValues.push_back(   31.25);
  defaultValues.push_back(   62.5);
  defaultValues.push_back(  125.0);
  defaultValues.push_back(  250.0);
  defaultValues.push_back(  500.0);
  defaultValues.push_back( 1000.0);
  defaultValues.push_back( 2000.0);
  defaultValues.push_back( 4000.0);
  defaultValues.push_back( 8000.0);
  defaultValues.push_back(16000.0);
  p->setDefaultValues(defaultValues);
  p->setValue(frequency, false, false);
  p->registerParameterObserver(this);
  setupManagers(p);
  frequencyParameters[channel].add(p);

  p = new Param("Gain", -48.0, 48.0, 0.0, Parameter::LINEAR_BIPOLAR, 0.1);
  p->setMutexToUse(lock);
  defaultValues.clear();
  defaultValues.push_back(-24.0);
  defaultValues.push_back(-18.0);
  defaultValues.push_back(-12.0);
  defaultValues.push_back( -6.0);
  defaultValues.push_back( RAPT::rsAmpToDb(sqrt(0.5))); // -3.01 dB - standard definition for cutoff frequencies
  defaultValues.push_back(  0.0);
  defaultValues.push_back( +6.0);
  defaultValues.push_back(+12.0);
  defaultValues.push_back(+18.0);
  defaultValues.push_back(+24.0);
  p->setDefaultValues(defaultValues);
  p->setValue(gain, false, false);
  p->registerParameterObserver(this);
  setupManagers(p);
  gainParameters[channel].add(p);

  double defaultBandwidth = 2.0*asinh(1.0/sqrt(2.0))/log(2.0);  // ca. 1.9, Q = sqrt(1/2), maximum steepness without overshoot for shelves
  p = new Param("Bandwidth", 0.25, 6.0, defaultBandwidth, Parameter::EXPONENTIAL, 0.01);
  p->setMutexToUse(lock);
  defaultValues.clear();
  defaultValues.push_back(1.0/3.0);
  defaultValues.push_back(0.5);
  defaultValues.push_back(1.0);
  defaultValues.push_back(defaultBandwidth);
  defaultValues.push_back(2.0);
  defaultValues.push_back(3.0);
  defaultValues.push_back(4.0);
  p->setDefaultValues(defaultValues);
  p->setValue(bandwidth, false, false);
  p->registerParameterObserver(this);
  setupManagers(p);
  bandwidthParameters[channel].add(p);

  assignCallbacksForDynamicParameters();

  if( selectNewBand == true )
  {
    selectedChannel = channel;
    selectedIndex   = newIndex;
  }

  if( suppressNotification != true )
    sendParameterSetChangeNotification(this);
}

bool EqualizerAudioModule::removeBand(int channel, int indexToRemove)
{
  ScopedLock scopedLock(*lock);
  if( indexToRemove < 0 || indexToRemove >= wrappedEqualizerStereo->getNumBands(channel) )
    return false;

  filterModeParameters[channel][indexToRemove]->deRegisterParameterObserver(this);
  frequencyParameters[channel][indexToRemove]->deRegisterParameterObserver(this);
  gainParameters[channel][indexToRemove]->deRegisterParameterObserver(this);
  bandwidthParameters[channel][indexToRemove]->deRegisterParameterObserver(this);

  filterModeParameters[channel].remove(indexToRemove, true);
  frequencyParameters[channel].remove( indexToRemove, true);
  gainParameters[channel].remove(      indexToRemove, true);
  bandwidthParameters[channel].remove( indexToRemove, true);

  // we can now safely remove the band from the wrapped Equalizer object because the parameters wich have held the callback pointers have
  // just been deleted:
  wrappedEqualizerStereo->removeBand(channel, indexToRemove);

  assignCallbacksForDynamicParameters();
  selectedIndex = -1;
  sendParameterSetChangeNotification(this);
  return true;
}

void EqualizerAudioModule::removeAllBands()
{
  ScopedLock scopedLock(*lock);

  for(int c=0; c<2; c++)
  {
    filterModeParameters[c].clear(true);
    frequencyParameters[c].clear(true);
    gainParameters[c].clear(true);
    bandwidthParameters[c].clear(true);
  }

  wrappedEqualizerStereo->removeAllBands();
  selectedIndex = -1;
  sendParameterSetChangeNotification(this);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

int EqualizerAudioModule::getNumChannelsToPlot()
{
  ScopedLock scopedLock(*lock);
  if(  wrappedEqualizerStereo->getStereoMode() == rosic::EqualizerStereo::STEREO_LINKED
    || wrappedEqualizerStereo->getStereoMode() == rosic::EqualizerStereo::MONO )
  {
    return 1;
  }
  else
    return 2;
}

int EqualizerAudioModule::getNumBands(int channel)
{
  ScopedLock scopedLock(*lock);
  int result = wrappedEqualizerStereo->getNumBands(channel);
  return result;
}

void EqualizerAudioModule::getMagnitudeResponse(int channel, double *frequencies, double *magnitudes, int numBins)
{
  ScopedLock scopedLock(*lock);
  wrappedEqualizerStereo->getMagnitudeResponse(channel, frequencies, magnitudes, numBins);
}

double EqualizerAudioModule::getDecibelsAt(int channel, double frequency)
{
  return 0; // preliminary
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void EqualizerAudioModule::reset()
{
  ScopedLock scopedLock(*lock);
  wrappedEqualizerStereo->reset();
}

void EqualizerAudioModule::assignCallbacksForDynamicParameters()
{
  ScopedLock scopedLock(*lock);
  int numChannels = 2;
  for(int c=0; c<numChannels; c++)
  {
    int numBands = filterModeParameters[c].size();
    for(int i=0; i<numBands; i++)
    {
      filterModeParameters[c][i]->clearValueChangeCallbacks();
      frequencyParameters[c][i]->clearValueChangeCallbacks();
      gainParameters[c][i]->clearValueChangeCallbacks();
      bandwidthParameters[c][i]->clearValueChangeCallbacks();
      // why do we need these call to clear - check and delete if not needed

      filterModeParameters[c][i]->setValueChangeCallback(&wrappedEqualizerStereo->equalizers[c].bands[i], &TwoPoleFilter::setMode);
      frequencyParameters[c][i]->setValueChangeCallback( &wrappedEqualizerStereo->equalizers[c].bands[i], &TwoPoleFilter::setFrequency);
      gainParameters[c][i]->setValueChangeCallback(      &wrappedEqualizerStereo->equalizers[c].bands[i], &TwoPoleFilter::setGain);
      bandwidthParameters[c][i]->setValueChangeCallback( &wrappedEqualizerStereo->equalizers[c].bands[i], &TwoPoleFilter::setBandwidth);
    }
  }
}

//=================================================================================================

EqualizerPlotEditor::EqualizerPlotEditor(CriticalSection *newPlugInLock, EqualizerAudioModule* newEqualizerModuleToEdit)
  : rsSpectrumPlot(juce::String("EqualizerEditor"))
{
  setDescription("Left: insert or grab band-handle, right: remove band");

  ParameterObserver::setIsGuiElement(true);

  plugInLock            = newPlugInLock;
  equalizerModuleToEdit = NULL;  // will be assigned later via call to setEqualizerModuleToEdit

  // set up the plot range:
  setAutoReRendering(false);
  setMaximumRange(15.625, 32000.0, -48.0, 48.0);
  setCurrentRange(15.625, 32000.0, -24.0, 24.0);
  setHorizontalCoarseGrid(6.0, true);
  setHorizontalFineGrid(  1.0, false);
  setVerticalCoarseGridVisible( true);
  setVerticalFineGridVisible(   false);
  rsPlot::setAxisValuesPositionX(rsPlotSettings::ABOVE_AXIS);
  rsPlot::setAxisValuesPositionY(rsPlotSettings::RIGHT_TO_AXIS);
  setAutoReRendering(true);

  currentMouseCursor = MouseCursor(MouseCursor::NormalCursor);
  setMouseCursor(currentMouseCursor);

  filterModeParameter = NULL;
  frequencyParameter  = NULL;
  gainParameter       = NULL;
  bandwidthParameter  = NULL;
  globalGainParameter = NULL;

  // this stuff will be (re-) assigned in resized():
  numBins       = 0;
  frequencies   = NULL;
  magnitudes1   = NULL;
  magnitudes2   = NULL;
  magnitudes[0] = magnitudes1;
  magnitudes[1] = magnitudes2;

  currentlyDraggedHandle = NONE;

  // activate automation for this ParameterObserver:
  ParameterObserver::setLocalAutomationSwitch(true);

  setEqualizerModuleToEdit(newEqualizerModuleToEdit);
}

EqualizerPlotEditor::~EqualizerPlotEditor(void)
{
  setEqualizerModuleToEdit(NULL); // to remove ourselves as ChangeListener
  deleteAndZero(frequencies);
  deleteAndZero(magnitudes1);
  deleteAndZero(magnitudes2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void EqualizerPlotEditor::setEqualizerModuleToEdit(EqualizerAudioModule* newEqualizerModuleToEdit)
{
  ScopedPointerLock spl(plugInLock);

  if( newEqualizerModuleToEdit == equalizerModuleToEdit )
    return;

  if( equalizerModuleToEdit != NULL )
    equalizerModuleToEdit->deRegisterParameterSetObserver(this);

  equalizerModuleToEdit = newEqualizerModuleToEdit;

  if( globalGainParameter != NULL )
    globalGainParameter->deRegisterParameterObserver(this);
  globalGainParameter = NULL;

  if( equalizerModuleToEdit != NULL )
  {
    equalizerModuleToEdit->registerParameterSetObserver(this);
    globalGainParameter = equalizerModuleToEdit->getParameterByName("GlobalGain");
    globalGainParameter->registerParameterObserver(this);
  }

  assignParametersToSelectedBand();
  updatePlot();
}

void EqualizerPlotEditor::unAssignParameters()
{
  ScopedPointerLock spl(plugInLock);

  if( filterModeParameter != NULL )
    filterModeParameter->deRegisterParameterObserver(this);
  filterModeParameter = NULL;

  if( frequencyParameter != NULL )
    frequencyParameter->deRegisterParameterObserver(this);
  frequencyParameter = NULL;

  if( gainParameter != NULL )
    gainParameter->deRegisterParameterObserver(this);
  gainParameter = NULL;

  if( bandwidthParameter != NULL )
    bandwidthParameter->deRegisterParameterObserver(this);
  bandwidthParameter = NULL;
}

void EqualizerPlotEditor::assignParametersToSelectedBand()
{
  ScopedPointerLock spl(plugInLock);

  unAssignParameters();

  if( equalizerModuleToEdit == NULL )
    return;

  int channel = equalizerModuleToEdit->selectedChannel;
  int index   = equalizerModuleToEdit->selectedIndex;

  if( index > -1 )
  {
    filterModeParameter = equalizerModuleToEdit->filterModeParameters[channel][index];
    filterModeParameter->registerParameterObserver(this);

    frequencyParameter = equalizerModuleToEdit->frequencyParameters[channel][index];
    frequencyParameter->registerParameterObserver(this);

    gainParameter = equalizerModuleToEdit->gainParameters[channel][index];
    gainParameter->registerParameterObserver(this);

    bandwidthParameter = equalizerModuleToEdit->bandwidthParameters[channel][index];
    bandwidthParameter->registerParameterObserver(this);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

int EqualizerPlotEditor::getBandIndexAtPixelPosition(int x, int y)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return -1;

  double globalGain = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();
  double dotRadius  = 4.0;
  double xd         = (double) x;
  double yd         = (double) y;
  for(int i=0; i<equalizerModuleToEdit->getNumBands(equalizerModuleToEdit->selectedChannel); i++)
  {
    double xi = equalizerModuleToEdit->wrappedEqualizerStereo->getBandFrequency(equalizerModuleToEdit->selectedChannel, i);
    double yi = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(equalizerModuleToEdit->selectedChannel, i);
    yi       += globalGain;
    toPixelCoordinates(xi, yi);
    double d = sqrt( (xi-xd)*(xi-xd) + (yi-yd)*(yi-yd) );
    if( d <= dotRadius )
      return i;

    if( i == equalizerModuleToEdit->selectedIndex )
    {
      // check if x,y is over the left bandwidth handle of the selected band:
      double xt = equalizerModuleToEdit->wrappedEqualizerStereo->getLowerBandedgeFrequency(
        equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
      double yt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(
        equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
      yt       += globalGain;
      toPixelCoordinates(xt, yt);
      if( RAPT::rsEuclideanDistance(xt, yt, xd, yd) < 3.0 )
        return i;

      // check if x,y is over the right bandwidth handle of the selected band:
      xt = equalizerModuleToEdit->wrappedEqualizerStereo->getUpperBandedgeFrequency(
        equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
      yt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(
        equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
      yt       += globalGain;
      toPixelCoordinates(xt, yt);
      if( RAPT::rsEuclideanDistance(xt, yt, xd, yd) < 3.0 )
        return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

/*
void EqualizerPlotEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
ScopedPointerLock spl(plugInLock);
assignParametersToSelectedBand();
updatePlot();
}
*/

void EqualizerPlotEditor::parameterSetChanged(ParameterSetHolder* parameterSetHolderThatHasChanged)
{
  ScopedPointerLock spl(plugInLock);
  assignParametersToSelectedBand();
  updatePlot();
}

void EqualizerPlotEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  ScopedPointerLock spl(plugInLock);
  updatePlot();
}

void EqualizerPlotEditor::parameterWillBeDeleted(Parameter* p)
{
  ScopedPointerLock spl(plugInLock);

  p->deRegisterParameterObserver(this);

  if( p == filterModeParameter )
    filterModeParameter = NULL;
  else if( p == frequencyParameter )
    frequencyParameter = NULL;
  else if( p == gainParameter )
    gainParameter = NULL;
  else if( p == bandwidthParameter )
    bandwidthParameter = NULL;

  updatePlot();
}

void EqualizerPlotEditor::mouseMove(const MouseEvent &e)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return;

  int index = getBandIndexAtPixelPosition(e.x, e.y);
  if( index != -1 )
  {
    //int    m    = equalizerModuleToEdit->wrappedEqualizerStereo->getBandMode(     equalizerModuleToEdit->selectedChannel, index);
    double f    = equalizerModuleToEdit->wrappedEqualizerStereo->getBandFrequency(equalizerModuleToEdit->selectedChannel, index);
    double g    = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(     equalizerModuleToEdit->selectedChannel, index);
    double b    = equalizerModuleToEdit->wrappedEqualizerStereo->getBandBandwidth(equalizerModuleToEdit->selectedChannel, index);

    juce::String mStr;
    int mode = equalizerModuleToEdit->wrappedEqualizerStereo->getBandMode(equalizerModuleToEdit->selectedChannel, index);
    switch( mode )
    {
    case rosic::TwoPoleFilter::PEAK:       mStr = juce::String("Peak/Dip");           break;
    case rosic::TwoPoleFilter::LOW_SHELF:  mStr = juce::String("Low Shelving");       break;
    case rosic::TwoPoleFilter::HIGH_SHELF: mStr = juce::String("High Shelving");      break;
    case rosic::TwoPoleFilter::LOWPASS6:   mStr = juce::String("Lowpass 6 dB/oct");   break;
    case rosic::TwoPoleFilter::LOWPASS12:  mStr = juce::String("Lowpass 12 dB/oct");  break;
    case rosic::TwoPoleFilter::HIGHPASS6:  mStr = juce::String("Highpass 6 dB/oct");  break;
    case rosic::TwoPoleFilter::HIGHPASS12: mStr = juce::String("Highpass 12 dB/oct"); break;
    case rosic::TwoPoleFilter::BANDREJECT: mStr = juce::String("Notch 2*6 dB/oct");   break;
    }
    juce::String fStr = juce::String("Frequency: ") + hertzToStringWithUnitTotal5(f) + juce::String(", ");
    juce::String gStr = juce::String("Gain: ")      + decibelsToStringWithUnit2(g)   + juce::String(", ");
    juce::String bStr = juce::String("Bandwidth: ") + octavesToStringWithUnit2(b);

    juce::String description = mStr + juce::String(", ") + fStr;
    if( equalizerModuleToEdit->wrappedEqualizerStereo->doesModeSupportGain(equalizerModuleToEdit->selectedChannel, index) )
      description += gStr;
    if( equalizerModuleToEdit->wrappedEqualizerStereo->doesModeSupportBandwidth(equalizerModuleToEdit->selectedChannel, index) )
      description += bStr;

    setDescription(description);
  }
  else
    setDescription(juce::String("Left: insert or grab band-handle, right: remove band"));

  int dragHandle = getDragHandleAt(e.x, e.y);
  if( dragHandle == NONE )
    setMouseCursor(MouseCursor::NormalCursor);
  else if( dragHandle == FREQUENCY_AND_GAIN )
    setMouseCursor(MouseCursor::PointingHandCursor);
  else if( dragHandle == BANDWIDTH_AND_GAIN_LEFT || dragHandle == BANDWIDTH_AND_GAIN_RIGHT || dragHandle == FREQUENCY )
    setMouseCursor(MouseCursor::LeftRightResizeCursor);
  else if( dragHandle == GLOBALGAIN_LINE || dragHandle == GAIN )
    setMouseCursor(MouseCursor::UpDownResizeCursor);

}

void EqualizerPlotEditor::mouseDown(const MouseEvent &e)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return;

  currentlyDraggedHandle = getDragHandleAt(e.x, e.y);

  if( equalizerModuleToEdit->selectedIndex != -1 )
  {
    // check first, if one of the left/right or up/down handles is grabbed:
    if(  currentlyDraggedHandle == BANDWIDTH_AND_GAIN_LEFT  || currentlyDraggedHandle == BANDWIDTH_AND_GAIN_RIGHT
      || currentlyDraggedHandle == FREQUENCY || currentlyDraggedHandle == GAIN   )
    {
      //sendChangeMessage();
      return;
    }
  }

  int tmpIndex = getBandIndexAtPixelPosition(e.x, e.y);
  if( tmpIndex == -1 )
  {
    if( e.mods.isLeftButtonDown() )
    {
      if( currentlyDraggedHandle != GLOBALGAIN_LINE )
      {
        // create a new band and mark it selected:
        double f = e.x;
        double g = e.y;
        xyToFrequencyAndGain(f, g);
        equalizerModuleToEdit->addBand(equalizerModuleToEdit->selectedChannel, rosic::TwoPoleFilter::PEAK, f, g);
        equalizerModuleToEdit->selectBand(equalizerModuleToEdit->selectedChannel,
          equalizerModuleToEdit->getNumBands(equalizerModuleToEdit->selectedChannel)-1);
        currentlyDraggedHandle = FREQUENCY_AND_GAIN;
      }
    }
    // maybe de-select on mouse-click
    /*
    else if( e.mods.isRightButtonDown() )
    {
    currentlyDraggedHandle = NONE;
    openRightClickPopupMenu(e.x, e.y);
    return;
    }
    */
  }
  else
  {
    if( e.mods.isRightButtonDown() )
    {
      /*
      // this stuff should not be necesarry, but let's try for debugging:
      if( frequencyParameter != NULL )
      frequencyParameter->deRegisterParameterObserver(this);
      frequencyParameter = NULL;
      if( gainParameter != NULL )
      gainParameter->deRegisterParameterObserver(this);
      gainParameter = NULL;
      // end debug
      */

      equalizerModuleToEdit->removeBand(equalizerModuleToEdit->selectedChannel, tmpIndex);
      currentlyDraggedHandle = NONE;
    }
    else
    {
      equalizerModuleToEdit->selectBand(equalizerModuleToEdit->selectedChannel, tmpIndex);
      currentlyDraggedHandle = FREQUENCY_AND_GAIN;
    }
  }

  updatePlot();

  //sendChangeMessage();
}

void EqualizerPlotEditor::mouseDrag(const juce::MouseEvent &e)
{
  if( e.mods.isRightButtonDown() || e.mouseWasClicked() )
    return;   // ignore right-drags because the band was just removed, alos ignore just-clicks

  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return;

  if( frequencyParameter == NULL || gainParameter == NULL || bandwidthParameter == NULL
    || globalGainParameter == NULL )
  {
    if( currentlyDraggedHandle != GLOBALGAIN_LINE )
      return;
  }

  double globalGain = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();

  // get the position of the event in components coordinates:
  double x = e.getMouseDownX() + e.getDistanceFromDragStartX();
  double y = e.getMouseDownY() + e.getDistanceFromDragStartY();

  x = RAPT::rsClip(x, 0.0, (double) getWidth());
  y = RAPT::rsClip(y, 0.0, (double) getHeight());
  fromPixelCoordinates(x, y);

  if( currentlyDraggedHandle == FREQUENCY_AND_GAIN )
  {
    frequencyParameter->setValue(x, true, true);
    gainParameter->setValue(y-globalGain, true, true);
  }
  else if( currentlyDraggedHandle == BANDWIDTH_AND_GAIN_LEFT )
  {
    double bw = rosic::TwoPoleFilter::lowerBandedgeFrequencyToBandwdith(x,
      frequencyParameter->getValue());
    bandwidthParameter->setValue(bw, true, true);
    gainParameter->setValue(y-globalGain, true, true);
  }
  else if( currentlyDraggedHandle == BANDWIDTH_AND_GAIN_RIGHT )
  {
    double bw = rosic::TwoPoleFilter::upperBandedgeFrequencyToBandwdith(x,
      frequencyParameter->getValue());
    bandwidthParameter->setValue(bw, true, true);
    gainParameter->setValue(y-globalGain, true, true);
  }
  else if( currentlyDraggedHandle == FREQUENCY )
  {
    frequencyParameter->setValue(x, true, true);
  }
  else if( currentlyDraggedHandle == GAIN )
  {
    gainParameter->setValue(y-globalGain, true, true);
  }
  else if( currentlyDraggedHandle == GLOBALGAIN_LINE )
  {
    globalGainParameter->setValue(y, true, true);
  }

  mouseMove(e);
  updatePlot();
}

void EqualizerPlotEditor::mouseUp(const juce::MouseEvent &e)
{
  currentlyDraggedHandle = NONE;
}

void EqualizerPlotEditor::mouseWheelMove(const MouseEvent& e, const MouseWheelDetails& wheel)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return;

  int channel = equalizerModuleToEdit->selectedChannel;
  int index   = getBandIndexAtPixelPosition(e.x, e.y);
  if( index != -1 )
  {
    Parameter *p = equalizerModuleToEdit->bandwidthParameters[channel][index];
    double b     = p->getValue();
    b  = RAPT::rsExpToLin(b, 0.25, 4.0, 0.0, 1.0);
    b += 0.0625*wheel.deltaY;
    b  = RAPT::rsLinToExp(b, 0.0, 1.0, 0.25, 4.0);
    p->setValue(b, true, true);
    updatePlot();
  }
}

int EqualizerPlotEditor::getDragHandleAt(int x, int y)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return NONE;

  double globalGain = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();
  double xd         = (double) x;
  double yd         = (double) y;
  double xt, yt;             // target coordinates for matching

  if( equalizerModuleToEdit->selectedIndex != -1 )
  {
    // check if x,y is over the left bandwidth handle of the selected band:
    xt = equalizerModuleToEdit->wrappedEqualizerStereo->getLowerBandedgeFrequency(
      equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
    yt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(
      equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex) + globalGain;
    toPixelCoordinates(xt, yt);
    if( RAPT::rsEuclideanDistance(xt, yt, xd, yd) < 3.0 )
      return BANDWIDTH_AND_GAIN_LEFT;

    // check if x,y is over the right bandwidth handle of the selected band:
    xt = equalizerModuleToEdit->wrappedEqualizerStereo->getUpperBandedgeFrequency(
      equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
    yt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(
      equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex) + globalGain;
    toPixelCoordinates(xt, yt);
    if( RAPT::rsEuclideanDistance(xt, yt, xd, yd) < 3.0 )
      return BANDWIDTH_AND_GAIN_RIGHT;
  }

  // check if x,y is over the freq/gain handle of some band:
  for(int i=0; i<equalizerModuleToEdit->getNumBands(equalizerModuleToEdit->selectedChannel); i++)
  {
    xt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandFrequency(
      equalizerModuleToEdit->selectedChannel, i);
    yt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(
      equalizerModuleToEdit->selectedChannel, i) + globalGain;
    toPixelCoordinates(xt, yt);
    if( RAPT::rsEuclideanDistance(xt, yt, xd, yd) < 4.0 )
      return FREQUENCY_AND_GAIN;
  }

  if( equalizerModuleToEdit->selectedIndex != -1 )
  {
    // check if y is over the horizontal gain line for the selected band:
    xt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandFrequency(
      equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex);
    yt = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(
      equalizerModuleToEdit->selectedChannel, equalizerModuleToEdit->selectedIndex) + globalGain;
    toPixelCoordinates(xt, yt);
    if( fabs(yt-yd) <= 2.0 )
      return GAIN;

    // check if x is over the vertical frequency line for the selected band:
    if( fabs(xt-xd) <= 2.0 )
      return FREQUENCY;
  }

  // check if y is over the global gain line:
  xt = 1000.0;                            // dummy;
  yt = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();
  toPixelCoordinates(xt, yt);
  if( fabs(y-yt) < 4.0 )
    return GLOBALGAIN_LINE;

  return NONE;
}

//-------------------------------------------------------------------------------------------------
// drawing:

void EqualizerPlotEditor::resized()
{
  rsSpectrumPlot::resized();

  // (re) allocate and fill the arrays for the magnitude plot
  numBins = getWidth();
  if( frequencies == NULL )
    delete[] frequencies;
  if( magnitudes1 == NULL )
    delete[] magnitudes1;
  if( magnitudes2 == NULL )
    delete[] magnitudes2;

  frequencies   = new double[numBins];
  magnitudes1   = new double[numBins];
  magnitudes2   = new double[numBins];
  magnitudes[0] = magnitudes1;
  magnitudes[1] = magnitudes2;


  getDisplayedFrequencies(frequencies, numBins);
  updatePlot();
}

void EqualizerPlotEditor::updatePlot()
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
  {
    fillWithZeros(magnitudes1, numBins);
    fillWithZeros(magnitudes2, numBins);
  }
  else
  {
    equalizerModuleToEdit->getMagnitudeResponse(0, frequencies, magnitudes1, numBins);
    equalizerModuleToEdit->getMagnitudeResponse(1, frequencies, magnitudes2, numBins);
  }
  setSpectra(numBins, 2, frequencies, magnitudes);
}

void EqualizerPlotEditor::plotCurveFamily(Graphics &g, juce::Image* targetImage,
  XmlElement *targetSVG)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return;

  rsDataPlot::setNumCurves( equalizerModuleToEdit->getNumChannelsToPlot() );
  rsDataPlot::plotCurveFamily(g, targetImage, targetSVG);

  double x, y;

  // draw the horizontal line for the reference gain:

  Colour curveColour = plotColourScheme.getCurveColourUniform(0);
  g.setColour(curveColour);
  y = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();
  x = 1000.0;
  toPixelCoordinates(x, y);
  const float dashLengths[2] = {4.f, 6.f};
  g.drawDashedLine(Line<float>(0.f, (float) y, (float) getWidth(), (float) y), dashLengths, 2, 2.f);


  int channel = equalizerModuleToEdit->selectedChannel;
  if( equalizerModuleToEdit->getNumChannelsToPlot() == 1 )
    channel = 0; // ignore selection in thsi case

  curveColour = getCurveColour(channel);
  g.setColour(curveColour);

  float  dotRadius = 4.f;
  for(int i=0; i<equalizerModuleToEdit->getNumBands(channel); i++)
  {
    double globalGain = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();
    x = equalizerModuleToEdit->wrappedEqualizerStereo->getBandFrequency(channel, i);
    y = equalizerModuleToEdit->wrappedEqualizerStereo->getBandGain(channel, i) + globalGain;
    toPixelCoordinates(x, y);

    g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius),
      (float) (2*dotRadius), (float) (2*dotRadius) );

    if( i == equalizerModuleToEdit->selectedIndex )
    {
      g.setColour(curveColour.withMultipliedAlpha(0.5f));
      //g.setColour(Colour(0xff404090).withMultipliedAlpha(0.5f));
      g.drawLine((float) x,        0.f, (float)          x, (float) getHeight(), 1.f);
      g.drawLine(      0.f,  (float) y, (float) getWidth(), (float) y          , 1.f);

      g.fillEllipse((float) (x-dotRadius-2), (float) (y-dotRadius-2),
        (float) (2*dotRadius+4), (float) (2*dotRadius+4) );

      g.setColour(curveColour.darker(0.5f));
      //g.setColour(Colour(0xff404090));

      g.fillEllipse((float) (x-dotRadius), (float) (y-dotRadius),
        (float) (2*dotRadius), (float) (2*dotRadius) );

      if( equalizerModuleToEdit->wrappedEqualizerStereo->doesModeSupportBandwidth(channel, i) )
      {
        double x1 = equalizerModuleToEdit->wrappedEqualizerStereo
          ->getLowerBandedgeFrequency(channel, i);
        double y1 = equalizerModuleToEdit->wrappedEqualizerStereo
          ->getBandGain(channel, i) + globalGain;
        toPixelCoordinates(x1, y1);

        double x2 = equalizerModuleToEdit->wrappedEqualizerStereo
          ->getUpperBandedgeFrequency(channel, i);
        double y2 = equalizerModuleToEdit->wrappedEqualizerStereo
          ->getBandGain(channel, i) + globalGain;
        toPixelCoordinates(x2, y2);

        g.drawLine((float) x1, (float) y1, (float) x2, (float) y1, 2.0f);
        g.drawLine((float) x1, (float) (y1-5.0), (float) x1, (float) (y1+5.0), 2.0f);
        g.drawLine((float) x2, (float) (y1-5.0), (float) x2, (float) (y1+5.0), 2.0f);
      }

      g.setColour(curveColour);
    }
  }
}

/*
void EqualizerPlotEditor::createRightClickPopupMenu(PopupMenu*& menu)
{
menu = new PopupMenu();
menu->addItem(1, juce::String(T("Insert node: Peak/Dip"))           );
menu->addItem(2, juce::String(T("Insert node: Low Shelving"))       );
menu->addItem(3, juce::String(T("Insert node: High Shelving"))      );
menu->addItem(4, juce::String(T("Insert node: Lowpass 6 dB/oct"))   );
menu->addItem(5, juce::String(T("Insert node: Lowpass 12 dB/oct"))  );
menu->addItem(6, juce::String(T("Insert node: Highpass 6 dB/oct"))  );
menu->addItem(7, juce::String(T("Insert node: Highpass 12 dB/oct")) );
menu->addItem(8, juce::String(T("Insert node: Notch 2*6 dB/oct"))   );
}

void EqualizerPlotEditor::handleRightClickPopupMenuResult(int result, int x, int y)
{
ScopedPointerLock spl(plugInLock);
if( equalizerModuleToEdit == NULL )
return;

double g = (double) y;
double f = (double) x;
xyToFrequencyAndGain(f, g);
switch( result )
{
case 1: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::PEAK,       f, g); break;
case 2: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::LOW_SHELF,  f, g); break;
case 3: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::HIGH_SHELF, f, g); break;
case 4: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::LOWPASS6,   f, g); break;
case 5: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::LOWPASS12,  f, g); break;
case 6: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::HIGHPASS6,  f, g); break;
case 7: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::HIGHPASS12, f, g); break;
case 8: equalizerModuleToEdit->addBand(rosic::TwoPoleFilter::BANDREJECT, f, g); break;
}

updatePlot();
}

void EqualizerPlotEditor::openRightClickPopupMenu(int x, int y)
{
PopupMenu* menu = NULL;
createRightClickPopupMenu(menu);
int result = menu->show();
if( menu != NULL )
delete menu;
handleRightClickPopupMenuResult(result, x, y);
}
*/

void EqualizerPlotEditor::xyToFrequencyAndGain(double &x, double &y)
{
  ScopedPointerLock spl(plugInLock);
  if( equalizerModuleToEdit == NULL )
    return;

  double globalGain = equalizerModuleToEdit->wrappedEqualizerStereo->getGlobalGain();
  double f          = x;
  double g          = y;
  fromPixelCoordinates(f, g);
  g                -= globalGain;
  x                 = f;
  y                 = g;
}

//=================================================================================================

rsEqualizerPlotEditor::rsEqualizerPlotEditor(EqualizerAudioModule* eqModule) 
  : equalizerModule(eqModule)
{
  freqRespPlot = new rsFunctionPlot;
  freqRespPlot->setupForDecibelsAgainstLogFrequency(15.625, 32000.0, -24.0, 24.0, 6);
  freqRespPlot->addFunction([this](double f)->double { return equalizerModule->getDecibelsAt(0, f); });
  freqRespPlot->addFunction([this](double f)->double { return equalizerModule->getDecibelsAt(1, f); });
  addPlot(freqRespPlot);

  nodeEditor = new rsNodeEditor;
  addWidget(nodeEditor);
}

void rsEqualizerPlotEditor::nodeWasAdded(rsNodeEditor* editor, int nodeIndex)
{

}

void rsEqualizerPlotEditor::nodeWillBeRemoved(rsNodeEditor* editor, int nodeIndex)
{

}

void rsEqualizerPlotEditor::nodeWasMoved(rsNodeEditor* editor, int nodeIndex)
{

}

//=================================================================================================

// construction/destruction:

EqualizerModuleEditor::EqualizerModuleEditor(CriticalSection *newPlugInLock,
  EqualizerAudioModule* newEqualizerAudioModule)
  : AudioModuleEditor(newEqualizerAudioModule)
{
  ScopedPointerLock spl(lock);
  equalizerModule = newEqualizerAudioModule;
  createWidgets();
  layout              = SLIDERS_RIGHT;
  useShortSliderNames = false;
  useSmallComboBox    = false;
  stateWidgetSet->addChangeListener(this);
  setEqualizerModuleToEdit(newEqualizerAudioModule);
  setSize(500, 234); // needs to be called before updateWidgetsAccordingToState or we may call resize during paint or soemthing (hits a jassert)
  updateWidgetsAccordingToState();
}

EqualizerModuleEditor::~EqualizerModuleEditor()
{
  if( plotEditor->equalizerModuleToEdit != NULL )
  {
    plotEditor->equalizerModuleToEdit->deRegisterParameterSetObserver(this);
    //plotEditor->equalizerModuleToEdit->removeStateWatcher(stateWidgetSet);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void EqualizerModuleEditor::setEqualizerModuleToEdit(EqualizerAudioModule* newEqualizerModuleToEdit)
{
  ScopedPointerLock spl(lock);
  //if( newEqualizerModuleToEdit == plotEditor->equalizerModuleToEdit )
  //  return;

  // unassign widgets from parameters of the old equalizer:
  bypassButton->assignParameter(      NULL);
  stereoModeComboBox->assignParameter(NULL);
  gainRangeComboBox->assignParameter( NULL);
  globalGainSlider->assignParameter(  NULL);
  filterModeComboBox->assignParameter(NULL);
  frequencySlider->assignParameter(   NULL);
  gainSlider->assignParameter(        NULL);
  bandwidthSlider->assignParameter(   NULL);

  if( plotEditor->equalizerModuleToEdit != NULL )
  {
    plotEditor->equalizerModuleToEdit->deRegisterParameterSetObserver(this);
    plotEditor->equalizerModuleToEdit->removeStateWatcher(stateWidgetSet);
  }

  plotEditor->setEqualizerModuleToEdit(newEqualizerModuleToEdit);

  if( plotEditor->equalizerModuleToEdit != NULL )
  {
    plotEditor->equalizerModuleToEdit->registerParameterSetObserver(this);
    plotEditor->equalizerModuleToEdit->addStateWatcher(stateWidgetSet);

    // assign widgets to parameters of the new equalizer (only static parameters - dynamic parameters will be assigned in updateSliders):
    bypassButton->assignParameter(
      plotEditor->equalizerModuleToEdit->getParameterByName("Bypass") );
    stereoModeComboBox->assignParameter(
      plotEditor->equalizerModuleToEdit->getParameterByName("StereoMode") );
    gainRangeComboBox->assignParameter(
      plotEditor->equalizerModuleToEdit->getParameterByName("GainRange") );
    globalGainSlider->assignParameter(
      plotEditor->equalizerModuleToEdit->getParameterByName("GlobalGain") );
  }

  updateWidgetsAccordingToState();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void EqualizerModuleEditor::setUseShortSliderNames(bool shouldBeShort)
{
  useShortSliderNames = shouldBeShort;
  updateWidgetAppearance();
}

void EqualizerModuleEditor::setUseSmallComboBox(bool shouldBeSmall)
{
  useSmallComboBox = shouldBeSmall;
  updateWidgetAppearance();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void EqualizerModuleEditor::rButtonClicked(RButton *buttonThatWasClicked)
{
  ScopedPointerLock spl(lock);
  if( plotEditor->equalizerModuleToEdit == NULL )
    return;

  if( buttonThatWasClicked == channelSelectButton1 || buttonThatWasClicked == channelSelectButton2 )
    plotEditor->equalizerModuleToEdit->selectChannel( (bool) channelSelectButton2->getToggleState() );
  else
    AudioModuleEditor::rButtonClicked(buttonThatWasClicked);
}

void EqualizerModuleEditor::rComboBoxChanged(RComboBox  *rComboBoxThatHasChanged)
{
  ScopedPointerLock spl(lock);
  updateWidgetVisibility();
  updateWidgetAppearance();
  updatePlotRange();
}

void EqualizerModuleEditor::changeListenerCallback(ChangeBroadcaster *objectThatHasChanged)
{
  ScopedPointerLock spl(lock);
  if( plotEditor->equalizerModuleToEdit == NULL )
    return;

  // the call must have been due to preset recall - deselect band in this case:
  plotEditor->equalizerModuleToEdit->deSelectBand();
  AudioModuleEditor::changeListenerCallback(objectThatHasChanged);

  /*
  if( objectThatHasChanged == plotEditor->equalizerModuleToEdit )
  {
  //plotEditor->equalizerModuleToEdit->markStateAsDirty(); // move this call into EqualizerAudioModule
  updateWidgetsAccordingToState();

  // maybe we can get rid of the branch now that we have parameterSetChanged?
  }
  else
  {
  // the call must have been due to preset recall - deselect band in this case:
  plotEditor->equalizerModuleToEdit->deSelectBand();
  AudioModuleEditor::changeListenerCallback(objectThatHasChanged);
  }
  */
}

void EqualizerModuleEditor::parameterSetChanged(ParameterSetHolder* parameterSetHolderThatHasChanged)
{
  ScopedPointerLock spl(lock);
  if( plotEditor->equalizerModuleToEdit == NULL )
    return;

  if( parameterSetHolderThatHasChanged == plotEditor->equalizerModuleToEdit )
    updateWidgetsAccordingToState();
}

void EqualizerModuleEditor::updateWidgetsAccordingToState()
{
  ScopedPointerLock spl(lock);

  AudioModuleEditor::updateWidgetsAccordingToState();
  if( plotEditor->equalizerModuleToEdit == NULL )
  {
    updateWidgetVisibility();
    updateWidgetAppearance();
    return;
  }

  // un-assign widgets (maybe we should do this first?):
  filterModeComboBox->assignParameter(nullptr);
  frequencySlider->assignParameter(nullptr);
  gainSlider->assignParameter(nullptr);
  bandwidthSlider->assignParameter(nullptr);

  // connect the sliders to the parameters of the band selected channel/band:
  int channel = plotEditor->equalizerModuleToEdit->selectedChannel;
  int index   = plotEditor->equalizerModuleToEdit->selectedIndex;
  if( index != -1 )
  {
    filterModeComboBox->assignParameter(plotEditor->equalizerModuleToEdit->filterModeParameters[channel][index]);
    frequencySlider->assignParameter(   plotEditor->equalizerModuleToEdit->frequencyParameters[channel][index]);
    gainSlider->assignParameter(        plotEditor->equalizerModuleToEdit->gainParameters[channel][index]);
    bandwidthSlider->assignParameter(   plotEditor->equalizerModuleToEdit->bandwidthParameters[channel][index]);
  }

  if( channel == 0 )
    channelSelectButton1->setToggleState(true, false);
  else
    channelSelectButton2->setToggleState(true, false);

  updatePlotRange();
  updateWidgetVisibility();
  updateWidgetAppearance();
  plotEditor->assignParametersToSelectedBand();
  plotEditor->updatePlot();
}

void EqualizerModuleEditor::paint(Graphics &g)
{
  AudioModuleEditor::paint(g);
  if(layout != SLIDERS_ABOVE)
  {
    fillRectWithBilinearGradient(g, rightSectionRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
      editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
    fillRectWithBilinearGradient(g, bottomSectionRectangle, editorColourScheme.topLeft, editorColourScheme.topRight,
      editorColourScheme.bottomLeft, editorColourScheme.bottomRight);

    g.setColour(editorColourScheme.outline);
    g.drawRect(rightSectionRectangle);
    g.drawRect(bottomSectionRectangle);
  }
  // ...this could actually be done in the baseclass via the guiRectangles array, i think

  /*
  // gray out the editor if it's disabled:
  if( !isEnabled() )
  g.fillAll(Colours::lightgrey.withAlpha(0.75f));
  */
}

void EqualizerModuleEditor::resized()
{
  AudioModuleEditor::resized();

  int x = 0;
  int y = getPresetSectionBottom();
  int w = getWidth();
  int h = getHeight();

  if( layout == SLIDERS_RIGHT )
  {
    // the size/position computations are a mess - improve!
    bandParametersLabel->setVisible(true);

    int rightSectionWidth   = 128;
    int bottomSectionHeight = 24;

    x = w-rightSectionWidth;
    h = infoField->getY() - y - bottomSectionHeight;

    rightSectionRectangle.setBounds(x, y+4, rightSectionWidth, h-4);
    x = rightSectionRectangle.getX();
    y = rightSectionRectangle.getY();
    w = rightSectionRectangle.getWidth();

    y = rightSectionRectangle.getY();

    bandParametersLabel->setBounds(x+4, y+4, w-8, 16);

    y += 24;

    filterModeComboBox->setNameLabelPosition(RNamedComboBox::ABOVE_BOX);
    filterModeComboBox->setBounds(x+4, y+4, w-8, 32);

    y += 40;

    frequencySlider->setLayout(RSlider::NAME_ABOVE);
    frequencySlider->setBounds(x+4, y+4, w-8, 32);

    y += 36;

    gainSlider->setLayout(RSlider::NAME_ABOVE);
    gainSlider->setBounds(x+4, y+4, w-8, 32);

    y += 36;

    bandwidthSlider->setLayout(RSlider::NAME_ABOVE);
    bandwidthSlider->setBounds(x+4, y+4, w-8, 32);

    y = rightSectionRectangle.getBottom()-2*16-8;

    //globalGainSlider->setLayout(RSlider::NAME_ABOVE);
    //globalGainSlider->setBounds(x+4, y+4, w-8, 32);

    x = 0;
    w = rightSectionRectangle.getX()-x;
    y = getPresetSectionBottom();
    h = infoField->getY() - y - bottomSectionHeight;
    plotEditor->setBounds(x, y+4, w+2, h-4);

    // new widgets for v10.04:
    y = plotEditor->getBottom()-2;
    bottomSectionRectangle.setBounds(0, y, getWidth(), bottomSectionHeight+4); 

    x = bottomSectionRectangle.getX();

    y += 1;
    bypassButton->setBounds(x+4, y+4, 60, 16);
    stereoModeComboBox->setBounds(bypassButton->getRight()+4, y+4, 180, 16);
    channelSelectButton1->setBounds(stereoModeComboBox->getRight()+4,   y+4, 16, 16);
    channelSelectButton2->setBounds(channelSelectButton1->getRight()-2, y+4, 16, 16);

    x = rightSectionRectangle.getX();
    w = rightSectionRectangle.getWidth();
    globalGainSlider->setBounds(x+4, y+4, w-8, 16);

    w = 92;
    gainRangeComboBox->setBounds(globalGainSlider->getX()-w-4, y+4, w, 16);
  }
  if( layout == SLIDERS_BELOW )
  {
    bandParametersLabel->setVisible(false); // correct?

    x = 0;
    w = getWidth();
    y = getPresetSectionBottom();
    h = getHeight()-40-y;
    plotEditor->setBounds(4, y+4, w-8, h-8);

    y = plotEditor->getBottom();
    w = getWidth()/3;
    frequencySlider->setBounds(x+4, y+4, w-4, 16);
    x = frequencySlider->getRight();
    gainSlider->setBounds(x+4, y+4, w-4, 16);
    x = gainSlider->getRight();
    w = getWidth()-x;
    bandwidthSlider->setBounds(x+4, y+4, w-8, 16);

    y += 20;

    filterModeComboBox->setNameLabelPosition(RNamedComboBox::LEFT_TO_BOX);
    filterModeComboBox->setBounds(4, y+4, 120, 16);
    //modeLabel->setBounds(4, y+4, 40, 16);
    //filterModeComboBox->setBounds(modeLabel->getRight(), y+4, 120, 16);

    x = filterModeComboBox->getRight();
    w = getWidth()-x;
    globalGainSlider->setBounds(x+4, y+4, w-8, 16);


    //useShortSliderNames = true;

    /*
    globalGainSlider->setSliderName(juce::String(T("GG")));
    frequencySlider->setSliderName(juce::String(T("F")));
    gainSlider->setSliderName(juce::String(T("G")));
    bandwidthSlider->setSliderName(juce::String(T("B")));
    */
  }
  if( layout == SLIDERS_ABOVE )
  {
    bandParametersLabel->setVisible(false);

    //x = getHeadlineRight()+40;
    x = getWidth()/2;
    w = getWidth()-x;

    //globalGainSlider->setBounds(x+4, y+4, w-8, 16);
    //stateWidgetSet->stateFileNameLabel->setVisible(false);
    stateWidgetSet->stateFileNameLabel->setVisible(true);

    //x  = stateWidgetSet->stateSaveButton->getX()-16-80;
    //x += stateWidgetSet->getX();
    //globalGainSlider->setBounds(x, 6, 80, 16);
    x = getHeadlineRight() + 4;
    w = stateWidgetSet->getX() - x - 8;
    globalGainSlider->setBounds(x, 6, w, 16);

    x = 0;
    y = getPresetSectionBottom();

    filterModeComboBox->setBounds(x+4, y+4, 80, 16);

    x = filterModeComboBox->getRight();
    w = (getWidth()-x) /4;
    frequencySlider->setBounds(x+4, y+4, w-4, 16);

    x = frequencySlider->getRight();
    gainSlider->setBounds(x+4, y+4, w-4, 16);

    x = gainSlider->getRight();
    bandwidthSlider->setBounds(x+4, y+4, w-4, 16);

    x = bandwidthSlider->getRight();
    w = getWidth() - x;
    globalGainSlider->setBounds(x+4, y+4, w-8, 16);

    y = gainSlider->getBottom()+4;
    w = getWidth();
    h = getHeight()-y;
    plotEditor->setBounds(0, y, w, h);

    /*
    globalGainSlider->setSliderName(juce::String(T("GG")));
    frequencySlider->setSliderName(juce::String(T("F")));
    gainSlider->setSliderName(juce::String(T("G")));
    bandwidthSlider->setSliderName(juce::String(T("B")));
    filterModeComboBox->setItemText(2, juce::String(T("Low Shelf")));
    filterModeComboBox->setItemText(3, juce::String(T("High Shelf")));
    filterModeComboBox->setItemText(4, juce::String(T("Lowpass 6")));
    filterModeComboBox->setItemText(5, juce::String(T("Lowpass 12")));
    filterModeComboBox->setItemText(6, juce::String(T("Highpass 6")));
    filterModeComboBox->setItemText(7, juce::String(T("Highpass 12")));
    filterModeComboBox->setItemText(8, juce::String(T("Notch")));
    */
  }

  updateWidgetVisibility();
  updateWidgetAppearance();
}

void EqualizerModuleEditor::createWidgets()
{
  typedef rsModulatableSlider Sld;
  typedef RNamedComboBox Box;
  typedef rsAutomatableButton Btn;
  Sld* s;
  //Box* c;
  Btn* b;

  addWidget( bypassButton = b = new Btn("Bypass") );
  b->setDescription("Bypass the whole equalizer");
  b->setDescriptionField(infoField);
  b->setClickingTogglesState(true);  //...? do we need this ? ...set this active by default in RButton

  addWidget( channelSelectButton1 = new RRadioButton("L") );
  channelSelectButton1->setDescription("Edit left channel curve");
  channelSelectButton1->setDescriptionField(infoField);
  channelSelectButton1->setClickingTogglesState(true);
  channelSelectButton1->addRButtonListener(this);
  //channelSelectButton1->setRadioGroupId(1);  // \todo: make this work again
  channelSelectButton1->addToRadioButtonGroup(&channelSelectRadioGroup);
  channelSelectButton1->setToggleState(true, false);

  addWidget( channelSelectButton2 = new RRadioButton("R") );
  channelSelectButton2->setDescription("Edit right channel curve");
  channelSelectButton2->setDescriptionField(infoField);
  channelSelectButton2->setClickingTogglesState(true);
  channelSelectButton2->addRButtonListener(this);
  //channelSelectButton2->setRadioGroupId(1); // \todo: make this work again
  channelSelectButton2->addToRadioButtonGroup(&channelSelectRadioGroup);
  channelSelectButton2->setToggleState(false, false);

  addWidget( stereoModeComboBox = new RNamedComboBox(juce::String("StereoModeComboBox"),
    juce::String("Stereo Mode:")) );
  stereoModeComboBox->setDescription("Select mode for processing stereo signals");
  stereoModeComboBox->setDescriptionField(infoField);
  stereoModeComboBox->registerComboBoxObserver(this);

  addWidget( gainRangeComboBox = new RNamedComboBox(juce::String("RangeComboBox"),
    juce::String("Range:")) );
  gainRangeComboBox->setDescription("Select the range for the plot");
  gainRangeComboBox->setDescriptionField(infoField);
  gainRangeComboBox->registerComboBoxObserver(this);

  addWidget( globalGainSlider = s = new Sld );
  s->setDescription("Global gain to compensate for loudness change");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( bandParametersLabel = new RTextField("Band Parameters") );
  bandParametersLabel->setDescription("Parameters for th selected band");
  bandParametersLabel->setDescriptionField(infoField);
  bandParametersLabel->setJustification(Justification::centred);
  bandParametersLabel->setNoBackgroundAndOutline(true);

  addWidget( filterModeComboBox = new RNamedComboBox(juce::String("FilterModeComboBox"),
    juce::String("Mode:")) );
  filterModeComboBox->setDescription("Filter mode of selected band");
  filterModeComboBox->setDescriptionField(infoField);
  filterModeComboBox->registerComboBoxObserver(this);

  addWidget( frequencySlider = new Sld );
  frequencySlider->setSliderName(juce::String("Frequency"));
  frequencySlider->setDescription(juce::String("Frequency of selected band"));
  frequencySlider->setDescriptionField(infoField);
  frequencySlider->setStringConversionFunction(&hertzToStringWithUnitTotal5);

  addWidget( gainSlider = new Sld );
  gainSlider->setSliderName(juce::String("Gain"));
  gainSlider->setDescription(juce::String("Gain of selected band"));
  gainSlider->setDescriptionField(infoField);
  gainSlider->setStringConversionFunction(&decibelsToStringWithUnit1);

  addWidget( bandwidthSlider = new Sld );
  bandwidthSlider->setSliderName(juce::String("Bandwidth"));
  bandwidthSlider->setDescription(juce::String("Bandwidth of selected band"));
  bandwidthSlider->setDescriptionField(infoField);
  bandwidthSlider->setStringConversionFunction(&octavesToStringWithUnit2);

  plotEditor = new EqualizerPlotEditor(lock, equalizerModule);
  plotEditor->setDescriptionField(infoField);
  //plotEditor->addChangeListener(this);
  //if( equalizerModuleToEdit != NULL )
  //  plotEditor->setEqualizerToEdit(equalizerModuleToEdit->wrappedEqualizer);
  addPlot( plotEditor );
}

void EqualizerModuleEditor::updateWidgetVisibility()
{
  ScopedPointerLock spl(lock);

  frequencySlider->setVisible(     false);
  gainSlider->setVisible(          false);
  bandwidthSlider->setVisible(     false);
  filterModeComboBox->setVisible(  false);
  channelSelectButton1->setVisible(false);
  channelSelectButton2->setVisible(false);
  globalGainSlider->setVisible(    false);
  stateWidgetSet->setVisible(      false);

  if( plotEditor->equalizerModuleToEdit == NULL )
    return;

  globalGainSlider->setVisible(    true);
  stateWidgetSet->setVisible(      true);

  int channel = plotEditor->equalizerModuleToEdit->selectedChannel;
  int index   = plotEditor->equalizerModuleToEdit->selectedIndex;
  if( index != -1 )
  {
    filterModeComboBox->setVisible(true);
    frequencySlider->setVisible(true);
    if( plotEditor->equalizerModuleToEdit->wrappedEqualizerStereo->doesModeSupportBandwidth(channel, index) )
      bandwidthSlider->setVisible(true);
    if( plotEditor->equalizerModuleToEdit->wrappedEqualizerStereo->doesModeSupportGain(channel, index) )
      gainSlider->setVisible(true);
  }
}


void EqualizerModuleEditor::updateWidgetAppearance()
{
  ScopedPointerLock spl(lock);

  if( plotEditor->equalizerModuleToEdit == NULL )
    return;

  if( plotEditor->equalizerModuleToEdit->wrappedEqualizerStereo->getStereoMode() == rosic::EqualizerStereo::LEFT_RIGHT )
  {
    channelSelectButton1->setVisible(true);
    channelSelectButton2->setVisible(true);
    channelSelectButton1->setButtonText("L");
    channelSelectButton2->setButtonText("R");
    channelSelectButton1->setDescription(juce::String("Edit curve for left channel"));
    channelSelectButton2->setDescription(juce::String("Edit curve for right channel"));
  }
  else if( plotEditor->equalizerModuleToEdit->wrappedEqualizerStereo->getStereoMode() == rosic::EqualizerStereo::MID_SIDE )
  {
    channelSelectButton1->setVisible(true);
    channelSelectButton2->setVisible(true);
    channelSelectButton1->setButtonText("M");
    channelSelectButton2->setButtonText("S");
    channelSelectButton1->setDescription(juce::String("Edit curve for mid channel"));
    channelSelectButton2->setDescription(juce::String("Edit curve for side channel"));
  }

  if( useShortSliderNames == true )
  {
    globalGainSlider->setSliderName("GG");
    frequencySlider->setSliderName( "F");
    gainSlider->setSliderName(      "G");
    bandwidthSlider->setSliderName( "B");

    // currently, when we use short names, we also want the compact layout:
    globalGainSlider->setLayout(RSlider::NAME_INSIDE);
    frequencySlider->setLayout( RSlider::NAME_INSIDE);
    gainSlider->setLayout(      RSlider::NAME_INSIDE);
    bandwidthSlider->setLayout( RSlider::NAME_INSIDE);

  }

  if( useSmallComboBox == true )
  {
    filterModeComboBox->setNameLabelWidth(0); // workaround because INVISIBLE doesn't work
    filterModeComboBox->setNameLabelPosition(RNamedComboBox::LEFT_TO_BOX);
    filterModeComboBox->setItemText(2, "Low Shelf");
    filterModeComboBox->setItemText(3, "High Shelf");
    filterModeComboBox->setItemText(4, "Lowpass 6");
    filterModeComboBox->setItemText(5, "Lowpass 12");
    filterModeComboBox->setItemText(6, "Highpass 6");
    filterModeComboBox->setItemText(7, "Highpass 12");
    filterModeComboBox->setItemText(8, "Notch");
  }
}

void EqualizerModuleEditor::updatePlotRange()
{
  int rangeIndex = (int) plotEditor->equalizerModuleToEdit->getParameterByName(
    juce::String("GainRange"))->getValue();
  switch( rangeIndex )
  {
  case 0:
  {
    plotEditor->setCurrentRangeY( -3.0, +3.0);
    plotEditor->setHorizontalCoarseGrid(1.0, true);
  } break;
  case 1:
  {
    plotEditor->setCurrentRangeY( -6.0, +6.0);
    plotEditor->setHorizontalCoarseGrid(2.0, true);
  }
  break;
  case 2:
  {
    plotEditor->setCurrentRangeY(-12.0, +12.0);
    plotEditor->setHorizontalCoarseGrid(3.0, true);
  }
  break;
  case 3:
  {
    plotEditor->setCurrentRangeY(-24.0, +24.0);
    plotEditor->setHorizontalCoarseGrid(6.0, true);
  }
  break;
  case 4:
  {
    plotEditor->setCurrentRangeY(-48.0, +48.0);
    plotEditor->setHorizontalCoarseGrid(12.0, true);
  }
  break;
  }
  plotEditor->updatePlot();
}

/*

maybe parametrize a biquad via:
-midFreq, midGain, midWidth (regular bell filter params)
-lowGain, highGain (would be 0 dB in bell EQ, -inf in BP, 0,-inf in LP etc)
-or maybe instead of lowGain/highGain have meanGain and gainDelta
-lowGain/highGain would apply to 0 and infinite freq

*/
