PhaseScope::PhaseScope(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);

  setModuleTypeName("Scope");
  phaseScopeBuffer = new RAPT::rsPhaseScopeBuffer<double, float, double>;

  pixelScale = 1.0;
  displayWidth  = 100;
  displayHeight = 100;
  bypassPixelDecay = false;
  //phaseScopeBuffer->setUseAlphaMask(false);
  phaseScopeBuffer->setMaxSizeWithoutReAllocation(1000, 1000);
  //phaseScopeBuffer->setMaxSizeWithoutReAllocation(3840, 2160); // UHD
  //phaseScopeBuffer->setMaxSizeWithoutReAllocation(4096, 4096);
  updateBufferSize();

  createParameters();
  reset();

  //colorMap.setDefaultMap(ColorMap::fire);
  //colorMap.setDefaultMap(ColorMap::ice);

  //juce::ColourGradient g;
  //g.addColour(0.0, Colour(  0,   0,   0));
  //g.addColour(0.2, Colour(  0,   0, 255));
  //g.addColour(0.4, Colour(  0, 255, 255));
  //g.addColour(0.6, Colour(  0, 255,   0));
  //g.addColour(0.8, Colour(255, 255,   0));
  //g.addColour(1.0, Colour(255,   0,   0)); // maybe add magenta as last
  //colorMap.setFromColourGradient(g);
}

PhaseScope::~PhaseScope()
{
  delete phaseScopeBuffer;
}

void PhaseScope::createParameters()
{
  ScopedLock scopedLock(*lock);

  Parameter* p;

  p = new Parameter(lock, "Brightness", 0.001, 100.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  //p = new Parameter(plugInLock, "Brightness", -100.0, 100.0, 0.0, 1.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setBrightness);

  p = new Parameter(lock, "AfterGlow", 0.001, 50.0, 0.0, 0.1, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setAfterGlow);

  p = new Parameter(lock, "PixelSpread", 0.0, 1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setPixelSpread);

  p = new Parameter(lock, "PixelScale", 1.0, 8.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setPixelScale);

  p = new Parameter(lock, "LineDensity", 0.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setLineDensity);

  p = new Parameter(lock, "DotLimit", 1.0, 500.0, 1.0, 500.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setDotLimit);

  p = new Parameter(lock, "FrameRate", 1.0, 100.0, 0.0, 25.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setFrameRate);

  p = new Parameter("DrawMode", 0.0, 1.0, 0.0, Parameter::STRING);
  p->addStringValue("Linear");
  p->addStringValue("Quadratic");
  p->addStringValue("Cubic"); // maybe have density-compensated and uncompensated modes
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setDrawMode);

  p = new Parameter(lock, "AntiAlias", 0.0, 1.0, 0.0, 1.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setAntiAlias);
  addObservedParameter(p);

  p = new Parameter(lock, "OneDimensional", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setOneDimensionalMode);
  addObservedParameter(p);

  p = new Parameter(lock, "ScanFrequency", 0.1, 10.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setScanningFrequency);

  p = new Parameter(lock, "NumCycles", 1.0, 10.0, 1.0, 2.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setNumCyclesShown);

  //p = new Parameter(lock, "Zoom", 0.1, 10.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  //addObservedParameter(p);
  //p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setZoom);

  p = new Parameter(lock, "Sync", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setSyncMode);
  addObservedParameter(p);




  // the geometric trafo parameters are of a different type because they were copy/pasted from
  // PrettyScope - eventually, they should all be the same type:

  typedef MetaControlledParameter Param; // shortcut for convenience
  //Param* p;

  // geometric transforms:
  p = new Param("ScaleX", 0.1, 10.0, 1.0, Parameter::EXPONENTIAL);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setScaleX);
  addObservedParameter(p);
  p = new Param("ScaleY", 0.1, 10.0, 1.0, Parameter::EXPONENTIAL);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setScaleY);
  addObservedParameter(p);
  // maybe instead of ScaleX/Y, we could have Scale and Ratio where the former scales both, x and y 
  // and the latter modfied the aspect ratio - this could be more convenient, especially in XY mode

  //p = new Param("ShearX", -8.0, 8.0, 0.0, Parameter::LINEAR);
  //p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setShearX);
  //addObservedParameter(p);
  //p = new Param("ShearY", -8.0, 8.0, 0.0, Parameter::LINEAR);
  //p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setShearY);
  //addObservedParameter(p);

  p = new Param("Rotation", -180, 180, 0.0, Parameter::LINEAR);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setRotation);
  addObservedParameter(p);

  //p = new Param("ShiftX", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setShiftX);
  //addObservedParameter(p);
  //p = new Param("ShiftY", -1.0, 1.0, 0.0, Parameter::LINEAR);
  //p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setShiftY);
  //addObservedParameter(p);

}

void PhaseScope::setDisplayPixelSize(int width, int height)
{
  ScopedLock scopedLock(*lock);
  displayWidth  = width;
  displayHeight = height; 
  updateBufferSize();
}

void PhaseScope::setBrightness(double newBrightness)
{
  phaseScopeBuffer->setBrightness((float)newBrightness);
}
void PhaseScope::setAfterGlow(double newGlow)
{
  phaseScopeBuffer->setDecayTime(newGlow);
  if(newGlow == 50.0)
    bypassPixelDecay = true;
  else
    bypassPixelDecay = false;
}
void PhaseScope::setLineDensity(double newDensity)
{
  phaseScopeBuffer->setLineDensity((float)newDensity);
}
void PhaseScope::setDotLimit(double newLimit)
{
  phaseScopeBuffer->setLineDensityLimit((int) ::round(newLimit));
}
void PhaseScope::setPixelSpread(double newSpread)
{
  phaseScopeBuffer->setPixelSpread((float)newSpread);
}
void PhaseScope::setPixelScale(double newFactor)
{
  jassert(newFactor >= 1.0);
  pixelScale = newFactor;
  updateBufferSize();
}
void PhaseScope::setDrawMode(int mode)
{
  phaseScopeBuffer->setDrawMode(mode);
}
void PhaseScope::setAntiAlias(bool shouldAntiAlias)
{
  phaseScopeBuffer->setAntiAlias(shouldAntiAlias);
}
void PhaseScope::setFrameRate(double newRate)
{
  phaseScopeBuffer->setFrameRate(newRate);
  updateRepaintInterval();
}

void PhaseScope::setScaleX(double newScale)
{
  phaseScopeBuffer->setScaleX(newScale);
}
void PhaseScope::setScaleY(double newScale)
{
  phaseScopeBuffer->setScaleY(newScale);
}
//void PhaseScope::setShearX(double newShear)
//{
//  phaseScopeBuffer->setShearX(newShear);
//}
//void PhaseScope::setShearY(double newShear)
//{
//  phaseScopeBuffer->setShearY(newShear);
//}
void PhaseScope::setRotation(double degrees)
{
  phaseScopeBuffer->setRotation(degrees);
}
//void PhaseScope::setShiftX(double newShift)
//{
//  phaseScopeBuffer->setShiftX(newShift);
//}
//void PhaseScope::setShiftY(double newShift)
//{
//  phaseScopeBuffer->setShiftY(newShift);
//}
void PhaseScope::setOneDimensionalMode(bool shouldBe1D)
{
  phaseScopeBuffer->setOneDimensionalMode(shouldBe1D);
}
void PhaseScope::setScanningFrequency(double newFrequency)
{
  phaseScopeBuffer->setScanningFrequency(newFrequency);
}

void PhaseScope::setNumCyclesShown(int newNumCycles)
{
  phaseScopeBuffer->screenScanner.setNumCyclesShown(newNumCycles);
}

//void PhaseScope::setZoom(double newZoom)
//{
//  phaseScopeBuffer->screenScanner.setZoom(newZoom);
//}

void PhaseScope::setSyncMode(bool shouldSync)
{
  phaseScopeBuffer->setSyncMode(shouldSync);
}

AudioModuleEditor* PhaseScope::createEditor(int type)
{
  return new PhaseScopeEditor(this);
}

void PhaseScope::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  for(int n = 0; n < numSamples; n++)
    processStereoFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void PhaseScope::processStereoFrame(double *left, double *right)
{
  phaseScopeBuffer->processSampleFrame((float)(*left), (float)(*right));
  repaintCounter++;
  if(repaintCounter > repaintIntervalInSamples)
  {
    updateScopeImage();
    sendImageUpdateNotification(&image);
    repaintCounter = 0;
  }
}

void PhaseScope::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*lock);
  phaseScopeBuffer->setSampleRate(newSampleRate);
  updateRepaintInterval();
}

void PhaseScope::reset()
{
  ScopedLock scopedLock(*lock);
  phaseScopeBuffer->reset();
  repaintCounter = 0;
}

void PhaseScope::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);
  XmlElement* xmlColorMap = xmlState.getChildByName("ColorMap");
  if(xmlColorMap != nullptr)
    colorMap.setFromXml(*xmlColorMap);
}

XmlElement* PhaseScope::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  XmlElement* xml = AudioModule::getStateAsXml(stateName, markAsClean);
  XmlElement* colorMapXml = colorMap.getAsXml();
  xml->addChildElement(colorMapXml);
  return xml;
}

void PhaseScope::updateBufferSize()
{
  ScopedLock scopedLock(*lock);
  int w = (int) ::round(displayWidth  / pixelScale);
  int h = (int) ::round(displayHeight / pixelScale);
  w = jmax(w, 1);
  h = jmax(h, 1);

  jassert(w <= phaseScopeBuffer->getImage()->getMaxWidth());
  jassert(h <= phaseScopeBuffer->getImage()->getMaxHeight());

  phaseScopeBuffer->setSize(w, h);
  image = juce::Image(juce::Image::ARGB, w, h, false);
}

void PhaseScope::updateScopeImage()
{
  normalizedDataToImage(phaseScopeBuffer->getImage()->getPixelPointer(0, 0), image, colorMap);
  if(!bypassPixelDecay)
    phaseScopeBuffer->applyPixelDecay();
}

void PhaseScope::updateRepaintInterval()
{
  repaintIntervalInSamples = 
    (int) ::round(phaseScopeBuffer->getSampleRate() / phaseScopeBuffer->getFrameRate());
  repaintCounter = 0;
}

//=================================================================================================

PhaseScopeDisplay::PhaseScopeDisplay(jura::PhaseScope *newPhaseScopeToEdit)
{
  phaseScope = newPhaseScopeToEdit;
  phaseScope->addImageUpdateListener(this);
  addChangeListener(this);
}

PhaseScopeDisplay::~PhaseScopeDisplay()
{
  phaseScope->removeImageUpdateListener(this);
}

void PhaseScopeDisplay::resized()
{
  phaseScope->setDisplayPixelSize(getWidth(), getHeight());
}

void PhaseScopeDisplay::paint(Graphics &g)
{
  //ScopedLock scopedLock(*(phaseScope->lock));

  g.setImageResamplingQuality(Graphics::lowResamplingQuality);
  //g.setImageResamplingQuality(Graphics::mediumResamplingQuality);
  //g.setImageResamplingQuality(Graphics::highResamplingQuality);
  g.drawImage(phaseScope->image, juce::Rectangle<float>(0.f, 0.f, (float) getWidth(), 
    (float) getHeight()));
}

void PhaseScopeDisplay::imageWasUpdated(juce::Image* image)
{
  sendChangeMessage();
  // We will receive the change message ourselves and when we do, we will trigger a repaint. We 
  // can't directly call repaint here because then we hit an assertion which says we should acquire
  // a MessageManagerLock - but when we do so, it becomes unresponsive.

  // todo: use repaintOnMessageThread

  // this doesn't work:
  //const MessageManagerLock mmLock;
  //repaint();
}

void PhaseScopeDisplay::changeListenerCallback(ChangeBroadcaster *source)
{
  repaint();
}

//=================================================================================================

PhaseScopeEditor::PhaseScopeEditor(jura::PhaseScope *newPhaseScopeToEdit)
  : AudioModuleEditor(newPhaseScopeToEdit)
  , display(newPhaseScopeToEdit)
{
  ScopedLock scopedLock(*lock);
  scope = newPhaseScopeToEdit;
  widgetMargin = 150; 

  addAndMakeVisible(display);
  createWidgets();

  colorMapLoader->setColorMapDirectory(getSupportDirectory() + File::getSeparatorString() 
    + "ColorMaps");

  int headerMargin = 26;  // this is the height we need for headline and preset-section
  int displaySize = 250;
  setSize(displaySize+widgetMargin, displaySize+headerMargin);

  // needs to be done to show/hid sliders depending on whether sync is on or off
  setLocalAutomationSwitch(true);
  setIsGuiElement(true);
  scope->getParameterByName("Sync")->registerParameterObserver(this);
  parameterChanged(scope->getParameterByName("Sync"));

  // todo: we should really use a ResizableImage in PhaseScopeBuffer, set up a maximum size and 
  // then set up the ComponentResizeBoundsConstraine accordingly - otherwise, we may get threading
  // problems due to reallocating memory in the gui thread while writing into the same memory
  // in the audio thread
}

PhaseScopeEditor::~PhaseScopeEditor()
{
  scope->getParameterByName("Sync")->deRegisterParameterObserver(this);
}

void PhaseScopeEditor::createWidgets()
{
  RButton *b;
  RSlider *s;

  addWidget( sliderBrightness = s = new RSlider("BrightnessSlider") );
  s->assignParameter( scope->getParameterByName("Brightness") );
  s->setSliderName("Brightness");
  s->setDescription("Brightness");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( sliderAfterglow = s = new RSlider("GlowSlider") );
  s->assignParameter( scope->getParameterByName("AfterGlow") );
  s->setSliderName("Glow");
  s->setDescription("Afterglow time in seconds");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&secondsToStringWithUnitTotal4);

  addWidget( sliderPixelSpread = s = new RSlider("SpreadSlider") );
  s->assignParameter( scope->getParameterByName("PixelSpread") );
  s->setSliderName("Spread");
  s->setDescription("Pixel spreading");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( sliderPixelScale = s = new RSlider("ScaleSlider") );
  s->assignParameter( scope->getParameterByName("PixelScale") );
  s->setSliderName("Scale");
  s->setDescription("Pixel size rescaling");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( sliderLineDensity = s = new RSlider("DensitySlider") );
  s->assignParameter( scope->getParameterByName("LineDensity") );
  s->setSliderName("Density");
  s->setDescription("Line Density");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( sliderDotLimit = s = new RSlider("DotLimitSlider") );
  s->assignParameter( scope->getParameterByName("DotLimit") );
  s->setSliderName("DotLimit");
  s->setDescription("Limit for number of dots per line");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( sliderFrameRate = s = new RSlider("FrameRateSlider") );
  s->assignParameter( scope->getParameterByName("FrameRate") );
  s->setSliderName("FrameRate");
  s->setDescription("Frame rate for redrawing");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( boxDrawMode = new RComboBox() );
  boxDrawMode->assignParameter( scope->getParameterByName("DrawMode") );
  boxDrawMode->setDescription("Select drawing mode (interpolation between incoming samples)");
  boxDrawMode->setDescriptionField(infoField);

  addWidget( buttonAntiAlias = b = new RButton("AntiAlias") );
  b->assignParameter( scope->getParameterByName("AntiAlias") );
  b->setDescription("Anti aliased drawing (bilinear deinterpolation)");
  b->setDescriptionField(infoField);
  //b->setButtonPainter(&buttonPainter); // temporary, for test

  addWidget( button1D = b = new RButton("1D") );
  b->assignParameter( scope->getParameterByName("OneDimensional") );
  b->setDescription("Replace x-input with sawtooth screen scanner");
  b->setDescriptionField(infoField);

  addWidget( buttonSync = b = new RButton("Sync") );
  b->assignParameter( scope->getParameterByName("Sync") );
  b->setDescription("Synchronize 1D scanning frequency to input");
  b->setDescriptionField(infoField);

  addWidget( s = sliderScanFreq = new rsAutomatableSlider );
  s->assignParameter( scope->getParameterByName("ScanFrequency") );
  s->setSliderName("ScanFreq");
  s->setDescription("Scanning frequency in 1D mode");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( s = sliderNumCycles = new rsAutomatableSlider );
  s->assignParameter( scope->getParameterByName("NumCycles") );
  s->setSliderName("NumCycles");
  s->setDescription("Number of cycles in synced 1D mode");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString0);




  // geometric transforms:
  addWidget(s = sliderScaleX = new rsAutomatableSlider);
  s->assignParameter(scope->getParameterByName("ScaleX"));
  s->setSliderName("ScaleX");
  s->setDescription("Scaling along x-direction");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget(s = sliderScaleY = new rsAutomatableSlider);
  s->assignParameter(scope->getParameterByName("ScaleY"));
  s->setSliderName("ScaleY");
  s->setDescription("Scaling along y-direction");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  //addWidget(s = sliderShearX = new rsAutomatableSlider);
  //s->assignParameter(scope->getParameterByName("ShearX"));
  //s->setSliderName("ShearX");
  //s->setDescription("Shearing along x-direction");
  //s->setDescriptionField(infoField);
  //s->setStringConversionFunction(&valueToString3);

  //addWidget(s = sliderShearY = new rsAutomatableSlider);
  //s->assignParameter(scope->getParameterByName("ShearY"));
  //s->setSliderName("ShearY");
  //s->setDescription("Shearing along y-direction");
  //s->setDescriptionField(infoField);
  //s->setStringConversionFunction(&valueToString3);

  addWidget(s = sliderRotation = new rsAutomatableSlider);
  s->assignParameter(scope->getParameterByName("Rotation"));
  s->setSliderName("Rotation");
  s->setDescription("Rotation around origin");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&degreesToStringWithUnit0);

  //addWidget(s = sliderShiftX = new rsAutomatableSlider);
  //s->assignParameter(scope->getParameterByName("ShiftX"));
  //s->setSliderName("ShiftX");
  //s->setDescription("Shifting along x-direction");
  //s->setDescriptionField(infoField);
  //s->setStringConversionFunction(&valueToString3);

  //addWidget(s = sliderShiftY = new rsAutomatableSlider);
  //s->assignParameter(scope->getParameterByName("ShiftY"));
  //s->setSliderName("ShiftY");
  //s->setDescription("Shifting along y-direction");
  //s->setDescriptionField(infoField);
  //s->setStringConversionFunction(&valueToString3);

  colorMapLoader = new ColorMapLoader(scope->getColorMapPointer());
  addWidgetSet(colorMapLoader);
}

void PhaseScopeEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();

  int w = getWidth();
  int h = getHeight();
  int y = getPresetSectionBottom() + 4;

  //display.setBounds(0, y, w-widgetMargin, h-y); // old - allows non-square display
  int ds = jmin(w-widgetMargin, h-y);  // display size - always square
  ds = jmax(ds, 1);                    // at least 1x1
  display.setBounds(0, y, ds, ds);

  // set up widgets:
  int x = display.getRight() + 4;
  w = getWidth() - x - 8;          // slider width
  h = 16;                          // slider height
  int dy = h-2;                    // vertical distance ("delta-y") between widgets

  sliderBrightness ->setBounds(x, y, w,   h); y += dy;
  sliderAfterglow  ->setBounds(x, y, w,   h); y += dy;
  sliderPixelSpread->setBounds(x, y, w,   h); y += dy;
  sliderPixelScale ->setBounds(x, y, w,   h); y += dy;
  sliderLineDensity->setBounds(x, y, w,   h); y += dy;
  sliderDotLimit   ->setBounds(x, y, w,   h); y += dy;
  sliderFrameRate  ->setBounds(x, y, w,   h); y += dy;
  boxDrawMode      ->setBounds(x,     y, w/2, h); //y += dy;
  buttonAntiAlias  ->setBounds(x+w/2, y, w/2, h); y += dy;

  // transform controls:
  y += 8;
  sliderScaleX  ->setBounds(x, y, w, h); y += dy;
  sliderScaleY  ->setBounds(x, y, w, h); y += dy;
  //sliderShearX  ->setBounds(x, y, w, h); y += dy;
  //sliderShearY  ->setBounds(x, y, w, h); y += dy;
  sliderRotation->setBounds(x, y, w, h); y += dy;
  //sliderShiftX  ->setBounds(x, y, w, h); y += dy;
  //sliderShiftY  ->setBounds(x, y, w, h); y += dy;
  sliderScanFreq  ->setBounds(x, y, w, h);
  sliderNumCycles ->setBounds(x, y, w, h); y += dy;
  button1D        ->setBounds(x,     y, w/2, h); 
  buttonSync      ->setBounds(x+w/2, y, w/2, h);
  y += dy;


  // preliminary:
  y += 8;
  colorMapLoader->setBounds(x, y, w, 48);
}

void PhaseScopeEditor::parameterChanged(Parameter* parameterThatHasChanged)
{
  bool sync = scope->getParameterByName("Sync")->getValue() >= 0.5;
  if(sync)
  {
    sliderScanFreq->setVisible(false);
    sliderNumCycles->setVisible(true);
  }
  else
  {
    sliderScanFreq->setVisible(true);
    sliderNumCycles->setVisible(false);
  }
}

//
//
////=================================================================================================
//// the artistically extended version of the PhaseScope:
//
//PhaseScope2::PhaseScope2(CriticalSection *lockToUse) : PhaseScope(lockToUse)
//{
//  ScopedLock scopedLock(*lock);
//  createParameters();  // creates the additional parameters
//
//  //phaseScopeBuffer->setUseAlphaMask(true);
//  // i think, either the dot brightness should scale inversely to the dotSize or the lineDensity
//}
//
//void PhaseScope2::setPixelDecayByValue(double newDecayByValue)
//{ 
//  phaseScopeBuffer->setPixelDecayByValue(newDecayByValue); 
//}
//void PhaseScope2::setPixelDecayByAverage(double newDecayByAverage)
//{
//  phaseScopeBuffer->setPixelDecayByAverage(newDecayByAverage);
//}
//void PhaseScope2::setDrawDots(bool shouldDraw)
//{
//  phaseScopeBuffer->setDrawDots(shouldDraw);
//}
//void PhaseScope2::setUseBigDot(bool shouldUseBigDot)
//{
//  phaseScopeBuffer->setUseAlphaMask(shouldUseBigDot);
//}
//void PhaseScope2::setDotSize(double newSize)
//{
//  phaseScopeBuffer->dotMask.setSize(newSize);
//}
//void PhaseScope2::setDotBlur(double newBlur)
//{
//  phaseScopeBuffer->dotMask.setTransitionWidth(newBlur);
//}
//void PhaseScope2::setDotInnerSlope(double newSlope)
//{
//  phaseScopeBuffer->dotMask.setInnerSlope(newSlope);
//}
//void PhaseScope2::setDotOuterSlope(double newSlope)
//{
//  phaseScopeBuffer->dotMask.setOuterSlope(newSlope);
//}
//void PhaseScope2::setDrawLines(bool shouldDraw)
//{
//  phaseScopeBuffer->setDrawLines(shouldDraw);
//}
//void PhaseScope2::setLineBrightness(double newBrightness)
//{
//  phaseScopeBuffer->setLineBrightness(newBrightness);
//}
//void PhaseScope2::setLineWidth(double newWidth)
//{
//  phaseScopeBuffer->setLineWidth(newWidth);
//}
//void PhaseScope2::setLineProfile(int newProfile)
//{
//  phaseScopeBuffer->setLineProfile(newProfile);
//}
//
//
//void PhaseScope2::createParameters()
//{
//  ScopedLock scopedLock(*lock);
//  Parameter* p;
//
//  p = new Parameter(lock, "DecayByValue", -20.0, +20.0, 0.0, 0.0, Parameter::LINEAR_BIPOLAR);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setPixelDecayByValue);
//
//  p = new Parameter(lock, "DecayByAverage", 0.0, +5.0, 0.0, 0.0, Parameter::LINEAR);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setPixelDecayByAverage);
//
//
//  p = new Parameter(lock, "DrawDots", 0.0, 1.0, 0.0, 1.0, Parameter::BOOLEAN);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDrawDots);
//  addObservedParameter(p);
//
//  p = new Parameter(lock, "UseBigDot", 0.0, 1.0, 0.0, 1.0, Parameter::BOOLEAN);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setUseBigDot);
//  addObservedParameter(p);
//
//  p = new Parameter(plugInLock, "DotSize", 1.0, 30.0, 0.0, 2.0, Parameter::LINEAR);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDotSize);
//
//  p = new Parameter(plugInLock, "DotBlur", 0.0, 1.0, 0.0, 0.5, Parameter::LINEAR);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDotBlur);
//
//  p = new Parameter(plugInLock, "DotInnerSlope", 0.0, 3.0, 0.0, 0.0, Parameter::LINEAR);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDotInnerSlope);
//
//  p = new Parameter(plugInLock, "DotOuterSlope", 0.0, 3.0, 0.0, 0.0, Parameter::LINEAR);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDotOuterSlope);
//
//
//  p = new Parameter(plugInLock, "DrawLines", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDrawLines);
//  addObservedParameter(p);
//
//  p = new Parameter(plugInLock, "LineBrightness", 0.001, 100.0, 0.0, 1.0, Parameter::EXPONENTIAL);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setLineBrightness);
//
//  p = new Parameter(plugInLock, "LineWidth", 1.0, 50.0, 0.0, 1.0, Parameter::EXPONENTIAL);
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setLineWidth);
//
//  p = new Parameter(plugInLock, "LineProfile", 0.0, 3.0, 1.0, 0.0, Parameter::STRING);
//  p->addStringValue("Flat");
//  p->addStringValue("Linear");
//  p->addStringValue("Parabolic");
//  p->addStringValue("Cubic");
//  addObservedParameter(p);
//  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setLineProfile);
//}
//
//AudioModuleEditor* PhaseScope2::createEditor()
//{
//  return new PhaseScopeEditor2(this);
//}
//
////-------------------------------------------------------------------------------------------------
//
//PhaseScopeEditor2::PhaseScopeEditor2(jura::PhaseScope2 *newPhaseScopeToEdit)
//  : PhaseScopeEditor(newPhaseScopeToEdit)
//{
//  createWidgets();
//
//  // init preview dot (maybe factor out):
//  //dotPreviewMask.setSize(widgetMargin-12.0); 
//  dotPreviewMask.setSize((widgetMargin-4.0)/2); 
//  dotPreviewImage = juce::Image(juce::Image::ARGB, 
//    dotPreviewMask.getWidth(), dotPreviewMask.getHeight(), false);
//  updatePreviewDot();
//
//  //resized();       // to position additional widgets
//  setSize(800, 600); // calls resized()
//}
//
//void PhaseScopeEditor2::createWidgets()
//{
//  RButton *b;
//  RSlider *s;
//  RComboBox *c;
//
//  addWidget( sliderDecayByValue = s = new RSlider("DecayByValueSlider") );
//  s->assignParameter( scope->getParameterByName("DecayByValue") );
//  s->setSliderName("DecayByValue");
//  s->setDescription("Dependency of pixel decay time on pixel brightness");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//
//  addWidget( sliderDecayByAverage = s = new RSlider("DecayByAverageSlider") );
//  s->assignParameter( scope->getParameterByName("DecayByAverage") );
//  s->setSliderName("DecayByAverage");
//  s->setDescription("Dependency of pixel decay time on average brightness of screen");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//
//  addWidget( buttonDrawDots = b = new RButton("Dots") );
//  b->assignParameter( scope->getParameterByName("DrawDots") );
//  b->setDescription("Switches dot drawing on/off");
//  b->setDescriptionField(infoField);
//
//  addWidget( buttonBigDot = b = new RButton("Big Dot") );
//  b->assignParameter( scope->getParameterByName("UseBigDot") );
//  b->setDescription("Switches to use the big expensive dot");
//  b->setDescriptionField(infoField);
//
//  addWidget( sliderDotSize = s = new RSlider("DotSizeSlider") );
//  s->assignParameter( scope->getParameterByName("DotSize") );
//  s->setSliderName("DotSize");
//  s->setDescription("Dot size in pixels");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//
//  addWidget( sliderDotBlur = s = new RSlider("DotBlurSlider") );
//  s->assignParameter( scope->getParameterByName("DotBlur") );
//  s->setSliderName("DotBlur");
//  s->setDescription("Dot blur from 0..1");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//  s->addListener(this);
//
//  addWidget( sliderDotInnerSlope = s = new RSlider("DotInnerSlopeSlider") );
//  s->assignParameter( scope->getParameterByName("DotInnerSlope") );
//  s->setSliderName("DotInnerSlope");
//  s->setDescription("Dot brightness slope at center");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//  s->addListener(this);
//
//  addWidget( sliderDotOuterSlope = s = new RSlider("DotOuterSlopeSlider") );
//  s->assignParameter( scope->getParameterByName("DotOuterSlope") );
//  s->setSliderName("DotOuterSlope");
//  s->setDescription("Dot brightness slope at border");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//  s->addListener(this);
//
//
//  addWidget( buttonDrawLines = b = new RButton("Lines") );
//  b->assignParameter( scope->getParameterByName("DrawLines") );
//  b->setDescription("Switches line drawing on/off");
//  b->setDescriptionField(infoField);
//
//  addWidget( sliderLineBrightness = s = new RSlider("LineBrightnessSlider") );
//  s->assignParameter( scope->getParameterByName("LineBrightness") );
//  s->setSliderName("LineBrightness");
//  s->setDescription("Brightness of lines");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//
//  addWidget( sliderLineWidth = s = new RSlider("LineWidthSlider") );
//  s->assignParameter( scope->getParameterByName("LineWidth") );
//  s->setSliderName("LineWidth");
//  s->setDescription("Width of lines");
//  s->setDescriptionField(infoField);
//  s->setStringConversionFunction(&valueToString3);
//  //s->addListener(this);  // may be needed when we have a line preview
//
//  addWidget( boxLineProfile = c = new RComboBox("LineProfileBox") );
//  c->assignParameter( scope->getParameterByName("LineProfile") );
//  c->setDescription("Select line profile");
//  c->setDescriptionField(infoField);
//}
//
//void PhaseScopeEditor2::resized()
//{
//  PhaseScopeEditor::resized();
//
//
//  // this is from the baseclass - but we should arrange the differently here
//  int x  = display.getRight() + 4;
//  int y  = getPresetSectionBottom() + 4;
//  int w  = getWidth() - x - 8;          // slider width
//  int h  = 16;                          // slider height
//  int dy = h-2;                           // vertical distance ("delta-y") between widgets
//  int w3 = (w-8)/3;
//  int x2 = x  + w3 + 4;
//  int x3 = x2 + w3 + 4;
//
//  // general scope controls:
//  sliderFrameRate     ->setBounds(x, y, w, h); y += dy;
//  sliderAfterglow     ->setBounds(x, y, w, h); y += dy;
//  sliderDecayByValue  ->setBounds(x, y, w, h); y += dy;
//  sliderDecayByAverage->setBounds(x, y, w, h); y += dy;
//  sliderPixelScale    ->setBounds(x, y, w, h); y += dy;
//
//  // dot controls:
//  y += 8;
//  buttonDrawDots ->setBounds(x,  y, w3, h); 
//  buttonBigDot   ->setBounds(x2, y, w3, h); 
//  buttonAntiAlias->setBounds(x3, y, w3, h); 
//  y += dy;
//  sliderBrightness   ->setBounds(x, y, w, h); y += dy;
//  sliderPixelSpread  ->setBounds(x, y, w, h); y += dy;
//  sliderLineDensity  ->setBounds(x, y, w, h); y += dy;
//  sliderDotLimit     ->setBounds(x, y, w, h); y += dy;
//  sliderDotSize      ->setBounds(x, y, w, h); y += dy;
//  sliderDotBlur      ->setBounds(x, y, w, h); y += dy;
//  sliderDotInnerSlope->setBounds(x, y, w, h); y += dy;
//  sliderDotOuterSlope->setBounds(x, y, w, h); y += dy;
//
//  // line controls:
//  y += 8;
//  buttonDrawLines     ->setBounds(x, y, w3, h); y += dy;
//  sliderLineBrightness->setBounds(x, y, w,  h); y += dy;
//  sliderLineWidth     ->setBounds(x, y, w,  h); y += dy;
//  boxLineProfile      ->setBounds(x, y, w,  h); y += dy;
//}
//
//void PhaseScopeEditor2::paint(Graphics& g)
//{
//  PhaseScopeEditor::paint(g);
//
//  RWidget *bottomWidget = boxLineProfile;
//  float x, y, w, h;
//  x = (float)bottomWidget->getX();
//  y = (float)bottomWidget->getBottom() + 8.f;
//  w = (float)dotPreviewImage.getWidth();
//  h = (float)dotPreviewImage.getHeight();
//  g.drawImage(dotPreviewImage, Rectangle<float>(x, y, w, h));
//}
//
//void PhaseScopeEditor2::rSliderValueChanged(RSlider* slider)
//{
//  updatePreviewDot();
//}
//
//void PhaseScopeEditor2::updatePreviewDot()
//{
//  dotPreviewMask.copyShapeParametersFrom(scope->phaseScopeBuffer->dotMask);
//  //normalizedDataToImage(dotPreviewMask.getPixelPointer(0, 0), dotPreviewImage); // gray values
//  normalizedDataToImage(dotPreviewMask.getPixelPointer(0, 0), dotPreviewImage, scope->colorMap);
//  repaint();
//}

