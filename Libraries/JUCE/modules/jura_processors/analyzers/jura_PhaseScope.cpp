PhaseScope::PhaseScope(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhaseScope";
  setActiveDirectory(getApplicationDirectory() + "/PhaseScopePresets");

  pixelScale = 1.0;
  displayWidth  = 100;
  displayHeight = 100; 
  updateBufferSize();

  createParameters();
  reset();

  //colorMap.setDefaultMap(ColorMap::fire);
  colorMap.setDefaultMap(ColorMap::ice);

  //juce::ColourGradient g;
  //g.addColour(0.0, Colour(  0,   0,   0));
  //g.addColour(0.2, Colour(  0,   0, 255));
  //g.addColour(0.4, Colour(  0, 255, 255));
  //g.addColour(0.6, Colour(  0, 255,   0));
  //g.addColour(0.8, Colour(255, 255,   0));
  //g.addColour(1.0, Colour(255,   0,   0)); // maybe add magenta as last
  //colorMap.setFromColourGradient(g);
}

void PhaseScope::createParameters()
{
  ScopedLock scopedLock(*plugInLock);

  Parameter* p;

  p = new Parameter(plugInLock, "Brightness", 0.1, 100.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setBrightness);

  p = new Parameter(plugInLock, "AfterGlow", 0.001, 50.0, 0.0, 0.1, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setAfterGlow);

  p = new Parameter(plugInLock, "PixelSpread", 0.0, 1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setPixelSpread);

  p = new Parameter(plugInLock, "PixelScale", 1.0, 8.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setPixelScale);

  p = new Parameter(plugInLock, "LineDensity", 0.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setLineDensity);

  p = new Parameter(plugInLock, "FrameRate", 1.0, 100.0, 0.0, 25.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setFrameRate);

  p = new Parameter(plugInLock, "AntiAlias", 0.0, 1.0, 0.0, 1.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setAntiAlias);
  addObservedParameter(p);
}

void PhaseScope::setDisplayPixelSize(int width, int height)
{
  ScopedLock scopedLock(*plugInLock);
  displayWidth  = width;
  displayHeight = height; 
  updateBufferSize();
}

void PhaseScope::setBrightness(double newBrightness)
{
  phaseScopeBuffer.setBrightness((float)newBrightness);
}
void PhaseScope::setAfterGlow(double newGlow)
{
  phaseScopeBuffer.setDecayTime(newGlow);
}
void PhaseScope::setLineDensity(double newDensity)
{
  phaseScopeBuffer.setLineDensity((float)newDensity);
}
void PhaseScope::setPixelSpread(double newSpread)
{
  phaseScopeBuffer.setPixelSpread((float)newSpread);
}
void PhaseScope::setPixelScale(double newFactor)
{
  jassert(newFactor >= 1.0);
  pixelScale = newFactor;
  updateBufferSize();
}
void PhaseScope::setAntiAlias(bool shouldAntiAlias)
{
  phaseScopeBuffer.setAntiAlias(shouldAntiAlias);
}
void PhaseScope::setFrameRate(double newRate)
{
  phaseScopeBuffer.setFrameRate(newRate);
  updateRepaintInterval();
}

AudioModuleEditor* PhaseScope::createEditor()
{
  return new PhaseScopeEditor(this);
}

void PhaseScope::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(numChannels == 2);

  for(int n = 0; n < numSamples; n++)
    phaseScopeBuffer.bufferSampleFrame(inOutBuffer[0][n], inOutBuffer[1][n]);

  repaintCounter += numSamples;
  while(repaintCounter > repaintIntervalInSamples)
  {
    updateScopeImage();
    sendImageUpdateNotification(&image);
    repaintCounter -= repaintIntervalInSamples;
  }
}

void PhaseScope::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);
  phaseScopeBuffer.setSampleRate(newSampleRate);
  updateRepaintInterval();
}

void PhaseScope::reset()
{
  ScopedLock scopedLock(*plugInLock);
  phaseScopeBuffer.reset();
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
  ScopedLock scopedLock(*plugInLock);
  int w = (int) round(displayWidth  / pixelScale);
  int h = (int) round(displayHeight / pixelScale);
  phaseScopeBuffer.setSize(w, h);
  image = juce::Image(juce::Image::ARGB, w, h, false);
}

void PhaseScope::updateScopeImage()
{
  normalizedDataToImage(phaseScopeBuffer.getImage()->getPixelPointer(0, 0), image, colorMap);
  phaseScopeBuffer.applyPixelDecay();
}

void PhaseScope::updateRepaintInterval()
{
  repaintIntervalInSamples = 
    (int) round(phaseScopeBuffer.getSampleRate() / phaseScopeBuffer.getFrameRate());
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
  //ScopedLock scopedLock(*(phaseScope->plugInLock));

  g.setImageResamplingQuality(Graphics::lowResamplingQuality);
  //g.setImageResamplingQuality(Graphics::mediumResamplingQuality);
  //g.setImageResamplingQuality(Graphics::highResamplingQuality);
  g.drawImage(phaseScope->image, Rectangle<float>(0.f, 0.f, (float) getWidth(), 
    (float) getHeight()));
}

void PhaseScopeDisplay::imageWasUpdated(juce::Image* image)
{
  sendChangeMessage();
  // We will receive the change message ourselves and when we do, we will trigger a repaint. We 
  // can't directly call repaint here because then we hit an assertion which says we should acquire
  // a MessageManagerLock - but when we do so, it becomes unresponsive.

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
  ScopedLock scopedLock(*plugInLock);
  scope = newPhaseScopeToEdit;
  widgetMargin = 150; 

  addAndMakeVisible(display);
  createWidgets();

  int headerMargin = 26;  // this is the height we need for headline and preset-section
  setSize(400+widgetMargin, 400+headerMargin);
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

  addWidget( sliderFrameRate = s = new RSlider("FrameRateSlider") );
  s->assignParameter( scope->getParameterByName("FrameRate") );
  s->setSliderName("FrameRate");
  s->setDescription("Frame rate for redrawing");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( buttonAntiAlias = b = new RButton("AntiAlias") );
  b->assignParameter( scope->getParameterByName("AntiAlias") );
  b->setDescription("Anti aliased drawing (bilinear deinterpolation)");
  b->setDescriptionField(infoField);
}

void PhaseScopeEditor::resized()
{
  ScopedLock scopedLock(*plugInLock);
  AudioModuleEditor::resized();

  int w = getWidth();
  int h = getHeight();
  int y = getPresetSectionBottom() + 4;

  display.setBounds(0, y, w-widgetMargin, h-y);

  // set up widgets:
  int x = display.getRight() + 4;
  w = getWidth() - x - 8;          // slider width
  h = 16;                          // slider height
  int dy = h+4;                    // vertical distance ("delta-y") between widgets

  sliderBrightness ->setBounds(x, y, w, h); y += dy;
  sliderAfterglow  ->setBounds(x, y, w, h); y += dy;
  sliderPixelSpread->setBounds(x, y, w, h); y += dy;
  sliderPixelScale ->setBounds(x, y, w, h); y += dy;
  sliderLineDensity->setBounds(x, y, w, h); y += dy;
  sliderFrameRate  ->setBounds(x, y, w, h); y += dy;
  buttonAntiAlias  ->setBounds(x, y, w, h); y += dy;
}


//=================================================================================================
// the artistically entended version of the PhaseScope:

PhaseScope2::PhaseScope2(CriticalSection *lockToUse) : PhaseScope(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  createParameters();  // creates the additional parameters
}

void PhaseScope2::setDotSize(double newSize)
{

}

void PhaseScope2::setDotBlur(double newBlur)
{

}

void PhaseScope2::createParameters()
{
  ScopedLock scopedLock(*plugInLock);
  Parameter* p;

  p = new Parameter(plugInLock, "DotSize", 1.0, 10.0, 0.0, 2.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDotSize);

  p = new Parameter(plugInLock, "DotBlur", 0.0, 1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope2>(this, &PhaseScope2::setDotBlur);
}

AudioModuleEditor* PhaseScope2::createEditor()
{
  return new PhaseScopeEditor2(this);
}

//-------------------------------------------------------------------------------------------------

PhaseScopeEditor2::PhaseScopeEditor2(jura::PhaseScope2 *newPhaseScopeToEdit)
  : PhaseScopeEditor(newPhaseScopeToEdit)
{
  createWidgets();
  resized();
}

void PhaseScopeEditor2::createWidgets()
{
  RSlider *s;

  addWidget( sliderDotSize = s = new RSlider("DotSizeSlider") );
  s->assignParameter( scope->getParameterByName("DotSize") );
  s->setSliderName("DotSize");
  s->setDescription("Dot size in pixels");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);

  addWidget( sliderDotBlur = s = new RSlider("DotBlurSlider") );
  s->assignParameter( scope->getParameterByName("DotBlur") );
  s->setSliderName("DotBlur");
  s->setDescription("Dot blur from 0..1");
  s->setDescriptionField(infoField);
  s->setStringConversionFunction(&valueToString3);
}

void PhaseScopeEditor2::resized()
{
  PhaseScopeEditor::resized();

  // set up additional widgets:
  int x, y, w, h, dy;
  x  = display.getRight() + 4;
  w  = buttonAntiAlias->getWidth();   // slider width
  h  = 16;                            // slider height
  dy = h+4;                           // vertical distance ("delta-y") between widgets
  y  = buttonAntiAlias->getBottom() + dy;
  
  sliderDotSize->setBounds(x, y, w, h); y += dy;
  sliderDotBlur->setBounds(x, y, w, h); y += dy;
}

