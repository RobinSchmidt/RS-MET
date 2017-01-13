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
}

void PhaseScope::createParameters()
{
  ScopedLock scopedLock(*plugInLock);

  Parameter* p;

  p = new Parameter(plugInLock, "Brightness", 0.1, 10.0, 0.0, 1.0, Parameter::EXPONENTIAL);
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

void PhaseScope::updateBufferSize()
{
  ScopedLock scopedLock(*plugInLock);
  int w = (int) round(displayWidth  / pixelScale);
  int h = (int) round(displayHeight / pixelScale);
  phaseScopeBuffer.setSize(w, h);
  image = juce::Image(juce::Image::ARGB, w, h, false);
}

// move to GraphicsTools
void normalizedDataToImageGrayScale(float *data, juce::Image &image)
{
  juce::Image::BitmapData bitmap(image, juce::Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);

  // indices for RGBA components in target image:
  int ri, gi, bi, ai;
  colorComponentIndices(image, ri, gi, bi, ai);

  uint8 *p = bitmap.getPixelPointer(0, 0);
  for(int i = 0; i < bitmap.height * bitmap.width; i++)
  {
    p[ri] = p[gi] = p[bi] = (uint8)(255 * data[i]);
    p[ai] = 255;  // full opacity
    p += 4;
  }
}
void normalizedDataToImage(float *data, juce::Image &image, juce::ColourGradient *colorMap = nullptr)
{
  if(colorMap == nullptr)
  {
    normalizedDataToImageGrayScale(data, image);
    return;
  }

  int w = image.getWidth();
  int h = image.getHeight();
  for(int y = 0; y < h; y++)
  {
    for(int x = 0; x < w; x++)
      image.setPixelAt(x, y, colorMap->getColourAtPosition(data[y*w+x]));
  }
}

void PhaseScope::updateScopeImage()
{
  // test color gradient (\todo: make this a member - maybe as pointer and have a class 
  // ColorGradientSelector or something):
  juce::ColourGradient gradient;
  //gradient.addColour(0.0, Colour(  0,   0,   0));
  //gradient.addColour(0.2, Colour(  0,   0, 255));
  //gradient.addColour(0.4, Colour(  0, 255, 255));
  //gradient.addColour(0.6, Colour(  0, 255,   0));
  //gradient.addColour(0.8, Colour(255, 255,   0));
  //gradient.addColour(1.0, Colour(255,   0,   0)); // maybe add magenta as last

  //gradient.addColour(0.0, Colour(  0,   0,   0));
  //gradient.addColour(0.2, Colour(  0,   0, 255));
  //gradient.addColour(0.4, Colour(  0, 127, 127));
  //gradient.addColour(0.6, Colour(  0, 255,   0));
  //gradient.addColour(0.8, Colour(255, 255,   0));
  //gradient.addColour(1.0, Colour(255, 255, 255));

  // "Sun"
  gradient.addColour(0.0, Colour(  0,   0,   0));  // black
  gradient.addColour(0.4, Colour(255,   0,   0));  // red
  gradient.addColour(0.6, Colour(255, 255,   0));  // yellow
  gradient.addColour(1.0, Colour(255, 255, 255));  // white

  //// "Ice"
  //gradient.addColour(0.0, Colour(  0,   0,   0));  // black
  //gradient.addColour(0.4, Colour(  0,   0, 255));  // blue
  //gradient.addColour(0.6, Colour(  0, 255, 255));  // cyan
  //gradient.addColour(1.0, Colour(255, 255, 255));  // white

  //// "Grass"
  //gradient.addColour(0.0, Colour(  0,   0,   0));  // black
  //gradient.addColour(0.4, Colour(  0, 255,   0));  // green
  //gradient.addColour(0.6, Colour(255, 255,   0));  // yellow
  //gradient.addColour(1.0, Colour(255, 255, 255));  // white

  // here's a good tutorial about colormaps:
  //http://www.research.ibm.com/people/l/lloydt/color/color.HTM



  //normalizedDataToImage(phaseScopeBuffer.getDataMatrix()[0], image);
  normalizedDataToImage(phaseScopeBuffer.getDataMatrix()[0], image, &gradient);
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
  widgetMargin = 120; 

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
