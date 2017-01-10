
PhaseScope::PhaseScope(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhaseScope";
  setActiveDirectory(getApplicationDirectory() + "/PhaseScopePresets");
  //needsPixelDecay = false;
  rainbow = false;

  pixelScale = 1.0;
  displayWidth  = 100;
  displayHeight = 100; 
  updateBufferSize();

  colorPeriod = 5.0;
  colorCounter = 0.0;
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

  p = new Parameter(plugInLock, "Rainbow", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setRainbowMode);
  addObservedParameter(p);

  //resetParametersToDefaultValues();
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
void PhaseScope::setRainbowMode(bool shouldUseRainbowColors)
{
  rainbow = shouldUseRainbowColors;
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

  //// this may be obsolete soon:
  //if(needsPixelDecay)
  //{
  //  phaseScopeBuffer.applyPixelDecay();
  //  needsPixelDecay = false;
  //}

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
  colorCounter = 0.0;
  repaintCounter = 0;
}

Colour PhaseScope::getAndUpdateColor()
{
  ScopedLock scopedLock(*plugInLock);

  if(!rainbow)
    return Colours::white;

  ColourAHSL colorAHSL((float)colorCounter, 1.f, 0.5f, 0.5f);

  double inc    = 1.0 / (colorPeriod * phaseScopeBuffer.getFrameRate()); // color increment
  colorCounter += inc;
  colorCounter  = fmod(colorCounter, 1.0);

  return colorAHSL.getAsJuceColour();
}

void PhaseScope::updateBufferSize()
{
  ScopedLock scopedLock(*plugInLock);
  int w = (int) round(displayWidth  / pixelScale);
  int h = (int) round(displayHeight / pixelScale);
  phaseScopeBuffer.setSize(w, h);
  image = juce::Image(juce::Image::ARGB, w, h, false);
}

// \todo move these helper functions to jura_GraphicsTools, maybe templatize so it can be used 
// for double also:
void dataMatrixToPixelBrightnessGray(float **data, uint8 *pixels, int width, int height)
{
  uint8 *p = pixels;
  for(int i = 0; i < height; i++)     // loop over lines
  {
    for(int j = 0; j < width; j++)    // loop over pixels
    {
      // we assume here, that the alpha channel comes last in the byte order of the pixels
      p[0] = p[1] = p[2] = (uint8) (255 * data[i][j]);  // data determines white-value
      p[3] = 255;                                       // set to full opacity ("alpha")
      p   += 4;                                         // jump to next pixel
    }
  }
}
void dataMatrixToPixelBrightness(float **data, uint8 *pixels, int width, int height,
  uint8 red = 255, uint8 green = 255, uint8 blue = 255)
{
  uint8 *p = pixels;
  for(int i = 0; i < height; i++)     // loop over lines
  {
    for(int j = 0; j < width; j++)    // loop over pixels
    {
      // we assume here, that the byte order of the pixels is BGRA (seems to be true on PC)
      p[0] = (uint8) (blue  * data[i][j]);
      p[1] = (uint8) (green * data[i][j]);
      p[2] = (uint8) (red   * data[i][j]);
      p[3] = 255; // full opacity ("alpha")
      p   += 4;
    }
  }
  // todo: write a version of this function that uses a colormap
}
void dataMatrixToImage(float **data, juce::Image &image, 
  uint8 red = 255, uint8 green = 255, uint8 blue = 255)
{
  // We assume here that the size of the image (width and height) matches the dimensions of the 
  // data matrix)
  juce::Image::BitmapData bitmap(image, juce::Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);
  uint8 *pixelPointer = bitmap.getPixelPointer(0, 0);
  if(red == green && green == blue)
    dataMatrixToPixelBrightnessGray(data, pixelPointer, bitmap.width, bitmap.height);
  else
    dataMatrixToPixelBrightness(data, pixelPointer, bitmap.width, bitmap.height, red, green, blue);
}

void PhaseScope::updateScopeImage()
{
  // old (without color change):
  //dataMatrixToImage(phaseScopeBuffer.getDataMatrix(), image);
  //dataMatrixToImage(phaseScopeBuffer.getDataMatrix(), image, 0, 0, 255);

  Colour c = getAndUpdateColor();
  dataMatrixToImage(phaseScopeBuffer.getDataMatrix(), image, 
    c.getRed(), c.getGreen(), c.getBlue());
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

  //setSize(400, 400);
  //startTimerHz(25);
  //phaseScope->phaseScopeBuffer.setFrameRate(25);
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
  //phaseScope->updateScopeImage();

  //ScopedLock scopedLock(*(phaseScope->plugInLock));

  g.setImageResamplingQuality(Graphics::lowResamplingQuality);
  //g.setImageResamplingQuality(Graphics::mediumResamplingQuality);
  //g.setImageResamplingQuality(Graphics::highResamplingQuality);
  g.drawImage(phaseScope->image, Rectangle<float>(0.f, 0.f, (float) getWidth(), 
    (float) getHeight()));

  // maybe all of this should be refactored in way, such that we need to write only
  // g.drawImageAt(phaseScope->getImage(), 0, 0);  here - we could take care of proper thread 
  // synchronization there also

  //// this is the old, direct and horribly inefficient way to do it:
  //g.fillAll(Colours::black);
  //for(int x = 0; x < getWidth(); x++)
  //{
  //  for(int y = 0; y < getHeight(); y++)
  //  {
  //    g.setColour(phaseScope->getColourAt(x, y));
  //    g.setPixel(x, y);
  //  }
  //}
}

//void PhaseScopeDisplay::timerCallback()
//{
//  repaint();
//  phaseScope->phaseScopeBuffer.applyPixelDecay();
//  //phaseScope->triggerPixelDecay();
//
//  //startTimerHz(phaseScope->getFrameRate());
//  startTimer((int) (1000.0 / phaseScope->getFrameRate())); // fps to ms
//
//  // there's a kind of bug - the display flickers and the flickering seems to be different if
//  // we apply the pixel decay in the GUI thread (by calling 
//  // phaseScope->phaseScopeBuffer.applyPixelDecay() here) or apply it in the audio thread
//  // (by calling phaseScope->triggerPixelDecay()). maybe we need to do some kind of thread
//  // synchronization
//}

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

  addWidget( buttonRainbow = b = new RButton("Rainbow") );
  b->assignParameter( scope->getParameterByName("Rainbow") );
  b->setDescription("Rainbow color rotation");
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

  buttonAntiAlias->setBounds(x, y,   w, h); y += dy;
  buttonRainbow  ->setBounds(x, y,   w, h); y += dy;
}
