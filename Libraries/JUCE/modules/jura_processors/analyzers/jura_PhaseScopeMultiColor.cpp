PhaseScopeMultiColor::PhaseScopeMultiColor(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhaseScopeMultiColor";
  setActiveDirectory(getApplicationDirectory() + "/PhaseScopeMultiColorPresets");
  rainbow = false;

  pixelScale = 1.0;
  displayWidth  = 100;
  displayHeight = 100; 
  updateBufferSize();

  colorPeriod = 0.3;
  colorCounter = 0.0;
  createParameters();

  reset();
}

void PhaseScopeMultiColor::createParameters()
{
  ScopedLock scopedLock(*plugInLock);

  Parameter* p;

  p = new Parameter(plugInLock, "Brightness", 0.1, 10.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setBrightness);

  p = new Parameter(plugInLock, "AfterGlow", 0.001, 50.0, 0.0, 0.1, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setAfterGlow);

  p = new Parameter(plugInLock, "ColorPeriod", 0.0, 3.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(
    this, &PhaseScopeMultiColor::setColorCyclePeriod);

  p = new Parameter(plugInLock, "PixelSpread", 0.0, 1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setPixelSpread);

  p = new Parameter(plugInLock, "PixelScale", 1.0, 8.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setPixelScale);

  p = new Parameter(plugInLock, "LineDensity", 0.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setLineDensity);

  p = new Parameter(plugInLock, "FrameRate", 1.0, 100.0, 0.0, 25.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setFrameRate);

  p = new Parameter(plugInLock, "AntiAlias", 0.0, 1.0, 0.0, 1.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setAntiAlias);
  addObservedParameter(p);

  //p = new Parameter(plugInLock, "Rainbow", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  //p->setValueChangeCallback<PhaseScopeMultiColor>(this, &PhaseScopeMultiColor::setRainbowMode);
  //addObservedParameter(p);

  //resetParametersToDefaultValues();
}

void PhaseScopeMultiColor::setDisplayPixelSize(int width, int height)
{
  ScopedLock scopedLock(*plugInLock);
  displayWidth  = width;
  displayHeight = height; 
  updateBufferSize();
}

void PhaseScopeMultiColor::setBrightness(double newBrightness)
{
  brightness = newBrightness;
  //phaseScopeBuffer.setBrightness((float)newBrightness);
}
void PhaseScopeMultiColor::setAfterGlow(double newGlow)
{
  phaseScopeBuffer.setDecayTime(newGlow);
}
void PhaseScopeMultiColor::setColorCyclePeriod(double newPeriod)
{
  colorPeriod = newPeriod;
}
void PhaseScopeMultiColor::setLineDensity(double newDensity)
{
  phaseScopeBuffer.setLineDensity((float)newDensity);
}
void PhaseScopeMultiColor::setPixelSpread(double newSpread)
{
  phaseScopeBuffer.setPixelSpread((float)newSpread);
}
void PhaseScopeMultiColor::setPixelScale(double newFactor)
{
  jassert(newFactor >= 1.0);
  pixelScale = newFactor;
  updateBufferSize();
}
void PhaseScopeMultiColor::setAntiAlias(bool shouldAntiAlias)
{
  phaseScopeBuffer.setAntiAlias(shouldAntiAlias);
}
//void PhaseScopeMultiColor::setRainbowMode(bool shouldUseRainbowColors)
//{
//  rainbow = shouldUseRainbowColors;
//}
void PhaseScopeMultiColor::setFrameRate(double newRate)
{
  phaseScopeBuffer.setFrameRate(newRate);
  updateRepaintInterval();
}

AudioModuleEditor* PhaseScopeMultiColor::createEditor()
{
  return new PhaseScopeMultiColorEditor(this);
}

void PhaseScopeMultiColor::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  ScopedLock scopedLock(*plugInLock);
  jassert(numChannels == 2);

  for(int n = 0; n < numSamples; n++)
  {
    // preliminary, very inefficient (optimize later and maybe factor out):
    float b       = float(brightness);
    if(colorPeriod != 0.0)
    {
      double inc    = 1.0 / (colorPeriod * phaseScopeBuffer.getSampleRate()); // color increment
      colorCounter += inc;
      colorCounter  = fmod(colorCounter, 1.0);
      float red     = float(0.5f + 0.5f*sin(2*PI*colorCounter));
      float green   = float(0.5f + 0.5f*sin(2*PI*colorCounter + PI*120/180));  // 120° phase shift
      float blue    = float(0.5f + 0.5f*sin(2*PI*colorCounter + PI*240/180));  // 240° phase shift

      phaseScopeBuffer.setBrightness(Float32x4(b*red, b*green, b*blue, 1.f));
    }
    else
      phaseScopeBuffer.setBrightness(Float32x4(b, b, b, 1.f));


    phaseScopeBuffer.bufferSampleFrame(inOutBuffer[0][n], inOutBuffer[1][n]);
  }

  repaintCounter += numSamples;
  while(repaintCounter > repaintIntervalInSamples)
  {
    updateScopeImage();
    sendImageUpdateNotification(&image);
    repaintCounter -= repaintIntervalInSamples;
  }
}

void PhaseScopeMultiColor::setSampleRate(double newSampleRate)
{
  ScopedLock scopedLock(*plugInLock);
  phaseScopeBuffer.setSampleRate(newSampleRate);
  updateRepaintInterval();
}

void PhaseScopeMultiColor::reset()
{
  ScopedLock scopedLock(*plugInLock);
  phaseScopeBuffer.reset();
  colorCounter = 0.0;
  repaintCounter = 0;
}

void PhaseScopeMultiColor::updateBufferSize()
{
  ScopedLock scopedLock(*plugInLock);
  int w = (int) round(displayWidth  / pixelScale);
  int h = (int) round(displayHeight / pixelScale);
  phaseScopeBuffer.setSize(w, h);
  image = juce::Image(juce::Image::ARGB, w, h, false);
}

void PhaseScopeMultiColor::updateScopeImage()
{
  //RAPT::Float32x4 *tmp = phaseScopeBuffer.getDataMatrix()[0]; // old
  RAPT::Float32x4 *tmp = phaseScopeBuffer.getImage()->getPixelPointer(0, 0);
  float *floatData = reinterpret_cast<float*>(tmp);
  dataToImageOpaqueFloat32x4(floatData, image);
  phaseScopeBuffer.applyPixelDecay();
}

void PhaseScopeMultiColor::updateRepaintInterval()
{
  repaintIntervalInSamples = 
    (int) round(phaseScopeBuffer.getSampleRate() / phaseScopeBuffer.getFrameRate());
  repaintCounter = 0;
}

//=================================================================================================

PhaseScopeMultiColorDisplay::PhaseScopeMultiColorDisplay(
  jura::PhaseScopeMultiColor *newPhaseScopeToEdit)
{
  phaseScope = newPhaseScopeToEdit;
  phaseScope->addImageUpdateListener(this);
  addChangeListener(this);
}

PhaseScopeMultiColorDisplay::~PhaseScopeMultiColorDisplay()
{
  phaseScope->removeImageUpdateListener(this);
}

void PhaseScopeMultiColorDisplay::resized()
{
  phaseScope->setDisplayPixelSize(getWidth(), getHeight());
}

void PhaseScopeMultiColorDisplay::paint(Graphics &g)
{
  //ScopedLock scopedLock(*(phaseScope->plugInLock));

  g.setImageResamplingQuality(Graphics::lowResamplingQuality);
  //g.setImageResamplingQuality(Graphics::mediumResamplingQuality);
  //g.setImageResamplingQuality(Graphics::highResamplingQuality);
  g.drawImage(phaseScope->image, Rectangle<float>(0.f, 0.f, (float) getWidth(), 
    (float) getHeight()));
}

void PhaseScopeMultiColorDisplay::imageWasUpdated(juce::Image* image)
{
  sendChangeMessage();
  // We will receive the change message ourselves and when we do, we will trigger a repaint. We 
  // can't directly call repaint here because then we hit an assertion which says we should acquire
  // a MessageManagerLock - but when we do so, it becomes unresponsive.

  // this doesn't work:
  //const MessageManagerLock mmLock;
  //repaint();
}

void PhaseScopeMultiColorDisplay::changeListenerCallback(ChangeBroadcaster *source)
{
  repaint();
}

//=================================================================================================

PhaseScopeMultiColorEditor::PhaseScopeMultiColorEditor(
  jura::PhaseScopeMultiColor *newPhaseScopeToEdit)
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

void PhaseScopeMultiColorEditor::createWidgets()
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

  addWidget( sliderColorPeriod = s = new RSlider("ColorPeriodSlider") );
  s->assignParameter( scope->getParameterByName("ColorPeriod") );
  s->setSliderName("ColorPeriod");
  s->setDescription("Periodicity of the hue cycling time in seconds");
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

  //addWidget( buttonRainbow = b = new RButton("Rainbow") );
  //b->assignParameter( scope->getParameterByName("Rainbow") );
  //b->setDescription("Rainbow color rotation");
  //b->setDescriptionField(infoField);
}

void PhaseScopeMultiColorEditor::resized()
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
  sliderColorPeriod->setBounds(x, y, w, h); y += dy;
  sliderPixelSpread->setBounds(x, y, w, h); y += dy;
  sliderPixelScale ->setBounds(x, y, w, h); y += dy;
  sliderLineDensity->setBounds(x, y, w, h); y += dy;
  sliderFrameRate  ->setBounds(x, y, w, h); y += dy;

  buttonAntiAlias->setBounds(x, y,   w, h); y += dy;
  //buttonRainbow  ->setBounds(x, y,   w, h); y += dy;
}
