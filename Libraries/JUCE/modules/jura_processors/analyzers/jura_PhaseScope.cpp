
PhaseScopeBuffer::PhaseScopeBuffer()
{
  antiAlias   = true;
  frameRate   = 25.0;
  decayTime   = 0.5;
  lineDensity = 0.0f;
  thickness   = 0.707f;
  brightness  = 1.0;
  updateDecayFactor();

  bufferFlat = nullptr;
  buffer     = nullptr;
  width      = 200;
  height     = 200;
  allocateBuffer();

  setSampleRate(44100.0);
  reset();
}

PhaseScopeBuffer::~PhaseScopeBuffer()
{
  freeBuffer();
}

void PhaseScopeBuffer::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  updateInsertFactor();
}

void PhaseScopeBuffer::setFrameRate(double newFrameRate)
{
  frameRate = newFrameRate;
  updateDecayFactor();
}

void PhaseScopeBuffer::setBrightness(float newBrightness)
{
  brightness = newBrightness;
  updateInsertFactor();
}

void PhaseScopeBuffer::setDecayTime(double newDecayTime)
{
  decayTime = newDecayTime;
  updateDecayFactor();
}

void PhaseScopeBuffer::setSize(int newWidth, int newHeight)
{
  if(newWidth != width || newHeight != height)
  {
    freeBuffer();
    width  = newWidth;
    height = newHeight;
    allocateBuffer();
    reset();           // avoid flickering at size-change
  }
}

void PhaseScopeBuffer::setAntiAlias(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
}

void PhaseScopeBuffer::setLineDensity(float newDensity)
{
  lineDensity = newDensity;
}

void PhaseScopeBuffer::setPixelSpread(float newSpread)
{
  thickness = newSpread;
}

void PhaseScopeBuffer::convertAmplitudesToMatrixIndices(double &x, double &y)
{
  x  = 0.5*(x+1);  // convert -1..+1 into 0..1
  y  = 0.5*(y+1);
  x *= width;
  y *= height;
  // maybe we should add 0.5 after multiplication by width/height?
}

void PhaseScopeBuffer::bufferSampleFrame(double x, double y)
{
  convertAmplitudesToMatrixIndices(x, y);
  addLineTo((float)x, (float)y);
}

void PhaseScopeBuffer::applyPixelDecay()
{
  ArrayTools::rsScale(bufferFlat, width*height, decayFactor);
}

void PhaseScopeBuffer::reset()
{
  ArrayTools::rsFillWithZeros(bufferFlat, width*height);

  // (xOld,yOld) = (0,0) - but in pixel coordinates:
  xOld = 0.5f*width;   
  yOld = 0.5f*height;
}

//// get rid of this function:
//float PhaseScopeBuffer::pixelDistance(float x1, float y1, float x2, float y2)
//{
//  float dx = x2-x1;
//  float dy = y2-y1;
//  return sqrt(dx*dx + dy*dy);   // Euclidean distance
//  //return max(fabs(dx), fabs(dy)); // maybe try Euclidean distance instead
//}

void PhaseScopeBuffer::addLineTo(float x, float y)
{
  if(lineDensity == 0.f)
  {
    addDot(x, y, insertFactor);
    return;
  }

  float dx = x-xOld;
  float dy = y-yOld;
  float pixelDistance = sqrt(dx*dx + dy*dy);
  int   numDots = max(1, (int)floor(lineDensity*pixelDistance));
  float intensity = (float) (insertFactor/numDots);
  float scaler = (float)(1.0 / numDots);
  float k;
  for(int i = 1; i <= numDots; i++)
  {
    k = scaler * i;  // == i / numDots
    //addDot((1-k)*xOld + k*x, (1-k)*yOld + k*y, intensity); // unoptimized
    addDot(xOld + k*dx, yOld + k*dy, intensity);
  }
  xOld = x;
  yOld = y;
}

void PhaseScopeBuffer::addDot(float x, float y, float intensity)
{
  if(!antiAlias)
  {
    addDotFast(x, y, intensity);
    return;
  }

  int j = (int)floor(x);  // integer part of x
  int i = (int)floor(y);  // integer part of y
  x -= j;                 // fractional part of x
  y -= i;                 // fractional part of y

  // compute weights for bilinear deinterpolation (maybe factor out):
  float a, b, c, d;
  d = x*y;
  c = y-d;
  b = x-d;
  a = 1+d-x-y;

  // compute values to accumulate into the 4 pixels at (i,j), (i+1,j), (i,j+1), (i+1,j+1):
  a *= intensity;  // (i,   j)
  b *= intensity;  // (i+1, j+1)
  c *= intensity;  // (i+1, j)
  d *= intensity;  // (i+1, j+1)

  // accumulate values into the matrix:
  if(j >= 0 && j < width-1 && i >= 0 && i < height-1)
  {
    accumulate(buffer[i]  [j],   a);
    accumulate(buffer[i]  [j+1], b);
    accumulate(buffer[i+1][j],   c);
    accumulate(buffer[i+1][j+1], d);
  }

  // apply thickness:
  if(thickness > 0.f && j >= 1 && j < width-2 && i >= 1 && i < height-2)
  {
    float t, s, sa, sb, sc, sd, ta, tb, tc, td;
    t = thickness;             // weight for direct neighbour pixels
    s = t * (float)SQRT2_INV;  // weight for diagonal neighbour pixels

    sa = s*a;
    sb = s*b;
    sc = s*c;
    sd = s*d;

    ta = t*a;
    tb = t*b;
    tc = t*c;
    td = t*d;

    accumulate(buffer[i-1][j-1], sa);
    accumulate(buffer[i-1][j],   ta+sb);
    accumulate(buffer[i-1][j+1], tb+sa);
    accumulate(buffer[i-1][j+2], sb);

    accumulate(buffer[i]  [j-1], ta+sc);
    accumulate(buffer[i]  [j],   sd+tb+tc);
    accumulate(buffer[i]  [j+1], sc+ta+td);
    accumulate(buffer[i]  [j+2], tb+sd);

    accumulate(buffer[i+1][j-1], tc+sa);
    accumulate(buffer[i+1][j],   sb+ta+td);
    accumulate(buffer[i+1][j+1], sa+tb+tc);
    accumulate(buffer[i+1][j+2], td+sb);

    accumulate(buffer[i+2][j-1], sc);
    accumulate(buffer[i+2][j],   tc+sd);
    accumulate(buffer[i+2][j+1], td+sc);
    accumulate(buffer[i+2][j+2], sd);
  }
}

void PhaseScopeBuffer::addDotFast(float x, float y, float intensity)
{
  int j = (int)round(x);
  int i = (int)round(y);

  if(j >= 0 && j < width && i >= 0 && i < height)
    accumulate(buffer[i][j], intensity);

  // apply thickness:
  if(thickness > 0.f && j >= 1 && j < width-1 && i >= 1 && i < height-1)
  {
    float a, ta, sa;
    a  = intensity;
    ta = a  * thickness;
    sa = ta * (float)SQRT2_INV;

    accumulate(buffer[i-1][j-1], sa);
    accumulate(buffer[i-1][j],   ta);
    accumulate(buffer[i-1][j+1], sa);

    accumulate(buffer[i]  [j-1], ta);
    accumulate(buffer[i]  [j+1], ta);

    accumulate(buffer[i+1][j-1], sa);
    accumulate(buffer[i+1][j],   ta);
    accumulate(buffer[i+1][j+1], sa);
  }
}

void PhaseScopeBuffer::allocateBuffer()
{
  bufferFlat = new float[width*height];
  buffer = new float*[height];
  for(int i = 0; i < height; i++)
    buffer[i] = &bufferFlat[i*width]; 

  //buffer = new float*[width];
  //for(int i = 0; i < width; i++)
  //  buffer[i] = &bufferFlat[i*height];
  
  // We have the second index running over the y-coordinate and the first index running over the 
  // x-coordinate here...this is actually the other way around than in images where we have the 
  // first index pointing to a particular horizontal line. Maybe we should use a memory layout here
  // that corresponds to images (x-coordinate as 2nd index) ..i have to verify/modify all formulas 
  // again...i think, we just need to replace i with j and vice versa in all the "accumulate" calls

  // BUG:
  // there's definitely still something wrong somewhere - currently, the lines wrap around when the 
  // display is not square - figure it out and fix it - maybe make the GUI resizable before that, 
  // so it's better for testing in various size settings
}

void PhaseScopeBuffer::freeBuffer()
{
  delete[] buffer;
  delete[] bufferFlat;
  buffer = nullptr;
  bufferFlat = nullptr;
}

void PhaseScopeBuffer::updateDecayFactor()
{
  decayFactor = (float) exp(-1 / (decayTime*frameRate));
}

void PhaseScopeBuffer::updateInsertFactor()
{
  insertFactor = (float) (5000*brightness / sampleRate);
    // The 5000 factor is totally ad-hoc - maybe come up with some more meaningful factor. 
    // However, the proportionality to the birghtness parameter and inverse proportionality to 
    // the sample rate seems to make sense.
}

//=================================================================================================

PhaseScope::PhaseScope(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhaseScope";
  setActiveDirectory(getApplicationDirectory() + "/PhaseScopePresets");
  needsPixelDecay = false;
  rainbow = false;
  colorPeriod = 5.0;
  colorCounter = 0.0;
  createParameters();
}

void PhaseScope::createParameters()
{
  ScopedLock scopedLock(*plugInLock);

  Parameter* p;

  p = new Parameter(plugInLock, "Brightness", 0.1, 2.0, 0.0, 1.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setBrightness);

  p = new Parameter(plugInLock, "AfterGlow", 0.001, 10.0, 0.0, 0.1, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setAfterGlow);

  p = new Parameter(plugInLock, "PixelSpread", 0.0, 1.0, 0.0, 0.5, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setPixelSpread);

  p = new Parameter(plugInLock, "LineDensity", 0.0, 1.0, 0.0, 0.0, Parameter::LINEAR);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setLineDensity);

  p = new Parameter(plugInLock, "FrameRate", 1.0, 50.0, 0.0, 25.0, Parameter::EXPONENTIAL);
  addObservedParameter(p);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setFrameRate);

  p = new Parameter(plugInLock, "AntiAlias", 0.0, 1.0, 0.0, 1.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setAntiAlias);
  addObservedParameter(p);

  p = new Parameter(plugInLock, "Rainbow", 0.0, 1.0, 0.0, 0.0, Parameter::BOOLEAN);
  p->setValueChangeCallback<PhaseScope>(this, &PhaseScope::setRainbowMode);
  addObservedParameter(p);
}

void PhaseScope::setPixelSize(int width, int height)
{
  phaseScopeBuffer.setSize(width, height);
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
}


AudioModuleEditor* PhaseScope::createEditor()
{
  return new PhaseScopeEditor(this);
}

void PhaseScope::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
  if(needsPixelDecay)
  {
    phaseScopeBuffer.applyPixelDecay();
    needsPixelDecay = false;
  }
  for(int n = 0; n < numSamples; n++)
    phaseScopeBuffer.bufferSampleFrame(inOutBuffer[0][n], inOutBuffer[1][n]);
}

void PhaseScope::setSampleRate(double newSampleRate)
{
  phaseScopeBuffer.setSampleRate(newSampleRate);
}

void PhaseScope::reset()
{
  phaseScopeBuffer.reset();
  colorCounter = 0.0;
}

Colour PhaseScope::getAndUpdateColor()
{
  if(!rainbow)
    return Colours::white;

  ColourAHSL colorAHSL((float)colorCounter, 1.f, 0.5f, 0.5f);

  double inc    = 1.0 / (colorPeriod * phaseScopeBuffer.getFrameRate()); // color increment
  colorCounter += inc;
  colorCounter  = fmod(colorCounter, 1.0);

  return colorAHSL.getAsJuceColour();
}

//=================================================================================================

PhaseScopeDisplay::PhaseScopeDisplay(jura::PhaseScope *newPhaseScopeToEdit)
{
  phaseScope = newPhaseScopeToEdit;
  setSize(400, 400);

  startTimerHz(25);
  phaseScope->phaseScopeBuffer.setFrameRate(25);
}

void PhaseScopeDisplay::resized()
{
  phaseScope->setPixelSize(getWidth(), getHeight());
  image = Image(Image::ARGB, getWidth(), getHeight(), false);
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
}
// todo: write a version of this function that uses a colormap

void dataMatrixToImage(float **data, Image &image, 
  uint8 red = 255, uint8 green = 255, uint8 blue = 255)
{
  // We assume here that the size of the image (width and height) matches the dimensions of the 
  // data matrix)
  Image::BitmapData bitmap(image, Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);
  uint8 *pixelPointer = bitmap.getPixelPointer(0, 0);
  if(red == green && green == blue)
    dataMatrixToPixelBrightnessGray(data, pixelPointer, bitmap.width, bitmap.height);
  else
    dataMatrixToPixelBrightness(data, pixelPointer, bitmap.width, bitmap.height, red, green, blue);
}

void PhaseScopeDisplay::paint(Graphics &g)
{
  //dataMatrixToImage(phaseScope->phaseScopeBuffer.getDataMatrix(), image);
  //dataMatrixToImage(phaseScope->phaseScopeBuffer.getDataMatrix(), image, 0, 0, 255);

  Colour c = phaseScope->getAndUpdateColor();
  dataMatrixToImage(phaseScope->phaseScopeBuffer.getDataMatrix(), image, c.getRed(), c.getGreen(), 
    c.getBlue());
  g.drawImageAt(image, 0, 0);

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

void PhaseScopeDisplay::timerCallback()
{
  repaint();
  //phaseScope->phaseScopeBuffer.applyPixelDecay();
  phaseScope->triggerPixelDecay();

  //startTimerHz(phaseScope->getFrameRate());
  startTimer((int) (1000.0 / phaseScope->getFrameRate())); // fps to ms

  // there's a kind of bug - the display flickers and the flickering seems to be different if
  // we apply the pixel decay in the GUI thread (by calling 
  // phaseScope->phaseScopeBuffer.applyPixelDecay() here) or apply it in the audio thread
  // (by calling phaseScope->triggerPixelDecay()). maybe we need to do some kind of thread
  // synchronization
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

  setSize(400+widgetMargin, 400);  // preliminary
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

//void PhaseScopeEditor::setDisplayPixelSize(int newWidth, int newHeight)
//{
//  ScopedLock scopedLock(*plugInLock);
//  setSize(newWidth+widgetMargin, newHeight);
//
//  // hmm...maybe, we should also allow for the headline, preset section, etc.
//}

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
  sliderLineDensity->setBounds(x, y, w, h); y += dy;
  sliderFrameRate  ->setBounds(x, y, w, h); y += dy;

  buttonAntiAlias->setBounds(x, y,   w, h); y += dy;
  buttonRainbow  ->setBounds(x, y,   w, h); y += dy;
}
