
PhaseScopeBuffer::PhaseScopeBuffer()
{
  sampleRate = 44100.0;
  frameRate  = 25.0;
  decayTime  = 0.1;
  updateDecayFactor();

  bufferFlat = nullptr;
  buffer     = nullptr;
  width      = 200;
  height     = 200;
  allocateBuffer();

  reset();
}

PhaseScopeBuffer::~PhaseScopeBuffer()
{
  freeBuffer();
}

void PhaseScopeBuffer::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  double brightness = 4000.0;        // make this a member variable later
  insertFactor = (float) (brightness/sampleRate);
}

void PhaseScopeBuffer::setFrameRate(double newFrameRate)
{
  frameRate = newFrameRate;
  updateDecayFactor();
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
  }
}

void PhaseScopeBuffer::bufferSampleFrame(double left, double right)
{
  double x, y;
  x = left;   // preliminary - later, rotate by 45 degrees
  y = right;

  // compute coordinates, where this new pixel should be placed
  x  = 0.5*(x+1);  // convert -1..+1 into 0..1
  y  = 0.5*(y+1);
  x *= width;
  y *= height;
  int i = (int)round(x);
  int j = (int)round(y);

  // place the new pixel:
  if(i >= 0 && i < width && j >=0 && j < height)
    buffer[i][j] = min(1.f, buffer[i][j]+insertFactor);

  // todo: we could refine the display by de-interpolating the value into the buffer (i.e. modify 
  // not one but 4 pixels according to the frcational parts of x and y)
}

void PhaseScopeBuffer::applyPixelDecay()
{
  ArrayTools::rsScale(bufferFlat, width*height, decayFactor);
}

void PhaseScopeBuffer::reset()
{
  ArrayTools::rsFillWithZeros(bufferFlat, width*height);
}

void PhaseScopeBuffer::allocateBuffer()
{
  bufferFlat = new float[width*height];
  buffer = new float*[height];
  for(int i = 0; i < height; i++)
    buffer[i] = &bufferFlat[i*width];
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

//=================================================================================================

PhaseScope::PhaseScope(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*plugInLock);
  moduleName = "PhaseScope";
  setActiveDirectory(getApplicationDirectory() + "/PhaseScopePresets");
}

void PhaseScope::setPixelSize(int width, int height)
{
  phaseScopeBuffer.setSize(width, height);
}

AudioModuleEditor* PhaseScope::createEditor()
{
  return new PhaseScopeDisplay(this);
}

void PhaseScope::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  jassert(numChannels == 2);
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
}

//=================================================================================================

PhaseScopeDisplay::PhaseScopeDisplay(jura::PhaseScope *newPhaseScopeToEdit)
  : AudioModuleEditor(newPhaseScopeToEdit)
{
  phaseScope = newPhaseScopeToEdit;
  setSize(200, 200);

  startTimerHz(25);
  phaseScope->phaseScopeBuffer.setFrameRate(25);
}

void PhaseScopeDisplay::resized()
{
  phaseScope->setPixelSize(getWidth(), getHeight());
  image = Image(Image::ARGB, getWidth(), getHeight(), false);
}

void PhaseScopeDisplay::paint(Graphics &g)
{
  // this is very preliminary and (probably) horribly inefficient (we should work with an 
  // internally bufferd image and somehow obtain a pixel pointer and manipulate this data 
  // directly):
  g.fillAll(Colours::black);
  for(int x = 0; x < getWidth(); x++)
  {
    for(int y = 0; y < getHeight(); y++)
    {
      g.setColour(phaseScope->getColourAt(x, y));
      g.setPixel(x, y);
    }
  }


  //g.drawImageAt(image, 0, 0);
}

void PhaseScopeDisplay::timerCallback()
{
  repaint();
  phaseScope->phaseScopeBuffer.applyPixelDecay();
}
