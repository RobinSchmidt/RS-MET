
PhaseScopeBuffer::PhaseScopeBuffer()
{
  antiAlias  = true;
  frameRate  = 25.0;
  decayTime  = 0.1;
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
  double brightness = 10000.0;        // make this a member variable later
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
    reset();           // avoid flickering at size-change
  }
}

void PhaseScopeBuffer::setAntiAlias(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
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
  addDot((float)x, (float)y);
}

void PhaseScopeBuffer::applyPixelDecay()
{
  ArrayTools::rsScale(bufferFlat, width*height, decayFactor);
}

void PhaseScopeBuffer::reset()
{
  ArrayTools::rsFillWithZeros(bufferFlat, width*height);
}

void PhaseScopeBuffer::addDot(float x, float y)
{
  if(!antiAlias)
  {
    addDotFast(x, y);
    return;
  }

  int i = (int)floor(x);  // integer part of x
  int j = (int)floor(y);  // integer part of y
  x -= i;                 // fractional part of x
  y -= j;                 // fractional part of y

  // maybe factor this computation out into a function bilinearDeInterpolationWeights:
  double a, b, c, d;
  d = x*y;
  b = x-d;
  c = y-d;
  a = d-x-y+1;

  if(i >= 0 && i < width-1 && j >= 0 && j < height-1)
  {
    buffer[i]  [j]   = min(1.f, buffer[i]  [j]   + (float)a * insertFactor);
    buffer[i+1][j]   = min(1.f, buffer[i+1][j]   + (float)b * insertFactor);
    buffer[i]  [j+1] = min(1.f, buffer[i]  [j+1] + (float)c * insertFactor);
    buffer[i+1][j+1] = min(1.f, buffer[i+1][j+1] + (float)d * insertFactor);
  }

  // \todo: for optimization, take the min-operation out here and let it be applied at frame-rate
  // whenever the matrix is being read out (we need to provide a function saturateMatrixValues or 
  // something which the user can call at framerate and which applies the min function to the whole 
  // matrix at once) - then we can simply do 
  // buffer[i][j] += (float)a * insertFactor);  etc. here in this innermost code
  // ....but: it might be more expensive because it must be applied to the whole matrix,...hmmm
}

void PhaseScopeBuffer::addDotFast(float x, float y)
{
  int i = (int)round(x);
  int j = (int)round(y);
  if(i >= 0 && i < width && j >= 0 && j < height)
    buffer[i][j] = min(1.f, buffer[i][j]+insertFactor);
  //buffer[i][j] += insertFactor;
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
    //phaseScopeBuffer.bufferSampleFrameAntiAliased(inOutBuffer[0][n], inOutBuffer[1][n]);
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
  setSize(300, 300);

  startTimerHz(25);
  phaseScope->phaseScopeBuffer.setFrameRate(25);
}

void PhaseScopeDisplay::resized()
{
  phaseScope->setPixelSize(getWidth(), getHeight());
  image = Image(Image::ARGB, getWidth(), getHeight(), false);
}


// \todo move these 2 helper functions to jura_GraphicsTools, maybe templatize so it can be used 
// for double also:
void dataMatrixToPixelBrightness(float **data, uint8 *pixels, int width, int height)
{
  uint8 *p = pixels;
  for(int y = 0; y < height; y++)     // loop over lines
  {
    for(int x = 0; x < width; x++)    // loop over pixels
    {
      // we assume here, that the alpha channel comes last in the byte order of the pixels
      p[0] = p[1] = p[2] = (uint8) (255 * data[y][x]);  // data determines white-value
      p[3] = 255;                                       // set to full opacity ("alpha")
      p   += 4;                                         // jump to next pixel
    }
  }
  // todo: write a version of this function that uses a colormap
}
void dataMatrixToImage(float **data, Image &image)
{
  // We assume here that the size of the image (width and height) matches the dimensions of the 
  // data matrix)
  Image::BitmapData bitmap(image, Image::BitmapData::writeOnly);
  jassert(bitmap.pixelStride == 4);
  uint8 *pixelPointer = bitmap.getPixelPointer(0, 0);
  dataMatrixToPixelBrightness(data, pixelPointer, bitmap.width, bitmap.height);
}


void PhaseScopeDisplay::paint(Graphics &g)
{
  dataMatrixToImage(phaseScope->phaseScopeBuffer.getDataMatrix(), image);
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
  phaseScope->phaseScopeBuffer.applyPixelDecay();
}
