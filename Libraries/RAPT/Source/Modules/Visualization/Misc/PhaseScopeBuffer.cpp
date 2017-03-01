template<class TSig, class TPix, class TPar>
PhaseScopeBuffer<TSig, TPix, TPar>::PhaseScopeBuffer()
  : painter(&image, nullptr)
{
  painter.setUseAlphaMask(false);

  frameRate      = 25.0;
  decayTime      = 0.5;
  lineDensity    = 0.0f;
  maxDotsPerLine = 100;
  thickness      = 0.707f;
  brightness     = 1.0;
  updateDecayFactor();

  setSampleRate(44100.0);
  reset();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateInsertFactor();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setFrameRate(TPar newFrameRate)
{
  frameRate = newFrameRate;
  updateDecayFactor();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setBrightness(TPix newBrightness)
{
  brightness = newBrightness;
  updateInsertFactor();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setDecayTime(TPar newDecayTime)
{
  decayTime = newDecayTime;
  updateDecayFactor();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setSize(int newWidth, int newHeight)
{
  image.setSize(newWidth, newHeight);
  painter.setImageToPaintOn(&image); // so it can update it's internal size variables, too
  reset();                           // avoid flickering at size-change
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setAntiAlias(bool shouldAntiAlias)
{
  painter.setAntiAlias(shouldAntiAlias);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setMaxSizeWithoutReAllocation(int newMaxWidth, int newMaxHeight)
{
  image.setMaxSize(newMaxWidth, newMaxHeight);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setLineDensity(TPar newDensity)
{
  lineDensity = newDensity;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setLineDensityLimit(int newMaxNumDotsPerLine)
{
  maxDotsPerLine = newMaxNumDotsPerLine;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setPixelSpread(TPar newSpread)
{
  thickness = newSpread;
  painter.setNeighbourWeightsForSimpleDot((TSig)thickness, thickness*thickness);
  //painter.setNeighbourWeightsForSimpleDot((TSig)thickness, 0.5f*thickness);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::toPixelCoordinates(TSig &x, TSig &y)
{
  x  = TSig(0.5) * (x+1);  // convert -1..+1 into 0..1
  y  = TSig(0.5) * (y+1);
  x *= image.getWidth();
  y *= image.getHeight();
  y  = image.getHeight()-1 - y;
  // maybe we should add 0.5 after multiplication by width/height?
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::processSampleFrame(TSig x, TSig y)
{
  toPixelCoordinates(x, y);
  addLineTo(x, y);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::applyPixelDecay()
{
  ArrayTools::rsScale(image.getPixelPointer(0, 0), image.getNumPixels(), decayFactor);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::reset()
{
  ArrayTools::rsFillWithZeros(image.getPixelPointer(0, 0), image.getNumPixels());

  // (xOld,yOld) = (0,0) - but in pixel coordinates:
  xOld = 0.5f * image.getWidth();
  yOld = 0.5f * image.getHeight();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::addLineTo(TSig x, TSig y)
{
  if(lineDensity == 0.f)
    painter.paintDot(x, y, (TPix) insertFactor);
  else
    painter.drawDottedLine((TSig)xOld, (TSig)yOld, (TSig)x, (TSig)y, (TPix) insertFactor, 
      (TSig)lineDensity, maxDotsPerLine);
  xOld = x;
  yOld = y;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::updateDecayFactor()
{
  decayFactor = (TPar) exp(-1 / (decayTime*frameRate));
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::updateInsertFactor()
{
  insertFactor = (TPix(10000) * brightness / TPix(sampleRate));
  // The factor is totally ad-hoc - maybe come up with some more meaningful factor. 
  // However, the proportionality to the birghtness parameter and inverse proportionality to 
  // the sample rate seems to make sense.
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPix, class TPar>
PhaseScopeBuffer2<TSig, TPix, TPar>::PhaseScopeBuffer2() 
: lineDrawer(&image)
{
  painter.setAlphaMaskForDot(&dotMask);
  dotMask.setMaxSize(20, 20);
  dotMask.setSize(5);
  dotMask.setTransitionWidth(0.5);

  lineDrawer.setBlendMode(ImageDrawer<TPix, TSig, TSig>::BLEND_ADD_SATURATE);
  lineDrawer.setRoundCaps(true); // preliminary
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setUseAlphaMask(bool shouldUseMask)
{
  painter.setUseAlphaMask(shouldUseMask);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setPixelDecayByValue(TPar newDecayByValue)
{
  decayByValue = newDecayByValue;
  updateDecayFactor();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setPixelDecayByAverage(TPar newDecayByAverage)
{
  decayByAverage = newDecayByAverage;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setLineBrightness(TPar newBrightness)
{
  lineBrightness = TPix(newBrightness);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setLineWidth(TPar newWidth)
{
  lineDrawer.setLineWidth(newWidth);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setLineProfile(int newProfile)
{
  lineDrawer.setLineProfile(newProfile);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::applyPixelDecay()
{
  int   N = image.getNumPixels();  
  TPix *p = image.getPixelPointer(0, 0);

  if(decayByValue == 0)
    PhaseScopeBuffer::applyPixelDecay();
  else
  {

    for(int n = 0; n < N; n++)
      p[n] *= p[n]*TPix(decayFactorAt1) + (1-p[n])*TPix(decayFactor);
    // i think, a linear interpolation between two different decay factors might not be ideal
    // maybe we should use a table instead
  }

  if(decayByAverage != 0)
  {
    // apply additional decay that depends on the global average value
    TPix mean   = ArrayTools::rsMean(p, N);
    //TPix scaler = 1 / (1 + decayByAverage*mean); // maybe try a formula with exp
    TPix scaler =  (TPix) (1-exp(-1 / (decayByAverage*mean)));
    ArrayTools::rsScale(p, N, scaler);
  }
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::updateDecayFactor()
{
  if(decayByValue == 0)
  {
    decayFactor = (TPar)exp(-1 / (decayTime*frameRate)); // taken from baseclass
    decayFactorAt1 = decayFactor;
  }
  else
  {
    TPar factor    = pow(2, decayByValue);
    TPar decayAt1  = decayTime * factor;
    TPar decayAt0  = decayTime;
    decayFactor    = (TPar)exp(-1 / (decayAt0*frameRate));
    decayFactorAt1 = (TPar)exp(-1 / (decayAt1*frameRate));

    //TPar factor    = pow(2, 0.5*decayByValue);
    //TPar decayAt1  = decayTime * factor;
    //TPar decayAt0  = decayTime / factor;
    //decayFactor    = (TPar)exp(-1 / (decayAt0*frameRate));
    //decayFactorAt1 = (TPar)exp(-1 / (decayAt1*frameRate));
  }
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::addLineTo(TSig x, TSig y)
{
  if(drawLines){
    TSig dx = x - xOld;
    TSig dy = y - yOld;
    TPar L  = sqrt(dx*dx + dy*dy);
    TPar scaler = 1;
    if(L > 1)
      scaler = 1/L;
    lineDrawer.setColor(lineBrightness * TPix(scaler));
    lineDrawer.drawLine(xOld, yOld, x, y); }
  if(drawDots)
    PhaseScopeBuffer::addLineTo(x, y); // draws dotted line, updates xOld, yOld
  else {
    xOld = x;
    yOld = y; }
}