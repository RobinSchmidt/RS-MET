template<class TSig, class TPix, class TPar>
PhaseScopeBuffer<TSig, TPix, TPar>::PhaseScopeBuffer()
  : painter(&image, nullptr)
{
  painter.setUseAlphaMask(false);

  //scanFreq       = 1.0;
  frameRate      = 25.0;
  decayTime      = 0.5;
  lineDensity    = 0.0f;
  maxDotsPerLine = 100;
  thickness      = 0.707f;
  brightness     = 1.0;
  useGradient    = false;
  updateDecayFactor();
  updateTransformCoeffs();

  setSampleRate(44100.0);
  reset();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateInsertFactor();
  screenScanner.setSampleRate(sampleRate);
  //updateScanIncrement();
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
void PhaseScopeBuffer<TSig, TPix, TPar>::setUseColorGradient(bool shouldUseGradient)
{
  useGradient = shouldUseGradient;
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
  painter.setNeighbourWeightsForSimpleDot((TSig)thickness, TSig(thickness*thickness));
  //painter.setNeighbourWeightsForSimpleDot((TSig)thickness, 0.5f*thickness);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setOneDimensionalMode(bool shouldBe1D)
{
  oneDimensonal = shouldBe1D;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setScanningFrequency(TPar newFrequency)
{
  screenScanner.setScanFreqNoSync(newFrequency);
  //scanFreq = newFrequency;
  //updateScanIncrement();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setScaleX(TSig newScale)
{
  scaleX = newScale;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setScaleY(TSig newScale)
{
  scaleY = newScale;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setShearX(TSig newShear)
{
  shearX = newShear;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setShearY(TSig newShear)
{
  shearY = newShear;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setRotation(TSig degrees)
{
  rotation = TSig(PI/180) * degrees;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setShiftX(TSig newShift)
{
  shiftX = newShift;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setShiftY(TSig newShift)
{
  shiftY = newShift;
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
  // apply affine transformation in normalized coordinates:
  TSig xt = x;  // temporary
  x = Axx * x  + Axy * y + shiftX;
  y = Ayx * xt + Ayy * y + shiftY;

  // replace x with sawtooth-scanner in 1D mode:
  if(oneDimensonal == true)
    x = getScannerSaw(x+y);

  // transform to pixel coordinates and draw line:
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
  scanPos = 0.0;

  // (xOld,yOld) = (0,0) - but in pixel coordinates:
  xOld = 0.5f * image.getWidth();
  yOld = 0.5f * image.getHeight();
  cOld = TPix(0);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::addLineTo(TSig x, TSig y)
{
  if(lineDensity == 0.f)
    painter.paintDot(x, y, (TPix) insertFactor);
  else
    drawDottedLine((TSig)xOld, (TSig)yOld, (TSig)x, (TSig)y, (TPix) insertFactor,
      (TSig)lineDensity, maxDotsPerLine, true);
  xOld = x;
  yOld = y;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::drawDottedLine(TSig x1, TSig y1, TSig x2, 
  TSig y2, TPix color, TPar density, int maxNumDots, bool scaleByNumDots, TPar minDotDistance)
{
  TSig dx = x2-x1;
  TSig dy = y2-y1;
  TSig pixelDistance = sqrt(dx*dx + dy*dy);
  int  numDots = rsMax(1, (int)floor(density*pixelDistance/minDotDistance));
  if(maxNumDots > 0)
    numDots = rsMin(numDots, maxNumDots);

  TPix c;
  if(useGradient)
  {
    TPix ct = color / (TPix)numDots;  // target color that would be used if we don't do gradients
    c = (2.f*ct) - cOld;              // desired enpoint color
    c = rsMax(c, ct);                 // c could come out negative, use ct as lower bound
    painter.drawLineDotted(x1, y1, x2, y2, cOld, c, numDots);
    cOld = c;
  }
  else
  {
    c = color;
    if(scaleByNumDots)
      c = color / (TPix)numDots;
    painter.drawLineDotted(x1, y1, x2, y2, c, c, numDots);
  }
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

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::updateTransformCoeffs()
{
  TSig s = sin(rotation);
  TSig c = cos(rotation);
  Axx = c*scaleX + s*shearY;
  Axy = s*scaleY + c*shearX;
  Ayx = c*shearY - s*scaleX;
  Ayy = c*scaleY - s*shearX;
}

//template<class TSig, class TPix, class TPar>
//void PhaseScopeBuffer<TSig, TPix, TPar>::updateScanIncrement()
//{
//  scanInc = scanFreq / sampleRate;
//}

template<class TSig, class TPix, class TPar>
TSig PhaseScopeBuffer<TSig, TPix, TPar>::getScannerSaw(TSig x)
{
  return 2*screenScanner.getSample(x) - 1;

  //// replace by usen the ScopeScreenScanner:
  //TSig scanVal = 2*scanPos - 1; 
  //scanPos += TSig(scanInc);
  //if(scanPos > 1)
  //  scanPos -= 1;
  //return scanVal;
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPix, class TPar>
PhaseScopeBuffer2<TSig, TPix, TPar>::PhaseScopeBuffer2()
: lineDrawer(&this->image)
{
  this->painter.setAlphaMaskForDot(&this->dotMask);
  dotMask.setMaxSize(20, 20);
  dotMask.setSize(5);
  dotMask.setTransitionWidth(0.5);

  lineDrawer.setBlendMode(ImageDrawer<TPix, TSig, TSig>::BLEND_ADD_SATURATE);
  lineDrawer.setRoundCaps(true); // preliminary
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer2<TSig, TPix, TPar>::setUseAlphaMask(bool shouldUseMask)
{
  this->painter.setUseAlphaMask(shouldUseMask);
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
  int   N = this->image.getNumPixels();
  TPix *p = this->image.getPixelPointer(0, 0);

  if(decayByValue == 0)
    PhaseScopeBuffer<TSig, TPix, TPar>::applyPixelDecay();
  else
  {
    for(int n = 0; n < N; n++)
      p[n] *= p[n]*TPix(decayFactorAt1) + (1-p[n])*TPix(this->decayFactor);
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
    this->decayFactor = (TPar)exp(-1 / (this->decayTime*this->frameRate)); // taken from baseclass
    decayFactorAt1 = this->decayFactor;
  }
  else
  {
    TPar factor       = pow(2, decayByValue);
    TPar decayAt1     = this->decayTime * factor;
    TPar decayAt0     = this->decayTime;
    this->decayFactor = (TPar)exp(-1 / (decayAt0*this->frameRate));
    decayFactorAt1    = (TPar)exp(-1 / (decayAt1*this->frameRate));

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
    TSig dx = x - this->xOld;
    TSig dy = y - this->yOld;
    TPar L  = sqrt(dx*dx + dy*dy);
    TPar scaler = 1;
    if(L > 1)
      scaler = 1/L;
    lineDrawer.setColor(lineBrightness * TPix(scaler));
    //lineDrawer.drawLine(xOld, yOld, x, y); // old - creates phantom dots at joints
    lineDrawer.lineTo(x, y);
  }
  if(drawDots)
    PhaseScopeBuffer<TSig, TPix, TPar>::addLineTo(x, y); // draws dotted line, updates xOld, yOld
  else {
    this->xOld = x;
    this->yOld = y; }
}
