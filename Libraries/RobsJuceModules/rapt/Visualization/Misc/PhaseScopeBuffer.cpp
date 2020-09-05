template<class TSig, class TPix, class TPar>
rsPhaseScopeBuffer<TSig, TPix, TPar>::rsPhaseScopeBuffer()
  : painter(&image, nullptr)
{
  painter.setUseAlphaMask(false);

  frameRate  = 25.0;
  decayTime  = 0.5;
  thickness  = 0.707f;
  brightness = 1.0;
  scanPos    = TSig(0);

  updateDecayFactor();
  updateTransformCoeffs();
  setSampleRate(44100.0);
  reset();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateInsertFactor();
  screenScanner.setSampleRate(TSig(sampleRate));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setFrameRate(TPar newFrameRate)
{
  frameRate = newFrameRate;
  updateDecayFactor();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setBrightness(TPix newBrightness)
{
  brightness = newBrightness;
  updateInsertFactor();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setDecayTime(TPar newDecayTime)
{
  decayTime = newDecayTime;
  updateDecayFactor();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setSize(int newWidth, int newHeight)
{
  image.setSize(newWidth, newHeight);
  painter.setImageToPaintOn(&image); // so it can update it's internal size variables, too
  updateDotBufferSizes();
  reset();                           // avoid flickering at size-change
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setAntiAlias(bool shouldAntiAlias)
{
  painter.setAntiAlias(shouldAntiAlias);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setMaxSizeWithoutReAllocation(int newMaxWidth, int newMaxHeight)
{
  image.setMaxSize(newMaxWidth, newMaxHeight);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setPixelSpread(TPar newSpread)
{
  thickness = newSpread;
  painter.setNeighbourWeightsForSimpleDot((TSig)thickness, TSig(thickness*thickness));
  //painter.setNeighbourWeightsForSimpleDot((TSig)thickness, 0.5f*thickness);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setOneDimensionalMode(bool shouldBe1D)
{
  oneDimensonal = shouldBe1D;
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setScanningFrequency(TPar newFrequency)
{
  screenScanner.setScanFreqNoSync(TSig(newFrequency));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setSyncMode(bool shouldSync)
{
  screenScanner.setSync(shouldSync);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setScaleX(TSig newScale)
{
  scaleX = newScale;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setScaleY(TSig newScale)
{
  scaleY = newScale;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setShearX(TSig newShear)
{
  shearX = newShear;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setShearY(TSig newShear)
{
  shearY = newShear;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setRotation(TSig degrees)
{
  rotation = TSig(PI/180) * degrees;
  updateTransformCoeffs();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setShiftX(TSig newShift)
{
  shiftX = newShift;
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setShiftY(TSig newShift)
{
  shiftY = newShift;
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::toPixelCoordinates(TSig &x, TSig &y)
{
  x  = TSig(0.5) * (x+1);  // convert -1..+1 into 0..1
  y  = TSig(0.5) * (y+1);
  x *= image.getWidth();
  y *= image.getHeight();
  y  = image.getHeight()-1 - y;
  // maybe we should add 0.5 after multiplication by width/height?
}
// maybe move to splineGen

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::processSampleFrame(TSig x, TSig y)
{
  // apply affine transformation in normalized coordinates:
  TSig xt = x;  // temporary
  x = Axx * x  + Axy * y + shiftX;
  y = Ayx * xt + Ayy * y + shiftY;

  // replace x with sawtooth-scanner in 1D mode:
  if(oneDimensonal == true)
  {
    xt = x+y;                         // use sum for analysis and display (*)
    x  = scaleX * getScannerSaw(xt);
    y  = xt;
    // (*) todo: use arbitrary linear combination of x,y maybe have convenience function to select
    // left (1,0), right(0,1), left+right(1,1), (left+right)/2 (.5,.5), 
    // left+right/sqrt(2) = mid, left-right/sqrt(2) = side
  }

  // transform to pixel coordinates and draw line:
  toPixelCoordinates(x, y);

  if(oneDimensonal && screenScanner.resetOccurred())
    moveTo(x, y);        // start a new segment without drawing a connection from old datapoint
  else
    addSegmentTo(x, y);  // connect old datapoint with current
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::applyPixelDecay()
{
  rsArrayTools::scale(image.getPixelPointer(0, 0), image.getNumPixels(), decayFactor);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::reset()
{
  clearImage();
  scanPos = 0.0;
  splineGen.reset(TSig(0.5*image.getWidth()), TSig(0.5*image.getHeight()));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::moveTo(TSig newX, TSig newY)
{
  splineGen.updatePointBuffers(newX, newY);
  splineGen.shiftPointBuffers(-TSig(image.getWidth()), TSig(0));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::addSegmentTo(TSig newX, TSig newY)
{
  int numDots = splineGen.updateDotBuffers(newX, newY);
  painter.drawDots(&dotsX[0], &dotsY[0], &dotsC[0], numDots);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateDecayFactor()
{
  decayFactor = (TPar) exp(-1 / (decayTime*frameRate));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateInsertFactor()
{
  TPix insertFactor = (TPix(10000) * brightness / TPix(sampleRate));
  //TPix insertFactor = (TPix(1000) * brightness / TPix(sampleRate));
  splineGen.setBrightness(insertFactor);

  // The factor is totally ad-hoc - maybe come up with some more meaningful factor.
  // However, the proportionality to the brightness parameter and inverse proportionality to
  // the sample rate seems to make sense. Maybe it should also be related to the frameRate?
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateTransformCoeffs()
{
  TSig s = sin(rotation);
  TSig c = cos(rotation);
  Axx = c*scaleX + s*shearY;
  Axy = s*scaleY + c*shearX;
  Ayx = c*shearY - s*scaleX;
  Ayy = c*scaleY - s*shearX;
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateDotBufferSizes()
{
  int size = 2048; // preliminary - maybe use k*sqrt(w^2 + h^2) for some k >= 1
                   // hmm...or maybe some other formula based on width and height
  dotsX.resize(size);
  dotsY.resize(size);
  dotsC.resize(size);
  splineGen.setDotBuffers(&dotsX[0], &dotsY[0], &dotsC[0], size);
  // ...because resizing may cause re-allocation, also, it needs to know the new size
}

//template<class TSig, class TPix, class TPar>
//void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateScanIncrement()
//{
//  scanInc = scanFreq / sampleRate;
//}

template<class TSig, class TPix, class TPar>
TSig rsPhaseScopeBuffer<TSig, TPix, TPar>::getScannerSaw(TSig x)
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
rsPhaseScopeBuffer2<TSig, TPix, TPar>::rsPhaseScopeBuffer2()
: lineDrawer(&this->image)
{
  this->painter.setAlphaMaskForDot(&this->dotMask);
  dotMask.setMaxSize(20, 20);
  dotMask.setSize(5);
  dotMask.setTransitionWidth(0.5);

  lineDrawer.setBlendMode(rsImageDrawer<TPix, TSig, TSig>::BLEND_ADD_SATURATE);
  lineDrawer.setRoundCaps(true); // preliminary
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::setUseAlphaMask(bool shouldUseMask)
{
  this->painter.setUseAlphaMask(shouldUseMask);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::setPixelDecayByValue(TPar newDecayByValue)
{
  decayByValue = newDecayByValue;
  updateDecayFactor();
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::setPixelDecayByAverage(TPar newDecayByAverage)
{
  decayByAverage = newDecayByAverage;
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::setLineBrightness(TPar newBrightness)
{
  lineBrightness = TPix(newBrightness);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::setLineWidth(TPar newWidth)
{
  lineDrawer.setLineWidth(newWidth);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::setLineProfile(int newProfile)
{
  lineDrawer.setLineProfile(newProfile);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::applyPixelDecay()
{
  int   N = this->image.getNumPixels();
  TPix *p = this->image.getPixelPointer(0, 0);

  if(decayByValue == 0)
    rsPhaseScopeBuffer<TSig, TPix, TPar>::applyPixelDecay();
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
    TPix mean = rsArrayTools::mean(p, N);
    //TPix scaler = 1 / (1 + decayByAverage*mean); // maybe try a formula with exp
    TPix scaler =  (TPix) (1-exp(-1 / (decayByAverage*mean)));
    rsArrayTools::scale(p, N, scaler);
  }
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::updateDecayFactor()
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
void rsPhaseScopeBuffer2<TSig, TPix, TPar>::addSegmentTo(TSig x, TSig y)
{
  rsAssert(false);
  // function needs update since we implemented the spline drawing in the baseclass

  /*
  if(drawLines){
    //TSig dx = x - this->xOld;
    //TSig dy = y - this->yOld;
    TSig dx = x - this->x[1];
    TSig dy = y - this->y[1];
    TPar L  = sqrt(dx*dx + dy*dy);
    TPar scaler = 1;
    if(L > 1)
      scaler = 1/L;
    lineDrawer.setColor(lineBrightness * TPix(scaler));
    //lineDrawer.drawLine(xOld, yOld, x, y); // old - creates phantom dots at joints
    lineDrawer.lineTo(x, y);
  }
  if(drawDots)
    rsPhaseScopeBuffer<TSig, TPix, TPar>::addSegmentTo(x, y); // draws dotted line, updates xOld, yOld
  else {
    //this->xOld = x;
    //this->yOld = y; 
    this->x[1] = x;
    this->y[1] = y; 
  }
  */
}
