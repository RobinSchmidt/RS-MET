template<class TSig, class TPix, class TPar>
rsPhaseScopeBuffer<TSig, TPix, TPar>::rsPhaseScopeBuffer()
  : painter(&image, nullptr)
{
  painter.setUseAlphaMask(false);

  frameRate      = 25.0;
  decayTime      = 0.5;
  lineDensity    = 0.0f;
  maxDotsPerLine = 100;
  thickness      = 0.707f;
  brightness     = 1.0;
  useGradient    = false;

  //xOld = yOld = TSig(0);
  cOld = TPix(0);
  scanPos = TSig(0);

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
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setUseColorGradient(bool shouldUseGradient)
{
  useGradient = shouldUseGradient;
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
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setLineDensity(TPar newDensity)
{
  lineDensity = newDensity;
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::setLineDensityLimit(int newMaxNumDotsPerLine)
{
  maxDotsPerLine = newMaxNumDotsPerLine;
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
  addSegmentTo(x, y);  // rename to addCurveTo, addSegmentTo
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::applyPixelDecay()
{
  rsArray::scale(image.getPixelPointer(0, 0), image.getNumPixels(), decayFactor);
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::reset()
{
  rsArray::fillWithZeros(image.getPixelPointer(0, 0), image.getNumPixels());
  scanPos = 0.0;
  splineGen.reset();

  // remove soon:
  cOld = TPix(0);
  rsArray::fillWithValue(x, 4, TSig(0.5 * image.getWidth()) );
  rsArray::fillWithValue(y, 4, TSig(0.5 * image.getHeight()));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::addSegmentTo(TSig newX, TSig newY)
{
  rsArray::pushFrontPopBack4(newX, x);
  rsArray::pushFrontPopBack4(newY, y);

  if(lineDensity == 0.f)
    painter.paintDot(newX, newY, (TPix) insertFactor);
  else
    drawDottedLine((TPix) insertFactor, (TSig)lineDensity, maxDotsPerLine, true);
}
// obsolete soon

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::drawDottedLine(TPix color, TPar density, int maxNumDots, 
  bool scaleByNumDots, TPar minDotDistance)
{
  // x2 was always the destination and x1 the source

  TSig dx = x[1]-x[2];
  TSig dy = y[1]-y[2];
  TSig pixelDistance = sqrt(dx*dx + dy*dy);
  int  numDots = rsMax(1, (int)floor(density*pixelDistance/minDotDistance));
  if(maxNumDots > 0)
    numDots = rsMin(numDots, maxNumDots);

  TPix c;
  if(useGradient)
  {
    TPix ct = color / (TPix)numDots;  // target color that would be used if we don't do gradients
    c = (2.f*ct) - cOld;              // desired endpoint color
    c = rsMax(c, ct);                 // c could come out negative, use ct as lower bound
    drawDottedSegment(cOld, c, numDots);
    cOld = c;
  }
  else
  {
    c = color;
    if(scaleByNumDots)
      c = color / (TPix)numDots;
    drawDottedSegment(c, c, numDots);
  }
}
// rename to drawConnection or drawCurve or something and dispatch between line-drawing, 
// spline-drawing and maybe other drawing modes (maybe non-dotted like bresenham, wu, etc.)
// obsolete soon

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::drawDottedSegment(TPix color1, TPix color2, int numDots)
{
  // maybe this function should not have a numDots parameter

  switch(drawMode)
  {
  case DOTTED_LINE: painter.drawLineDotted(x[2], y[2], x[1], y[1], color1, color2, numDots); break;
  case DOTTED_SPLINE:
  {
    // compute (2nd order) derivative estimates at inner points x[1], x[2]:
    TSig dx1, dx2, dy1, dy2;
    dx1 = TSig(0.5)*(x[0]-x[2]);    // estimated dx/dt at x[1]
    dy1 = TSig(0.5)*(y[0]-y[2]);    // estimated dy/dt at y[1]
    dx2 = TSig(0.5)*(x[1]-x[3]);    // estimated dx/dt at x[2]
    dy2 = TSig(0.5)*(y[1]-y[3]);    // estimated dy/dt at y[2]

    // draw hermite spline with given derivative values:
    painter.drawDottedSpline(
      x[2], dx2, y[2], dy2,      // older inner point comes first
      x[1], dx1, y[1], dy1,      // newer inner point comes second
      color1, color2, numDots);

    // todo: maybe scale the number of dots according to the length of the spline segment 
    // (currently, the number is determined by the length of the line connecting the inner points
    // which serves as a crude approximation of the arc-length)

    // there's a color mismatch at the joints (draw the test lissajous figure with 40 data points
    // to see it). it's especially apparent in regions of high curvature...has it to do with the 
    // number of dots? does it happen when a low-curvature segment joins a high-curvature segment?
    // that would make sense bcs the high-curv segment spreads its dots along a longer arc

    // -> figure out, if there's an abrupt change in dot-color/brightness or if the apparent color
    // mismatch is indeed due to a density-change of the dots...maybe print out dot-color and 
    // dot-distance between subsequent dots...or maybe try to plot it
    // or: use smaller density to plot non-overlapping dots

    // i think, it has to do with the dots inside a single segment not being equidistant
    // a segment may start with low dot density, end with high dot density and join to a segment 
    // that again starts with low density. ideas:
    // -try to figure out a better sequence of t-values, making the dots more equidistant
    // -or: compensate varying dot-density by varying brightness in the opposite way
    //  -keep track of distance d between successive dots and scale brightness by d/da where da
    //   is the assumed average distance
    //  -maybe instead of keeping track of the actual distance (requiring a square-root for each 
    //   dot), we could use something based on the 2nd derivative (which is easy to compute), maybe
    //   brightness *= 1 / (1+cx*cx+cy*cy), cx = d2x/dt2, cy = d2y/dt2 - 2nd derivatives of x,y 
    //   with respect to t

  } break;

  default: { } break; // don't draw anything if drawMode has invalid value

    // todo: dispatch according to drawMode, store xOld, yOld here instead of in addLineTo
  }

}



template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateDecayFactor()
{
  decayFactor = (TPar) exp(-1 / (decayTime*frameRate));
}

template<class TSig, class TPix, class TPar>
void rsPhaseScopeBuffer<TSig, TPix, TPar>::updateInsertFactor()
{
  insertFactor = (TPix(10000) * brightness / TPix(sampleRate));
  // The factor is totally ad-hoc - maybe come up with some more meaningful factor.
  // However, the proportionality to the brightness parameter and inverse proportionality to
  // the sample rate seems to make sense.
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
    TPix mean = rsArray::mean(p, N);
    //TPix scaler = 1 / (1 + decayByAverage*mean); // maybe try a formula with exp
    TPix scaler =  (TPix) (1-exp(-1 / (decayByAverage*mean)));
    rsArray::scale(p, N, scaler);
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
  // function needs update since we implemented the spline drawing in the baseclass

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
}
