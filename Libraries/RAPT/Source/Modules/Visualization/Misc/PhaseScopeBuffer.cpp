template<class TSig, class TPix, class TPar>
PhaseScopeBuffer<TSig, TPix, TPar>::PhaseScopeBuffer()
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

template<class TSig, class TPix, class TPar>
PhaseScopeBuffer<TSig, TPix, TPar>::~PhaseScopeBuffer()
{
  freeBuffer();
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
  if(newWidth != width || newHeight != height)
  {
    freeBuffer();
    width  = newWidth;
    height = newHeight;
    allocateBuffer();
    reset();           // avoid flickering at size-change
  }
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setAntiAlias(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setLineDensity(TPar newDensity)
{
  lineDensity = newDensity;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::setPixelSpread(TPar newSpread)
{
  thickness = newSpread;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::convertAmplitudesToMatrixIndices(TSig &x, TSig &y)
{
  x  = TSig(0.5) * (x+1);  // convert -1..+1 into 0..1
  y  = TSig(0.5) * (y+1);
  x *= width;
  y *= height;
  // maybe we should add 0.5 after multiplication by width/height?
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::bufferSampleFrame(TSig x, TSig y)
{
  convertAmplitudesToMatrixIndices(x, y);
  addLineTo(x, y);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::applyPixelDecay()
{
  ArrayTools::rsScale(bufferFlat, width*height, decayFactor);
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::reset()
{
  ArrayTools::rsFillWithZeros(bufferFlat, width*height);

  // (xOld,yOld) = (0,0) - but in pixel coordinates:
  xOld = 0.5f*width;   
  yOld = 0.5f*height;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::addLineTo(TSig x, TSig y)
{
  if(lineDensity == 0.f)
  {
    addDot(x, y, insertFactor);
    return;
  }

  TSig dx = x-xOld;
  TSig dy = y-yOld;
  TSig pixelDistance = sqrt(dx*dx + dy*dy);
  int  numDots = rsMax(1, (int)floor(lineDensity*pixelDistance));
  //int   numDots = max(1, (int)floor(lineDensity*pixelDistance));
  TPix intensity = (TPix) (insertFactor/numDots);
  TSig scaler = (TSig)(1.0 / numDots);
  TSig k;
  for(int i = 1; i <= numDots; i++)
  {
    k = scaler * i;  // == i / numDots
                     //addDot((1-k)*xOld + k*x, (1-k)*yOld + k*y, intensity); // unoptimized
    addDot(xOld + k*dx, yOld + k*dy, intensity);
  }
  xOld = x;
  yOld = y;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::addDot(TSig x, TSig y, TPix intensity)
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
  TPix a, b, c, d;
  d = TPix(x*y);
  c = TPix(y)-d;
  b = TPix(x)-d;
  a = 1+d-TPix(x+y);

  // compute values to accumulate into the 4 pixels at (i,j), (i+1,j), (i,j+1), (i+1,j+1):
  a *= intensity;  // (i,   j)
  b *= intensity;  // (i,   j+1)
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
    TPix t, s, sa, sb, sc, sd, ta, tb, tc, td;
    t = (TPix)thickness;      // weight for direct neighbour pixels
    s = t * (TPix)SQRT2_INV;  // weight for diagonal neighbour pixels

    s = t*t;  // test

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

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::addDotFast(TSig x, TSig y, TPix intensity)
{
  int j = (int)round(x);
  int i = (int)round(y);

  if(j >= 0 && j < width && i >= 0 && i < height)
    accumulate(buffer[i][j], intensity);

  // apply thickness:
  if(thickness > 0.f && j >= 1 && j < width-1 && i >= 1 && i < height-1)
  {
    TPix a, ta, sa;
    a  = intensity;
    ta = a  * (TPix)thickness;
    sa = ta * (TPix)SQRT2_INV;

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

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::allocateBuffer()
{
  bufferFlat = new TPix[width*height];
  buffer = new TPix*[height];
  for(int i = 0; i < height; i++)
    buffer[i] = &bufferFlat[i*width]; 
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::freeBuffer()
{
  delete[] buffer;
  delete[] bufferFlat;
  buffer = nullptr;
  bufferFlat = nullptr;
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::updateDecayFactor()
{
  decayFactor = (TPar) exp(-1 / (decayTime*frameRate));
}

template<class TSig, class TPix, class TPar>
void PhaseScopeBuffer<TSig, TPix, TPar>::updateInsertFactor()
{
  insertFactor = (TPix) (10000*brightness / sampleRate);
  // The factor is totally ad-hoc - maybe come up with some more meaningful factor. 
  // However, the proportionality to the birghtness parameter and inverse proportionality to 
  // the sample rate seems to make sense.
}
