template<class TPix, class TWgt, class TCor>
ImageDrawer<TPix, TWgt, TCor>::ImageDrawer(Image<TPix> *imageToDrawOn)
{
  setImageToDrawOn(imageToDrawOn);
  setColor(TPix(1));  // typically white
  setBlendMode(0);
}

// setup:

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::setImageToDrawOn(Image<TPix> *imageToDrawOn)
{
  image = imageToDrawOn;
}

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::setBlendMode(int newMode)
{
  rsAssert(newMode >= 0);
  rsAssert(newMode < NUM_BLEND_MODES);
  blendMode = newMode;
  switch(blendMode)
  {
  case BLEND_LINEAR:       blendFunction = linearBlend;    break;
  case BLEND_ADD_CLIP:     blendFunction = addAndClip;     break;
  case BLEND_ADD_SATURATE: blendFunction = addAndSaturate; break;
  }
}

// blend functions:

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::linearBlend(TPix &pixel, TPix color, TWgt blend)
{
  pixel = (1-TPix(blend)) * pixel + TPix(blend) * color;
}

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::addAndClip(TPix &pixel, TPix color, TWgt blend)
{
  pixel = rsMin(TPix(1), pixel + TPix(blend)*color);
}

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::addAndSaturate(TPix &pixel, TPix color, TWgt blend)
{
  color *= TPix(blend);
  pixel  = (pixel + color) / (TPix(1) + color);
}

//=================================================================================================

template<class TPix, class TWgt, class TCor>
LineDrawer<TPix, TWgt, TCor>::LineDrawer(Image<TPix> *imageToDrawOn) 
  : ImageDrawer(imageToDrawOn) 
{
  setLineProfile(0);
  setLineWidth(1);
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::setLineProfile(int newProfile)
{
  rsAssert(newProfile >= 0);
  rsAssert(newProfile < NUM_LINE_PROFILES);
  profileIndex = newProfile;
  switch(profileIndex)
  {
  case PROFILE_FLAT:      lineProfile = profileFlat;      break;
  case PROFILE_LINEAR:    lineProfile = profileLinear;    break;
  case PROFILE_PARABOLIC: lineProfile = profileParabolic; break;
  case PROFILE_CUBIC:     lineProfile = profileCubic;     break;
  }
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::setLineWidth(TCor newWidth)
{
  rsAssert(newWidth > 0);
  w2 = 0.5*(newWidth+1);   // +1 is a hack, otherwise lines are 1 pixel too narrow
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::setRoundCaps(bool shouldBeRound)
{
  roundCaps = shouldBeRound;
}

// drawing:

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::drawLine(TCor x0, TCor y0, TCor x1, TCor y1)
{
  TCor dx, dy, a, b, yf, dp, d, L, A, B, C0, C1, AxBy;
  TWgt sc;
  int xMax, yMax, xs, xe, ys, ye, x, y, dvy;
  bool steep;

  xMax  = image->getWidth()-1;   // guard for x index
  yMax  = image->getHeight()-1;  // guard for y index
  dx    = x1 - x0;               // x-distance
  dy    = y1 - y0;               // y-distance
  L     = sqrt(dx*dx + dy*dy);   // length of the line
  steep = abs(dy) > abs(dx); 
  if(steep){                     // swap roles of x and y for steep lines
    swap(dx, dy);
    swap(x0, y0);
    swap(x1, y1);
    swap(xMax, yMax); }
  if(x0 > x1){                   // swap roles of start and end for leftward lines
    swap(x0, x1);
    swap(y0, y1); 
    dx = -dx;
    dy = -dy; }

  // slope/intercept equation y = a*x + b coefficients for main line:
  a = dy / dx;
  b = y0 - a*x0;

  // implicit line equation A*x + B*y + C = 0 coeffs for perpendicular lines through the endpoints:
  A  = dx / L;         // A and B are the same for both endpoints, we also have A^2 + B^2 = 1, so
  B  = dy / L;         // A*x + B*y + C = d gives the signed distance of any x,y to the line
  C0 = -(A*x0 + B*y0); // C coeff for left endpoint
  C1 = -(A*x1 + B*y1); // C coeff for right endpoint

  // end-cap extension:
  d = w2;   
  if(!roundCaps)
    d *= (abs(dx)+abs(dy))/L;

  // main loop:
  xs  = rsLimit((int)floor(x0-d), 0, xMax);   // start x-index 
  xe  = rsLimit((int)ceil( x1+d), 0, xMax);   // end x-index
  dvy = (int)ceil(w2/A);                      // maximum vertical pixel distance from line
  for(x = xs; x <= xe; x++){                  // outer loop over x
    yf = a*x + b;                             // ideal y (float)
    y  = roundToInt(yf);                      // rounded y
    ys = rsMax(y-dvy, 0);                     // scanline start
    ye = rsMin(y+dvy, yMax);                  // scanline end
    for(y = ys; y <= ye; y++){                // inner loop over y-scanline
      dp = A * abs(yf-y);                     // perpendicuar pixel distance from line
      sc = lineProfile(dp, w2);               // intensity/color scaler
      AxBy = A*x + B*y;
      if((d = -AxBy - C0) > 0.f){              // left end cap
        if(roundCaps){
          d  = sqrt(dp*dp+d*d);
          sc = lineProfile(d, w2); }
        else
          sc *= lineProfile(d, w2); }
      if((d = AxBy + C1) > 0.f){              // right end cap
        if(roundCaps){
          d = sqrt(dp*dp+d*d);
          sc = lineProfile(d, w2);  }
        else
          sc *= lineProfile(d, w2); }
      plot(x, y, sc, steep);                  // plot pixel (may swap x,y according to "steep") 
    }// for y
  }// for x  
}


// profile functions:

template<class TPix, class TWgt, class TCor>
TWgt LineDrawer<TPix, TWgt, TCor>::profileFlat(TCor d, TCor w2)
{
  if(d <= w2-1)
    return 1;
  if(d <= w2)
    return w2-d;
  return 0;
}

template<class TPix, class TWgt, class TCor>
TWgt LineDrawer<TPix, TWgt, TCor>::profileLinear(TCor d, TCor w2)
{
  if(d > w2)
    return 0;
  return (w2-d)/w2;
}

template<class TPix, class TWgt, class TCor>
TWgt LineDrawer<TPix, TWgt, TCor>::profileParabolic(TCor d, TCor w2)
{
  TCor x = d/w2;
  if(abs(x) > 1)
    return 0;
  return 1 - x*x;
}

template<class TPix, class TWgt, class TCor>
TWgt LineDrawer<TPix, TWgt, TCor>::profileCubic(TCor d, TCor w2)
{
  TCor x = d/w2;
  return rsPositiveBellFunctions<TCor>::cubic(abs(x));
}