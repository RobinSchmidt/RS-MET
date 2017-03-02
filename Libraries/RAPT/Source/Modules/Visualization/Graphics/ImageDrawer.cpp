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

  //pixel  = (pixel + color) / (TPix(1) + abs(color)); // must be able to handle negative values
  // maybe try pixel = (pixel + color) / (TPix(1) + color*color); as alternative
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
  w2 = 0.5f*(newWidth+1);   // +1 is a hack, otherwise lines are 1 pixel too narrow
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
  setupAlgorithmVariables(x0, y0, x1, y1);

  // hack/workaround to deal with zero-length lines (which result in division by zero):
  if(dx == 0)
    return;  // preliminary - we should draw a circle or rectangle (depending on cap setting)

  // draw segments:
  drawMiddleSection();
  drawLeftCap();
  drawRightCap();

  // store line endpoint to be used as startpoint for subsequent lineTo calls
  xOld = x1; 
  xOld = y1;
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::lineTo(TCor x, TCor y)
{
  setupAlgorithmVariables(xOld, yOld, x, y);
  if(dx == 0)
    return;   // todo: maybe draw circle

  drawMiddleSection();

  // draw caps, some logic is needed to figure out which cap is to be drawn wnd which isn't
  // because of potential swapping of start/end and/or x/y:
  if(!back){                                        // forward
    if(steep){                                      // forward and steep
      drawCapForJoint(xs, xel, yOld, xOld);
      drawCapForJoint(xsr, xe, y, x);        }
    else{                                           // forward and flat 
      drawCapForJoint(xs, xel, xOld, yOld);
      drawCapForJoint(xsr, xe, xOld, yOld);  }}
  else{                                             // backward
    if(steep){                                      // backward and steep 
      drawCapForJoint(xs, xel, yOld, xOld);
      drawCapForJoint(xsr, xe, y,    x);     }
    else{                                           // backward and flat
      drawCapForJoint(xs, xel, xOld, yOld);
      drawCapForJoint(xsr, xe, xOld, yOld);  }}

  // it still sometimes doesn't work (check with random lines) - maybe we have also take into 
  // account whether or not the previous line was forward or backward and/or steep or flat?


  xOld = x; 
  yOld = y;
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

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::setupAlgorithmVariables(TCor x0, TCor y0, TCor x1, TCor y1)
{
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
  back = x0 > x1;
  if(back){                      // swap roles of start and end for leftward lines
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

  // compute loop limits:
  d = w2;                                     // end-cap extension
  if(!roundCaps)
    d *= (abs(dx)+abs(dy))/L;
  xs  = rsLimit((int)floor(x0-d), 0, xMax);   // start of left cap (and overall line)
  xel = rsLimit((int)ceil( x0+d), 0, xMax);   // end of left cap
  xsr = rsLimit((int)floor(x1-d), 0, xMax);   // start of right cap
  xe  = rsLimit((int)ceil( x1+d), 0, xMax);   // end of right cap (and overall line)
  dvy = (int)ceil(w2/A);                      // maximum vertical pixel distance from line  

  //// store line endpoint to be used as startpoint (for lineTo calls):
  //this->x0 = x1; 
  //this->y0 = y1;
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::drawMiddleSection()
{
  for(x = xel+1; x <= xsr-1; x++){            // outer loop over x
    yf = a*x + b;                             // ideal y (float)
    y  = roundToInt(yf);                      // rounded y
    ys = rsMax(y-dvy, 0);                     // scanline start
    ye = rsMin(y+dvy, yMax);                  // scanline end
    for(y = ys; y <= ye; y++){                // inner loop over y-scanline
      dp = A * abs(yf-y);                     // perpendicuar pixel distance from line
      sc = lineProfile(dp, w2);               // intensity/color scaler
      plot(x, y, sc, steep);                  // plot pixel (may swap x,y according to "steep") 
    }// for y
  }// for x 
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::drawCap(int start, int end)
{
  for(x = start; x <= end; x++){
    yf = a*x + b;
    y  = roundToInt(yf);
    ys = rsMax(y-dvy, 0);
    ye = rsMin(y+dvy, yMax);
    for(y = ys; y <= ye; y++){
      dp = A * abs(yf-y);
      sc = lineProfile(dp, w2);
      AxBy = A*x + B*y;
      if((d = -AxBy - C0) > 0.f){             // pixel is left to left trimline (in left cap)
        if(roundCaps){
          d  = sqrt(dp*dp+d*d);
          sc = lineProfile(d, w2); }
        else
          sc *= lineProfile(d, w2); }
      if((d = AxBy + C1) > 0.f){              // right to right trimline (right end cap) - a pixel
        if(roundCaps){                        // may be in both, so we need to check both
          d = sqrt(dp*dp+d*d);                // conditions in both cap-drawing loops :-(
          sc = lineProfile(d, w2);  }
        else
          sc *= lineProfile(d, w2); }
      plot(x, y, sc, steep);
    }// for y
  }// for x  
}

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::drawCapForJoint(int start, int end, TCor xj, TCor yj)
{
  //TCor r2;
  //TCor r0, r1;
  TCor dj; // distance from join (we may actually reuse the regular d for that);
  for(x = start; x <= end; x++){
    yf = a*x + b;
    y  = roundToInt(yf);
    ys = rsMax(y-dvy, 0);
    ye = rsMin(y+dvy, yMax);
    for(y = ys; y <= ye; y++){
      dp = A * abs(yf-y);
      sc = lineProfile(dp, w2);
      AxBy = A*x + B*y;

      // regular end cap code:
      if((d = -AxBy - C0) > 0.f){             // pixel is left to left trimline (in left cap)
        if(roundCaps){
          d  = sqrt(dp*dp+d*d);
          sc = lineProfile(d, w2); }
        else
          sc *= lineProfile(d, w2); }
      if((d = AxBy + C1) > 0.f){              // right to right trimline (right end cap) - a pixel
        if(roundCaps){                        // may be in both, so we need to check both
          d = sqrt(dp*dp+d*d);                // conditions in both cap-drawing loops :-(
          sc = lineProfile(d, w2);  }
        else
          sc *= lineProfile(d, w2); }

      // additional joining code:
      dj = (x-xj)*(x-xj) + (y-yj)*(y-yj);
      if(dj < w2*w2)
      {
        sc *= 1-profileFlat(sqrt(dj), w2);
        //sc *= 0.0; // use anti-aliasing - call (1-profileFlat(sqrt(dj)) 
      }

      //// old:
      //r0 = w2*w2;
      //if(back)
      //  r0 = (x-xOld)*(x-xOld) + (y-yOld)*(y-yOld);
      //else
      //  r0 = (x-x1)*(x-x1) + (y-y1)*(y-y1);
      //if(r0 < w2*w2)
      //  sc *= 0.0;


      plot(x, y, sc, steep);
    }// for y
  }// for x  
}
