template<class TPix, class TWgt, class TCor>
rsImageDrawer<TPix, TWgt, TCor>::rsImageDrawer(rsImage<TPix> *imageToDrawOn)
{
  setImageToDrawOn(imageToDrawOn);
  setColor(TPix(1));  // typically white
  setBlendMode(0);
}

// setup:

template<class TPix, class TWgt, class TCor>
void rsImageDrawer<TPix, TWgt, TCor>::setImageToDrawOn(rsImage<TPix> *imageToDrawOn)
{
  image = imageToDrawOn;
}

template<class TPix, class TWgt, class TCor>
void rsImageDrawer<TPix, TWgt, TCor>::setBlendMode(int newMode)
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
void rsImageDrawer<TPix, TWgt, TCor>::linearBlend(TPix &pixel, TPix color, TWgt blend)
{
  pixel = (1-TPix(blend)) * pixel + TPix(blend) * color;
}

template<class TPix, class TWgt, class TCor>
void rsImageDrawer<TPix, TWgt, TCor>::addAndClip(TPix &pixel, TPix color, TWgt blend)
{
  pixel = rsMin(TPix(1), pixel + TPix(blend)*color);
}

template<class TPix, class TWgt, class TCor>
void rsImageDrawer<TPix, TWgt, TCor>::addAndSaturate(TPix &pixel, TPix color, TWgt blend)
{
  color *= TPix(blend);
  pixel  = (pixel + color) / (TPix(1) + color);

  //pixel  = (pixel + color) / (TPix(1) + abs(color)); // must be able to handle negative values
  // maybe try pixel = (pixel + color) / (TPix(1) + color*color); as alternative
}

//=================================================================================================

template<class TPix, class TWgt, class TCor>
rsLineDrawer<TPix, TWgt, TCor>::rsLineDrawer(rsImage<TPix> *imageToDrawOn)
  : rsImageDrawer<TPix, TWgt, TCor>(imageToDrawOn)
{
  setLineProfile(0);
  setLineWidth(1);
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::setLineProfile(int newProfile)
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
void rsLineDrawer<TPix, TWgt, TCor>::setLineWidth(TCor newWidth)
{
  rsAssert(newWidth > 0);
  w2 = 0.5f*(newWidth+1);   // +1 is a hack, otherwise lines are 1 pixel too narrow
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::setRoundCaps(bool shouldBeRound)
{
  roundCaps = shouldBeRound;
}

// drawing:

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::drawLine(TCor x0, TCor y0, TCor x1, TCor y1,
  bool joinableStart, bool joinableEnd)
{
  setupAlgorithmVariables(x0, y0, x1, y1);

  // hack/workaround to deal with zero-length lines (which result in division by zero):
  if(dx == 0)
    return;  // preliminary - we should draw a circle or rectangle (depending on cap setting)

  // draw segments:
  drawMiddleSection();
  drawLeftCap(joinableStart);
  drawRightCap(joinableEnd);

  // store line endpoint to be used as startpoint for subsequent lineTo calls:
  xOld = x1;
  yOld = y1;
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::lineTo(TCor x, TCor y, bool uniformColor)
{
  if(!uniformColor){
    drawLine(xOld, yOld, x, y, true, true);
    return; }

  setupAlgorithmVariables(xOld, yOld, x, y);
  if(dx == 0)
    return;   // todo: draw circle

  drawMiddleSection();
  if(steep){  // steep line
    drawCapForJointUniformColor(xs, xel, yOld, xOld);
    drawCapForJointUniformColor(xsr, xe, yOld, xOld); }
  else{       // flat line
    drawCapForJointUniformColor(xs, xel, xOld, yOld);
    drawCapForJointUniformColor(xsr, xe, xOld, yOld); }

  // \todo: make a version of it that always draws both caps (start and end) but with half
  // brightness/alpha - because the current startegy is only suitable for polylines in which each
  // adjacent line segment has the same color. if the color is changed between segments, artifacts
  // appear

  xOld = x;
  yOld = y;
}

// profile functions:

template<class TPix, class TWgt, class TCor>
TWgt rsLineDrawer<TPix, TWgt, TCor>::profileFlat(TCor d, TCor w2)
{
  if(d <= w2-1)
    return 1;
  if(d <= w2)
    return w2-d;
  return 0;
}

template<class TPix, class TWgt, class TCor>
TWgt rsLineDrawer<TPix, TWgt, TCor>::profileLinear(TCor d, TCor w2)
{
  if(d > w2)
    return 0;
  return (w2-d)/w2;
}

template<class TPix, class TWgt, class TCor>
TWgt rsLineDrawer<TPix, TWgt, TCor>::profileParabolic(TCor d, TCor w2)
{
  TCor x = d/w2;
  if(std::abs(x) > 1)
    return 0;
  return 1 - x*x;
}

template<class TPix, class TWgt, class TCor>
TWgt rsLineDrawer<TPix, TWgt, TCor>::profileCubic(TCor d, TCor w2)
{
  TCor x = d/w2;
  return rsPositiveBellFunctions<TCor>::cubic(std::abs(x));
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::setupAlgorithmVariables(TCor x0, TCor y0, TCor x1, TCor y1)
{
  xMax  = this->image->getWidth()-1;   // guard for x index
  yMax  = this->image->getHeight()-1;  // guard for y index
  dx    = x1 - x0;                     // x-distance
  dy    = y1 - y0;                     // y-distance
  L     = sqrt(dx*dx + dy*dy);         // length of the line
  steep = std::abs(dy) > std::abs(dx);
  if(steep){                           // swap roles of x and y for steep lines
    rsSwap(dx, dy);
    rsSwap(x0, y0);
    rsSwap(x1, y1);
    rsSwap(xMax, yMax); }
  back = x0 > x1;
  if(back){                            // swap roles of start and end for leftward lines
    rsSwap(x0, x1);
    rsSwap(y0, y1);
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
    d *= (std::abs(dx)+std::abs(dy))/L;
  xs  = rsClip((int)floor(x0-d), 0, xMax);    // start of left cap (and overall line)
  xel = rsClip((int)ceil( x0+d), 0, xMax);    // end of left cap
  xsr = rsClip((int)floor(x1-d), 0, xMax);    // start of right cap
  xe  = rsClip((int)ceil( x1+d), 0, xMax);    // end of right cap (and overall line)
  dvy = (int)ceil(w2/A);                      // maximum vertical pixel distance from line
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::drawMiddleSection()
{
  for(x = xel+1; x <= xsr-1; x++){            // outer loop over x
    yf = a*x + b;                             // ideal y (float)
    y  = rsRoundToInt(yf);                    // rounded y
    ys = rsMax(y-dvy, 0);                     // scanline start
    ye = rsMin(y+dvy, yMax);                  // scanline end
    for(y = ys; y <= ye; y++){                // inner loop over y-scanline
      dp = A * std::abs(yf-y);                     // perpendicuar pixel distance from line
      sc = lineProfile(dp, w2);               // intensity/color scaler
      plot(x, y, sc, steep);                  // plot pixel (may swap x,y according to "steep")
    }// for y
  }// for x
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::drawCap(int start, int end, bool joinable)
{
  TWgt js = 1;  // joint scaler
  if(joinable)
    js = 0.5;
  for(x = start; x <= end; x++){
    yf = a*x + b;
    y  = rsRoundToInt(yf);
    ys = rsMax(y-dvy, 0);
    ye = rsMin(y+dvy, yMax);
    for(y = ys; y <= ye; y++){
      dp = A * std::abs(yf-y);
      sc = lineProfile(dp, w2);
      AxBy = A*x + B*y;
      if((d = -AxBy - C0) > 0.f){             // pixel is left to left trimline (in left cap)
        if(roundCaps){
          d  = sqrt(dp*dp+d*d);
          sc = js*lineProfile(d, w2); }
        else
          sc *= lineProfile(d, w2); }
      if((d = AxBy + C1) > 0.f){              // right to right trimline (right end cap) - a pixel
        if(roundCaps){                        // may be in both, so we need to check both
          d = sqrt(dp*dp+d*d);                // conditions in both cap-drawing loops :-(
          sc = js*lineProfile(d, w2);  }
        else
          sc *= lineProfile(d, w2); }
      plot(x, y, sc, steep);
    }// for y
  }// for x
}

template<class TPix, class TWgt, class TCor>
void rsLineDrawer<TPix, TWgt, TCor>::drawCapForJointUniformColor(int start, int end,
  TCor xj, TCor yj)
{
  for(x = start; x <= end; x++){
    yf = a*x + b;
    y  = rsRoundToInt(yf);
    ys = rsMax(y-dvy, 0);
    ye = rsMin(y+dvy, yMax);
    for(y = ys; y <= ye; y++){
      dp = A * std::abs(yf-y);
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

      // additional joining code (optimize!):
      d = (x-xj)*(x-xj) + (y-yj)*(y-yj);
      if(d < w2*w2)
        sc *= 1-profileFlat(sqrt(d), w2);

      plot(x, y, sc, steep);
    }// for y
  }// for x
}

// rasterization of filled triangles:
// https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
// https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-stage

// https://www.youtube.com/watch?v=9A5TVh6kPLA
// https://github.com/planetchili/3D_Fundamentals