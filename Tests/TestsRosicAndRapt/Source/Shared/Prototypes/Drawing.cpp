//-------------------------------------------------------------------------------------------------
// Lines

// Sources:
// Anti Aliasing in general:
// http://hinjang.com/articles/01.html

// Gupta Sproull Algorithm:
// http://www.vis.uky.edu/~ryang/Teaching/cs535-2012spr/Lectures/18-anti-aliasing.pdf
// https://courses.engr.illinois.edu/ece390/archive/archive-f2000/mp/mp4/anti.html // 1-pixel

// Graphics Gems:
// http://www.graphicsgems.org/
// http://www.realtimerendering.com/resources/GraphicsGems/category.html#2D Rendering_link


// This is the Wu line drawing algorithm, directly translated from the pseudocode here:
// https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm
// todo: move to a Prototypes.h/cpp pair of files in the Shared folder
inline int   ipart(float x)           { return (int) x;              }
inline int   roundToInt(float x)      { return ipart(x + 0.5f);      }
inline float fpart(float x)           { return x - ipart(x);         }
inline float rfpart(float x)          { return 1 - fpart(x);         }
inline void  swap(float& x, float& y) { float t = x; x = y; y = t;   }
//inline float min(float x, float y)    { return x < y ? x : y; }
//inline void  plot(ImageF& im, int x, int y, float c){ im(x, y) += c; }
//inline void  plot(ImageF& im, int x, int y, float c){ im(x,y) = min(1.f, im(x,y)+c); }
inline void  plot(rsImageF& im, int x, int y, float c){ im(x,y) = (im(x,y)+c)/(1+c) ; }
void drawLineWuPrototype(rsImageF& img, float x0, float y0, float x1, float y1, float color)
{
  bool steep = abs(y1 - y0) > abs(x1 - x0);

  if(steep){
    swap(x0, y0);
    swap(x1, y1); }
  if(x0 > x1){
    swap(x0, x1);
    swap(y0, y1); }

  float dx = x1 - x0;
  float dy = y1 - y0;
  float gradient = dy / dx;
  if(dx == 0.0)
    gradient = 1.0;

  // handle first endpoint:
  int   xend  = roundToInt(x0);                     
  float yend  = y0 + gradient * (xend - x0);
  float xgap  = rfpart(x0 + 0.5f);
  int   xpxl1 = xend;                  // will be used in the main loop
  int   ypxl1 = ipart(yend);
  if(steep){
    plot(img, ypxl1,   xpxl1, rfpart(yend) * xgap * color);
    plot(img, ypxl1+1, xpxl1,  fpart(yend) * xgap * color); } 
  else {
    plot(img, xpxl1, ypxl1,   rfpart(yend) * xgap * color);
    plot(img, xpxl1, ypxl1+1,  fpart(yend) * xgap * color); }
  float intery = yend + gradient;      // first y-intersection for the main loop

  // handle second endpoint:  
  xend = roundToInt(x1);
  yend = y1 + gradient * (xend - x1);
  xgap = fpart(x1 + 0.5f);
  int xpxl2 = xend;                    // will be used in the main loop
  int ypxl2 = ipart(yend);
  if(steep){
    plot(img, ypxl2,   xpxl2, rfpart(yend) * xgap * color);
    plot(img, ypxl2+1, xpxl2,  fpart(yend) * xgap * color); }
  else {
    plot(img, xpxl2, ypxl2,   rfpart(yend) * xgap * color);
    plot(img, xpxl2, ypxl2+1,  fpart(yend) * xgap * color); }

  // main loop:
  if(steep){
    for(int x = xpxl1+1; x <= xpxl2-1; x++){
      plot(img, ipart(intery),   x, rfpart(intery) * color);
      plot(img, ipart(intery)+1, x,  fpart(intery) * color);
      intery = intery + gradient; }}
  else{
    for(int x = xpxl1+1; x <= xpxl2-1; x++){
      plot(img, x, ipart(intery),  rfpart(intery) * color);
      plot(img, x, ipart(intery)+1, fpart(intery) * color);
      intery = intery + gradient; }}
}

// Bresenham line drawing algorithm (integer arithmetic version):
void drawLineBresenham(rsImageF& img, int x0, int y0, int x1, int y1, float color)
{
  bool steep = abs(y1 - y0) > abs(x1 - x0);
  if(steep){
    swap(x0, y0);
    swap(x1, y1); }
  if(x0 > x1){
    swap(x0, x1);
    swap(y0, y1); }
  int ystep;
  if(y0 < y1)
    ystep = 1;
  else
    ystep = -1;

  int deltaX = x1 - x0;
  int deltaY = abs(y1 - y0);
  int error  = deltaX / 2;
  int y = y0;

  for(int x = x0; x <= x1; x++){
    if(steep)
      plot(img, y, x, color);
    else
      plot(img, x, y, color);
    error -= deltaY;
    if(error < 0){
      y += ystep;
      error += deltaX; }}
}
// https://en.wikipedia.org/wiki/Bresenham's_line_algorithm
// http://graphics.idav.ucdavis.edu/education/GraphicsNotes/Bresenhams-Algorithm.pdf

//-------------------------------------------------------------------------------------------------
// Thick Lines

// http://kt8216.unixcab.org/murphy/index.html
// http://www.zoo.co.uk/murphy/thickline/  // parallel line perpendicular or parallel to primary line
// http://ect.bell-labs.com/who/hobby/87_2-04.pdf // actually curves, no anti alias
// http://e-collection.library.ethz.ch/view/eth:3201
// http://e-collection.library.ethz.ch/eserv/eth:3201/eth-3201-01.pdf
// Ch3: geometric approach, Ch4: brush trajectory approach
// other ideas: scanlines - like Wu/Bresenham/Gupta but with a nested orthogonal loop over a
// a scanline orthogonal to the major direction, for example along y when line is x-dominant
// the interior of the scanline can just be filled, only at the ends we must compute coverages
// http://ect.bell-labs.com/who/hobby/thesis.pdf // Digitized Brush Trajectories (John Hobby)

void drawThickLineViaWu(rsImageF& img, float x0, float y0, float x1, float y1, float color, 
  float thickness)
{
  // This function tries to draw a thick line by drawing several parallel 1 pixel wide Wu lines. It 
  // doesn't work. There are artifacts from overdrawing pixels.

  if(thickness <= 1)
    drawLineWuPrototype(img, x0, y0, x1, y1, thickness*color); 

  float dx, dy, s, ax, ay, t2, x0p, y0p, x0m, y0m, x1p, y1p, x1m, y1m, xs, ys;

  // compute 4 corners (x0p,y0p), (x0m,y0m), (x1p,y1p), (x1m, y1m) of the rectangle reprenting the
  // thick line:
  dx  = x1 - x0;                 // (dx,dy) is a vector in the direction of our line
  dy  = y1 - y0;                 //
  s   = 1 / sqrt(dx*dx + dy*dy); // 1 / length(dx,dy)
  ax  = -dy*s;                   // (ax,ay) is a unit vector in a direction perpendicular to the
  ay  =  dx*s;                   // direction of our line
  t2  = 0.5f*(thickness-1);      // why -1? it works, but why?
  x0p = x0 + t2 * ax;
  y0p = y0 + t2 * ay;
  x0m = x0 - t2 * ax;
  y0m = y0 - t2 * ay;
  x1p = x1 + t2 * ax;
  y1p = y1 + t2 * ay;
  x1m = x1 - t2 * ax;
  y1m = y1 - t2 * ay;

  // draw lines:
  //drawLineWuPrototype(img, x0p, y0p, x1p, y1p, color); // line with max positive offset
  //drawLineWuPrototype(img, x0m, y0m, x1m, y1m, color); // line with max negative offset
  //drawLineWuPrototype(img, x0,  y0,  x1,  y1,  color); // center line
  // preliminary - we have to call this inside a loop several times to draw several parallel
  // Wu lines as explained for the case of Bresenham lines here: 
  // http://www.zoo.co.uk/murphy/thickline/. We use the same general principle her (the 2nd 
  // version, drawing lines parallel and stepping perpendicularly), with the only difference
  // that each line is a Wu line instead of a bresenham line. Maybe this can be further optimized
  // by drawing indeed Bresenham lines for the inner 1-pixel lines and using Wu lines only for
  // the 2 outermost lines?

  // draw lines:
  int numLines = (int)ceil(thickness);
  xs = (x0p-x0m) / (numLines); // step in x-direction
  ys = (y0p-y0m) / (numLines); // step in y-direction
  for(int i = 0; i < numLines; i++){
    dx = i*xs;
    dy = i*ys;
    drawLineWuPrototype(img, x0m+dx, y0m+dy, x1m+dx, y1m+dy, color); }
  // todo: scale color by line's intensity profile ...but it doesn't work properly yet
  // ther are artifacts we see a strange intensity profile pattern - maybe try drawing a dotted 
  // line instead of a Wu line? If that doesn't work either, i think, we need a scanline approach
  // that visits each pixel once
  // or maybe we are still using wrong stepsizes? -> experiment a bit 
  // xs = (x0p-x0m) / (numLines) is different from xs = (x0p-x0m) / (numLines-1) etc.
}

void plotLineWidth(rsImageF& img, int x0, int y0, int x1, int y1, float wd)
{ 
  // Adapted from http://members.chello.at/~easyfilter/bresenham.c. The setPixelColor calls had to
  // be modified, the original code was kept as comment.

  // plot an anti-aliased line of width wd
  int dx = abs(x1-x0), sx = x0 < x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0 < y1 ? 1 : -1;
  int err = dx-dy, e2, x2, y2;                                   // error value e_xy
  float ed = dx+dy == 0 ? 1 : sqrt((float)dx*dx+(float)dy*dy);

  for (wd = (wd+1)/2; ; ) {                                      // pixel loop
    plot(img, x0, y0, 1-rsMax(0.f, abs(err-dx+dy)/ed-wd+1));     //setPixelColor(x0, y0, max(0,255*(abs(err-dx+dy)/ed-wd+1)));
    e2 = err; x2 = x0;
    if (2*e2 >= -dx) {                                            // x step
      for (e2 += dy, y2 = y0; e2 < ed*wd && (y1 != y2 || dx > dy); e2 += dx)
        plot(img, x0, y2 += sy, 1-rsMax(0.f, abs(e2)/ed-wd+1));   //setPixelColor(x0, y2 += sy, max(0,255*(abs(e2)/ed-wd+1)));
      if (x0 == x1) break;
      e2 = err; err -= dy; x0 += sx;
    }
    if (2*e2 <= dy) {                                             // y step
      for (e2 = dx-e2; e2 < ed*wd && (x1 != x2 || dx < dy); e2 += dy)
        plot(img, x2 += sx, y0, 1-rsMax(0.f, abs(e2)/ed-wd+1));   //setPixelColor(x2 += sx, y0, max(0,255*(abs(e2)/ed-wd+1)));
      if (y0 == y1) break;
      err += dx; y0 += sy;
    }
  }

  // The left endpoint looks wrong.
}

float lineIntensity1(float d, float t2)
{
  // Computes intensity for a solid color line as function of the (perpendicular) distance d of a 
  // pixel from a line where t2 = thickness/2.
  if(d <= t2-1)
    return 1;
  if(d <= t2)
    return t2-d;
  return 0;
}
float lineIntensity2(float d, float t2)
{
  // Can be used alternatively to lineIntensity1 - instad of a solid color, the line has an 
  // intensity profile that decreases linearly from the center to the periphery of the line.
  if(d > t2)
    return 0;
  return (t2-d)/t2;
}
float lineIntensity3(float d, float t2)
{
  float x = d/t2;
  if(fabs(x) > 1)
    return 0;
  return 1 - x*x;
}
float lineIntensity4(float d, float t2)
{
  float x = d/t2;
  return rsPositiveBellFunctions<float>::cubic(fabs(x));
}
inline void plot(rsImageF& img, int x, int y, float color, bool swapXY)
{
  if(swapXY)
    plot(img, y, x, color);
  else
    plot(img, x, y, color);
}
void drawThickLine2(rsImageF& img, float x0, float y0, float x1, float y1, float color,
  float thickness, int endCaps)
{
  // ...Under construction...
  // Draws a thick line using a Bresenham stepper along the major axis in an outer loop and for 
  // each Bresenham pixel, it draws a scanline along the minor axis.

  // ToDo: 
  // -endCaps: 0: no special handling (hard cutoff of the main loop), 1: flat, 
  //  2: round (half circular)
  // -handle end caps properly (draw half circles) - to do that, we need to figure out, if the 
  //  current pixel belongs to the left (or right) end cap and if so, use the distance from the
  //  endpoint to the pixel (instead of the pixel-line distance). Even better would be to not use
  //  the lineProfile function but a corresponding dotProfile function (which, i think, should be
  //  the derivative of the lineProfile function)
  // -test with all possible cases

  thickness += 1.f; // hack, because the line seems one pixel too narrow

                    // adjustments for steep lines and negative slopes:
  bool steep = abs(y1 - y0) > abs(x1 - x0);
  if(steep){         // swap roles of x and y
    swap(x0, y0);
    swap(x1, y1); }
  if(x0 > x1){       // swap roles of start and end
    swap(x0, x1);
    swap(y0, y1); }
  int ystep;
  if(y0 < y1)        
    ystep = 1;
  else
    ystep = -1;

  // From http://graphics.idav.ucdavis.edu/education/GraphicsNotes/Bresenhams-Algorithm.pdf
  // page 9, deals with arbitrary (non-integer) endpoints. We also use a different convention for
  // the error - instead of having it between -1..0, we have it between -0.5..+0.5

  // variables for original Bresenham algo:
  float dx  =     x1 - x0;  // x distance, x1 >= x0 is already ensured
  float dy  = abs(y1 - y0); // y distance
  float s   = dy / dx;      // slope
  int   i0  = (int)x0;      // 1st x-index
  int   i1  = (int)x1;      // last x-index
  int   jb  = (int)y0;      // y-index in Bresenham algo
  float e   = -(1-(y0-jb)-s*(1-(x0-i0)))+0.5f; // Bresenham y-error, different from pdf by +0.5
                                               // -> check, if this formula is right

  // additional variables for thickness: 
  float L   = sqrt(dx*dx + dy*dy); // length of the line
  float sp  = dx / L;              // conversion factor between vertical and perpendicular distance
  float t2  = 0.5f*thickness;      // half-thickness
  int dj    = (int)ceil(t2/sp);    // maximum vertical pixel distance from line
  int jMax;
  if(steep)
    jMax = img.getWidth()-1;
  else
    jMax = img.getHeight()-1;

  // variables for end cap handling:
  float A  = dx / L;  // = sp
  float B  = dy / L;
  float C0 = -(A*x0 + B*y0);
  float C1 = -(A*x1 + B*y1);
  //tmp = A*A + B*B;  // for check - should be 1
  //int dummy = 0;
  // To handle the left end-cap, we take a line going through x0,y0 which is perpendicular to the 
  // main line, express this line with the implicit line equation A*x + B*y + C = 0. When the 
  // coefficients are normalized such that A^2 + B^2 = 1 (which is the case for the formulas 
  // above), then the right hand side gives the signed distance of any point x,y from the line. If 
  // this distance is negative for a given x,y, the point is left to the perpendicular line through
  // x0,y0 and belongs to the left cap so we need a different formula to determine its brightness. 
  // A similar procedure is used for right end cap, just that x1,y1 is used as point on the
  // perpendicular line and the rhs must be positive (the point must be to the right of the right
  // border line). For the right end cap, only the C coefficient is different, A and B are equal, 
  // so we have two C coeffs C0 for the left and C1 for the right endpoint.
  // ToDo - to get this right, we should extend the line by at most t2...but later...
  // Maybe we should wrap these formulas into a function 
  // lineTwoPointsToImplicit(T x0, T y0, T x1, T y2, T& A, T& B, T& C)

  // variables use in the loop:
  int j0, j1;
  float jr, dp, sc, d;


  // main loop to draw the line, stepping through the major axis (x):
  for(int i = i0; i <= i1; i++)
  {
    // Regular Bresenham algo would plot a pixel at (i,jb) for non-steep or (jb,i) for steep lines 
    // here. Instead, we plot a whole scanline along the minor j-direction, extending from jb-dj
    // to jb+dj:
    j0 = rsMax(jb-dj, 0);           // scanline start
    j1 = rsMin(jb+dj, jMax);        // scanline end
    jr = jb+ystep*e;                // reference minor coordinate to which y-distance is taken
    for(int j = j0; j <= j1; j++)   // loop over scanline
    {    
      // end cap handling:
      float AiBj = A*i + B*j;
      if((d = AiBj + C0) < 0.f){ // left end cap
                                 // something to do...
        continue;
      }
      if((d = AiBj + C1) > 0.f){ // right end cap
                                 // something to do
        continue;
      }

      // no cap, use perpendicular pixel/line distance:
      dp = sp * abs(jr-j);               // perpendicuar pixel distance from line
      sc = lineIntensity3(dp, t2);       // intensity/color scaler
      plot(img, i, j, sc*color, steep);  // color pixel (may swap i,j according to "steep") 
    }          

    if(e >= 0.5f){   // different from pdf by +0.5                                                
      jb += ystep;   // conditional Bresenham step along minor axis (y)
      e  -= 1; }     // error update
    e += s; 
  }

  // Here are some sources:
  // http://members.chello.at/~easyfilter/bresenham.html --> GOOD, complete with 100 page pdf and..
  // http://members.chello.at/~easyfilter/Bresenham.pdf  ...C sourcecode, demo program, etc.
  // The general idea there is to use the implicit equation of a curve F(x,y) = 0 as a distance
  // function. F(x,y) = d gives the distance of any point from the curve. This distance could
  // directly be used as input to a lineIntensityProfile function.

  // https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation
  // https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-stage
  // https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-practical-implementation
  // A pretty cool website about computer graphics, 2nd link explains the edge-function which may 
  // be useful for the end caps 

  // https://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html
}

void drawThickLine(rsImageF& img, float x0, float y0, float x1, float y1, float color,
  float thickness, bool roundCaps)
{
  float (*lineProfile)(float, float) = lineIntensity1;

  thickness += 1.f; // hack, because line are one pixel too narrow (why?)

  float dx, dy, a, b, t2, yf, dp, sc, d, L, A, B, C0, C1, AxBy;
  int xMax, yMax, xs, xe, ys, ye, x, y, dvy;
  bool steep;

  t2   = 0.5f*thickness;        // half thickness
  xMax = img.getWidth()-1;      // guard for x index
  yMax = img.getHeight()-1;     // guard for y index
  dx   = x1 - x0;               // x-distance
  dy   = y1 - y0;               // y-distance
  L    = sqrt(dx*dx + dy*dy);   // length of the line
  steep = abs(dy) > abs(dx); 
  if(steep){                    // swap roles of x and y for steep lines
    swap(dx, dy);
    swap(x0, y0);
    swap(x1, y1);
    swap(xMax, yMax); }
  if(x0 > x1){                  // swap roles of start and end for leftward lines
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
  d = t2;   
  if(!roundCaps)
    d *= (abs(dx)+abs(dy))/L;

  // main loop:
  xs  = rsClip((int)floor(x0-d), 0, xMax);    // start x-index 
  xe  = rsClip((int)ceil( x1+d), 0, xMax);    // end x-index
  dvy = (int)ceil(t2/A);                      // maximum vertical pixel distance from line
  for(x = xs; x <= xe; x++){                  // outer loop over x
    yf = a*x + b;                             // ideal y (float)
    y  = roundToInt(yf);                      // rounded y
    ys = rsMax(y-dvy, 0);                     // scanline start
    ye = rsMin(y+dvy, yMax);                  // scanline end
    for(y = ys; y <= ye; y++){                // inner loop over y-scanline
      dp = A * abs(yf-y);                     // perpendicuar pixel distance from line
      sc = lineProfile(dp, t2);               // intensity/color scaler
      AxBy = A*x + B*y;
      if((d = -AxBy - C0) > 0.f){              // left end cap
                                               //d = -d; 
        if(roundCaps){
          d  = sqrt(dp*dp+d*d);
          sc = lineProfile(d, t2); }
        else
          sc *= lineProfile(d, t2); }
      if((d = AxBy + C1) > 0.f){              // right end cap
        if(roundCaps){
          d = sqrt(dp*dp+d*d);
          sc = lineProfile(d, t2);  }
        else
          sc *= lineProfile(d, t2); }
      plot(img, x, y, sc*color, steep);       // color pixel (may swap x,y according to "steep") 
    }// for y
  }// for x    
}

//-------------------------------------------------------------------------------------------------
// Triangles

// code based on:  https://www.youtube.com/watch?v=9A5TVh6kPLA
// v0 should be left to v1 and v2 below the line connecting v0 and v1
void drawTriangleFlatTop(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color)
{
  static const float d = 0.5f;                // offset of a pixel's center from its index/coord
  float m0 = (v2.x - v0.x) / (v2.y - v0.y);   // inverse of slope of line from v0 to v2
  float m1 = (v2.x - v1.x) / (v2.y - v1.y);   // inverse of slope of line from v1 to v2
  for(int y = (int)ceil(v0.y-d); y < (int)ceil(v2.y-d); y++) { // loop over scanlines
    float px0 = m0 * (float(y) - v0.y + d) + v0.x;             // start x-coord
    float px1 = m1 * (float(y) - v1.y + d) + v1.x;             // end x-coord
    for(int x = (int)ceil(px0-d); x < (int)ceil(px1-d); x++)   // loop over pixels in scanline
      drw.plot(x, y, color);
  }
  // the ceil-function together with the offset of 0.5 amounts to the top-left rule
}
// maybe compute px0, px1 incrementally, i.e. init to v0.x, v1.x and incerement by dx0, dx1 in each
// iteration, where dx0 =
// but not in the prototype

// v1 should be left to v2 and v0 above the line connecting v1 and v2
void drawTriangleFlatBottom(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color)
{
  float m0 = (v1.x - v0.x) / (v1.y - v0.y);   // inverse of slope of line from v0 to v1
  float m1 = (v2.x - v0.x) / (v2.y - v0.y);   // inverse of slope of line from v0 to v2
  int ys = (int) ceil(v0.y - 0.5f);           // y-coord of first scanline
  int ye = (int) ceil(v2.y - 0.5f);           // y-coord of scanline after the last line drawn
  for( int y = ys; y < ye; y++ ) {            // loop over scanlines
    float px0 = m0 * (float(y) + 0.5f - v0.y) + v0.x; // start x-coord
    float px1 = m1 * (float(y) + 0.5f - v0.y) + v0.x; // end x-coord
    int xs = (int) ceil(px0 - 0.5f);                  // start pixel
    int xe = (int) ceil(px1 - 0.5f);                  // end pixel (after the last pixel drawn)
    for(int x = xs; x < xe; x++)                      // loop over pixels in current scanline
      drw.plot(x, y, color);
  }
}
// compactify this further (get rid of ys,ye,xs,xe, use d=0.5)

// i think, the anti-aliased version should let loop indices start at = floor(...) and end at
// <= ceil(...)

// maybe wrap into class rsPolygonDrawer

void drawTriangle(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color)
{
  // use pointers so we can swap (for sorting purposes)
  typedef rsVector2DF Vec2; // for convenience
  const Vec2* pv0 = &v0;
  const Vec2* pv1 = &v1;
  const Vec2* pv2 = &v2;
  // todo: use pointers as arguments - maybe provide convenience function that takes const 
  // references

  // sort vertices by y:
  if(pv1->y < pv0->y) std::swap(pv0, pv1);
  if(pv2->y < pv1->y) std::swap(pv1, pv2);
  if(pv1->y < pv0->y) std::swap(pv0, pv1);

  if(pv0->y == pv1->y) {        // triangle is flat top
    if(pv1->x < pv0->x) std::swap(pv0, pv1); // sort top vertices by x
    drawTriangleFlatTop(drw, *pv0, *pv1, *pv2, color);
  }
  else if(pv1->y == pv2->y) {   // triangle is flat bottom
    if(pv2->x < pv1->x) std::swap( pv1,pv2 ); // sort bottom vertices by x
    drawTriangleFlatBottom(drw, *pv0, *pv1, *pv2, color);
  }
  else {
    // split general triangle into flat-top and flat-bottom:
    const float alpha = (pv1->y - pv0->y) / (pv2->y - pv0->y);
    const Vec2 vi = *pv0 + alpha * (*pv2 - *pv0);    // splitting vertex by linear interpolation between v0 and v2
    if(pv1->x < vi.x) { // long side is on the right (major right)
      drawTriangleFlatBottom(drw, *pv0, *pv1,   vi, color);
      drawTriangleFlatTop(   drw, *pv1,   vi, *pv2, color);
    }
    else {              // long side is on the left (major left)
      drawTriangleFlatBottom(drw, *pv0,   vi, *pv1, color);
      drawTriangleFlatTop(   drw,   vi, *pv1, *pv2, color);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Polygons:

bool insideEdge(const rsVector2DF& p, const rsVector2DF& v1, const rsVector2DF& v2)
{
  // from v1, v2, obtain implicit line equation for the edge defined by v1, v2 and plug point
  // p into that line equation - sign of the result determines if the point is considered inside

  return true; // preliminary
}

rsVector2DF intersection(const rsVector2DF& p1, const rsVector2DF& p2,
  const rsVector2DF& q1, const rsVector2DF& q2)
{
  // finds intersection point of the two lines through p1,p2 and q1,q2

  return rsVector2DF();  // preliminary
}

// p: general polygon to be clipped, c: convex clipping polygon
std::vector<rsVector2DF> clipConvexPolygons(const std::vector<rsVector2DF>& p, 
  const std::vector<rsVector2DF>& c)
{
  std::vector<rsVector2DF> r; // result

  // something to do...

  return r;
}
// Sutherland-Hodgman algorithm (Foley, page 124ff)
// can p really be non-convex? in this case the output may have to be a set of polygons (one 
// non-convex could split into many polygons), see page 125
