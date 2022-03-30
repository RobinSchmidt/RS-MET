typedef rsVector2DF Vec2;
typedef std::vector<Vec2> ArrVec2;

//-------------------------------------------------------------------------------------------------
// Utilities

float pixelCoverage(float x, float y, Vec2 a, Vec2 b, Vec2 c)
{
  ArrVec2 triangle = { a, b, c };
  ArrVec2 square   = { Vec2(x, y), Vec2(x, y+1), Vec2(x+1, y+1), Vec2(x+1, y) };
  ArrVec2 polygon  = clipPolygon(triangle, square);
  return abs(polygonSignedArea(polygon));
}
float pixelCoverage(int x, int y, const rsVector2DF& a, const rsVector2DF& b, 
  const rsVector2DF& c)
{
  return pixelCoverage((float) x, (float) y, a, b, c);
}
// make a simplified version of the polygon clipping algorithm to clip a triangle against a pixel
// take advantage of the simplicity of the shapes to optimize away unnecessary operations
// (some of the vector elements become 0 or 1 -> allows to remove the respective additions and 
// multiplications) ...maybe it can be based on the implicit line equations - maybe that would make
// it even simpler? ...more work to do...i started doing this in Polygon.cpp (not yet finished)

// todo: use center of mass of the polygon and de-interpolate it (maybe we need to take care of the
// convention for pixel coords (center or corner)

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
  float thickness, int /*endCaps*/)
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
  // the ceil-function together with the offset of d=0.5 amounts to the top-left rule
  // maybe replace the constat d = 0.5 with two variables dx, dy - this lets the user decide, 
  // where inside the pixel the sampling point is
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
// for production code, instead of just passing a color value, pass a std::function object
// that computes a color - it should implement the () operator which should return the color
// -> allows for all kinds of shading algorithms (texture, lighting, bump, etc.)

void drawTriangleAntiAliasedProto(rsImageDrawerFFF& drw,
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, float color)
{
  rsImageF* img = drw.getImageToDrawOn();
  for(int y = 0; y < img->getHeight(); y++) {
    for(int x = 0; x < img->getWidth(); x++) {
      float coverage = pixelCoverage(x, y, a, b, c);
      drw.plot(x, y, coverage*color);
    }
  }
}

void drawTriangleAntiAliased(rsImageDrawerFFF& drw,
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, float color)
{
  int w = drw.getImageToDrawOn()->getWidth();
  int h = drw.getImageToDrawOn()->getHeight();
  int min = rsMin(w, h);
  int thresh = 10;  // currently arbitrarily chosen - todo: make performance tests and pick a good value
  if(min < thresh)
    drawTriangleAntiAliasedBoxBased(drw, a, b, c, color);
  else
    drawTriangleAntiAliasedSpanBased(drw, a, b, c, color);
}

void drawTriangleAntiAliasedBoxBased(rsImageDrawerFFF& drw,
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, float color)
{
  int w = drw.getImageToDrawOn()->getWidth();
  int h = drw.getImageToDrawOn()->getHeight();
  int xMin = rsMax(0,   rsMin(rsFloorInt(a.x), rsFloorInt(b.x), rsFloorInt(c.x)));
  int xMax = rsMin(w-1, rsMax(rsCeilInt( a.x), rsCeilInt( b.x), rsCeilInt( c.x)));  // one too much?
  int yMin = rsMax(0,   rsMin(rsFloorInt(a.y), rsFloorInt(b.y), rsFloorInt(c.y)));
  int yMax = rsMin(h-1, rsMax( rsCeilInt(a.y), rsCeilInt( b.y), rsCeilInt( c.y)));  // one too much?
  for(int y = yMin; y <= yMax; y++) {
    for(int x = xMin; x <= xMax; x++) {
      float coverage = pixelCoverage(x, y, a, b, c);
      drw.plot(x, y, coverage*color);
    }
  }
}

bool areTriangleVerticesOrdered(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c)
{
  bool r = true;
  r &= (a.y <= b.y && a.y <= c.y);
  r &= (b.x < c.x);  // what about degenerate cases - for example, all x coords the same?
  return r;
}
void orderTriangleVertices(rsVector2DF& a, rsVector2DF& b, rsVector2DF& c)
{
  if(a.y > b.y)
    RAPT::rsSwap(a, b);
  if(a.y > c.y)
    RAPT::rsSwap(a, c);
  if(b.x > c.x)
    RAPT::rsSwap(b, c);
  // is that all?

  rsAssert(areTriangleVerticesOrdered(a, b, c), "Vertex ordering failed");
  // ...under construction...
}
// for production code, maybe use a coordinate-comparison function pointer instead of the
// operators to make it work for pixel and world coordinates (pixel coordinates use a downward 
// y-axis) - these functions above work with pixel coordinates - but maybe we don't need this
// function for world coordinates?

// this does not work yet!!!
void drawTriangleScanlineSpans(int y, float sHere, float sNext, float eHere, float eNext,
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c,
  float color, rsImageDrawerFFF& drw, int w)
{
  // threshold for distance between left and right triangle edge below which we use a simple loop 
  // to compute coverages for all pixels in scanline
  float thresh = 2.f; // in pixels
  int x, xMin, xMax;
  float dxHere = eHere - sHere;
  float dxNext = eNext - sNext;
  if(rsMin(dxHere, dxNext) < thresh)
  {
    xMin  = rsClip((int) floor(sNext),    0, w-1);
    xMax  = rsClip((int) ceil(eNext) - 1, 0, w-1);
    for(x = xMin; x <= xMax; x++) 
      drw.plot(x, y, color*pixelCoverage(x, y, a, b, c));
    return;
  }
  // maybe with these distance thresholds, we can get rid of the separate top- and bottom scanline
  // code (if we treat them like all others, the lines will fall into this case)



  // 1st loop: partially covered pixels on the left side of the scanline - compute coverages:
  if(sNext < sHere)  {  // a top-left side
    xMin = (int)floor(sNext);
    xMax = (int)ceil( sHere) - 1;
  }
  else {                // a bottom-left side
    xMin = (int)floor(sHere);
    xMax = (int)ceil( sNext) - 1;
  }
  xMin = rsClip(xMin, 0, w-1);   // rsMax(xMin, 0)   is not enough
  xMax = rsClip(xMax, 0, w-1);   // rsMin(xMax, w-1) is not enough
  for(x = xMin; x <= xMax ; x++) 
    drw.plot(x, y, color*pixelCoverage(x, y, a, b, c));

  // 2nd loop: fully covered pixels in the middle of the scanline - use color as is:
  xMin = xMax+1;
  xMax = (int) floor(eHere) - 1;
  xMax = rsClip(xMax, 0, w-1);
  for(x = xMin; x <= xMax ; x++) 
    drw.plot(x, y, color);  // replace color variable by shading function-call

  // 3rd loop: partially covered pixels on the right side of the scanline - compute coverages:
  //xMin = xMax+1; // is this sometimes wrong? it seems, i may fail when xMin comes out > xMax which
  //               // indicates an empty span
  if(xMax >= xMin)
    xMin = xMax+1;
  else
  {
    return; // is that correct?
    // maybe something else to do? recompute xMin?
  }

  if(eNext > eHere)
    xMax = (int)ceil(eNext) - 1;  // a top-right side
  else
    xMax = (int)ceil(eHere) - 1;  // a bottom-right side
  xMax = rsClip(xMax, 0, w-1);
  for(x = xMin; x <= xMax ; x++) 
    drw.plot(x, y, color*pixelCoverage(x, y, a, b, c));
}

void drawTriangleAntiAliasedSpanBased(rsImageDrawerFFF& drw,
  const rsVector2DF& aIn, const rsVector2DF& bIn, const rsVector2DF& cIn, float color)
{
  typedef rsVector2DF Vec;
  Vec a = aIn, b = bIn, c = cIn;
  orderTriangleVertices(a, b, c);
  int w = drw.getImageToDrawOn()->getWidth();
  int h = drw.getImageToDrawOn()->getHeight();
  int yMin = rsMax(0,   rsMin(rsFloorInt(a.y), rsFloorInt(b.y), rsFloorInt(c.y)));
  int yMax = rsMin(h-1, rsMax( rsCeilInt(a.y), rsCeilInt( b.y), rsCeilInt( c.y)) - 1);

  // left edge (from a to b):
  Vec leftEdgeStart = a;
  Vec leftEdgeEnd   = b;

  // right edge (from a to c):
  Vec rightEdgeStart = a;
  Vec rightEdgeEnd   = c;

  // somwhere in the middle of the loop over the scanlines, we must switch either the left or right
  // triangle edge - here we figure out where and which:
  int breakLine;
  bool breakLeft;
  if(c.y > b.y) {
    breakLeft = true;             // left edge is broken in two pieces
    breakLine = (int)floor(b.y);  // between floor(b.y) and ceil(b.y)
  }
  else {
    breakLeft = false;            // right edge is broken in two pieces
    breakLine = (int)floor(c.y);  // between floor(c.y) and ceil(c.y)
  }

  // span-range variables:
  Vec sHere, eHere, sNext, eNext;
  int x, xMin, xMax, y = yMin;
  float yf = (float)y;

  // top scanline - compute coverages for every pixel because the only case where we could have
  // fully covered pixels in the top-row would be a flat-top triangle with an integer y-coordinate
  // for the flat top (but this may be not so unlikely - maybe treat this special case in a further
  // optimized way - pixel-aligned flat-top rectangles (and therefore also triangles) are actually 
  // a common case for handling windows):
  sNext = lineIntersection(Vec(0, yf+1), Vec(1, yf+1), a, b);  // intersection of next scanline with left edge
  eNext = lineIntersection(Vec(0, yf+1), Vec(1, yf+1), a, c);  // intersection of next scanline with right edge
  xMin  = rsClip((int) floor(sNext.x),    0, w-1);
  xMax  = rsClip((int) ceil(eNext.x) - 1, 0, w-1);
  for(x = xMin; x <= xMax; x++) 
    drw.plot(x, y, color*pixelCoverage(x, y, a, b, c));

  // top triangle area:
  yf = (float) (yMin+1);
  sHere = sNext;
  eHere = eNext;
  for(y = yMin+1; y < breakLine; y++)
  {
    yf = (float) y;
    sNext = lineIntersection(Vec(0, yf+1), Vec(1, yf+1), a, b);
    eNext = lineIntersection(Vec(0, yf+1), Vec(1, yf+1), a, c); 
    drawTriangleScanlineSpans(y, sHere.x, sNext.x, eHere.x, eNext.x, a, b, c, color, drw, w);
    sHere = sNext;
    eHere = eNext;
  }

  // middle scanline:
  yf = (float) (breakLine+1);
  if(breakLeft) {
    leftEdgeStart = b;
    leftEdgeEnd   = c;
    sHere = b;
    sNext = lineIntersection(Vec(0, yf), Vec(1, yf), leftEdgeStart, leftEdgeEnd);
    eHere = eNext;
    eNext = lineIntersection(Vec(0, yf), Vec(1, yf), rightEdgeStart, rightEdgeEnd); 
  }
  else {
    rightEdgeStart = c;
    rightEdgeEnd   = b;
    eHere = c;
    eNext = lineIntersection(Vec(0, yf), Vec(1, yf), rightEdgeStart, rightEdgeEnd);
    sHere = sNext;
    sNext = lineIntersection(Vec(0, yf), Vec(1, yf), leftEdgeStart, leftEdgeEnd);
  }
  drawTriangleScanlineSpans(y, sHere.x, sNext.x, eHere.x, eNext.x, a, b, c, color, drw, w);

  // bottom triangle area:
  yf = (float) (breakLine+1);
  sHere = sNext;
  eHere = eNext;
  for(y = breakLine+1; y < yMax; y++)
  {
    yf = (float) y;
    sNext = lineIntersection(Vec(0, yf+1), Vec(1, yf+1), leftEdgeStart,  leftEdgeEnd);
    eNext = lineIntersection(Vec(0, yf+1), Vec(1, yf+1), rightEdgeStart, rightEdgeEnd); 
    drawTriangleScanlineSpans(y, sHere.x, sNext.x, eHere.x, eNext.x, a, b, c, color, drw, w);
    sHere = sNext;
    eHere = eNext;
  }
 
  // bottom scanline:
  xMin = rsClip((int) floor(sHere.x),    0, w-1);
  xMax = rsClip((int) ceil(eHere.x) - 1, 0, w-1);
  for(x = xMin; x <= xMax; x++) 
    drw.plot(x, y, color*pixelCoverage(x, y, a, b, c));
}



int RectangleF::getRegion(const rsVector2DF& p)
{
  if(p.x < xMin) {    // left regions
    if(p.y < yMin)
      return BOTTOM_LEFT;
    if(p.y > yMax)
      return TOP_LEFT;
    return CENTER_LEFT;
  }
  if(p.x > xMax) {    // right regions
    if(p.y < yMin)
      return BOTTOM_RIGHT;
    if(p.y > yMax)
      return TOP_RIGHT;
    return CENTER_RIGHT;
  }
  if(p.y < yMin)
    return BOTTOM_CENTER;
  if(p.y > yMax)
    return TOP_CENTER;
  return INSIDE;
}

int clipTriangleToUnitSquare(const rsVector2DF& a, const rsVector2DF& /*b*/, const rsVector2DF& /*c*/,
  rsVector2DF* p)
{
  // assume a is top vertex and left in a flat-top triangle (with a downward y-axis)
  // so a has minimum y-coordinate)

  int nv = 0; // number of vertice in output polygon
  if(a.y > 1) 
    return 0;

  RectangleF r;
  int ra = r.getRegion(a);
  switch(ra)
  {
  case r.INSIDE: {      // may leave
    p[nv] = a; nv++; 
  } break;
  case r.TOP_RIGHT: {   // may enter from right or top edge
    // ...
  } break;
  case r.TOP_CENTER: {  // may enter from from top edge
    //...
  } break;
  case r.TOP_LEFT: {  // may enter from from top or left edge
    //...
  } break;

  }



  return nv;
}



// maybe the best way would be to make a simplified sutherland-hodgman version
// https://www.codeguru.com/cpp/misc/misc/graphics/article.php/c8965/Polygon-Clipping.htm
// that clips against all 4 edges separately with simplified tests, simplified intersection
// computations and avoiding dynamic memory allocations, here's a skeleton:
int clipAgainstTop(rsVector2DF* in, int N, rsVector2DF* out) // maybe rename to clipAgainstY1, Y1 means y=1
{
  int nv = 0;                   // number of output vertices
  rsVector2DF S = in[0], E;     // current start- and end vertex
  for(int i = 0; i < N; i++)
  {
    if(S.y > 1)             // outside top boundary
    {

    }
    else
    {
      out[nv] = S;
      nv++;

    }
  }
  
  return nv;
}

int clipTriangleToUnitSquare2(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c,
  rsVector2DF* out)
{
  // init:
  rsVector2DF p[7], q[7];  
  int np = 3, nq = 0;
  p[0] = a; p[1] = b; p[2] = c;     // p contains triangle with 3 vertices, q is empty

  // 7 vertices occur for example with a = (0.5,1.25), b = (-0.25,-0.25), c = (1.25,0.5)

  // clip against the 4 edges:
  nq = clipAgainstTop(   p, np, q);  // q contains partially clipped polygon
  //np = clipAgainstLeft(  q, nq, p);  // p contains partially clipped polygon
  //nq = clipAgainstBottom(p, np, q);  // q contains partially clipped polygon
  //np = clipAgainstRight( q, nq, p);  // p contains partially clipped polygon

  // copy result to output (for production code, get rid of local p, use out directly (rename to 
  // p)):
  for(int i = 0; i < np; i++)
    out[i] = p[i];
  return np;
}

//=================================================================================================

// Maybe move into a class rsPixelClassifier, maybe it could have the img as member

// maybe get rid of the repititive "Pixel" in the names, maybe make them members of the image class
bool isInteriorPixel(int i, int j, const rsImageF& img)
{
  return i > 0 && j > 0 && i < img.getWidth()-1 && j < img.getHeight()-1;
}
bool isLeftEdgePixel(int i, int j, const rsImageF& img)
{
  return i == 0;
}
bool isRightEdgePixel(int i, int j, const rsImageF& img)
{
  return i == img.getWidth()-1;
}
bool isTopEdgePixel(int i, int j, const rsImageF& img)
{
  return j == 0;
}
bool isBottomEdgePixel(int i, int j, const rsImageF& img)
{
  return j == img.getHeight()-1;
}
bool isCornerPixel(int i, int j, const rsImageF& img)
{
  return i == 0                && j == 0                   // top-left
    ||   i == 0                && j == img.getHeight()-1   // bottom-left
    ||   i == img.getWidth()-1 && j == 0                   // top-right
    ||   i == img.getWidth()-1 && j == img.getHeight()-1;  // bottom-right
}

// Checks, if the given predicate is true for any of the neighboring pixels of pixel i,j. The 
// predicate should take the center pixel's value as first argument and the neighbor pixel's value
// as second argument and return true, if the predicate holds for this pair of pixels.
template<class P> // P: predicate
bool isTrueForAnyNeighbor_I(int i, int j, const rsImageF& img, P pred)
{
  rsAssert(isInteriorPixel(i, j, img), "Function is made only for interior pixels."); 
  // Trying to use it for boundary pixels will lead to an access violation. For boundary pixels,
  // separate implementations exist and higher level code is supposed to use these.

  float p = img(i, j);  // pixel value

  // Check against direct neighbors (maybe factor out):
  if( pred(p, img(i-1,j)) ) return true;
  if( pred(p, img(i+1,j)) ) return true;
  if( pred(p, img(i,j-1)) ) return true;
  if( pred(p, img(i,j+1)) ) return true;

  // Check against diagonal neighbors (maybe factor out):
  if( pred(p, img(i-1,j-1)) ) return true;
  if( pred(p, img(i-1,j+1)) ) return true;
  if( pred(p, img(i+1,j-1)) ) return true;
  if( pred(p, img(i+1,j+1)) ) return true;

  // Predicate holds for all neighbor pixels in 3x3 neighborhood:
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_T(int i, int j, const rsImageF& img, P pred)
{
  // Code is the same as in isTrueForAllNeighbors_I but all lines involving j-1 have been deleted 
  // because in the top row, j-1 is not a valid y-coordinate
  rsAssert(isTopEdgePixel(i, j, img), "Made for top-edge pixels");
  rsAssert(!isCornerPixel(i, j, img), "Not made for corner pixels");
  float p = img(i, j); 
  if( pred(p,img(i-1,j  )) ) return true;
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i-1,j+1)) ) return true;
  if( pred(p,img(i+1,j+1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_B(int i, int j, const rsImageF& img, P pred)
{
  // Lines with j+1 have been removed:
  rsAssert(isBottomEdgePixel(i, j, img), "Made for bottom-edge pixels");
  rsAssert(!isCornerPixel(i, j, img), "Not made for corner pixels");
  float p = img(i, j);
  if( pred(p,img(i-1,j  )) ) return true;
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j-1)) ) return true;
  if( pred(p,img(i-1,j-1)) ) return true;
  if( pred(p,img(i+1,j-1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_L(int i, int j, const rsImageF& img, P pred)
{
  // Lines with i-1 have been removed:
  rsAssert(isLeftEdgePixel(i, j, img), "Made for left-edge pixels");
  rsAssert(!isCornerPixel(i, j, img), "Not made for corner pixels");
  float p = img(i, j);
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j-1)) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i+1,j-1)) ) return true;
  if( pred(p,img(i+1,j+1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_R(int i, int j, const rsImageF& img, P pred)
{
  // Lines with i+1 have been removed:
  rsAssert(isRightEdgePixel(i, j, img), "Made for right-edge pixels");
  rsAssert(!isCornerPixel(i, j, img), "Not made for corner pixels");
  float p = img(i, j);
  if( pred(p, img(i-1,j  )) ) return true;
  if( pred(p, img(i,  j-1)) ) return true;
  if( pred(p, img(i,  j+1)) ) return true;
  if( pred(p, img(i-1,j-1)) ) return true;
  if( pred(p, img(i-1,j+1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_TL(const rsImageF& img, P pred)
{
  // Lines with i-1 and j-1 have been removed:
  int i = 0;
  int j = 0;
  float p = img(i, j);
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i+1,j+1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_TR(const rsImageF& img, P pred)
{
  // Lines with i+1 and j-1 have been removed:
  int i = 0;
  int j = img.getWidth()-1;  // maybe use img.getRight()
  float p = img(i, j);
  if( pred(p,img(i-1,j  )) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i-1,j+1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_BL(const rsImageF& img, P pred)
{
  // Lines with i+1 and j-1 have been removed:
  int i = img.getHeight()-1;  // maybe use img.getBottom
  int j = 0;
  float p = img(i, j);
  if( pred(p, img(i-1,j  )) ) return true;
  if( pred(p, img(i,  j+1)) ) return true;
  if( pred(p, img(i-1,j+1)) ) return true;
  return false;
}

template<class P> 
bool isTrueForAnyNeighbor_BR(const rsImageF& img, P pred)
{
  // Lines with i+1 and j+1 have been removed:
  int i = img.getHeight()-1; 
  int j = img.getWidth()-1; 
  float p = img(i, j);
  if( pred(p, img(i-1,j  )) ) return true;
  if( pred(p, img(i,  j-1)) ) return true;
  if( pred(p, img(i-1,j-1)) ) return true;
  return false;
}

// Classifies the pixels as belongnig to class c where the given two-pixel predicate holds for the
// given pixel at coordinates (i,j) with respect to any of its neighbor pixels. That means, if for 
// a pixel at coordinates i,j the predicate holds for the pixel and any of its neighbors, The 
// corresponding C(i,j) element is set to c.
template<class P> 
void classifyWhenTrueForAnyNeighbor(const rsImageF& img, rsImage<char>& C, char c, P pred)
{
  rsAssert(C.hasSameShapeAs(img));
  int w = img.getWidth();
  int h = img.getHeight();

  // Classify interior pixels:
  for(int j = 1; j < h-1; j++) 
    for(int i = 1; i < w-1; i++) 
      if(isTrueForAnyNeighbor_I(i, j, img, pred)) 
        C(i, j) = c;

  // Classify edge pixels (excluding corners):
  for(int i = 1; i < w-1; i++) { 
    if(isTrueForAnyNeighbor_T(   i,   0,   img, pred)) C(i, 0  ) = c;    // top row
    if(isTrueForAnyNeighbor_B(   i,   h-1, img, pred)) C(i, h-1) = c; }  // bottom row
  for(int j = 1; j < h-1; j++) {    
    if(isTrueForAnyNeighbor_L(   0,   j,   img, pred)) C(0,   j) = c;    // left column
    if(isTrueForAnyNeighbor_R(   w-1, j,   img, pred)) C(w-1, j) = c; }  // right column

  // Classify corner pixels:
  if(isTrueForAnyNeighbor_TL(img, pred)) C(0,   0  ) = c;
  if(isTrueForAnyNeighbor_TR(img, pred)) C(w-1, 0  ) = c;
  if(isTrueForAnyNeighbor_BL(img, pred)) C(0,   h-1) = c;
  if(isTrueForAnyNeighbor_BR(img, pred)) C(w-1, h-1) = c;
}
// I think, we may need two versions: classifyWhenTrue, classifyWhenFalse and in the wheFalse 
// version we should negate the return values from isTrueFor... This is not the same as when the 
// caller just passes the negated predicate! Maybe we also could need WhenFalseForAny, 
// WhenTrueForAll, WhenFalseForAll...i think, we really need a unit test for this function




 // ToDo:
 // -Refactor also all the edge-case variations this function below to take a predicate, maybe 
 //  append _I, _T, _L, _B, _R, _TL, _TR, _BL, _BR for interior, top, left, etc
bool isFlatInterior3x3(int i, int j, const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_I(i, j, img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}
bool isFlatTop3x3(int i, int j, const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_T(i, j, img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}
bool isFlatBottom3x3(int i, int j, const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_B(i, j, img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}
bool isFlatLeft3x3(int i, int j, const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_L(i, j, img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}
bool isFlatRight3x3(int i, int j, const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_R(i, j, img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}
bool isFlatTopLeft3x3(const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_TL(img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}

bool isFlatTopRight3x3(const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_TR(img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}

bool isFlatBottomLeft3x3(const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_BL(img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}

bool isFlatBottomRight3x3(const rsImageF& img, float tol = 0.f)
{
  return !isTrueForAnyNeighbor_BR(img, [=](float p, float n){ return rsAbs(p-n) > tol; } );
}


void classifyFlatPixels3x3(const rsImageF& img, rsImage<char>& C, char F, float tol)
{
  // nope! that doesn't work:
  //classifyPixels(img, C, F,  [=](float p, float n){ return rsAbs(p-n) > tol; } ); return;
  //classifyPixels(img, C, F,  [=](float p, float n){ return !rsAbs(p-n) > tol; } ); return;


  rsAssert(C.hasSameShapeAs(img));
  int w = img.getWidth();
  int h = img.getHeight();

  // Classify interior pixels:
  for(int j = 1; j < h-1; j++) 
    for(int i = 1; i < w-1; i++) 
      if(isFlatInterior3x3(i, j, img, tol)) 
        C(i, j) = F;

  /*
  // The code below seems to make gradientifyFlatRegions fail, but the PixelClasses.ppm looks 
  // actually right, so the code itself seems right but using it seems to break something on a
  // higher level:
  // Classify edge pixels (excluding corners):
  for(int i = 1; i < w-1; i++) { 
    if(isFlatTop3x3(   i,   0,   img, tol))  C(i, 0  ) = F;    // top row
    if(isFlatBottom3x3(i,   h-1, img, tol))  C(i, h-1) = F; }  // bottom row
  for(int j = 1; j < h-1; j++) {    
    if(isFlatLeft3x3(  0,   j,   img, tol))  C(0,   j) = F;    // left column
    if(isFlatRight3x3( w-1, j,   img, tol))  C(w-1, j) = F; }  // right column

  // Classify corner pixels:
  if(isFlatTopLeft3x3(    img, tol))  C(0,   0  ) = F;
  if(isFlatTopRight3x3(   img, tol))  C(w-1, 0  ) = F;
  if(isFlatBottomLeft3x3( img, tol))  C(0,   h-1) = F;
  if(isFlatBottomRight3x3(img, tol))  C(w-1, h-1) = F;
  */
}

std::vector<rsVector2D<int>> findAll(const rsImage<char>& C, char c)
{
  std::vector<rsVector2D<int>> r; // result
  for(int j = 0; j < C.getHeight(); j++)   
    for(int i = 0; i < C.getWidth(); i++) 
      if(C(i,j) == c)   
        r.push_back(rsVector2D<int>(i,j));
  return r;
}



int gradientifyFlatRegions(const rsImageF& in, rsImageF& out, int numPasses) 
{
  using Vec2D = rsVector2D<int>;

  int maxIts = 200;   // make parameter, maybe return the number of iterations taken
  float tol  = 1.e-5f;
  int w = in.getWidth();
  int h = in.getHeight();
  int i, j, k;                   // loop iteration indices
  out.copyPixelDataFrom(in);     // initialize output - do we need this?
  //writeImageToFilePPM(out, "AfterInit.ppm");  // for debug

  //...............................................................................................
  // Step 1: Extract the coordinates of pixels in flat-color regions and on boundaries between such
  // regions. We record this information in two ways: (1) as arrays of pixel coordinates to 
  // facilitate iterating over the subsets and (2) as a matrix of char values that stores for each
  // pixel its class (encoded as symbolic constants defined below), to facilitate to access the 
  // classification in O(1) during the iteration. Pixels that do not belong into either of these 
  // classes are of no interest and classified as "rest" and we don't record their coordinates. We
  // only take pixels into account that are not at the image's edges because we need to work with
  // the neighbors of the pixels.
  std::vector<Vec2D> F, B;           // sets of (F)lat, (B)oundary
  rsImage<char> C(w, h);             // pixel classification matrix
  static const char rest     = 0;    // symbolic constants used in the code below, the values are
  static const char flat     = 100;  // ..also used as grayscale to encode the classes in an image 
  static const char boundary = 175;  // ..that can be written to disk for debug purposes
  C.fillAll(rest);                   // initially, all are "rest"

  // Identify pixels that belong to flat color regions:
  classifyFlatPixels3x3(in, C, flat);  // C(i,j) == flat iff in(i,j) belongs to flat region
  F = findAll(C, flat);                // F stores coordinates of all pixels in flat regions
  auto isFlat = [&](int i, int j) { return C(i,j) == flat; }; // for convenience
  // hmm..this seems to give a different result than in the old implementation - but the old one 
  // was crap anyway. still - this part should actually have worked the same

  // In a second pass, we identify the pixels at the boundary of flat regions;
  auto hasFlatNeighbor = [&](int i, int j)
  {
    if(isFlat(i-1, j  ) || isFlat(i+1, j  )) return true;
    if(isFlat(i,   j-1) || isFlat(i,   j+1)) return true;
    if(isFlat(i-1, j-1) || isFlat(i-1, j+1)) return true;
    if(isFlat(i+1, j-1) || isFlat(i+1, j+1)) return true;
    return false;
    // Maybe factor out into isAdjacentTo(int i, int j, const rsImage<char>& C, char c) or something
  };
  auto isAtBoundarySlow = [&](int i, int j, const rsImageF& img)
  {
    return (!isFlat(i, j)) && hasFlatNeighbor(i, j);
  };
  for(j = 1; j < h-1; j++) {
    for(i = 1; i < w-1; i++) {
      if(isAtBoundarySlow(i, j, in)) {
        B.push_back(Vec2D(i,j));
        C(i,j) = boundary;  }}} 
  auto isAtBoundary = [&](int i, int j) { return C(i,j) == boundary; }; // now faster!
  writeImageToFilePPM(C, "PixelClasses.ppm");  
  // For debug: pixels at the image edges are not correctly classified when they belong to a 
  // boundary. It's because we loop hetre only over the interior pixels

  //...............................................................................................
  // Step 2: For all boundary pixels: replace them by the average of those of their neighbors which
  // are also boundary pixels. The pixel itself is also included in that average:
  auto isRelevant = [&](int i, int j) 
  { 
    return isAtBoundary(i, j);
    //return !isFlat(i, j);  // test
    //return true;  // may also be useful. maybe provide different modes
  }; 
  for(k = 0; k < (int)B.size(); k++)
  {
    i = B[k].x;
    j = B[k].y;

    // Accumulate sum of the relevant neighbors:
    int   n = 0;    // number of relevant neighbors
    float s = 0.f;  // sum of colors of relevant neighbors
    if(isRelevant(i-1, j  )) { s += in(i-1, j  ); n += 1; }
    if(isRelevant(i+1, j  )) { s += in(i+1, j  ); n += 1; }
    if(isRelevant(i,   j-1)) { s += in(i,   j-1); n += 1; }
    if(isRelevant(i,   j+1)) { s += in(i,   j+1); n += 1; }
    if(isRelevant(i-1, j-1)) { s += in(i-1, j-1); n += 1; }
    if(isRelevant(i-1, j+1)) { s += in(i-1, j+1); n += 1; }
    if(isRelevant(i+1, j-1)) { s += in(i+1, j-1); n += 1; }
    if(isRelevant(i+1, j+1)) { s += in(i+1, j+1); n += 1; }

    // Compute the average, including the pixel at (i,j), and assign it to output image pixel:
    s += in(i, j);
    n += 1;
    float a = s / float(n);  // average
    out(i, j) = a;
  }
  writeImageToFilePPM(out, "AfterStep2.ppm");

  //...............................................................................................
  // Step 3: Alternatingly do to the flat-region pixels and boundary pixels: iteratively replace 
  // their values by the average of their neighbors until convergence:

  // Helper function. Computes difference of pixel value with respect to neighborhood average and
  // updates it to get get closer to that average. Returns the computed difference:
  auto applyFilter = [](const rsImageF& in, rsImageF& out, int i, int j, float amount = 1.f)
  {
    float avg;
    avg  = in(i,   j-1) + in(i,   j+1) + in(i-1, j  ) + in(i+1, j  );
    avg += in(i-1, j-1) + in(i-1, j+1) + in(i+1, j-1) + in(i+1, j+1);
    avg *= 1.f/8.f;
    float d = in(i,j) - avg;
    out(i,j) = in(i,j) - amount * d;
    return d;
    // Factor out
  };

  int maxItsTaken = 0;
  float step = 1.f;    // may not be needed, maybe get rid - but first, let's experiment with it
                       // a little bit to see if it can be used to accelerate convergence

  for(int i = 1; i <= numPasses; i++)
  {
    int its;
    float dMax, d;

    for(its = 0; its < maxIts; its++)                 // iteration over flat-region
    {
      dMax = 0.f;                                     // maximum delta applied
      for(k = 0; k < F.size(); k++) {
        d = applyFilter(out, out, F[k].x, F[k].y, step);
        dMax = rsMax(d, dMax);   }
      if(dMax <= tol)                                 // Check convergence criterion
        break;
    }
    maxItsTaken = rsMax(maxItsTaken, its);

    for(its = 0; its < maxIts; its++)                 // iteration over boundary
    {
      dMax = 0.f;
      for(k = 0; k < B.size(); k++) {
        d = applyFilter(out, out, B[k].x, B[k].y, step);
        dMax = rsMax(d, dMax);   }
      if(dMax <= tol)
        break;
    }
    maxItsTaken = rsMax(maxItsTaken, its);
  }
  writeImageToFilePPM(out, "AfterStep3.ppm");

  return maxItsTaken;

  // -Can we speed up the convergence? maybe it's actually not such a good idea to work in place?
  //  Try to compute all the updates first and then do all the upates at once. compare convergence
  //  to what we do now (updating every pixel immediately after the update was computed, such that 
  //  into computation of the next pixel, the current pixel enters already with updated color)
  // -It seems to have problems when the boundaries of the flat regions are anti-aliased. When 
  //  applied to the image generated by "contours", it doesn't seem to do much. The flat regions 
  //  are correctly classified. But i think, we have problems with the boundaries because the 
  //  different flat color regions do not border each other directly. There's a thin (1-pixel wide) 
  //  transition which is either a flat intermediate color (when anti-alias is turned off) or an
  //  actual gradienty thing from proper anti-aliasing. Even when anti-alias is turned off we do 
  //  some sort of crude, improper transition. Both of them trip up the algo.
  //  -The classification with and without AA looks almost the same, but there a very few pixels
  //   that get classified differently in both cases (i spotted just 1 in high zoom - they are 
  //   *really* rare)
  //  -I think the iteration doesn't really do very much in these cases because the flat region is
  //   already at the right color...but wait...no...this makes no sense
  //  -maybe we need a 3rd class of pixels: transition pixels. these are those which are neither in
  //   a flat region nor directly at its boundary but within the 1-pixel wide transition zone

  // Other idea:
  // -Instead of identifying boundary pixels and keeping them fixed during the iterations of the
  //  flat regions, identify a sort of central/neutral fiber within each flat region and keep 
  //  *that* fixed during iteration.
  // -This will also solve the problem with anti-aliased boundaries.
  // -It may also render the alternation between the flat and boundary iterations obsolete.
  // -To find that neutral fiber, we could proceed as follows:
  //  -Initialize a matrix of (half of squared) distances D to the boundary with zeros (or -1, 
  //   just some special unused value to encode: not yet computed/assigned, maybe INT_MAX is 
  //   convenient, because the algo may replace the current value with a min between current and 
  //   some new computed value)
  //  -Initialize a temp array with F, call it F' - it represents the flat-region pixels that are
  //   still left to be processed. In the process, this will shrink
  //  -For each boundary pixel with coords (i,j), do:
  //   -If the boundary is left/right/top/bottom (i.e straight), set D(i,j) = 1, else (boundary 
  //    is diagonally), set D(i,j) = 2
  //  -For each pixel in the flat region F' with coordinates i,j, do:
  //   -If any of its neighbors has D assigned, do
  //    -find minimum M of D among the neighbors
  //    -assign D(i,j) = M + 1 if min-neighbor is straight or M + 2 if min-neighbor is diag
  //    -remove pixel from F', maybe add it to another array of finished flat pixels, maybe
  //     give it a 3rd coordinate z representing D(i,j)
  //  -when done, go back to start (to "For each pixel in F'")
  // -Now we have for each pixel in the flat region a distance value. The neutral fiber within each 
  //  flat region is a path through that region that traverses the max-distance pixels.
  // -but damn! this is an O(N^2) algorithm - this is not practical!
  // -It could be improved, if we would not have to iterate over F' again and again. The problem is
  //  that in each iteration, most pixels will not have neighbors with assigned D. Can we make sure
  //  that we visit the pixels in an appropriate order such that we never (or rarely) encounter a 
  //  pixel with unassigned D(i,j) for its neighbors? But even then, the random-access removal from
  //  F' will be O(N) (within an O(N) loop)...maybe that removal can be scrapped, too? Or at least 
  //  turned into O(1) removal, i.e. removal from the back, or maybe a linked list is more 
  //  appropriate here?
  // -Maybe we should visit the pixels by going along the boundaries, thereby recording 2nd order
  //  boundaries. Then go along the 2nd order boundaries bulding up 3rd order boundaries and so on.
  //  This should avoid the O(N^2) scaling behavior. Maybe it needs two passes over each "layer" 
  //  because in the first pass, we may overestimate the true distance because some of the neighbor
  //  distances that should enter the min-computation are yet unassigned
  // -But it may not necessarily be a "fiber". It could also be just a single spot
  // -Maybe define an "edgeness" feature for each pixel as the maximum of the absolute differences
  //  between the pixel and its neighbors. It should be zero in flat regions. It could be used as
  //  a scaler for pixel updates in the iteration (or rather 1-edgeness then edgeness itself). 
  //  Pixels in flat regions would take an update step of 1, pixels in non-flat regions would 
  //  take smaller steps.

  // ToDo:
  // -Implement another variant of this algo that uses only 2 classes: flat regions and the rest.
  //  Here, the flat pixels are identified by having the same color as their 5x5 neighborhood 
  //  (maybe without the corners) rather than a 3x3 neighborhood. The "heat" equation iteration 
  //  is done only for the flat pixels without an special handling of boundary pixels. Any non-flat
  //  pixels are kept as is. Rationale: the 3 classes compicate the algorithm and don't seem to
  //  be beneficial...at least not for gradientifying the flat regions in the Newton fractal.
  // -Maybe for the fractal, we should just do the iteration for pixels that are either in flat 
  //  regions or at a boundary but should not handle these two classes seperately but as one class
  //  "flatOrBoundary"? Maybe then we could still get away with 3x3 neighborhoods? But maybe that 
  //  class would thene includ all pixels in the whole image? Maybe not. ..we want to maintain 
  //  the boundaries between regions that converge to a different attractor but smear the 
  //  boundaries between regions that converge to the same attractor but with different speed.
  //  Maybe we should just form pixel classes based on the attractor (directly in the fractal 
  //  rendering) and do the heat-equation iteration separately with each such region?
  // -Maybe before running the iteration, extend the image by repeating boundary pixels and after
  //  the iteration, crop back to the original size.

}

//=================================================================================================

void rsConvertImage(
  const rsImage<float>& R, const rsImage<float>& G, const rsImage<float>& B, bool clip,
  rsImage<rsPixelRGB>& img)
{
  using uchar = unsigned char;
  int w = R.getWidth();
  int h = R.getHeight();
  rsAssert(  G.hasShape(w, h));
  rsAssert(  B.hasShape(w, h));
  rsAssert(img.hasShape(w, h));
  if(clip) {
    for(int j = 0; j < h; j++) {
      for(int i = 0; i < w; i++) {
        img(i, j).r = (uchar)(255.f * rsClip(R(i, j), 0.f, 1.f));
        img(i, j).g = (uchar)(255.f * rsClip(G(i, j), 0.f, 1.f));
        img(i, j).b = (uchar)(255.f * rsClip(B(i, j), 0.f, 1.f)); }}}
  else {
    for(int j = 0; j < h; j++) {
      for(int i = 0; i < w; i++) {
        img(i, j).r = (uchar)(255.f * R(i, j));
        img(i, j).g = (uchar)(255.f * G(i, j));
        img(i, j).b = (uchar)(255.f * B(i, j)); }}}
}


unsigned char float2uchar(float x, bool clip)
{
  if(clip) return (unsigned char) (255.f * rsClip(x, 0.f, 1.f));
  else     return (unsigned char) (255.f * x);
}

void rsConvertImage(const rsImage<rsFloat32x4>& in, rsImage<rsPixelRGB>& out, bool clip)
{
  using uchar = unsigned char;
  int w = in.getWidth();
  int h = in.getHeight();
  rsAssert(out.hasShape(w, h));
  rsAssert(out.hasShape(w, h));
  for(int j = 0; j < h; j++) {
    for(int i = 0; i < w; i++) {
      out(i, j).r = float2uchar(in(i, j)[0], clip);
      out(i, j).g = float2uchar(in(i, j)[1], clip);
      out(i, j).b = float2uchar(in(i, j)[2], clip); }}
}

rsImage<rsPixelRGB> rsConvertImage(const rsImage<rsFloat32x4>& in, bool clip)
{
  rsImage<rsPixelRGB> out(in.getWidth(), in.getHeight());
  rsConvertImage(in, out, clip);
  return out;
}

rsImage<rsPixelRGB> rsConvertImage(
  const rsImage<float>& R, const rsImage<float>& G, const rsImage<float>& B, bool clip)
{
  int w = R.getWidth();
  int h = R.getHeight();
  rsImage<rsPixelRGB> img(w, h);
  rsConvertImage(R, G, B, clip, img);
  return img;
}
//template class rsImage<rsPixelRGB>; 

// todo: 
// -make a convertImage function that takes an rsImage<rsFloat32x4> as input. maybe we need
//  different variants for interpresting the 4 floats in different ways (RGBA, HSLA, etc.)
// -maybe have functions exctractColoChannels and combineColorChannels that convert between 4
//  images of float and one image rsFloat32x4
// -maybe wrap all this functions into a class rsImageConverter
// -function above could then be static members with names: floatToChar
