using namespace RAPT;

void lineDrawing()
{
  // Compares different line drawing algorithms. We draw lines of different directions.
  // todo:
  // Move the prototype implementations to another shared file, so they can also be used in unit
  // tests
  // drag the prototypes to ImagePainter, write unit tests, do optimzations, write performance
  // comparison between Bresenham, Wu and Dotted line

  // user parameters:
  int imageWidth   = 400;
  int imageHeight  = 400;
  int numLines     = 10;    // number of lines per direction (horizontal'ish and vertical'ish)
  float margin     = 10.f;
  float brightness = 0.5f;

  // create objects:
  rsImageF image(imageWidth, imageHeight);
  rsImagePainterFFF painter(&image, nullptr);
  float x0 = margin;
  float y0 = margin;
  int i;

  // create arrays for line endpoints:
  vector<float> x1, y1;
  for(i = 0; i < numLines; i++){    // flat, horizontal'ish
    x1.push_back(imageWidth - margin);
    y1.push_back(margin + i * (imageHeight - margin) / numLines); }
  x1.push_back(imageWidth -margin); // 45° diagonal
  y1.push_back(imageHeight-margin);
  for(i = 0; i < numLines; i++){    // steep, vertical'ish
    x1.push_back(margin + i * (imageWidth - margin) / numLines);
    y1.push_back(imageHeight - margin); }

  //// dotted algorithm:
  //for(i = 0; i < x1.size(); i++)
  //  painter.drawLineDotted(x0, y0, x1[i], y1[i], brightness, brightness, numDots); // we need to pass a number of dots
  //writeImageToFilePPM(image, "LinesDotted.ppm");

  // Wu algorithm:
  image.clear();
  for(i = 0; i < x1.size(); i++)
    painter.drawLineWu(x0, y0, x1[i], y1[i], brightness);
    //drawLineWuPrototype(image, x0, y0, x1[i], y1[i], brightness);
  writeImageToFilePPM(image, "LinesWu.ppm");

  // Bresenham algorithm:
  image.clear();
  for(i = 0; i < x1.size(); i++)
    drawLineBresenham(image, rsRoundToInt(x0), rsRoundToInt(y0), 
      rsRoundToInt(x1[i]), rsRoundToInt(y1[i]), brightness);
  writeImageToFilePPM(image, "LinesBresenham.ppm");
}

void lineDrawingThick()
{
  // user parameters:
  int imageWidth   = 800;
  int imageHeight  = 800;
  int numAngles     = 7;
  float brightness = 0.75f;
  float thickness  = 50.f;

  // create objects:
  rsImageF image(imageWidth, imageHeight);
  rsLineDrawerFFF drawer(&image);
  drawer.setColor(brightness);
  drawer.setBlendMode(rsImageDrawerFFF::BLEND_ADD_SATURATE);
  drawer.setLineWidth(thickness);
  drawer.setLineProfile(rsLineDrawerFFF::PROFILE_LINEAR);
  drawer.setRoundCaps(false);

  // create endpoint arrays:
  float margin = 2*thickness;
  int numLines = 2*numAngles - 2;
  vector<float> x0(numLines), y0(numLines), x1(numLines), y1(numLines);
  int i, j;
  for(i = 0; i < numAngles; i++){
    x0[i] = margin + i * (imageWidth - 2*margin) / (numAngles-1);
    y0[i] = margin;
    x1[i] = imageWidth - x0[i];
    y1[i] = imageHeight - margin; }
  for(i = 0; i < numAngles-2; i++){
    j = numAngles + i;
    y0[j] = margin + (i+1) * (imageHeight - 2*margin) / (numAngles-1);
    x0[j] = margin;
    y1[j] = imageHeight - y0[j];
    x1[j] = imageWidth - margin; }

  // draw the lines and save file:
  for(i = 0; i < numLines; i++)
  {
    drawer.drawLine(x0[i], y0[i], x1[i], y1[i]);
    //drawThickLine(image, x0[i], y0[i], x1[i], y1[i], brightness, thickness, true);
  }
  writeImageToFilePPM(image, "LinesThick.ppm");
}
void lineDrawingThick2()
{
  // user parameters:
  int imageWidth   = 100;
  int imageHeight  = 100;
  float brightness = 0.5f;
  float thickness  = 10.f;
  float x0 = 10.3f, y0 = 10.6f, x1 = 90.2f, y1 = 40.4f;
  //float x0 = 10, y0 = 10, x1 = 90, y1 = 40;

  rsImageF image(imageWidth, imageHeight);
  //drawThickLine(image, 10, 10, 70, 30, 1.f, 15.f); // dx > dy, x0 < x1, base case
  //drawThickLine(image, 20, 20, 80, 80, 1.f, 15.f); // dx = dy, x0 < x1, 45° diagonal
  //drawThickLine(image, 10, 10, 30, 70, 1.f, 15.f); // dx < dy, x0 < x1, steep case
  //drawThickLine(image, 70, 10, 10, 30, 1.f, 15.f); // dx > dy, x0 > x1, x-swap case
  //drawThickLine(image, 10, 30, 70, 10, 1.f, 15.f);
  //drawThickLine(image, 30, 10, 10, 70, 1.f, 15.f);
  //drawThickLine(image, 10, 10, 25, 15, 1.f,  8.f, true);
  //drawThickLine(image, 20, 20, 50, 30, 1.f, 16.f, false);
  drawThickLine(image, 20, 20, 50, 50, 1.f, 20.f, true);
  drawThickLine(image, 20, 50, 50, 80, 1.f, 20.f, false);
  //drawThickLine(image, 20, 70, 50, 70, 1.f, 5.f, false);
                                            
  //drawThickLine(image, 10, 10, 50, 90, 1.f, 15.f);
  //drawThickLine(image, x0, y0, x1, y1, brightness, thickness);
  //plotLineWidth(image, (int)x0, (int)y0, (int)x1, (int)y1, thickness);

  writeImageToFilePPM(image, "ThickLineScanlineTest.ppm");
}

void lineJoints()
{
  // user parameters:
  int imageWidth   = 800;
  int imageHeight  = 800;
  int numAngles    = 10;
  float brightness = 0.5f;
  float thickness  = 20.f;

  // create objects:
  rsImageF image(imageWidth, imageHeight);
  rsLineDrawerFFF drawer(&image);
  drawer.setBlendMode(rsImageDrawerFFF::BLEND_ADD_SATURATE);
  //drawer.setLineProfile(rsLineDrawerFFF::PROFILE_LINEAR);
  drawer.setLineProfile(rsLineDrawerFFF::PROFILE_FLAT);
  //drawer.setLineProfile(rsLineDrawerFFF::PROFILE_CUBIC);
  drawer.setLineWidth(thickness);
  drawer.setColor(brightness);

  // create endpoint arrays:
  // draw lines:
  float margin = 2*thickness;
  //vector<float> x0, y0, x1, y1;
  float x0, y0, x1, y1;
  float w2 = 0.5f*imageWidth;
  float h2 = 0.5f*imageHeight;
  float dy = (float)margin;
  float offset;
  for(int i = 0; i < numAngles; i++)
  {
    offset = i*margin;
    x0 = margin;
    y0 = margin + offset;
    x1 = w2;
    y1 = margin + i*dy + offset;
    drawer.drawLine(x0, y0, x1, y1);
    //drawThickLine(image, x0, y0, x1, y1, brightness, thickness, true);
    x0 = x1;
    y0 = y1;
    x1 = imageWidth - margin;
    y1 = margin + offset;
    drawer.drawLine(x0, y0, x1, y1);
    //drawThickLine(image, x0, y0, x1, y1, brightness, thickness, true);
  }

  writeImageToFilePPM(image, "LineJoints.ppm");
}

void lineTo()
{
  // Test the lineTo function by drawing 4 lines using lineTo for the 4 combinations of back/steep

  // user parameters:
  int size = 800;         // image width and height
  float brightness = 0.5f;
  float thickness  = 20.f;

  // create objects:
  rsImageF image(size, size);
  rsLineDrawerFFF drawer(&image);
  drawer.setBlendMode(rsImageDrawerFFF::BLEND_ADD_SATURATE);
  //drawer.setLineProfile(rsLineDrawerFFF::PROFILE_LINEAR);
  drawer.setLineProfile(rsLineDrawerFFF::PROFILE_FLAT);
  //drawer.setLineProfile(rsLineDrawerFFF::PROFILE_CUBIC);
  drawer.setLineWidth(thickness);
  drawer.setColor(brightness);

  float margin = 2*thickness;

  // flat, forward:
  drawer.initPolyLine(     margin, 0.4f*size);
  drawer.lineTo(  size-margin, 0.6f*size);

  // flat, backward:
  drawer.initPolyLine(size-margin, 0.4f*size);
  drawer.lineTo(       margin, 0.6f*size);

  // steep, forward:
  drawer.initPolyLine(0.4f*size,      margin);
  drawer.lineTo(  0.6f*size, size-margin);

  // steep, backward:
  drawer.initPolyLine(0.4f*size, size-margin);
  drawer.lineTo(  0.6f*size,      margin);

  writeImageToFilePPM(image, "LineTo.ppm");
}

void polyLineRandom()
{
  // user parameters:
  int imageWidth      = 800;
  int imageHeight     = 800;
  int numLines        = 50;
  float minBrightness = 0.125f;
  float maxBrightness = 1.0f;
  float thickness     = 20.f;

  // create objects:
  rsImageF image(imageWidth, imageHeight);
  rsLineDrawerFFF drawer(&image);
  //drawer.setBlendMode(rsImageDrawerFFF::BLEND_ADD_SATURATE);
  drawer.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);
  //drawer.setLineProfile(rsLineDrawerFFF::PROFILE_LINEAR);
  drawer.setLineProfile(rsLineDrawerFFF::PROFILE_FLAT);
  //drawer.setLineProfile(rsLineDrawerFFF::PROFILE_CUBIC);
  drawer.setLineWidth(thickness);
  //drawer.setColor(brightness);
  drawer.setRoundCaps(true);

  float margin = 2*thickness;
  float xMin, yMin, xMax, yMax;
  float br;
  xMin = yMin = margin;
  xMax = imageWidth  - margin;
  yMax = imageHeight - margin;
  //drawer.drawLine(xMin, yMin, xMax, yMax);
  rsRandomUniform(0.0, 1.0, 1);
  float x0, y0, x1, y1;
  x0 = (float) rsRandomUniform(xMin, xMax);  
  y0 = (float) rsRandomUniform(yMin, yMax);
  drawer.initPolyLine(x0, y0);
  for(int i = 2; i <= numLines; i++)
  {
    x1 = (float) rsRandomUniform(xMin, xMax);
    y1 = (float) rsRandomUniform(yMin, yMax);
    br = (float) rsRandomUniform(minBrightness, maxBrightness);
    //drawer.drawLine(x0, y0, x1, y1);
    drawer.setColor(br);
    //drawer.lineTo(x1, y1);
    drawer.lineTo(x1, y1, true); // true: line joining code for uniform color polylines
    x0 = x1;
    y0 = y1;
  }

  writeImageToFilePPM(image, "PolyLineRandom.ppm");
}

void phaseScopeLissajous()
{
  // We create a PhaseScope image of a Lissajous figure to test the drawing code.
  // x(t) = sin(2*pi*a*t), y(t) = sin(2*pi*b*t)

  // input signal parameters:
  static const int N = 11;  // number of data points per cycle
  //static const int N = 80;
  //static const int N = 35;
  //static const int N = 23;
  float a = 2.f;
  float b = 3.f;
  float scale = 0.9f;
  //static const int numCycles = 1;

  typedef RAPT::rsRealTimeSpline<double, float> SG; // for spline-generator


  // create and set up rsPhaseScopeBuffer object:
  rsPhaseScopeBufferFFD psb;
  psb.setSampleRate(N);
  psb.setAntiAlias(true);
  //psb.setBrightness(300.f);  // wtf? - why did this work formerly?
  psb.setBrightness(0.1f); 
  psb.setLineDensity(1.f);
  psb.setPixelSpread(0.3f);
  //psb.setDrawMode(SG::LINEAR);
  //psb.setDrawMode(SG::CUBIC_HERMITE);
  psb.setDrawMode(SG::QUADRATIC);
  psb.setUseColorGradient(true);
  psb.setDensityNormalization(true);
  psb.setSize(400, 400);

  //// settings for testing color discontinuities in spline drawing (remove when problems are fixed):
  //psb.setAntiAlias(false);
  //psb.setUseColorGradient(false);
  //psb.setLineDensity(0.3f);
  //psb.setPixelSpread(0.0f);
  //psb.setSize(800, 800);
  //psb.setOneDimensionalMode(true); // for debug

  // create image:
  psb.reset();
  float x[N], y[N];
  float s = float(2*PI) / N;
  for(int n = 0; n < N; n++) {
    x[n] = scale*sin(s*a*n); y[n] = scale*sin(s*b*n); psb.processSampleFrame(x[n], y[n]);
  }
  // first run through the figure was only for warm-up (to avoid initial artifacts):
  psb.clearImage();
  for(int n = 0; n < N; n++) {
    x[n] = scale*sin(s*a*n); y[n] = scale*sin(s*b*n); psb.processSampleFrame(x[n], y[n]);
  } // factor out this loop into a function drawLissajous(...)

  // retrieve and save image:
  psb.getImage();
  writeImageToFilePPM(*psb.getImage(), "PhaseScopeLissajous.ppm");

  //// plot (for reference):
  //GNUPlotter plt;
  //plt.addDataArrays(N, x, y);
  //plt.plot();
}

void splineArc()
{
  // Tests the spline arc drawing - especially the numeric integration of the arc-length

  //int numDots = 100;

  float density = 0.125;
  int N = 100;        // number of spline evaluation samples
  int width   = 400;
  int height  = 400;

  //float x1 = 0, y1 = 0, x1s = 0, y1s =  1; // start at (0,0), pointing upward
  //float x2 = 1, y2 = 0, x2s = 0, y2s = -1; // end   at (1,0), pointing downward
  float x1 =  10, y1 = 10, x1s = 0, y1s =  1500; // start center left, pointing upward
  float x2 = 390, y2 = 10, x2s = 0, y2s = -1500; // end center right, pointing downward
  float distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); // distance between the two points

  //int numDots = density*distance

  // maybe to test the density, have all points on one lines

  // compute polynomial coefficients:
  float a[4], b[4];
  cubicSplineArcCoeffs2D(x1, x1s, y1, y1s, x2, x2s, y2, y2s, a, b);

  typedef rsPolynomial<float> PL;

  // create points via simple algorithm (without density compensation):
  std::vector<float> t(N), x(N), y(N), s(N);
  float dx, dy, scaler = (1.f / (N-1));
  float length = 0.f;  // accumulated arc-length estimate
  t[0] = 0;
  x[0] = PL::evaluate(0, a, 3); // use optimized evaluateCubic (maybe inlined)
  y[0] = PL::evaluate(0, b, 3);
  for(int n = 1; n < N; n++) {
    t[n] = scaler * n;  // == n / N
    x[n] = PL::evaluate(t[n], a, 3);
    y[n] = PL::evaluate(t[n], b, 3);
    dx   = x[n]-x[n-1];
    dy   = y[n]-y[n-1];
    length += sqrt(dx*dx+dy*dy);
  }

  cubicSplineArcLength2D(&a[0], &b[0], &t[0], &s[0], N);
  // actually, we could try to use a different (less dense) t-array here, the s-array must then 
  // also be shorter - this number should be a second parameter

  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0], &y[0]); // the actual spline curve
  plt.addDataArrays(N, &t[0], &s[0]);   // arc-length s as function of parameter t
  //plt.plot();


  float splineLength = s[N-1]; // last value in s is total length: s(t=1)
  int numSplineDots = std::max(1, (int)round(splineLength * density));


  rsImageF image(width, height);
  rsImagePainter<float, float, float> painter(&image);
  //painter.setAntiAlias(false);
  //painter.drawDottedSpline1(a, b, 1.0, 1.0, numSplineDots); // no density compensation
  //painter.drawDottedSpline2(a, b, 1.0, 1.0); // removed
   

  writeImageToFilePPM(*painter.getImage(), "CubicSpline.ppm");
  int dummy = 0;
}

void triangles()
{
  // create and set up objects and parameters:
  typedef rsVector2DF Vec2;    // for convenience
  typedef rsVector2DF V;       // even shorter - for inlined constructor calls
  float c = 1.0f;              // color (gray value)
  rsImageF img(35, 20);        // image to draw on
  rsImageDrawerFFF drw(&img);  // drawer object
  drw.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);

  void (*pDrawTriangle)(rsImageDrawerFFF&, const Vec2&, const Vec2&, const Vec2&, float);
  pDrawTriangle = &drawTriangle;
  //pDrawTriangle = &drawTriangleAntiAliasedProto;

  //void drawTriangle(rsImageDrawerFFF& drw, 
  //  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color);

  // parallelogram made from a flat-top and flat-bottom triangle:
  Vec2 
    p1(  2.f,  1.f),  // V(2,1) 
    p2(  7.f,  1.f),  // V(7,1)
    p3(  4.f,  5.f),  // V(4,5)
    p4(  9.f,  5.f);  // V(9,5)
  pDrawTriangle(drw, p1, p2, p3, c);  // flat-top
  pDrawTriangle(drw, p2, p4, p3, c);  // flat-bottom
  //drawTriangleFlatTop(   drw, p1, p2, p3, c);  // p1,p2 above p3, clockwise
  //drawTriangleFlatBottom(drw, p2, p3, p4, c);  // p2 above p3,p4, counterclockwise

  // a general triangle:
  pDrawTriangle(drw, V(25,1), V(31,10), V(17,5), c);

  // polygon from 3 general triangles:
  pDrawTriangle(drw, V( 6, 7), V(15, 8), V( 9,10), c);
  pDrawTriangle(drw, V(15, 8), V(14,13), V( 9,10), c);
  pDrawTriangle(drw, V( 6, 7), V( 9,10), V( 3,11), c);

  // polygon from 7 triangles:
  pDrawTriangle(drw, V(17,12), V(15,17), V(17,16), c);
  pDrawTriangle(drw, V(17,12), V(17,16), V(20,17), c);
  pDrawTriangle(drw, V(19,11), V(17,12), V(20,17), c);
  pDrawTriangle(drw, V(19,11), V(20,17), V(23,12), c);
  pDrawTriangle(drw, V(23,12), V(20,17), V(23,15), c);
  pDrawTriangle(drw, V(23,12), V(23,15), V(26,14), c);
  pDrawTriangle(drw, V(26,14), V(23,15), V(24,17), c);

  // polygon from 5 triangles:
  pDrawTriangle(drw, V( 1,13), V(3,17), V( 7,13), c);
  pDrawTriangle(drw, V( 7,13), V(3,17), V(11,13), c);
  pDrawTriangle(drw, V(11,13), V(3,17), V( 9,16), c);
  pDrawTriangle(drw, V(11,13), V(9,16), V(11,16), c);
  pDrawTriangle(drw, V( 3,17), V(7,19), V( 9,16), c);

  // half-pixel mini triangles:
  pDrawTriangle(drw, V(27,12), V(27,13), V(28,12), c); // top left
  pDrawTriangle(drw, V(29,12), V(29,13), V(30,13), c); // bottom left
  pDrawTriangle(drw, V(31,12), V(32,13), V(32,12), c); // top right
  pDrawTriangle(drw, V(34,12), V(33,13), V(34,13), c); // bottom right

  // small "arrow-head":
  pDrawTriangle(drw, V(29,15), V(30,17), V(30,16), c);
  pDrawTriangle(drw, V(29,15), V(30,16), V(31,16), c);

  // for copy/paste:
  //pDrawTriangle(drw, V(,), V(,), V(,), c);

  // when assuming a downward y-axis, alway start at the top vertex (in case of a flat-top, at
  // the top-left) and traverse vertices in counterclockwise order

  // save to file:
  //writeImageToFilePPM(img, "PolygonsViaTriangles.ppm");
  writeScaledImageToFilePPM(img, "PolygonsViaTriangles.ppm", 16);
    // todo: write a function that includes a magnification factor (or maybe two, for x and y
    // separately)
}

void pixelCoverage()
{
  typedef rsVector2DF Vec2;    // for convenience
  typedef rsVector2DF V;       // even shorter - for inlined constructor calls

  Vec2 a(1.5f, 1.5f), b(-1.0f, 0.25f), c(2.0, -0.75);

  float x = 0, y = 0;

  float cov = pixelCoverage2(x, y, a, b, c);

  cov = pixelCoverage2(x, y, V(2,1), V(-1,0.25f), V(2.5f,-0.25f));

  int dummy = 0;
}

// idea for contour drawing:
// -input: image the function values, array of levels
// -output: image with the level lines / contours
//


void getContourSubPixelPosition1(float z00, float z01, float z10, float z11, float c,
  float* x, float* y)
{
  // average value along pixel borders:
  float L = (z00 + z01) * 0.5f;  // left
  float R = (z10 + z11) * 0.5f;  // right
  float T = (z00 + z10) * 0.5f;  // top
  float B = (z01 + z11) * 0.5f;  // bottom

  // x-offset is determined by where the level would fall between L and R:
  //*x = ((L+R) * 0.5 - c) / (R-L);  // is that correct? what about div-by-zero
  // sometimes |x| > 1 - that should not happen - actually it should be always in the 
  // range 0...1 - so that formula is not yet correct

  //if(L < R)
  //  *x = rsLinToLin(c, L, R, 0.f, 1.f);
  //else
  //  *x = 1.f - rsLinToLin(c, R, L, 0.f, 1.f);  // is that correct?

  *x = rsLinToLin(c, L, R, 0.f, 1.f);
  *y = rsLinToLin(c, T, B, 0.f, 1.f); // top comes before bottom - y increases downward

  // output may go negative - maybe we have to do cases if(L < R) .. else ...
}
// this doesn't work

void getContourSubPixelPosition2(float z00, float z01, float z10, float z11, float c,
  float min, float max, float* x, float* y)
{
  z00 -= c;
  z01 -= c;
  z10 -= c;
  z11 -= c;

  //float zAvg = 0.25f * (z00 + z01 + z10 + z11); 
    // optimize: delete the -= and do zAvg = 0.25f * (z00 + z01 + z10 + z11 - 4*c);

  float L = (z00 + z01) * 0.5f;  // left
  float R = (z10 + z11) * 0.5f;  // right
  float T = (z00 + z10) * 0.5f;  // top
  float B = (z01 + z11) * 0.5f;  // bottom

  *x = rsLinToLin(c, L, R, 0.f, 1.f);
  *y = rsLinToLin(c, T, B, 0.f, 1.f);
}
// not yet tested

/** Given values at a pixel z00 and its neigbours, this function computes the coeffs of a 
parametric line equation x(t) = x0 + t*(x1-x0), y(t) = y0 + t*(y1-y0). The line is supposed to 
approximate a segment of a contour of a function/image with level given by c. This is how the 
function values are distributed over the pixels: z00 is the value at the pixel under investigation 
and z01, z10, z11 are bottom, right and diagonal neighbours (in that order):
   z00--z10
    |    |
   z01--z11  
 so the first index refers to the x-coordinate, the second to the y-coordinate and y increases 
 downward as is common in image processing. The function also returns an integer indicating which 
 branch we ended up in: 0: top-left, 1: top-right, 2: bottom-left, 3: bottom-right, 
 4: horizontalish, 5: verticalish - where "top-left" etc. means that the segment "cuts off" the 
 top-left corner of the pixel. */
int contourSegmentCoeffs(float z00, float z01, float z10, float z11, float c,
  float& x0, float& y0, float& x1, float& y1)
{
  int branch;
  if((z00 < c && z01 >= c) || (z00 >= c && z01 < c)) {        // segment goes through left border
    x0 = 0.f;
    y0 = rsLinToLin(c, z00, z01, 0.f, 1.f);
    if((z00 < c && z10 >= c) || (z00 >= c && z10 < c)) {      // segment goes through top border
      branch = 0;                                             //   -> top-left
      x1 = rsLinToLin(c, z00, z10, 0.f, 1.f);
      y1 = 0.f; }
    else if((z01 < c && z11 >= c) || (z01 >= c && z11 < c)) { // segment goes through bottom border
      branch = 2;                                             //   -> bottom-left
      x1 = rsLinToLin(c, z01, z11, 0.f, 1.f);
      y1 = 1.f; }
    else {                                                    // segment goes through right border
      branch = 4;                                             //   -> horizontalish
      x1 = 1.f;                                               
      y1 = rsLinToLin(c, z10, z11, 0.f, 1.f); }}
  else {                                                   // doesn't go through left border
    if((z00 < c && z10 >= c) || (z00 >= c && z10 < c)) {   // goes through top border
      x0 = rsLinToLin(c, z00, z10, 0.f, 1.f);
      y0 = 0.f;
      if((z10 < c && z11 >= c) || (z10 >= c && z11 < c)) { // goes through right border
        branch = 1;                                        //   -> top-right
        x1 = 1.f;
        y1 = rsLinToLin(c, z10, z11, 0.f, 1.f); }
      else  {                                              // goes through bottom border
        branch = 5;                                        //   -> verticalish
        x1 = rsLinToLin(c, z01, z11, 0.f, 1.f);            
        y1 = 1.f; }}
    else  {                                                // doesn't go through top border 
      branch = 3;                                          //   -> bottom-right
      x0 = rsLinToLin(c, z01, z11, 0.f, 1.f);
      y0 = 1.f;
      x1 = 1.f;
      y1 = rsLinToLin(c, z10, z11, 0.f, 1.f); }}
  return branch;
}
// optimize the calls to rsLinToLin to get rid of divisions where possible - keep this code as 
// prototype for unit testing the optimized code - i think, it's not possible, but we may get rid 
// of some of the multiplications because outMax-outMin = 1 - make a function rsLinTo01, have a 
// similar rs01ToLin
// note that the order of (x0,y0),(x1,y1) can't be changed without breaking contourPixelCoverage

// maybe the logical statements can be simplified by checking things like 
// if (z00-c)*(z01-c) < 0,  >= 0 instead of the complicated and-or statements - but keep this 
// version for unit tests

template<class T>
T triangleArea(T x1, T y1, T x2, T y2, T x3, T y3)
{
  return T(0.5) * rsAbs(x1*y2 + x3*y1 + x2*y3 - x3*y2 - x1*y3 - x2*y1);
}
// needs test
//         1    |x1 y1 1|
// A = +- --- * |x2 y2 1|
//         2    |x3 y3 1|
// https://www.onlinemathlearning.com/area-triangle.html

// we also need a formula for the area of a quadrangle - i think, i have once implemented a general
// polygon-area function - the quadrangle can be obtained as special case

template<class T>
T contourPixelCoverage(T z00, T z01, T z10, T z11, T c)
{
  T x0, x1, y0, y1;
  T A = 0.f;  // covered area
  T h(0.5);   // half
  T I(1.0);   // one
  int branch = contourSegmentCoeffs(z00, z01, z10, z11, c, x0, y0, x1, y1);
  switch(branch) {
  case 0: { A = h *    x1  * y0;     if(z00 >= c) A = I-A; } break; // top-left
  case 1: { A = h * (I-x0) * y1;     if(z10 >= c) A = I-A; } break; // top-right
  case 2: { A = h *    x1  * (I-y0); if(z01 >= c) A = I-A; } break; // bottom-left
  case 3: { A = h * (I-x0) * (I-y1); if(z11 >= c) A = I-A; } break; // bottom-right
  case 4: {                              // horizontalish
    if(y0 < y1)  A = y1 + h * (y0-y1);   //   going up
    else         A = y0 + h * (y1-y0);   //   going down
    if(z00 >= c || z10 >= c)  
      A = I-A; } break;
  case 5: {                              // verticalish
    if(x0 < x1)  A = x1 + h * (x0-x1);   //   leaning left
    else         A = x0 + h * (x1-x0);   //   leaning right
    if(z00 >= c || z01 >= c) 
      A = I-A; } break; }
  return A;
}
// these simplified formulas work only because we know in which order contourSegmentCoeffs 
// returns the coeffs. maybe we should make it swappable whethr to use >= or < - sometimes we may 
// want to invert the result - when drawing the bin-fills, we sometimes want to fill with the 
// inverted weight ..i think - figure out - if so, maybe use a boolean and or let the user pass a 
// comparison function cmp(z00, c), etc... or call it like inside(z00, c) or outside(z00, c)



void contourSubPixelPosition(float z00, float z01, float z10, float z11, float c,
  float* x, float* y, float* weight)
{
  float x0, x1, y0, y1;
  contourSegmentCoeffs(z00, z01, z10, z11, c, x0, y0, x1, y1);

  // We have our line equation - evaluate line equation at midpoint to get the center of the 
  // contour segment. The weight is given by the length divided by sqrt(2) such that diagonals get 
  // weight 1.0
  float dx = (x1-x0);
  float dy = (y1-y0);
  *x = x0 + 0.5f * dx;
  *y = y0 + 0.5f * dy;

  //*weight = max(dx, dy);  // nope - this looks worse - screw-effect stronger and there are holes  
  *weight = sqrt(dx*dx + dy*dy) / sqrt(2.f);  // full weight only for diagonals
    // optimize, maybe use max(dx, dy)
}







// simpler idea:
// -compute z0 = (z00 + z01) / 2, z1 = (z10 + z11) / 2
//  z0 is the average value on the left, z1 on the right
// -the x-value/offset is determined by how much this average is above/below the target
//  level
// -similar for y
// -or is it the other way around?
// -might be even better than the center of the line 

template<class TPix, class TLvl>
void drawContour(const rsImage<TLvl>& z, TLvl level, rsImage<TPix>& target, TPix color, 
  bool antiAlias)
{
  rsImagePainter<TPix, TLvl, TLvl> painter(&target);
  painter.setDeTwist(true);  
  for(int i = 0; i < z.getWidth()-1; i++) {
    for(int j = 0; j < z.getWidth()-1; j++) {
      TLvl z00 = z(i,   j  );
      TLvl z01 = z(i,   j+1);
      TLvl z10 = z(i+1, j  );
      TLvl z11 = z(i+1, j+1);
      TLvl min = rsMin(z00, z01, z10, z11);
      TLvl max = rsMax(z00, z01, z10, z11);
      if(min < level && max >= level) {
        TLvl x(0), y(0), w(1); // x,y offsets and weight for color
        if(antiAlias)
          contourSubPixelPosition(z00, z01, z10, z11, level, &x, &y, &w);
        painter.paintDot(TLvl(i) + x, TLvl(j) + y, w * color); }}}
}
// if we do not anti-alias, we need not to call the expensive paintDot and can use the cheaper 
// painter.plot instead ...i think

// maybe don't loop over all pixels and follow the contours instead - but then there's no guarantee that 
// nothing is missed


template<class TPix, class TLvl>
void fillBetweenContours(const rsImage<TLvl>& z, TLvl lo, TLvl hi, rsImage<TPix>& target,
  TPix fillColor, bool antiAlias = false)
{
  rsImagePainter<TPix, TLvl, TLvl> painter(&target);
  for(int i = 0; i < z.getWidth()-1; i++) {
    for(int j = 0; j < z.getWidth()-1; j++) {
      TLvl z00 = z(i, j);
      TLvl z01 = z(i, j+1);
      TLvl z10 = z(i+1, j);
      TLvl z11 = z(i+1, j+1);
      TLvl min = rsMin(z00, z01, z10, z11);
      TLvl max = rsMax(z00, z01, z10, z11);
      if(min >= lo && max < hi)     // this seems to be artifact-free
        painter.plot(i, j, fillColor);
      else  {
        // we are either on the contour or totally outside the drawing area
        if(!antiAlias)  {                  // ok - this looks right
          if(min < lo && max >= lo)                      // on low contour
            painter.plot(i, j, fillColor * TPix(0.5));
          else if(min < hi && max >= hi)                 // on hi contour
            painter.plot(i, j, fillColor * TPix(0.5)); }
        else {
          TLvl c; // coverage
          if(min < lo && max >= lo) {
            c = contourPixelCoverage(z00, z01, z10, z11, lo);
            painter.plot(i, j, TPix(TLvl(1)-c)*fillColor); }
          else if(min < hi && max >= hi) {
            c = contourPixelCoverage(z00, z01, z10, z11, hi);
            painter.plot(i, j, TPix(c)*fillColor);  }}}}}
}
// if instead of using:
//   if(min >= lo && max < hi) 
// we would use
//   if(min > lo && max < hi)   -> leaves extra pixels blank (test with circle)
//   if(min >= lo && max <= hi) -> colors extra pixels in
//   if(min > lo && max <= hi)  -> no etra blank or colored pixels but ugly jaggies
// so the chosen variant seems best. this can be tested using the circles (and maybe commenting)
// out the code that handles the contour lines - i think it was set to somewhere around 11 or 12 
// levels...not sure anymore

template<class TPix, class TWgt>
TPix blend(TPix c1, TPix c2, TWgt w)
{
  return TPix((TWgt(1)-w))*c1 + TPix(w)*c2;
}

template<class TLvl, class TPix>
rsImage<TPix> getContours(const rsImage<TPix>& z, const std::vector<TLvl>& levels, 
  const std::vector<TPix>& colors, bool antiAlias,
  const std::vector<TPix>& fillColors = std::vector<TPix>() )
{
  rsImageF c(z.getWidth(), z.getHeight());
  for(size_t i = 0; i < levels.size(); i++)
    drawContour(z, levels[i], c, colors[i % colors.size()], antiAlias);
  return c;
}
// get rid of fillColors - it's not used (i think - if it is, it shouldn't be)

template<class TLvl, class TPix>
rsImage<TPix> getBinFills(
  const rsImage<TPix>& z, 
  const std::vector<TLvl>& levels,
  const std::vector<TPix>& colors, 
  bool antiAlias)
{
  rsImageF imgBins(z.getWidth(), z.getHeight());  // fills
  size_t j = 0; // color index
  size_t nc = colors.size();
  for(size_t i = 0; i < levels.size()-1; i++) {
    fillBetweenContours(z, levels[i], levels[i+1], imgBins, colors[j % nc], antiAlias);
    j++; }
  return imgBins;
}


template<class T>
void normalizeFast(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  int N = img.getNumPixels();
  T min = rsArrayTools::minValue(p, N);
  T max = rsArrayTools::maxValue(p, N);
  T scl = 1.f / (max-min);
  for(int i = 0; i < N; i++)
    p[i] = scl * (p[i] - min);
}
// maybe make member

// this *may* be better numerically (less prone to roundoff errors) - needs test
template<class T>
void normalize(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  int N = img.getNumPixels();
  T min = rsArrayTools::minValue(p, N);
  for(int i = 0; i < N; i++)
    p[i] -= min;
  T max = rsArrayTools::maxValue(p, N);
  for(int i = 0; i < N; i++)
    p[i] /= max;
}

/** Joint normalization of two images */
template<class T>
void normalizeJointly(rsImage<T>& img1, rsImage<T>& img2)
{
  //rsAssert(img2.hasSameShapeAs(img1));  // activate later
  using AT = rsArrayTools;
  int N = img1.getNumPixels();
  T* p1 = img1.getPixelPointer(0, 0);
  T* p2 = img2.getPixelPointer(0, 0);
  T min = rsMin(AT::minValue(p1, N), AT::minValue(p2, N));
  for(int i = 0; i < N; i++) {
    p1[i] -= min;
    p2[i] -= min; }
  T max = rsMax(AT::maxValue(p1, N), AT::maxValue(p2, N));
  for(int i = 0; i < N; i++) {
    p1[i] /= max;
    p2[i] /= max;
  }
}
// maybe have a version for three images as well - can this be generalized with variadic 
// templates?



bool testContourSubPixelStuff()
{
  bool r = true;

  // test - turn into unit-test
  //int b; // branch
  float x, y, w;
  contourSubPixelPosition(2.f, 8.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.25f  && y == 0.25f;
  contourSubPixelPosition(8.f, 2.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.25f  && y == 0.75f;
  contourSubPixelPosition(8.f, 8.f, 2.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.75f  && y == 0.25f;
  contourSubPixelPosition(8.f, 8.f, 8.f, 2.f, 5.f, &x, &y, &w); r &= x == 0.75f  && y == 0.75f;
  contourSubPixelPosition(2.f, 2.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.5f   && y == 0.5f;
  contourSubPixelPosition(2.f, 8.f, 2.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.5f   && y == 0.5f;
  contourSubPixelPosition(2.f, 4.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.375f && y == 0.5f;
  contourSubPixelPosition(2.f, 8.f, 4.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.5f   && y == 0.375f;

  // todo: test coverage compuation
  float c;
  c = contourPixelCoverage(2.f, 8.f, 8.f, 8.f, 5.f); r &= c == 0.125;
  c = contourPixelCoverage(8.f, 2.f, 8.f, 8.f, 5.f); r &= c == 0.125;
  c = contourPixelCoverage(8.f, 8.f, 2.f, 8.f, 5.f); r &= c == 0.125;
  c = contourPixelCoverage(8.f, 8.f, 8.f, 2.f, 5.f); r &= c == 0.125;

  // horizontalish lines:
  //c = contourPixelCoverage(2.f, 8.f, 2.f, 8.f, 5.f);
  //c = contourPixelCoverage(2.f, 2.f, 8.f, 8.f, 5.f);
  c = contourPixelCoverage(2.f, 8.f, 4.f, 8.f, 5.f); r &= c == 0.375;
  c = contourPixelCoverage(4.f, 8.f, 2.f, 8.f, 5.f); r &= c == 0.375;
  c = contourPixelCoverage(8.f, 2.f, 8.f, 4.f, 5.f); r &= c == 0.375;
  c = contourPixelCoverage(8.f, 4.f, 8.f, 2.f, 5.f); r &= c == 0.375;

  // verticalish lines:
  c = contourPixelCoverage(2.f, 4.f, 8.f, 8.f, 5.f); r &= c == 0.375;
  c = contourPixelCoverage(4.f, 2.f, 8.f, 8.f, 5.f); r &= c == 0.375;
  c = contourPixelCoverage(8.f, 8.f, 2.f, 4.f, 5.f); r &= c == 0.375;
  c = contourPixelCoverage(8.f, 8.f, 4.f, 2.f, 5.f); r &= c == 0.375;
  // wait! shouldn't the latter two give 0.625? hmm..no - maybe try cases where the inversion
  // kicks in - oh - it actually does

  return r;
}
// move to unit tests

/** Generates an image from a function by looping through the pixels, computing the pixel 
coordinates by transforming pixel index pairs (i,j) to the ranges xMin..xMax, yMin..yMax and 
evaluating the function f at the resulting (x,y) location. TPix is a separate type for the pixels, 
so you may use functions operating on double and having float for the pixels, for example. */
template<class T, class TPix> 
void generateFunctionImage(const function<T(T, T)>& f, T xMin, T xMax, T yMin, T yMax,
  rsImage<TPix>& img)
{
  for(int i = 0; i < img.getWidth(); i++) {
    for(int j = 0; j < img.getHeight(); j++) {
      T x = xMin + i * (xMax-xMin) / (img.getWidth()  - 1);
      T y = yMin + j * (yMax-yMin) / (img.getHeight() - 1);
      T z = f(x, y);
      img.setPixelColor(i, j, TPix(z)); }}
}
// todo:
// -maybe take two optional coordinate transformation functions c1(x,y), c2(x,y) where by default,
//  c1(x,y) = x and c2(x,y) = y - for polar coordinates, we would use c1(x,y) = sqrt(x*x + y*y), 
//  c2(x,y) = atan2(y,x)

template<class T, class TPix> 
void generateFunctionImageReIm(const function<complex<T>(complex<T>)>& f, 
  T xMin, T xMax, T yMin, T yMax, rsImage<TPix>& imgRe, rsImage<TPix>& imgIm)
{
  function<T(T, T)> fr, fi;
  fr = [=](T x, T y) { return f(complex<T>(x, y)).real(); };
  fi = [=](T x, T y) { return f(complex<T>(x, y)).imag(); };
  generateFunctionImage(fr, xMin, xMax, yMin, yMax, imgRe);
  generateFunctionImage(fi, xMin, xMax, yMin, yMax, imgIm);
}

void contours()
{
  // We plot the 2D function z = f(x,y) = x^2 - y^2 into an image where the height translates
  // to the pixel brightness

  rsAssert(testContourSubPixelStuff());

  int w = 129;               // width in pixels
  int h = 129;               // height in pixels

  w = h = 513;
  //w = h = 1025;

  float r = 18;
  int numLevels = 20;
  int numColors = numLevels + 1;


  float xMin = -r;
  float xMax = +r;
  float yMin = -r;
  float yMax = +r;


  std::function<float(float, float)> f;

  // choose your function here:
  //f = [&] (float x, float y) { return x*x + y*y; };
  //f = [&] (float x, float y) { return x*x - y*y; };
  //f = [&] (float x, float y) { return x*x - y*y + 2.f*x*y; };
  //f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y; };
  //f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y + 0.1f*x*x - 0.1f*y*y; };
  //f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y + 0.1f*x*x - 0.1f*x - 0.1f*y*y + 0.1f*y; };
    // try exchanging sin and cos an combining
  f = [&] (float x, float y) { return (float) spiralRidge(x, y); };


  // create image with function values:
  rsImageF imgFunc(w, h);
  generateFunctionImage(f, xMin, xMax, yMin, yMax, imgFunc);;
  normalize(imgFunc);

  // create images with contours:
  std::vector<float> levels = rsRangeLinear(0.f, 1.f, numLevels);
  rsImageF imgCont = getContours(imgFunc, levels, { 1.0f }, true);
  // with anti-aliasing, we need to use about twice as much brightness to get the same visual 
  // brightness

  // create images with bin-fills:
  std::vector<float> colors = rsRangeLinear(0.f, 1.f, numColors);
  rsImageF imgFills = getBinFills(imgFunc, levels, colors, true);

  // write images to files:
  writeScaledImageToFilePPM(imgFunc,  "Function.ppm", 1);
  //writeScaledImageToFilePPM(imgCont,  "Contours.ppm", 1);
  //writeScaledImageToFilePPM(imgFills, "BinFills.ppm", 1);

  // the right column and bottom row has no countour values - no surprise - the loop only goes up 
  // to w-1,h-1
  // maybe use powers of two +1 for the size and cut off bottom-row and right-column aftewards

  // todo: 
  // -make a composited image with function values and contours - maybe have a compose function 
  //  that takes a combiner-function for pixel values
  // -it could make sense to superimpose the contour lines on the original function - but how - 
  //  just add? that might be bad in bright areas
  // -maybe use the color-channels to plot more than one function 
  //  -plot complex functions - real -> red, imag -> green or blue
  // -could this be used as a drawing primitive? it can draw circles, for example - but it's 
  //  very inefficient - can this be optimized into a reasonable implicit curve drawing algo?

  // maybe try to overlay images with multiple settings for the number/positions of the level-lines
  // -could create intersting patterns
}

void complexContours()
{
  // -make a function w = f(z) = z^2
  // -plot real part into red channel and imaginary part into blue channel
  //  -maybe plot the contour lines into green channel

  // setup:
  int w           = 513;               // width in pixels
  int h           = 513;               // height in pixels
  int numLevels   = 20;                // number of contour levels
  int numColors   = numLevels + 1;     // number of contour colors
  double r        = 5;                 // range x = -r..+r, y = -r...+r
  bool antiAlias  = true;              // switch anti-alias on/off
  float levelPow  = 1.f;               // power/exponent for nonlinear spacing of contour levels
                                       // ...not useful - we need to squish toward the center

  double xMin = -r;
  double xMax = +r;
  double yMin = -r;
  double yMax = +r;

  //numColors = 10;

  using Complex = complex<double>;
  function<Complex(Complex)> f;

  // pick complex function to plot
  //f = [=](Complex z) { return z; };
  //f = [=](Complex z) { return z + 0.5*z*z + (1./6)*z*z*z; };
  f = [=](Complex z) { return z*z; };
  //f = [=](Complex z) { return z*z*z; };
  //f = [=](Complex z) { return z*z*z*z; };
  //f = [=](Complex z) { return 1./(1. + z); };


  // render images of function values for real and imaginary part:
  rsImageF funcRe(w, h), funcIm(w, h), empty(w, h);
  generateFunctionImageReIm(f, xMin, xMax, yMin, yMax, funcRe, funcIm);
  normalizeJointly(funcRe, funcIm); // joint normalization preserves re,im relationship (right?)
  //normalize(funcRe);
  //normalize(funcIm);
  // what about absolute value and phase? ca we do something useful with them, too?

  // get contour lines:
  std::vector<float> levels = rsRangeLinear(0.f, 1.f, numLevels);
  for(int i = 0; i < numLevels; i++)
    levels[i] = pow(levels[i], levelPow);
  rsImageF contRe = getContours(funcRe, levels, { 1.0f }, antiAlias);
  rsImageF contIm = getContours(funcIm, levels, { 1.0f }, antiAlias);
  // rename to getContourLines

  // get countour fills:
  std::vector<float> colors = rsRangeLinear(0.f, 1.f, numColors);
  rsImageF fillsRe = getBinFills(funcRe, levels, colors, antiAlias);
  rsImageF fillsIm = getBinFills(funcIm, levels, colors, antiAlias);
  // maybe rename to getContourFills, getCountourSteps

  // -maybe try to mix the contourFills with the orignal function - should give some indication of 
  //  lines but still smooth'ish progression
  // -combine contour lines of re and im
  // -overlay contour lines with func
  // -try a nonlinear spacing of the contour levels - when the function is z^3, i think, the 
  //  contour levels should also increase with x^3 (or x^(1/3) ?) - this effect can be realized by
  //  raising the elements of the levels array to a power - this makes the contour levels roughly
  //  equally spaced - for other functions, we may want an exponential relationship rather than
  //  a power rule - can we somehow estimate a good spacing from teh image data? maybe by fitting
  //  some sort of curve to the statistics of the data - if there are many low values, we want 
  //  the line-spacing to be denser in the low range - we wnat that the different bins occupy 
  //  roughly equal areas - maybe this should be done separately for real and imaginary part
  //  -maybe take the cumulative distribution of height-values and select the levels according to
  //   that
  //  -oh - for functions with odd exponents, we actually want to squish the level line density
  //   toward the center --or actually that applies to even exponents, too - just using an exponent
  //   is not good enough!
  // -maybe we need even more template types: double for the function, float for the heights and
  //  pixel colors may be some RGBA type (rsFloat32x4)
  // -make a class rsContourPlotter in rosic
  //   -have rendering options with presets draft (fastest), presentation (highest quality) and 
  //    maybe intermediate settings

  writeImageToFilePPM(funcRe, "FunctionRe.ppm");
  writeImageToFilePPM(funcIm, "FunctionIm.ppm");

  writeImageToFilePPM(contRe, "ContourLinesRe.ppm");
  writeImageToFilePPM(contIm, "ContourLinesIm.ppm");

  writeImageToFilePPM(fillsRe, "ContourFillsRe.ppm");
  writeImageToFilePPM(fillsIm, "ContourFillsIm.ppm");

  //writeImageToFilePPM(funcRe, funcIm, empty, "ComplexFunctionRG.ppm");
  writeImageToFilePPM(funcRe,  empty, funcIm,  "FunctionRB.ppm");
  writeImageToFilePPM(fillsRe, empty, fillsIm, "ContourFillsRB.ppm");
  writeImageToFilePPM(fillsRe, fillsIm, empty, "ContourFillsRG.ppm");
}


// implicit curve drawing algo: 
// input: f(x,y) = c, f as functor, c as value, one solution x0,y0 that solves the equation
//
// (1) draw/paint pixel at x0,y0
// (2) find neighbouring pixel to draw next among the 8 neighbours
// (3) fix one coordinate and solve for the other by 1D root-finding
// (3) check, if we are back at the same pixel where we started or reached the image boundary 
//     -if yes: return - we are done with the curve
// (4) back to 1
//
// -maybe fill the shape by coloring all pixels for which f(x,y) < c
// -which should we increment or decrement and for which should we solve? x0 or y0? maybe that 
//  should also be decided based on a condition
//  -maybe y should be use a fixed increment - because we may later extend it to taking two 
//   solutions (x0,y0),(x1,y2) and advancing y0,y1 in a loop and compute x0,x1 to get a span which
//   can then be filled - but the boundaries of the span should be only partially colored for 
//   anti-aliasing - that should give a reasonable anti-aliased ellipse-drawing algo
// -maybe make these algos available for use in python with the image processing library

/** Draws the curve defined by f(x,y) = c onto the image. It needs one solution x0,y0 for which
f(x0,y0) = c holds as starting point. */
template<class T, class TPix> 
void drawImplicitCurve(const function<T(T, T)>& f, T xMin, T xMax, T yMin, T yMax, T c, T x0, T y0,
  rsImage<TPix>& img, TPix color)
{
  rsAssert(f(x0, y0) == c, "x0,y0 should solve f(x0,y0) = c" );  // todo: use tolerance

  // maybe pass the painter object - this painte then also already should have the image assigned
  // so we don't need to pass it as additional parameter - this is similar to juce's Graphics 
  // object - we would need to have inquiry functions like getMaxPixelCoordinateX/Y
  rsImagePainter<TPix, T, T> painter(&img);
  painter.setDeTwist(true);

  // figure out start pixel:
  T xMaxPixel = T(img.getWidth()  - 1);   // maximum x-coordinate in pixel coordinates
  T yMaxPixel = T(img.getHeight() - 1);   // same for y-coordinate
  T x  = x0;
  T y  = y0;

  T sx = xMaxPixel / (xMax-xMin);   // one x-pixel in world coordinates
  T sy = yMaxPixel / (yMax-yMin);

  T sxi = (xMax-xMin) / xMaxPixel;
  T syi = (yMax-yMin) / yMaxPixel;

  int iterations = 0;
  while(true)
  {
    // Convert (x,y) to pixel coordinates and draw:
    T px = rsLinToLin(x, xMin, xMax, T(0), xMaxPixel);
    T py = rsLinToLin(y, yMin, yMax, T(0), yMaxPixel);
    painter.paintDot(px, py, color);

    // Figure out gradient (dx,dy) and contour direction (rx,ry) which is perpendicular to the 
    // gradient:
    T h  = 1.e-8;  // ad-hoc - make parameter
    T dx = (f(x+h, y) - f(x-h, y)) / (T(2)*h);  // x-component of gradient
    T dy = (f(x, y+h) - f(x, y-h)) / (T(2)*h);  // y-component of gradient
    T rx = -dy;  // (rx,ry) = (-dy,dx) - gradient, rotated by 90° counterclockwise..
    T ry =  dx;  // ..this is a direction along the contour (approximately)

    // Check, if the current segment is horizontalish/flat or verticalish/steep. In the flat case, 
    // advance x by one pixel and y by a distance derived from the direction vector - in the steep 
    // case, the other way around:
    bool flat = rsAbs(rx*sx) > rsAbs(ry*sy);   // curve sgement is horizontalish
    if(flat) {
      dx = rsSign(rx) / sx;    // this step should translate to 1 pixel left or right -> check this!
      dy = dx * ry/rx;         // the y-step is proportional to the x-step - is thsi formula the best we can do?
      x += dx;                 // walk one pixel left or right
      y += dy; }
    else {              // verticalish
      dy = rsSign(ry) / sy;
      dx = dy * rx/ry;
      x += dx;
      y += dy; }

    // In the step just taken, we may have drifted off the contour line due to approximation 
    // errors, so we fix this by refining x or y such that we land on the contour again. We use 1D 
    // Netwon iteration with numeric derivatives
    // if our direction is horizontalish (i.e.
    // rx*sclX > ry*sclY), we change y, otherwise, we change x - we do this by 1D Newton iteration
    // using numeric derivatives (is this a good idea? what about convergence problems?)
    // this does not yet work well - it seems to sometimes grind to a halt
    T err = f(x,y) - c;
    T tol = 1.e-12;
    if(!flat) {               // y-step is larger (steep) -> refine x
      while(rsAbs(err) > tol)  {
        dx    = (f(x+h, y) - f(x-h, y)) / (T(2)*h);
        x     = x - err / dx;
        T old = err;
        err   = f(x,y) - c;
        if(rsAbs(old) <= rsAbs(err)) {
          x += old / dx;  // restore old value bcs old error was better
          break; }}}
    else {
      while(rsAbs(err) > tol)  {
        dy = (f(x, y+h) - f(x, y-h)) / (T(2)*h);
        y   = y - err / dy;
        T old = err;
        err = f(x,y) - c;  
        if(rsAbs(old) <= rsAbs(err)) {
          y += old / dy;
          break; }}}

    if( rsAbs(x-x0) < sxi && rsAbs(y-y0) < syi ) // maybe && iterations >= 2 so we don't spuriously
      break;                                     // break in the very first iteration?
    // there's a gap sometimes - the last pixel is not drawn -unit circle with -2..+2 and 129x129
    // shows this


    iterations++;
    if(iterations > 7000)  // preliminary
      break;  // use condition later
    // possible stopping criteria: we are close to the starting point x0,y0 (within one pixel
    // distance?) or outside the image boundaries - maybe we should also have a maximum number of
    // iterations
  }

  //int i = (int) round(px);
  //int j = (int) round(py);

  int dummy = 0;
}

void implicitCurve()
{
  double width  = 300;
  double height = 300;

  double xMin   = -2.0;
  double xMax   = +2.0;
  double yMin   = -2.0;
  double yMax   = +2.0;


  function<double(double, double)> f;
  f = [=](double x, double y) { return x*x + y*y; };  // unit circle

  //f = [=](double x, double y) { return x*x + 1.5*y*y; }; 
  // we need one starting point - maybe the function should figure it out itself

  rsImageF imgCurve(width, height);
  drawImplicitCurve(f, xMin, xMax, yMin, yMax, 1.0, 1.0, 0.0, imgCurve, 1.f);


  writeScaledImageToFilePPM(imgCurve, "ImplicitCurve.ppm", 1);
}



// maybe make animations with
// http://www.softpedia.com/get/Multimedia/Graphic/Graphic-Others/APNG-Anime-Maker.shtml
