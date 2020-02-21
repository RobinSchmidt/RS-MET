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


/*
// these were other ideas for implementing rsImageContourPlotter::contourSubPixelPosition - they 
// have been discarded - but maybe it's worth to keep the code, so we may later continue to 
// experiment with these ideas:

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
*/


template<class T>
T triangleArea(T x1, T y1, T x2, T y2, T x3, T y3)
{
  return T(0.5) * rsAbs(x1*y2 + x3*y1 + x2*y3 - x3*y2 - x1*y3 - x2*y1);
}
// needs test
// Formulas for triangle areas:
// https://www.onlinemathlearning.com/area-triangle.html
// -base b and height h: A = b*h/2
// -three sides a,b,c: s = (a+b+c)/2; A = sqrt(s*(s-a)*(s-b)*(s-c))
// -two sides a,b and enclosed angle p: A = a*b*sin(p) / 2
// -equilateral with side-length s: A = s^2 * sqrt(3) / 4
// -1 vertex, 2 vectors v,w: A = norm(cross(v,w)) / 2
// -3 vertices:
//            1    |x1 y1 1|
//    A = +- --- * |x2 y2 1|
//            2    |x3 y3 1|

// we also need a formula for the area of a quadrangle - i think, i have once implemented a general
// polygon-area function - the quadrangle can be obtained as special case

template<class TPix, class TWgt>
TPix blend(TPix c1, TPix c2, TWgt w)
{
  return TPix((TWgt(1)-w))*c1 + TPix(w)*c2;
}

bool testContourSubPixelStuff()
{
  bool r = true;


  rsImageContourPlotter<float, float> cp;

  // test - turn into unit-test
  //int b; // branch
  float x, y, w;
  cp.contourSubPixelPosition(2.f, 8.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.25f  && y == 0.25f;
  cp.contourSubPixelPosition(8.f, 2.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.25f  && y == 0.75f;
  cp.contourSubPixelPosition(8.f, 8.f, 2.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.75f  && y == 0.25f;
  cp.contourSubPixelPosition(8.f, 8.f, 8.f, 2.f, 5.f, &x, &y, &w); r &= x == 0.75f  && y == 0.75f;
  cp.contourSubPixelPosition(2.f, 2.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.5f   && y == 0.5f;
  cp.contourSubPixelPosition(2.f, 8.f, 2.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.5f   && y == 0.5f;
  cp.contourSubPixelPosition(2.f, 4.f, 8.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.375f && y == 0.5f;
  cp.contourSubPixelPosition(2.f, 8.f, 4.f, 8.f, 5.f, &x, &y, &w); r &= x == 0.5f   && y == 0.375f;

  // todo: test coverage compuation
  float c;
  c = cp.contourPixelCoverage(2.f, 8.f, 8.f, 8.f, 5.f); r &= c == 0.125;
  c = cp.contourPixelCoverage(8.f, 2.f, 8.f, 8.f, 5.f); r &= c == 0.125;
  c = cp.contourPixelCoverage(8.f, 8.f, 2.f, 8.f, 5.f); r &= c == 0.125;
  c = cp.contourPixelCoverage(8.f, 8.f, 8.f, 2.f, 5.f); r &= c == 0.125;

  // horizontalish lines:
  //c = contourPixelCoverage(2.f, 8.f, 2.f, 8.f, 5.f);
  //c = contourPixelCoverage(2.f, 2.f, 8.f, 8.f, 5.f);
  c = cp.contourPixelCoverage(2.f, 8.f, 4.f, 8.f, 5.f); r &= c == 0.375;
  c = cp.contourPixelCoverage(4.f, 8.f, 2.f, 8.f, 5.f); r &= c == 0.375;
  c = cp.contourPixelCoverage(8.f, 2.f, 8.f, 4.f, 5.f); r &= c == 0.375;
  c = cp.contourPixelCoverage(8.f, 4.f, 8.f, 2.f, 5.f); r &= c == 0.375;

  // verticalish lines:
  c = cp.contourPixelCoverage(2.f, 4.f, 8.f, 8.f, 5.f); r &= c == 0.375;
  c = cp.contourPixelCoverage(4.f, 2.f, 8.f, 8.f, 5.f); r &= c == 0.375;
  c = cp.contourPixelCoverage(8.f, 8.f, 2.f, 4.f, 5.f); r &= c == 0.375;
  c = cp.contourPixelCoverage(8.f, 8.f, 4.f, 2.f, 5.f); r &= c == 0.375;
  // wait! shouldn't the latter two give 0.625? hmm..no - maybe try cases where the inversion
  // kicks in - oh - it actually does

  rsAssert(r);
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
      //T y = yMin + j * (yMax-yMin) / (img.getHeight() - 1);  // wrong! we need to flip vertically
      T y = yMax - j * (yMax-yMin) / (img.getHeight() - 1);
      T z = f(x, y);
      img.setPixelColor(i, j, TPix(z)); }}
}
// move to rsImageGenerator
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

  //w = h = 100;
  //w = h = 513;
  w = h = 500;
  //w = h = 1025;
  //w = h = 800;

  float r = 18;
  int numLevels = 20;
  int numColors = numLevels + 1;

  //numColors = 3;  // test

  bool antiAlias = true;

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
  f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y + 0.1f*x*x - 0.1f*x - 0.1f*y*y + 0.1f*y; };
    // try exchanging sin and cos an combining
  //f = [&] (float x, float y) { return (float) pow(spiralRidge(x, y, 0.25), 3.0); };
   // exponent 3 makes for good balance between black and white - but middle gray is 
   // underrepresented - todo: apply expansion of middle gray and compression of black/white values

  rsImageContourPlotter<float, float> cp;

  using IP = rsImageProcessor<float>;

  // create image with function values:
  rsImageF imgFunc(w, h);
  generateFunctionImage(f, xMin, xMax, yMin, yMax, imgFunc);;
  IP::normalize(imgFunc);

  // create images with contours:
  std::vector<float> levels = rsRangeLinear(0.f, 1.f, numLevels);
  rsImageF imgCont = cp.getContourLines(imgFunc, levels, { 1.0f }, antiAlias);
  // with anti-aliasing, we need to use about twice as much brightness to get the same visual 
  // brightness

  // create images with bin-fills:
  std::vector<float> colors = rsRangeLinear(0.f, 1.f, numColors);
  rsImageF imgFills = cp.getContourFills(imgFunc, levels, colors, antiAlias);
  // the highest levels are not white but gray - ah: it was because the painter used the saturating
  // mode - saturating mode should *NOT* be used for filling contours!!!

  // write images to files:
  writeScaledImageToFilePPM(imgFunc,  "Function.ppm", 1);
  writeScaledImageToFilePPM(imgCont,  "ContourLines.ppm", 1);
  writeScaledImageToFilePPM(imgFills, "ContourFills.ppm", 1);

  // the right column and bottom row has no contour values - no surprise - the loop only goes up 
  // to w-1,h-1
  // maybe use powers of two +1 for the size and cut off bottom-row and right-column aftewards

  // -with the spirals, we sometimes get spurios, short contour segments that are not connected to 
  //  any contour that continues through the whole image...it must be that the pixels around the 
  //  spurious segemnts happen to be wholly above or below the contour level...due to sampling the
  //  continuous function - remedy: create an image of boolean values that just keep the info, if
  //  there's a contour through the pixel or not and, as post processing, keep only those contour
  //  pixels, whose neighbours are also on a contour - remove all  pixels form the contour for 
  //  which none of the 8 neighbours is on a contour.. we need a function isValidContourPixel
  //  that follows the contour using neighbour-contour info - return true, if the contour comes 
  //  back to the start point or goes out of the image

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
  int w           = 400;               // width in pixels
  int h           = 400;               // height in pixels
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
  //f = [=](Complex z) { return z*z; };
  f = [=](Complex z) { return z*z*z; };
  //f = [=](Complex z) { return z*z*z*z; };
  //f = [=](Complex z) { return 1./(1. + z); };


  // render images of function values for real and imaginary part:
  using IP = rsImageProcessor<float>;
  rsImageF funcRe(w, h), funcIm(w, h), empty(w, h);
  generateFunctionImageReIm(f, xMin, xMax, yMin, yMax, funcRe, funcIm);
  IP::normalizeJointly(funcRe, funcIm); // joint normalization preserves re,im relationship (right?)
  //normalize(funcRe);
  //normalize(funcIm);
  // what about absolute value and phase? ca we do something useful with them, too?

  // get contour lines:
  rsImageContourPlotter<float, float> cp;
  std::vector<float> levels = rsRangeLinear(0.f, 1.f, numLevels);
  for(int i = 0; i < numLevels; i++)
    levels[i] = pow(levels[i], levelPow);
  rsImageF contRe = cp.getContourLines(funcRe, levels, { 1.0f }, antiAlias);
  rsImageF contIm = cp.getContourLines(funcIm, levels, { 1.0f }, antiAlias);

  // get countour fills:
  std::vector<float> colors = rsRangeLinear(0.f, 1.f, numColors);
  rsImageF fillsRe = cp.getContourFills(funcRe, levels, colors, antiAlias);
  rsImageF fillsIm = cp.getContourFills(funcIm, levels, colors, antiAlias);
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

  // todo: overlay the contour lines with the contour fills and/or the original function 
  // separately for re and im
}

template<class T, class TPix> 
void drawConicSection(T A, T B, T C, T D, T E, T F, T xMin, T xMax, T yMin, T yMax,
  rsImage<TPix>& img, TPix color)
{
  // under construction

  function<T(T, T)> f;
  f = [=](T x, T y) { return A*x*x + B*x*y + C*y*y + D*x + E*y + F; };

  // todo: figure out a point x0,y0 on the conic section - maybe set x = 0 and solve for y? or
  // set x to some value between xMin..xMax...but if the resulting value for y is outside the range
  // yMin..yMax, we need to choose a different x

  // we should perhaps figure out, if the conic section is an ellipse - if so, it means it's closed
  // and we need only draw one curve...except when the ellipse is too large to fit on the image - in 
  // such a case, we may have to draw at most 4 separate segments

  // actually, all the newton-iteration stuff is not needed - we can solve explicitly for y in terms
  // of x or vice versa, so a drawConicsection function doe not need to be based on 
  // drawImplicitCurve (but can be - but that's wasteful)

  T x = 0.5*(xMax+xMin);
  //T y0 = ;  
  // todo: solve conic equation for y - we may get two solutions for the two hyperbolas

  // find coeffs for the quadratic equation y = a0 + a1*y + a2*y^2 for that chosen x:
  T a0 = A*x*x + D*x + F;
  T a1 = B*x + E;
  T a2 = C;

  // solve quadratic equation:
  //std::complex<T> y1, y2;
  //rsPolynomial<T>::rootsQuadraticComplex(complex<T>(a0), complex<T>(a1), complex<T>(a2), &y1, &y2);
  // todo: have a function that takes real coeffs and returns complex results

 
  int dummy = 0;

  //drawImplicitCurve(f, xMin, xMax, yMin, yMax, 0.0, x, y1, img, color);
}
// A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0.

void implicitCurve()
{
  int width  = 800;
  int height = 800;
  double range = 2.1;


  using IP = rsImageProcessor<float>;
  rsImageF imgCurve(width, height);
  function<double(double, double)> f;


  rsImagePlotter<float, double> ig;
  ig.setRange(-range, range, -range, range);
  ig.painter.setDeTwist(false);  // should be only used for single pixel lines
  ig.painter.setNeighbourWeightsForSimpleDot(0.375, 0.375*sqrt(0.5));


  // test:
  //drawConicSection(1.0, 0.0, 1.0, 0.0, 0.0, -1.0, xMin, xMax, yMin, yMax, imgCurve, 1.f);

  //f = [=](double x, double y) { return x*x + 1.5*y*y; }; 
  // we need one starting point - maybe the function should figure it out itself

  float color = 0.375f;
  // looks good with the saturating accumulation in rsImagePainter


  f = [=](double x, double y) { return x*x + y*y; };  // unit circle
  ig.plotImplicitCurve(f, 1.0, 1.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return x*x - y*y; };  // unit hyperbola - opens to right
  ig.plotImplicitCurve(f, 1.0, 1.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return x*x - y*y; };  // unit hyperbola - opens to left
  ig.plotImplicitCurve(f, 1.0, -1.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return y*y - x*x; };  // unit hyperbola - opens to top
  ig.plotImplicitCurve(f, 1.0, 0.0, 1.0, imgCurve, color);

  f = [=](double x, double y) { return y*y - x*x; };  // unit hyperbola - opens to bottom
  ig.plotImplicitCurve(f, 1.0, 0.0, -1.0, imgCurve, color);


  f = [=](double x, double y) { return (x*x)/4 + y*y; };  // ellipse with width 2 and height 1
  ig.plotImplicitCurve(f, 1.0, 2.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return x*x + (y*y)/4; };  // ellipse with width 1 and height 2
  ig.plotImplicitCurve(f, 1.0, 0.0, 2.0, imgCurve, color);


  f = [=](double x, double y) { return y - x*x; };  // unit parabola - opens to top
  ig.plotImplicitCurve(f, 0.0, 0.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return y + x*x; };  // unit parabola - opens to bottom
  ig.plotImplicitCurve(f, 0.0, 0.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return x - y*y; };  // unit parabola - opens to right
  ig.plotImplicitCurve(f, 0.0, 0.0, 0.0, imgCurve, color);

  f = [=](double x, double y) { return x + y*y; };  // unit parabola - opens to left
  ig.plotImplicitCurve(f, 0.0, 0.0, 0.0, imgCurve, color);
  // make the imgCurve and color members of the plotter object such that they don't have to be 
  // passed over and over again - this has also other benefits such as being able to set up the
  // world-to-pixel/pixel-to-world coordinate transformations once (in setPlotImage and 
  // setRange) - use coordinate-transformer objects there and allow also for logarithmic scaling
  // of the x- and/or y-axis

  // maybe draw everything except the circle again but rotated by 45°


  // make higher level code: drawCircle(cx, cy, r), drawEllipse(cx, cy, width, height, rotation)
  // drawConicSection(a,b,c,d,e,f)
  // ...what about filling the inside of a curve, i.e points for which f(x,y) <= c

  IP::normalize(imgCurve);
  writeScaledImageToFilePPM(imgCurve, "ImplicitCurves.ppm", 1);

  // -with step = 4 (parameter inside the function), the ellipses are drawn heavier than with 
  //  step = 3 - why? with step = 5, they look as expected again, with 6 they get denser again
  //  (but not as dense as with 4) and here, also the circles have denser dots
  // -also, with step > larger than 1, we begin to see a doubl-drawing of the last pixel again
  //  (maybe we need to scale the weight of the final pixel by 1/step)
}
// other curves to try:
// https://en.wikipedia.org/wiki/Pedal_curve
// https://en.wikipedia.org/wiki/Roulette_(curve)
// https://en.wikipedia.org/wiki/Cyclocycloid
// https://en.wikipedia.org/wiki/Epicycloid
// https://en.wikipedia.org/wiki/Hypotrochoid
// https://en.wikipedia.org/wiki/Epitrochoid


template<class TPix, class TVal>
void plotParametricCurve(const std::function<TVal(TVal)>& fx, const std::function<TVal(TVal)>& fy,
  const std::vector<TVal>& t, rsImage<TPix>& img, TPix color, 
  TVal xMin, TVal xMax, TVal yMin, TVal yMax)
{
  rsImagePainter<TPix, TVal, TVal> painter;
  painter.setImageToPaintOn(&img);
  painter.setDeTwist(true);

  for(size_t i = 0; i < t.size(); i++) {
    TVal x = fx(t[i]);
    TVal y = fy(t[i]);
    x = rsLinToLin(x, xMin, xMax, TVal(0),  TVal(img.getWidth()-1));
    y = rsLinToLin(y, yMin, yMax, TVal(img.getHeight()-1), TVal(0));
    painter.paintDot(x, y, color); }
}


/** For parametric curve given by x(t) = fx(t), y(t) = fy(t) and a given set of (strictly 
increasing) sample-points t, this function computes the arc-length function and stores the values
in s. */
template<class T>
void arcLengthFunction(
  const std::function<T(T)>& fx, const std::function<T(T)>& fy,
  const std::vector<T>& t, std::vector<T>& s)
{
  rsAssert(s.size() == t.size()); // or maybe just resize s - or use raw arrays
  T x0 = fx(t[0]);
  T y0 = fy(t[0]);
  s[0] = T(0);
  for(size_t n = 1; n < t.size(); n++) {
    T x1 = fx(t[n]); T y1 = fy(t[n]);
    T dx = x1 - x0;  T dy = y1 - y0;
    T ds = sqrt(dx*dx + dy*dy);
    s[n] = s[n-1] + ds; 
    x0 = x1; y0 = y1; }
}
// todo: maybe first fill the s-array with the ds-values and then call the trapezoidal integration
// routine - the code above computes a Riemann sum which is less accurate


template<class TPix, class TVal>
void plotParametricCurve(const std::function<TVal(TVal)>& fx, const std::function<TVal(TVal)>& fy,
  TVal t0, TVal t1, int N, rsImage<TPix>& img, TPix color,
  TVal xMin, TVal xMax, TVal yMin, TVal yMax)
{
  using Vec = std::vector<TVal>;
  using AT  = rsArrayTools;

  // Avoid painting the last point twice in case of closed curves:
  bool avoidDoubleDraw = true; // make user parameter
  if(avoidDoubleDraw)
    t1 = t0 + TVal(N-1)*(t1-t0) / TVal(N);

  // Create array of time-stamps:
  Vec t(N);  
  AT::fillWithRangeLinear(&t[0], N, t0, t1);

  // If natural parameterization (i.e. parametrization by arc-length) is desired, warp 
  // t-array accordingly:
  bool natural = true;  // make user parameter
  if(natural) {
    Vec s(N);  // arc-length as function of raw parameter t
    Vec r(N);  // equally spaced values from 0..arcLength
    Vec t2(N); // new, modified array of time-stamps
    arcLengthFunction(fx, fy, t, s);
    AT::fillWithRangeLinear(&r[0], N, TVal(0), s[N-1]); // s[N-1] is total length of the curve
    resampleNonUniformLinear(&s[0], &t[0], N, &r[0], &t2[0], N);
    //rsPlotVectorsXY(t, s);   // arc-length as function of raw parameter t
    //rsPlotVectorsXY(r, t2);  // mapped parameter t as function of arc-length s
    plotParametricCurve(fx, fy, t2, img, color, xMin, xMax, yMin, yMax); }
  else
    plotParametricCurve(fx, fy, t, img, color, xMin, xMax, yMin, yMax);
}

void parametricCurve()
{
  // Plot parameters:
  int width    = 500;
  int height   = 500;
  int numDots  = 2500;
  //int numDots  = 500;
  double range = 1.5;
  bool natural = false;  // natural parametrization (by arc-length s)

  // Lissajous curve parameters:
  double a  = 2.0;
  double b  = 3.0;
  double p  = PI;
  double t0 = 0.0;
  double t1 = 2*PI * 1.0;

  // define functions x = fx(t), y = fy(t)
  std::function<double(double)> fx, fy;
  fx = [&](double t) { return cos(a*t - 0.5*p); };
  fy = [&](double t) { return sin(b*t + 0.5*p); };

  // plot curve:
  rsImageF img(width, height);
  float color = 0.375f;
  plotParametricCurve(fx, fy, t0, t1, numDots, img, color, -range, range, -range, range);
  rsImageProcessor<float>::normalize(img);
  writeScaledImageToFilePPM(img, "ParametricCurve.ppm", 1);

  // -when using numDots such that discrete dots become visible (like 500), those that land on a 
  //  pixel appear brighter than those which don't
  //  -> using de-twisting improves this
}

template<class TPix, class TVal>
void plotFunction(const std::function<TVal(TVal)>& f, rsImage<TPix>& img, TPix color,
  TVal xMin, TVal xMax, TVal yMin, TVal yMax)
{
  // loop over the x-coordinates, compute the next y-coordinate, compute length of segment,
  // compute number of dots to use for this segment, draw the segment

}


void plotFunction()
{
  // test function plotting with some example functions

}



//-------------------------------------------------------------------------------------------------
// image effects:

template<class T>
void logisticCompression(rsImage<T>& img, T k)
{
  T* p = img.getPixelPointer(0, 0);
  for(int i = 0; i < img.getNumPixels(); i++)
    p[i] = T(1) / (T(1) + exp(-k*(p[i]-0.5)));  // 0.5 put middle gray at the center of the sigmoid
}
// move to rsImageProcessor, maybe generalize to:
// https://en.wikipedia.org/wiki/Generalised_logistic_function
// https://en.wikipedia.org/wiki/Gompertz_function
// and/or let the center (the -0.5) be another paraneter

/** Applies an elliptic frame to the image, i.e. darkens the parts of the image that are outside 
the ellipse. The ellipse is centered at (cx,cy) and has x,y radii rx,ry. The "steepness" parameter 
controls, how steep/sharp the transition between the unaffected area inside the ellipse and the 
affected area outside the ellipse is, "shape" controls a sort of crossfade between rectangular and
elliptic shape of the frame, "invert" may be used to switch the roles inside/outside, i.e. darken 
the pixels inside the ellipse. */
template<class T>
void frameElliptic(rsImage<T>& img, T cx, T cy, T rx, T ry, T steepness = T(1), T shape = T(1),
  T amount = T(1), bool invert = false)
{
  for(int j = 0; j < img.getHeight(); j++) {
    for(int i = 0; i < img.getWidth(); i++) {
      T x = T(i); T dx = x - cx; dx /= rx;       // scaled x-distance
      T y = T(j); T dy = y - cy; dy /= ry;       // scaled y-distance
      T dE = dx*dx + dy*dy;                      // elliptic distance
      T dR = rsMax(rsAbs(dx), rsAbs(dy));        // rectangular distance ..rsMin instead of max has also
      T d  = (T(1)-shape) * dR  + shape * dE;    // rectangle vs ellipse
      T s  = T(1) / (T(1) + pow(d, steepness));  // "Butterworth" function
      s    = (1-amount) + amount*s;              // "low-shelf" vs "lowpass"
      if(invert) s = T(1) - s;                   // keep only stuff outside the frame
      img(i, j) *= s; }}                         // apply framing
}
// maybe allow for a rotation (apply rotation matrix between subtracting the center and scaling
//  ...maybe wrap into function weightedDistance(x, y, cx, cy, rx, ry, angle, &dx, &dy)
// maybe raise s to a power (before applying amount) - "apply filter multiple times" - should make
// the ellipse more gaussian - but also smaller - maybe compensate by scaling rx,ry
// we loop over the lines in the outer loop to be more cache-friendly
// maybe call it ellipticFrame and keep only the part in the ellipse - if "invert" is true, use 
// 1-s
// but: should we somehow also normalize the steepness with respect to the image size? ...such 
// that an image with twice the size (using twice rx,, ry) when bein downsampled looks the same
// as the image with normal size? ..make the algorithm invariant with respect to oversampling
// steepness /= oversampling?

// maybe have also a rectangularFrame function (using d = min(dx,dy)?)
// maybe continuously fade between elliptic and rectangular using a shape parameter?
// ...i think, this is the wrong way to try to get a rectangular frame - for this, we should 
// perhaps use one multiplier for the x-distance and another for the y-distance, like
// s  = T(1) / (T(1) + pow(dx, steepnessX));
// s *= T(1) / (T(1) + pow(dy, steepnessY));

// maybe make a function that only keep stuff within an elliptic annulus: apply two frames, one
// of them inverted


// todo: gammaCorrection, rationalMapUnipolar, rationalMapBipolar, blur (bi-direction IIR filter)

void testImageEffectFrame()
{
  int w = 500;
  int h = 300;

  // position and size:
  float cx = 200.f;
  float cy = 100.f;
  float rx = 100.f;
  float ry =  50.f;

  // parameters:
  float steepness = 0.5f; // at 100, it looks like an anti-aliased filled ellipse
  bool invert = false;
  float amount = 1.0f;
  float shape  = 1.0f;

  // apply framing effect to fully white image:
  rsImageF img(w, h);
  img.clear(1.f);      
  frameElliptic(img, cx, cy, rx, ry, steepness, shape, amount, invert);
  writeImageToFilePPM(img, "FrameTest.ppm");
  int dummy = 0;

  // todo: try to plot a contour at level 0.5 - it should be independent from the steepness
}

//-------------------------------------------------------------------------------------------------

void testDistanceMap()
{
  int w = 500;   // image witdh
  int h = 500;   // image height
  int N = 500;   // number of sample points on the curve

  // Liassajou curve parameters
  double a = 2.0;
  double b = 3.0;
  double p = PI;    // phase offset

  double xMin = -1.5;
  double xMax = +1.5;
  double yMin = -1.5;
  double yMax = +1.5;

  // create curve:
  std::vector<double> x(N), y(N);
  for(int n = 0; n < N; n++)
  {
    double t = double(2*PI*n) / double(N);  // div by N gives more symmetry than N-1
    x[n] = cos(a*t - 0.5f*p);
    y[n] = sin(b*t + 0.5f*p);

    // convet to pixel coordinates:
    x[n] = rsLinToLin(x[n], xMin, xMax, 0.0, double(w-1));
    y[n] = rsLinToLin(y[n], yMin, yMax, double(h-1), 0.0);
  }

  // craete images for curve and distance map:
  rsImageF imgDist(w, h), imgCurve(w, h);

  using IP = rsImageProcessor<float>;

  rsImagePlotter<float, double> ig;
  ig.setRange(xMin, xMax, yMin, yMax);


  ig.plotDistanceMap(imgDist, &x[0], &y[0], N);
  IP::normalize(imgDist);
  IP::invert(imgDist);
  IP::sineShape(imgDist);
  IP::gammaCorrection(imgDist, 16.f);

  //gammaCorrection(imgDist);

  // maybe apply an exponent/gamma and/or rational mapping or sinusoidal mapping

  writeImageToFilePPM(imgDist, "DistanceMap.ppm");
  //writeScaledImageToFilePPM(imgDist, "DistanceMapScaled.ppm", 2);

  // -when using too small N, there are artifacts that show up like plotting the curve dotted 
  //  -with N = 2000, w=h=500, these disappear - but the computation takes long

  // -try using different N for the 3 color channels - should be small'ish like 20..50 such that 
  //  the points are clearly visible
  // -llok a bit like a voronoi tesselation - which makes sense
  // -a=3,b=5,N=100 gives a square

  // -maybe try different distance measures: euclidean, manhattan
}

//-------------------------------------------------------------------------------------------------

void plotSpiralHeightProfile()
{
  // Plots the height of the spiralRidge function as function of the radius for a given angle using 
  // a spiral with growth-factor 2.

  // Settings:
  int    N      = 511;  // number of samples
  double rMin   = 1./4; // minimum radius
  double rMax   = 4;    // maximum radius
  double angle  = 0.0;  // angle along which we walk when increasing the radius
  double p      = 0.0;  // phase of the spiral
  double shrink = 2;    // maybe call it grow
  double k      = 1;    // exponent - default: 1

  // Generate data:
  std::vector<double> r(N), h0(N), h1(N), h2(N), h3(N);  // radius and heights
  rsArrayTools::fillWithRangeExponential(&r[0], N, rMin, rMax);
  rsImagePlotter<float, double> ig;
  double a = log(shrink) / (2*PI);
  for(int i = 0; i < N; i++) {
    double x = r[i] * cos(angle);
    double y = r[i] * sin(angle);             // when k == 1:
    h0[i] = ig.spiralRidge(x, y, a, p, 1., 0, k); // triangular
    h1[i] = ig.spiralRidge(x, y, a, p, 1., 1, k); // smooth sinusoid (sin-shaped from triangular)
    h2[i] = ig.spiralRidge(x, y, a, p, 1., 2, k); // rectified sine (original)
    h3[i] = ig.spiralRidge(x, y, a, p, 1., 3, k); // inverted rectified sine
  }

  // Plot results:
  GNUPlotter plt;
  plt.addDataArrays(N, &r[0], &h0[0], &h1[0], &h2[0], &h3[0]);
  plt.setLogScale("x", 2.0);
  plt.setRange(rMin, rMax);
  plt.addCommand("set xtics (0.125,0.25,0.5,1.0,2.0,4.0,8.0)"); // wrap into function setTicsX(vector<double>)
  plt.plot();
}

void testSpiralHeightProfile()
{
  // Create the height profile (as function of radius) via the spiralRidge function and compare it 
  // to a half-period of the sine wave - they do indeed match perfectly.

  int N = 511;  // number of samples
  double g = 2.0; // growth factor
  std::vector<double> r(N), R(N), h(N), s(N);                // radius and height and sine
  double a = log(g) / (2*PI);                                // shrink/grow factor is 2
  rsArrayTools::fillWithRangeExponential(&r[0], N, 1.0, g);  // interval 1...g, log-scaled
  rsArrayTools::fillWithRangeLinear(     &R[0], N, 0.0, PI); // interval 0..PI, lin-scaled
  rsImagePlotter<float, double> ig;
  for(int i = 0; i < N; i++) {
    h[i] = ig.spiralRidge(r[i], 0, a, 0., 1., 2);  // x=r, y=0
    s[i] = sin(R[i]); }
  std::vector<double> err = s - h;   // is numerically zero
  rsPlotVectorsXY(R, h, s, err);     // error is of the order of machine epsilon

  // can we prove/derive that the height-profile is a rectified sine? idea:
  // -consider the case shrink = 2 and the look at the interval 1..2 on the x-axis (i.e. set y=0)
  // -plug x into height-function and simplify
  // -the results generalize to arbitrary angles by taking a different interval for the exponential 
  //  range - it widens from 1..2 to 2..4 during one revolution (i think, after 180° it's at 
  //  sqrt(2)..2*sqrt(2) = 1.414..2.828 - verify)
  // -they also generalize to different growth factors g by not taking 1..2 but 1..g as the basic 
  //  interval on the x-axis
}

/** Structure to represent parameter for an image generation/processing algorithm which may have 
different values for the 3 color channels. */
template<class T>
struct rsParamRGB
{
  rsParamRGB(T _r, T _g, T _b) : r(_r), g(_g), b(_b) {}
  T r = 0, g = 0, b = 0;
};

/** Structure to represent paremeters for generating a spiral image. */
struct rsSpiralParams
{
  rsParamRGB<double> phase, grow, sign, gamma, exponent;
  rsParamRGB<int>    profile;
  double range;
};

void generateSpiralImage(const rsSpiralParams& p, rsImageF& R, rsImageF& G, rsImageF& B)
{

  /*
  std::function<double(double, double)> f;
  f = [&](double x, double y) 
  { 
    double s = pow(spiralRidge(x, y, density, phase, sign, profile), power); 
    return s;
  };
  */

}



void spirals()
{
  //plotSpiralHeightProfile();
  //testSpiralHeightProfile();
  //testImageEffectFrame(); return;
  testDistanceMap(); return;

  //int size = 1000;
  int w = 1200;
  int h = 800;

  double range    = 1.5;

  //double phase    = 120.0;
  //double phaseInc = 120;   // 120 is a nice defailt


  double power = 1.5;   // 3 lead to balance between black/white - lower value give mor white
    // this seems to be 1/gamma - try gamma-correction with irfan view

  double phaseR = 0;
  double phaseG = 120;
  double phaseB = 240;

  //double density  = 0.15;  // smaller values -> denser arms
  //double densInc  = 0.5;   // relative density increment

  double shrinkR = 3 / 1.4;
  double shrinkG = 3;
  double shrinkB = 3 * 1.4;

  double signR = +1.0;
  double signG = -1.0;
  double signB = +1.0;

  int profileR = 0;
  int profileG = 0;
  int profileB = 0;

  // maybe it would make sense to have a struct rsParamRGB { double r; double g; double b; }
  // that can hold the 3 parameter values for the 3 channels


  //double alt      = +1;    // -1: alternate directions, +1: don't alternate


  // test:
  //double density = log(shrinkR) / (2*PI);   // rename;
  double density, phase, sign; 
  int profile;
  // leads to shrinking of 1/2 per revolution - maybe the user parameter should be the shrink-factor
  // and the algo parameter a = log(shrinkFactor) / (2*PI) - nad maybe instead of incrementing
  // a linearly, we should have a "spread" factor that's used like 
  // shrinkRed = shrink/shrinkSpread, shrinkGreen = shrink, shrinkBlue = shrink*shrinkSpread
  // or better: have shrinkR, shrinkG, shrinkB and maybe a global shrink that applies to all
  // maybe instead of 2, use the golden ratio



  //// convert user to algo params
  //densInc *= density;  // increment should be relative
  //phase = rsDegreeToRadiant(phase);
  //phaseInc = rsDegreeToRadiant(phaseInc);


  rsImagePlotter<float, double> ig;
  std::function<double(double, double)> f;
  f = [&](double x, double y) 
  { 
    double s = pow(ig.spiralRidge(x, y, density, phase, sign, profile), power); 
    //double s = pow(spiralRidge2(x, y, density, phase, sign), power); 
    return s;

    // compression of black/white, expansion of gray:
    //double k = 1.0;  // amount of compression
    //s = 1 / (1 + exp(-k*s));
    //return s;
    // nope - compression needs to be applied *after* normalization
  };




  // create image with function values:
  using IP = rsImageProcessor<float>;
  rsImageF red(w, h), green(w, h), blue(w, h);


  density = log(shrinkR) / (2*PI);
  phase   = phaseR;
  sign    = signR;
  profile = profileR;
  generateFunctionImage(f, -range, range, -range, range, red);
  IP::normalize(red); 

  density = log(shrinkG) / (2*PI);
  phase   = phaseG;
  sign    = signG;
  profile = profileG;
  generateFunctionImage(f, -range, range, -range, range, green);
  IP::normalize(green); 

  density = log(shrinkB) / (2*PI);
  phase   = phaseB;
  sign    = signB;
  profile = profileB;
  generateFunctionImage(f, -range, range, -range, range, blue);
  IP::normalize(blue);  

  //// test:
  //float  comp = 4.0;
  //logisticCompression(red,   comp);
  //logisticCompression(green, comp);
  //logisticCompression(blue,  comp);
  //// hmm ..doesn't really look good - we need a different compression function

  writeImageToFilePPM(red,   "SpiralsR.ppm");
  writeImageToFilePPM(green, "SpiralsG.ppm");
  writeImageToFilePPM(blue,  "SpiralsB.ppm");



  writeImageToFilePPM(red,  green, blue, "Spirals.ppm");


  rosic::writeToMonoWaveFile("Spirals.wav", red.getPixelPointer(0,0), red.getNumPixels(), 44100);
  // todo: make stereo file from two of the channels
  // -maybe make a funtion writeImageToFileWav (this should include the mapping 0..1 -> -1..+1)

  // -when the parameter t in the parametric equation:
  //    x(t) = exp(a*t)*cos(t), y(t) = exp(a*t)*sin(t)
  //  increases by 2*pi, we make one full revolution around the spiral, but the radius will be 
  //  multplied exp(a*2*pi)
  // -this means that a pic that results by a given range r such that xMin=yMin=-r, xMax=yMax=r 
  //  will give the same image as the range r * exp(a*2*pi)
  // -this can be used for an infinite-zoom animation as follows:
  // -set range1 = ..., range2 = range1 * exp(a*2*pi)
  // -create images with ranges = rsLinToExp(rsRangeLinear(0,1), 0,1, range1, range2)
  // -the last image will look the same as the first, so we may loop the animation between frame1
  //  and last frame (but skip the duplicate frame!)
  // -but this infinite zoom will work only when the 3 spirals all have the same density parameter
  //  -if not, the periodicity will be larger (i think, it will be given by the gcd of ar,ag,ab 
  //   which are the a-paremeters of r,g,b channels, maybe scaled up to become all integers) - so 
  //   we need gcd(ar,ag,ab) revolutions instead of just one - i think


  // -fun: let one of the colors rotate in the other direction
  // -maybe combine it with a spiral that rotates in the other direction
  // -try color-inversion (per channel)
  // -nice: use densities: 0.25,0.30,0.35
  // -move the whole algo into a function that takes a bunch of parameters and returns the image

  // -try to interpret the color-channels not as r,g,b but as c,m,y where c=g+b, m=r+b, y=r+g
  //  -> solve for r,g,b - make function rgb2cmy, cmy2rgb

  // -maybe wrap this into a python function that we may use in jupyter - so we can qucikly vary 
  //  the parameters with sliders
  // -make animations with sweeping the phase from 0 to 360 - then loop

  // todo: 

  // -try a = log(2) / (2*pi) (or something) - this should lead to a shrink/grow factor of 2 per 
//    revolution:  t=0: (x,y)=(1,0), t=(2*pi): (x,y)=(.5,0), t=(4*pi): (x,y)=(.25,0), ...
  // -plot distance as function of radius r for a given angle phi to see how it oscillates - maybe
  //  with some sort of waveshaping we can make this oscillation sinuosidal
  //  ->try to use a rational mapping of the heights to expand the middle-gray range
}
// exponent 3 makes for good balance between black and white - but middle gray is 
// underrepresented - todo: apply expansion of middle gray and compression of black/white values

template<class T>
int numIterationsToDivergence(std::function<T(T, T)> fx, std::function<T(T, T)> fy, T x0, T y0,
  T thresh, int maxNumIterations)
{
  T x = x0;
  T y = y0;
  for(int i = 0; i <= maxNumIterations; i++) {
    T r2 = x*x + y*y;   // squared radius
    if(r2 > thresh)     // should it be >? maybe - the threshold itself may be the stability limit
      return i;
    r2 = fx(x, y);      // r2 used as temporary for x
    y  = fy(x, y);
    x  = r2; }
  return -1;  // -1 indicates that threshold was not exceeded in any iteration
}
// todo: maybe instead of returning just the number of iterations until divergence is assured, 
// return some more elaborate information about the trajectory - for example: 
// -the total (squared?) distance taken in all the steps

void mandelbrot(rsImage<float>& img, int maxIterations, 
  double xMin, double xMax, double yMin, double yMax)
{
  double cx, cy;
  std::function<double(double, double)> fx, fy;
  fx = [&](double x, double y) { return x*x - y*y + cx; };
  fy = [&](double x, double y) { return 2*x*y     + cy; };
  for(int j = 0; j < img.getHeight(); j++) {
    for(int i = 0; i < img.getWidth(); i++) {
      cx = rsLinToLin((double)i, 0.0,  double(img.getWidth()-1), xMin, xMax);
      cy = rsLinToLin((double)j, double(img.getHeight()-1), 0.0, yMin, yMax);
      int its = numIterationsToDivergence(fx, fy, 0.0, 0.0, 4.0, maxIterations);
      img(i, j) = float(its); }}
      // maybe do: if(its >= 0) -> black, else white
}

void fractal()
{
  int w = 400;
  int h = 400;

  int numIterations = 1000;

  double xMin = -1.5;
  double xMax = +0.5;
  double yMin = -1.0;
  double yMax = +1.0;

  using IP = rsImageProcessor<float>;
  rsImageF img(w, h);
  mandelbrot(img, numIterations, xMin, xMax, yMin, yMax);
  IP::normalize(img);
  // maybe apply gamma, contrast, etc.

  writeImageToFilePPM(img, "Mandelbrot.ppm");
}



// maybe make animations with
// http://www.softpedia.com/get/Multimedia/Graphic/Graphic-Others/APNG-Anime-Maker.shtml
// https://stackoverflow.com/questions/3191978/how-to-use-glut-opengl-to-render-to-a-file/14324292#14324292

// or:
// https://www.codeproject.com/Articles/4169/A-simple-interface-to-the-Video-for-Windows-API-fo
// https://www.learnopencv.com/read-write-and-display-a-video-using-opencv-cpp-python/


// http://www.adp-gmbh.ch/win/programming/avi/index.html
// http://www.adp-gmbh.ch/win/programming/avi/avi.html
// http://www.adp-gmbh.ch/win/programming/graphic/bitmap.html#bitmap_h


// see also:
// https://en.wikipedia.org/wiki/Comparison_of_video_container_formats

// https://www.quora.com/What-are-the-two-most-common-uncompressed-video-formats


// try this!
// https://rosettacode.org/wiki/Mandelbrot_set#C
// https://rosettacode.org/wiki/Dragon_curve#C
// https://rosettacode.org/wiki/Sierpinski_carpet#C


// https://rosettacode.org/wiki/Bitmap/B%C3%A9zier_curves/Quadratic#C
// https://rosettacode.org/wiki/Bitmap/B%C3%A9zier_curves/Cubic#C


// https://rosettacode.org/wiki/Window_creation#C

// https://rosettacode.org/wiki/Conway%27s_Game_of_Life#C