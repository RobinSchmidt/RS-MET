using namespace RAPT;

void colorGradientHSL()
{
  int  w = 200;             // pixel width
  int  h = 100;             // pixel height

  using Real  = double;
  using RGB   = rsColorRGB<Real>;
  using HSL   = rsColorHSL<Real>;
  using Color = rsColor<Real>;

  // Create an image with a bilinear gradient in HSL colorspace. We pick a fixed hue and let the
  // lightness gradient go from left to right and the saturation gradient from top to bottom:
  Real H = 210.0 / 360.0;   // hue
  HSL tl(H, 0.0, 0.0);      // top-left        use generic Color
  HSL tr(H, 0.0, 1.0);      // top-right
  HSL bl(H, 1.0, 0.0);      // bottom-left
  HSL br(H, 1.0, 1.0);      // bottom-right
  rsImageF imgR(w, h), imgG(w, h), imgB(w, h);
  for(int iy = 0; iy < h; iy++)
  {  
    Real  y  = Real(iy) / Real(h-1);
    HSL cl = (1.0-y)*tl + y*bl;      // color left
    HSL cr = (1.0-y)*tr + y*br;      // color right
    for(int ix = 0; ix < w; ix++)
    {
      Real x   = Real(ix) / Real(w-1);
      HSL  hsl = (1.0-x)*cl + x*cr;
      RGB  rgb;
      Color::hsl2rgb(hsl.x, hsl.y, hsl.z, &rgb.x, &rgb.y, &rgb.z);
      imgR(ix, iy) = float(rgb.x);
      imgG(ix, iy) = float(rgb.y);
      imgB(ix, iy) = float(rgb.z);
    }
  }
  writeImageToFilePPM(imgR, imgG, imgB, "GradientLightnessSaturationHSL.ppm");


  Color lch(0.8, 0.6, 0.2);
  Color rgb;
  Color::lch2rgb(lch.x, lch.y, lch.z, &rgb.x, &rgb.y, &rgb.z);

  // Now do a similar thing in LCH (lightness, chroma (~saturation), hue)
  tl.x = 0.f; tl.y = 0.f; tl.z = H;
  tr.x = 1.f; tr.y = 0.f; tr.z = H;
  bl.x = 0.f; bl.y = 1.f; bl.z = H;
  br.x = 1.f; br.y = 1.f; br.z = H;
  for(int iy = 0; iy < h; iy++)
  {  
    Real  y  = Real(iy) / Real(h-1);
    Color cl = (1.0-y)*tl + y*bl;      // color left
    Color cr = (1.0-y)*tr + y*br;      // color right
    for(int ix = 0; ix < w; ix++)
    {
      Real  x   = Real(ix) / Real(w-1);
      Color lch = (1.0-x)*cl + x*cr;
      RGB rgb;
      Color::lch2rgb(lch.x, lch.y, lch.z, &rgb.x, &rgb.y, &rgb.z);
      imgR(ix, iy) = float(rgb.x);
      imgG(ix, iy) = float(rgb.y);
      imgB(ix, iy) = float(rgb.z);
    }
  }
  writeImageToFilePPM(imgR, imgG, imgB, "GradientLightnessChromaLCH.ppm");




  //tl.x = 0.f; tl.y = 0.f; tl.z = H;

  // https://de.wikipedia.org/wiki/LCh-Farbraum




  int dummy = 0;

  // Observations:
  // -top pixel row is grayscale gradient
  // -left column is black, right column is white
  // -pure colors occur at center bottom
  // -hues (in degrees): 0: red, 30: orange, 60: yellow, 90: greenyellow, 120: green, 
  //  150: mint, 180: cyan, 210: sky, 240: blue, 270: purple, 300: magenta, 330: pink
  //
  // ToDo:
  // -create a gradient with fixed luminance of 0.5 - with full saturation, this should give the
  //  pure colors. hue should go from left to right, saturation from top to bottom (top should be
  //  gray/unsaturated)
  // -Implement unit test -> should test roundtrips and a couple of conversion examples with known
  //  target results


  // Ideas:
  // -other color space: 
  //    saturation := max(r,g,b) - min(r,g,b)
  //    luminance  := sqrt(wr*r^2 + wg*g^2 + wb*b^2)     // with weights wr,wg,wb
  //    hue        := 
}

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
  //int i;

  // create arrays for line endpoints:
  vector<float> x1, y1;
  for(int i = 0; i < numLines; i++){    // flat, horizontal'ish
    x1.push_back(imageWidth - margin);
    y1.push_back(margin + i * (imageHeight - margin) / numLines); }
  x1.push_back(imageWidth -margin); // 45° diagonal
  y1.push_back(imageHeight-margin);
  for(int i = 0; i < numLines; i++){    // steep, vertical'ish
    x1.push_back(margin + i * (imageWidth - margin) / numLines);
    y1.push_back(imageHeight - margin); }

  //// dotted algorithm:
  //for(i = 0; i < x1.size(); i++)
  //  painter.drawLineDotted(x0, y0, x1[i], y1[i], brightness, brightness, numDots); // we need to pass a number of dots
  //writeImageToFilePPM(image, "LinesDotted.ppm");

  // Wu algorithm:
  image.clear();
  for(size_t i = 0; i < x1.size(); i++)
    painter.drawLineWu(x0, y0, x1[i], y1[i], brightness);
    //drawLineWuPrototype(image, x0, y0, x1[i], y1[i], brightness);
  writeImageToFilePPM(image, "LinesWu.ppm");

  // Bresenham algorithm:
  image.clear();
  for(size_t i = 0; i < x1.size(); i++)
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
  // Tests rasterization based drawing of filled polygons via triangles.

  // Setup:
  bool antiAlias = true;

  // create and set up objects and parameters:
  typedef rsVector2DF Vec2;    // for convenience
  typedef rsVector2DF V;       // even shorter - for inlined constructor calls
  float c = 1.0f;              // color (gray value)
  rsImageF img(35, 20);        // image to draw on
  rsImageDrawerFFF drw(&img);  // drawer object
  drw.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);


  void (*pDrawTriangle)(rsImageDrawerFFF&, const Vec2&, const Vec2&, const Vec2&, float);
  if(antiAlias) pDrawTriangle = &drawTriangleAntiAliasedProto;
  else          pDrawTriangle = &drawTriangle;

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

  // ToDo:
  // -check, if all possible cases are covered
  // -make it efficient (compare to prototype code in unit test)
  // -move to library
  // -compare results of anti-aliased drawing to non-antialiased, oversampled drawing with 
  //  subsequent downsampling

  // see also:
  // https://bisqwit.iki.fi/jutut/kuvat/programming_examples/polytut/
  // https://www.youtube.com/watch?v=PahbNFypubE

  // raylib:
  // https://gamefromscratch.com/raylib-4-released/
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

/*
int gradientifyFlatRegions(const rsImageF& in, rsImageF& out, int numTrips = 1) 
{
  // maybe rename to smoothContours or gradientifyFlatRegions

  using Vec2D = rsVector2D<int>;

  int maxIts = 100;   // make parameter, maybe return the number of iterations taken
  float tol  = 1.e-5;
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
  //  the neighbors of the pixels.
  std::vector<Vec2D> F, B;           // sets of (F)lat, (B)oundary
  rsImage<char> C(w, h);             // pixel classification matrix
  static const char rest     = 0;    // symbolic constants used in the code below
  static const char flat     = 200;
  static const char boundary = 100;
  C.fillAll(rest);                   // initially, all are "rest"

  // In a first pass, we identify the flat regions:
  auto isFlatSlow = [](int i, int j, const rsImageF& img)
  {
    // This is the slow version of the function that needs to be used in the first pass. Later,
    // we can use the C-matrix to retrieve that information faster.
    float p = img(i, j);  // pixel value
    if(p != img(i-1,j) || p != img(i+1, j) || p != img(i,j-1) || p != img(i,j+1))
      return false;
    if(p != img(i-1,j-1) || p != img(i-1, j+1) || p != img(i+1,j-1) || p != img(i+1,j+1))
      return false;
    return true;
    // todo: allow tolerance
  };
  for(j = 1; j < h-1; j++) {
    for(i = 1; i < w-1; i++) {
      if(isFlatSlow(i, j, in)) {
        F.push_back(Vec2D(i,j));
        C(i,j) = flat;  }}} 
  auto isFlat = [&](int i, int j) { return C(i,j) == flat; }; // now faster!
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
  writeImageToFilePPM(C, "PixelClasses.ppm");  // for debug - looks good!

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
  };

  int maxItsTaken = 0;
  float step = 1.f;    // may not be needed, maybe get rid - but first, let's experiment with it
                       // a little bit to see if it can be used to accelerate convergence

  for(int i = 1; i <= numTrips; i++)
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

  return maxItsTaken;

  // -Can we speed up the convergence? maybe it's actually not such a good idea to work in place?
  //  Try to compute all the updates first and then do all the upates at once. compare convergence
  //  to what we do now (updating every pixel immediately after the update was computed, such that 
  //  into computation of the next pixel, the current pixel enters already with updated color)
  // -It seems to have problems when the boundaries of the flat regions are anti-aliased. When 
  //  applied to the image generated by contours, it doesn't seem to do much. The flat regions are
  //  correctly classified. But i think, we have problems with the boundaries because the different
  //  flat color regions do not border each other directly. There's a thin (1-pixel wide) 
  //  transition which is either a falt intermediate color (when anti-alias is turned off) or an
  //  actual gradienty thing that looks like proper anti-aliasing. Even when anti-alias is turned
  //  off we do some sort of crude, improper transition. Both of them trip up the algo.
  //  -The classification with and without AA looks almost the same, but there a very few pixels
  //   that get classified differently in both cases (i spotted just 1 in high zoom - they are 
  //   *really* rare)
  //  -I think the iteration doesn't really do very much in these cases because the flat region is
  //   already at the right color...but wait...no...this makes no sense
  //  -maybe we need a 3rd class of pixels: transition pixels. these are those which are neither in
  //   a flat region nor directly at its boundary but within the 1-pixel wide transition zone
}
*/
// move to prototypes -  done

void fillRectangle(rsImageF& img, int x0, int y0, int x1, int y1, float color)
{
  rsImageDrawerFFF drawer(&img);
  drawer.setColor(color);
  for(int j = y0; j <= y1; j++)
    for(int i = x0; i <= x1; i++)
      drawer.plot(i, j, 1.f);
}
// move to prototypes/drawing

/** A test image with 3 flat color regions...tbc... */
rsImageF testImg3Regions(int w, int h)
{
  rsImageF img(w, h);
  fillRectangle(img, w/2, 0,   w-1,   h-1,   1.0f); // white region on the right
  fillRectangle(img, w/4, h/4, 3*w/4, 3*h/4, 0.5f); // gray region in the center
  return img;
  // ToDo:
  // -Let the user specify the 3 colors, maybe templatize to allow also RGB colors
  //  -Then we also explicity need to fill the left region that is now implicitly "filled" black by
  //   leaving it as initialized.
}

rsImageF testImgVerticalStripes(int w, int h, int numStripes)
{
  rsImageF img(w, h);
  for(int i = 0; i < numStripes; i++)
  {
    int   L = ( i   *w) / numStripes;          // left
    int   R = ((i+1)*w) / numStripes - 1;      // right
    float c = float(i)  / float(numStripes-1); // color, gray value
    fillRectangle(img, L, 0, R,  h-1, c); 
  }
  return img;
}

void gradientify()
{
  // Tests the "gradientify" algorithm that turns flat regions in an image into smooth gradients.

  // User parameters for the test image:
  int s = 3;      // scaler to control size conveniently
  int w = s*80;   // width in pixels
  int h = s*60;   // height in pixels


  // Helper function to compute output of gradientification algo wnad write the result to disk:
  auto computeResult = [&](const rsImageF& imgIn, int numIterations)
  {
    rsImageF imgOut(imgIn.getWidth(), imgIn.getHeight());
    int maxItsTaken = gradientifyFlatRegions(imgIn, imgOut, numIterations);
    std::cout << "Max iterations taken: " + to_string(maxItsTaken) + "\n";
    std::string name = "gradientifyOut" + to_string(numIterations) + ".ppm";
    writeImageToFilePPM(imgOut, name.c_str());
    rsPlotArray(imgOut.getPixelPointer(0, h/2), w);
    return maxItsTaken;
  };

  // Create the test input image and write it to disk:

  //rsImageF imgIn = testImg3Regions(w, h);

  w = 1500; h = 50;
  rsImageF imgIn = testImgVerticalStripes(w, h, 5);
  //fillRectangle(imgIn, 0, h/2-h/6, w-1, h/2+h/6, 0.5f); // add horizontal stripe of gray


  writeImageToFilePPM(imgIn, "gradientifyIn.ppm");

  // Create outputs of the gradientify algorithm with a different setting for the number of 
  // iterations and write the results to disk (Shlemiel strikes again):
  int maxIts;
  //computeResult(imgIn, 1);
  //computeResult(imgIn, 2);
  //computeResult(imgIn, 3);
  //computeResult(imgIn, 4);
  //computeResult(imgIn, 8);
  maxIts = computeResult(imgIn, 25);

  // Report completion to console:
  //std::cout << "Max iterations taken: " + to_string(maxIts);
  rsPrintLine("gradientify() done");


  // Observations:
  // -The larger we choose s, the more iterations are needed to get a smooth gradient without 
  //  artifacts. With s=1, the artifacts are resaonably removed already after 3 iterations. With 
  //  s=3 and 3 iterations, the original class boundaries are quite smeared out but still easily
  //  visible - at least with the test image with gray area in the middle. For the test with the
  //  vertical stripes, it actually seems to converge quickly also for larger images. But as soon 
  //  as we add the additional horizonatl stripe of gray, it again converges slowly. Maybe it's
  //  problematic for convergence when 3 (or more) color regions meet in a single point? That seems
  //  plausible because also in the 3-region test, the artifacts are most severe at the 3-color 
  //  meeting points.
  //  -Maybe the algorithm can be optimized by first running it on a downsampled version and using
  //   the result as a sort of "initial guess" for a higher resolution version. That would be 
  //   similar to multigrid methods for solving PDEs.
  //  -When using 3 vertical stripes and a size of 1200x180, we can clearly see the effect: the 
  //   transition between black and gray is fast (around 120 pixels wide), then we have a flat 
  //   region of gray and then again a fast transition between gray and white. Maybe we should try
  //   it with a 1D version of the algo to figure out, how it behaves. Define a piecewise 1D 
  //   function with 3 levels and do throw sort of 1D heat equation solver at it..or maybe just fix
  //   the boundaries at two positions (to 0 and 1, say) and do the diffusion iteration and look at
  //   the final shape. ...Yes! Plotting the middle row of the gradientified vertical-stripes image
  //   as 1D function clearly reveals this "smoothed stairstep" like shape. With w=1500, h=100,
  //   numStripes=5, it's nicely visible. We see inflection points at 300,600,900,1200 and saddles
  //   in between. Reducing w to 500 clearly reduces the "saddleness" or "stairsteppyness" feature.
  //   The result looks more like a straight linear transition with some slight wavy/wobblyness. 
  //   It's apparently the absoluate size in pixels of the flat regions that determines whether or
  //   they will be sufficiently gradientified or remain (almost) flat. Maybe the tolerance it too
  //   high?
  //  -Bumping up maxIts = 2000; and reducing the tolerance to tol = 1.e-7f; clearly reduces the
  //   saddleness a bit (visible only in direct comparison of the 1D plots). It seems to be a 
  //   problem of painfully slow convergence and we indeed stop too early.
  //  -Maybe we need indeed a "multigrid" approach. I think, for downsampling, we should use a 
  //   naive algo without averaging because the averaging is like an anti-aliasing that may thwart
  //   the classification of pixels on the boundaries between flat regions.
  //
  // ToDo:
  // -Try to figure out why the convergence is so slow when the flat regions have a large absolute
  //  size. Try to speed it up. One idea could be to try to do a kind of "multigrid" method - solve
  //  the problem at increasing resolutions where a suitably interpolated lower-resolution result 
  //  is used as initial guess for the next higher resolution.
  //  -Maybe use such a multigrid method only for the flat-region pixels but not for the boundary
  //   pixels. 
  // -To the helper function applyFilter, we we pass "out" for both parameters "in" and "out", i.e. 
  //  we directly overwrite in-place. Maybe that is not a good idea and we should use separate 
  //  arrays?
  // -Try to make a 1D implementation and check, if it has the same convergence problems. If so, we
  //  may try to find tricks speed up the convergence for 1D first and hope that it will be 
  //  applicable to the 2D case, too.
  // -Test it with more complex input images. It seems to work well with this particular test image
  //  but not so well with the Newton fractal rendering. Try some test inputs with a complexity 
  //  somewhere in between to figure out what the problem may be.
  // -Test with using striped patterns, maybe make a function that creates these images
  // -Use a pattern with vertical stripes
  // -Producing images with various numbers of iterations is actually a "Shlemiel the painter"
  //  algorithm. Maybe optimize that.
  // -Try it with the center section being black or white, too. Also, try to use RGB in various
  //  combinations for the 3 regions
  // -Try the vertical sripes with some horizontal stripes overlaid
  // -The boundary pixels remain unmodified by the algo but one pixel in, we already see a nice 
  //  smoothing effect so it seems workable to deal with boundary pixels by creating a temporary
  //  image that adds a 1-pixel wide frame to the original image which repeats the boundary pixel, 
  //  do the computations on that image and afterwards crop back to the original size, i.e. remove
  //  the frame. This strategy would also make it easier to try a variant of the algo that uese 5x5
  //  neighborhoods (instead of 3x3) to identify flat regions because then we would just need to 
  //  add a 2 pixel wide frame and run the classifier only over interior pixels saving us from the 
  //  mess of treating all the various edge cases differently.

}

void contours()
{
  //testDeContourize();  // preliminary


  // We plot the 2D function z = f(x,y) = x^2 - y^2 into an image where the height translates
  // to the pixel brightness

  rsAssert(testContourSubPixelStuff()); // maybe move to unit test?

  int w = 129;               // width in pixels
  int h = 129;               // height in pixels

  //w = h = 100;
  //w = h = 513;
  w = h = 500;
  //w = h = 1025;
  //w = h = 800;

  float r = 18;              // range for x and y: x = -r..+r, y = -r..+r
  int numLevels = 20;
  bool antiAlias = false;


  std::function<float(float, float)> f;

  // choose your function here:
  //f = [&] (float x, float y) { return x*x + y*y; };    // Circles
  //f = [&] (float x, float y) { return x*x - y*y; };    // Hyperbolas

  // Elliptic curves:
  //f = [&] (float x, float y) { return x*x*x + y*y; }; r = 1.5; numLevels = 40;
  //f = [&] (float x, float y) { return x*x + y*y*y; }; r = 1.0; numLevels = 40; // Elliptic curve

  // Cassini curves:
  f = [&] (float x, float y) { return (x*x+y*y)*(x*x+y*y) - 2*(x*x-y*y) + 1; }; r = 1.0; numLevels = 21;
  // Drawing range is not yet optimal
  // would perhaps be better to have an x-range from -2..+2 and a y-range from -1..+1
  // we should make sure that 1 is mong the levels in order to see the lemniskate
  // https://de.wikipedia.org/wiki/Cassinische_Kurve
  // https://en.wikipedia.org/wiki/Cassini_oval

  // A landscape derived from taking the reciprocal of the Cassini curve landscape - has poles 
  // where the original has zeros:
  //f = [&] (float x, float y) 
  //{ 
  //  float z = (x*x+y*y)*(x*x+y*y) - 2*(x*x-y*y) + 1; 
  //  return rsMin(pow(1.f/z, 0.3f), 5.f);
  //}; 
  //r = 2.0; numLevels = 21;
  // hmmm..the contouring algo finds a black (i.e. zero valued) contour where we expect maximum 
  // values. Is there some wrap-around going on? Maybe write a function isNormalized that checks
  // if the minimum value is zero and the maximum is 1. maybe 
  // rsArrayTools::spansRange(*x, N, min, max) - returns true, if the min and max of the array is
  // equal to the desired min/max values - then use an assert to figure out if the max has wrapped 
  // around


  //f = [&] (float x, float y) { return x*x - y*y + 2.f*x*y; };
  //f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y; };
  //f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y + 0.1f*x*x - 0.1f*y*y; };
  //f = [&] (float x, float y) { return x*sin(y) + y*cos(x) + 0.1f*x*y + 0.1f*x*x - 0.1f*x - 0.1f*y*y + 0.1f*y; };
    // try exchanging sin and cos an combining
  //f = [&] (float x, float y) { return (float) pow(spiralRidge(x, y, 0.25), 3.0); };
   // exponent 3 makes for good balance between black and white - but middle gray is 
   // underrepresented - todo: apply expansion of middle gray and compression of black/white values


  // A nice rational function which shows a lot of features that may be interesting for
  // demo/visualization plots (minima, maxima, saddles), especially in the context of gradient 
  // descent:
  r = 4;
  f = [&](float x, float y)
  {
    float px = x * (x+1) * (x-1);  // p(x): polynomial in x with roots at -1,0,+1
    float py = y * (y+1) * (y-1);  // same in y
    float num = 1*px + 2*py + 1*x*y;
    float den = 1 + x*x*x*x + y*y*y*y;  
    //den = 1 + x*x + y*y; den *= den;  // a variation
    return num/den;
  };
  // The function doesn't feature any poles though. That's appropriate for some sort of plausible 
  // "terrain"


  float xMin = -r;
  float xMax = +r;
  float yMin = -r;
  float yMax = +r;
  int numColors = numLevels + 1;
  //numColors = 3;  // test


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
  writeScaledImageToFilePPM(imgFunc,  "ContourInput.ppm",  1);
  writeScaledImageToFilePPM(imgCont,  "ContourLines.ppm",  1);
  writeScaledImageToFilePPM(imgFills, "ContourFills.ppm",  1);



  /*
  // This takes long and does not work very well yet:

  // Try to undo the contouring by an experimental decontourize algorithm:
  rsImageF grad1(w, h);
  gradientifyFlatRegions(imgFills, grad1, 1);
  rsImageF grad25(w, h);
  gradientifyFlatRegions(imgFills, grad25, 25);

  // write images to files:
  writeScaledImageToFilePPM(grad1,    "ContourGrad1.ppm",  1);
  writeScaledImageToFilePPM(grad25,   "ContourGrad25.ppm", 1);
  */

  rsPrintLine("contours() done.");

  int dummy = 0;


  // Notes:
  // -The (re)gradientification works only with antiAlias turned off because otherwise the 
  //  detection of boundary pixels within gradientifyFlatRegions fails (i think). Oh - but it still
  //  doesn't work well because even without anti-aliasing, we do this strange thing with the 1 
  //  pixel wide boundary line of the intermediate color in the contour plot. Maybe make it 
  //  optional to use this. If the pixel can't be clearly classified to beloongig to one of the 
  //  2 neighboring levels, just pick one (maybe the nearest or something). It's the rudimentary
  //  anti-aliasing in rsImageContourPlotter::fillBetweenContours in the if(!antiAlias) branch. The
  //  coloring of pixels that lie *on* a contour is weird. I tried this with the Cassini curve 
  //  100x100. With 500x500, we see additional artifacts at the image boundaries. They are not 
  //  properly handled an the artifacts bleed into the interior. But perhaps that's due to the fact
  //  that the conturization algo itself produces an output with a black frame. Maybe that needs to
  //  operate on an extended image, too. The "bleeding edge" artifacts are correllated with the 
  //  PixelClasses.ppm image: where the black classes bleed into the image, there's no bleed in the
  //  result and vice versa. But apart from these artifacts, the "decontourized" version actually
  //  looks quite good. We want ContourGrad25,ppm to be close to ContourInput.ppm which it indeed 
  //  is aside from the artifacts.


  // Bugs:
  // -The antiAlias flag seems to affect only the outlines. The fills look always antialiased.
  //  ...ah - no - not really, without anti-alias, they have a 1-pixel thick boundary line with
  //  a mixed-color. this looks a bit smoothish, too - and it trips up the gradientifyFlatRegions
  //  just as much as proper anti-aliasing


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
  int w           = 800;               // width in pixels
  int h           = 800;               // height in pixels
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
  //f = [=](Complex z) { return z*z*z; };
  //f = [=](Complex z) { return z*z*z*z; };
  //f = [=](Complex z) { return 1./(1. + z); };

  // Some strange function i came up with when trying to invent non-orthogonal coordinate systems
  // for practicing calculations with metric tensors - might be called spiral-coordinates or 
  // something?:
  xMin = -1.5; xMax = +1.5; yMin = -1; yMax = +1;
  f = [=](Complex z) 
  { 
    double a = real(z);
    double t = imag(z);
    double x = exp(-a*t) * cos(PI*t);  // maybe use exp(-a*t) * cos(t)
    double y = exp(-a*t) * sin(PI*t);
    return Complex(x, y);
  };
  // ...maybe move this to somewhere else


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
  rsImageF contReIm(contRe);
  contReIm.blendWith(contIm, 1.f, 1.f);

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

  writeImageToFilePPM(contRe,   "ContourLinesRe.ppm");
  writeImageToFilePPM(contIm,   "ContourLinesIm.ppm");
  writeImageToFilePPM(contReIm, "ContourLinesReIm.ppm");

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

void implicitCurvesConic()
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
  writeScaledImageToFilePPM(imgCurve, "ImplicitCurvesConic.ppm", 1);

  // -with step = 4 (parameter inside the function), the ellipses are drawn heavier than with 
  //  step = 3 - why? with step = 5, they look as expected again, with 6 they get denser again
  //  (but not as dense as with 4) and here, also the circles have denser dots
  // -also, with step > larger than 1, we begin to see a doubl-drawing of the last pixel again
  //  (maybe we need to scale the weight of the final pixel by 1/step)

  // See also:
  // https://www.youtube.com/watch?v=6oMZb3yP_H8
  // https://www.youtube.com/watch?v=yOgIncKp0BE
}
// other curves to try:
// https://en.wikipedia.org/wiki/Pedal_curve
// https://en.wikipedia.org/wiki/Roulette_(curve)
// https://en.wikipedia.org/wiki/Cyclocycloid
// https://en.wikipedia.org/wiki/Epicycloid
// https://en.wikipedia.org/wiki/Hypotrochoid
// https://en.wikipedia.org/wiki/Epitrochoid


void implicitCurvesElliptic()
{
  int width  = 800;
  int height = 800;
  double range = 4.1;


  rsImageF imgCurve(width, height);
  function<double(double, double)> f;

  rsImagePlotter<float, double> ig;
  ig.setRange(-range, range, -range, range);
  ig.painter.setDeTwist(false);  // should be only used for single pixel lines
  ig.painter.setNeighbourWeightsForSimpleDot(0.375, 0.375*sqrt(0.5));

  float color = 0.375f;
  // looks good with the saturating accumulation in rsImagePainter

  f = [=](double x, double y) { return x*x*x + y*y; }; 

  ig.plotImplicitCurve(f,  -1.0, -1.0, 0.0, imgCurve, color); 
  // maybe there's a 2nd component somewhere?

  //ig.plotImplicitCurve(f, 0.0, 0.0, 0.0, imgCurve, color); // infinite loop
  ig.plotImplicitCurve(f,  +1.0, +1.0, 0.0, imgCurve, color);
  ig.plotImplicitCurve(f,  +8.0, +2.0, 0.0, imgCurve, color);
  ig.plotImplicitCurve(f, +27.0, +3.0, 0.0, imgCurve, color);
  ig.plotImplicitCurve(f, +64.0, +4.0, 0.0, imgCurve, color);
  //ig.plotImplicitCurve(f, 2.0, 1.0, 1.0, imgCurve, color);
  //ig.plotImplicitCurve(f, 5.0, 1.0, 2.0, imgCurve, color);
  // draw more ..maybe vary the contour line...yeah..maybe we should define this as a 
  // contour-plotting problem anyway...but we can do both


  using IP = rsImageProcessor<float>;
  IP::normalize(imgCurve);
  writeScaledImageToFilePPM(imgCurve, "ImplicitCurvesElliptic.ppm", 1);
}

void implicitCurves()
{
  implicitCurvesConic();
  implicitCurvesElliptic();
}

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
  //bool natural = false;  // natural parametrization (by arc-length s)

  // Lissajous curve parameters:
  double a  = 2.0;
  double b  = 3.0;
  double p  = PI;
  double t0 = 0.0;
  double t1 = 2*PI * 1.0;

  // define functions x = fx(t), y = fy(t)
  std::function<double(double)> fx, fy;

  // Lissajous figure:
  //a = 2; b = 3; t0 = 0, t1 = 2*PI;
  //fx = [&](double t) { return cos(a*t - 0.5*p); };
  //fy = [&](double t) { return sin(b*t + 0.5*p); };

  // Some nice parametric shapes: 
  a = 0.3, b = 1.5; t0 = 0, t1 = 2*PI; // a: "eggness", b: "ellipticity"
  fx = [&](double t) { return sin(t - a*sin(t)); };
  fy = [&](double t) { return b * cos(t);        };
  // https://www.desmos.com/calculator/jte3k1amlp
  // Circle:   a=0,n=1
  // Elllipse: a=0,b!=1
  // Egg:      a=0.3,b=1.5; a=0.4,b=2 (longer)
  // Drop:     a=0.85,b=2;  a=1,b=3 (longer, with cusp)
  // Fish:     a=2.3,b=1.5
  //
  // Variation: x(t) = sin(t - a*sin(t)), y(t) = b*cos(t - c*cos(t))
  // using b=1 and a=c=1.1, we get some sort of arrow-head:
  // https://www.desmos.com/calculator/v4crjw8cdo
  //
  // https://www.youtube.com/watch?v=tjyFw1BX4eM  another egg-equation


  // plot curve:
  rsImageF img(width, height);
  float color = 0.375f;
  plotParametricCurve(fx, fy, t0, t1, numDots, img, color, -range, range, -range, range);
  rsImageProcessor<float>::normalize(img);
  writeScaledImageToFilePPM(img, "ParametricCurve.ppm", 1);

  // -when using numDots such that discrete dots become visible (like 500), those that land on a 
  //  pixel appear brighter than those which don't
  //  -> using de-twisting improves this

  // Ideas for shapes:
  // -Egg: x(t) = sin(t - a*sin(t)), y(t) = b * cos(t), a: eggness, b: length, 
//    example: a=0.3,b=1.5 https://www.desmos.com/calculator/jte3k1amlp
}

template<class TPix, class TVal>
void plotFunction(const std::function<TVal(TVal)>& f, rsImage<TPix>& img, TPix color,
  TVal xMin, TVal xMax, TVal yMin, TVal yMax)
{
  // loop over the x-coordinates, compute the next y-coordinate, compute length of segment,
  // compute number of dots to use for this segment, draw the segment

}

/*
void plotFunction()
{
  // test function plotting with some example functions

}
*/


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

// move to unit tests:
bool testColrBHS()
{
  bool t = true;          // unit test result
  double r, g, b, B, H, S, r2, g2, b2;

  rsColorBHS<double> cs;  // colorspace
  //cs.setWeights(0.3, 0.5, 0.2);
  //cs.setWeights(0.3, 0.55, 0.15);
  //cs.setWeights(0.25, 0.6, 0.15);
  //cs.setWeights(0.3, 0.6, 0.1);
  cs.setWeights(0.3, 0.4, 0.3);
  double tol = 1.e-14;

  // between red and green, more red:
  r = 0.4, g = 0.3, b = 0.2;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // between red and green, more green:
  r = 0.3, g = 0.4, b = 0.2;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // between green and blue, more green:
  r = 0.2, g = 0.4, b = 0.3;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // between green and blue, more blue:
  r = 0.2, g = 0.3, b = 0.4;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // between blue and red, more blue:
  r = 0.3, g = 0.2, b = 0.4;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // between blue and red, more red:
  r = 0.4, g = 0.2, b = 0.3;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // white:
  r = 1.0, g = 1.0, b = 1.0;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // gray:
  r = 0.5, g = 0.5, b = 0.5;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // black:
  r = 0.0, g = 0.0, b = 0.0;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);

  // 50% yellow:
  r = 0.5, g = 0.5, b = 0.0;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);
  t &= rsIsCloseTo(S, 1.0,  tol);  // is considered fully saturated
  t &= rsIsCloseTo(H, 1./6, tol);  // 0 is red 1/3 is green - yellow is halfway, so it must be 1/6

  // 100% yellow:
  r = 1.0, g = 1.0, b = 0.0;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);
  t &= rsIsCloseTo(S, 1.0,  tol);
  t &= rsIsCloseTo(H, 1./6, tol);

  // 100% magenta:
  r = 1.0, g = 0.0, b = 1.0;
  cs.rgb2bhs(r, g, b, &B, &H, &S);
  cs.bhs2rgb(B, H, S, &r2, &g2, &b2);
  t &= rsIsCloseTo(r, r2, tol) && rsIsCloseTo(g, g2, tol) && rsIsCloseTo(b, b2, tol);
  t &= rsIsCloseTo(S, 1.0,  tol);
  t &= rsIsCloseTo(H, 5./6, tol);


  // todo: try some more random colors, try colors that are exactly on the border of the 
  // conditionals in bhs2rgb, test white, gray, black, pure red, green, blue, ..magenta,cyan,yellow

  // todo: fill an image with colors of a chosen brightness, hue goes into x-coordinate saturation 
  // into the y coordinate:
  B = 0.80;
  int w = 400, h = 200;
  rsImageF red(w, h), green(w, h), blue(w, h);
  for(int j = 0; j < h; j++) {
    for(int i = 0; i < w; i++) {
      double x = rsLinToLin(double(i), 0.0, w-1.0, 0.0, 1.0);
      double y = rsLinToLin(double(j), h-1.0, 0.0, 0.0, 1.0);
      H = x;
      S = y;
      cs.bhs2rgb(B, H, S, &r, &g, &b);

      //rsAssert(r <= 1 && g <= 1 && b <= 1); // this triggers: bhs2rgb does not ensure that r,g,b values are in 0..1
      r = rsMin(r, 1.0);
      g = rsMin(g, 1.0);
      b = rsMin(b, 1.0);

      red(  i, j) = float(r);
      green(i, j) = float(g);
      blue( i, j) = float(b); }}
  writeImageToFilePPM(red, green, blue, "ColorsBHS.ppm");



  rsAssert(t);
  return t;
}

void spirals()
{
  // move somewhere else:
  //plotSpiralHeightProfile();
  //testSpiralHeightProfile();
  //testImageEffectFrame(); return;
  //testDistanceMap(); return;
  //testColrBHS(); return;

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

/** Returns the number of iterations it takes for the iterative application of:
      (x,y) <- (fx(x,y), fy(x,y))
to reach point at which r^2 := x^2 + y^2 > thresh, i.e. the point at which the squared length of 
the vector (x,y) exceeds a given threshold. If the threshold is never exceeded within 
maxNumIterations, -1 is returned. The iteration starts at the initial point (x0,y0). */
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

  // ToDo: 
  // -Maybe instead of returning just the number of iterations until divergence is assured, 
  //  return some more elaborate information about the trajectory - for example: 
  //  -the total (squared?) distance taken in all the steps
  //  -the total x- and/or y-direction travelled
  //  -this data could also be used to shade the inside
  // -Maybe return a float instead of an int where the fractional part tries to somehow "refine" 
  //  the number in order to achieve continuous gradients in the resulting image. One idea could be
  //  to take the ratio of the (squared) radii in the previous and current iteration: rOld/rNew or
  //  something like (rNew-rOld) / ((rNew+rOld)/2). 
}

/** A variant that returns a floating point number whoe fractional part tries to somehow provide for
a smoother image ....tbc... */
template<class T>
float numIterationsToDivergence2(std::function<T(T, T)> fx, std::function<T(T, T)> fy, T x0, T y0,
  T thresh, int maxNumIterations)
{
  T xOld = x0;
  T yOld = y0;
  for(int i = 0; i <= maxNumIterations; i++)
  {
    float xNew  = fx(xOld, yOld);
    float yNew  = fy(xOld, yOld);
    float r2Old = xOld*xOld + yOld*yOld;
    float r2New = xNew*xNew + yNew*yNew;
    if(r2New > thresh)
    {
      //float r = r2Old / r2New;
      float r = 2.f * (r2New - r2Old) / (r2New + r2Old); // difference normalized by average

      float f = r;  // fractional part - preliminary - we may need some function f = func(r)
      //f = sqrt(r);

      f = pow(r, 8.0f);


      return (float) i + f;


      //return (float) i;  // preliminary, no fractional part yet
    }
    xOld = xNew;
    yOld = yNew;
  }
  return -1;
}


void mandelbrot(rsImage<float>& img, int maxIterations, 
  double xMin, double xMax, double yMin, double yMax)
{
  double cx, cy;
  std::function<double(double, double)> fx, fy;
  fx = [&](double x, double y) { return x*x - y*y + cx; };
  fy = [&](double x, double y) { return 2*x*y     + cy; };
  float scl = 1.f / float(maxIterations);
  for(int j = 0; j < img.getHeight(); j++) {
    for(int i = 0; i < img.getWidth(); i++) {
      cx = rsLinToLin((double)i, 0.0,  double(img.getWidth()-1), xMin, xMax);
      cy = rsLinToLin((double)j, double(img.getHeight()-1), 0.0, yMin, yMax);

      //float its = (float) numIterationsToDivergence(fx, fy, 0.0, 0.0, 4.0, maxIterations);
      float its = (float) numIterationsToDivergence2(fx, fy, 0.0, 0.0, 4.0, maxIterations);

      // simple black/white:
      //if(its >= 0.f) img(i, j) = 0.f; // iteration diverged -> pixel is outside the set
      //else           img(i, j) = 1.f; // not diverged -> pixel is inside the set

      // white inside, shaded outside:
      //if(its >= 0.f) img(i, j) = scl * float(its);
      //else           img(i, j) = 1.f;

      // test:
      img(i, j) = its;
    }}
    // maybe do: if(its >= 0) -> black, else white

}
// todo: maybe do not use the numIterationsToDivergence function - directly implement all the 
// formulas of the inner loop here - avoid overhead of using std::function...or maybe make an 
// optimized function int mandelbrot(x, y, maxNumIts)

void renderMandelbrot(int w, int h, int numIts = 300)
{
  double xMin = -1.5;
  double xMax = +0.5;
  double yMin = -1.0;
  double yMax = +1.0;

  //xMin = -2, xMax = +2, yMin = -2, yMax = +2; // wider range useful for low number of iterations

  using IP = rsImageProcessor<float>;
  rsImageF img(w, h);
  mandelbrot(img, numIts, xMin, xMax, yMin, yMax);
  IP::normalize(img);
  IP::gammaCorrection(img, 0.2f);
  // maybe apply gamma, contrast, etc.

  writeImageToFilePPM(img, "Mandelbrot.ppm");
  rsPrintLine("Done");


  // Observations:
  // -with a smaller number of iterations, some points that actually do not belong to the set are
  //  falsely considered to be within the set - so increasing the number of iterations 
  //  progressively removes spurious white* points (compare 100 vs 200 to see the effect) 
  //  (* if we draw points inside the set white and points outside the set black)
  // -when brightness is proportional to the number of iterations, points near the boundary of the
  //  set ten to get colored brighter that points further outside the set

  // -todo: implement color-maps (maybe move code from jura to rapt - at least partially)
  //  a colormap should have two template-parameters: the input value and the pixel-color - for
  //  example float and rsFloat32x2 - maybe use rsNodeBasedFunction 
}
// https://blogs.scientificamerican.com/roots-of-unity/a-few-of-my-favorite-spaces-the-mandelbrot-set/
// https://math.stackexchange.com/questions/1099/mandelbrot-like-sets-for-functions-other-than-fz-z2c
// https://math.stackexchange.com/questions/1398218/determine-coordinates-for-mandelbrot-set-zoom/1398356


// Returns the maximum lightness value of the given image, assumed to be in HSLA format. That means
// it returns the maximum value of the 3rd element of the simd-vectors.
float getMaxLightnessHSLA(const rsImage<rsFloat32x4>& img)
{
  float maxL = RS_MIN(float);
  rsFloat32x4 *p = img.getPixelPointer(0, 0);
  for(int i = 0; i < img.getNumPixels(); i++)
  {
    float L = p[i][2];    // L-value sits in 3rd element (index 2)
    if(L > maxL)
      maxL = L;
  }
  return maxL;
}
// hmm...that's inconvenient. maybe we should split a HSLA image into 4 seperate images for the 4 
// channels, process them seperately, then recombine

void splitChannels(const rsImage<rsFloat32x4>& img, 
  rsImage<float>& ch1, rsImage<float>& ch2, rsImage<float>& ch3, rsImage<float>& ch4)
{
  int w = img.getWidth();
  int h = img.getHeight();
  rsAssert(ch1.hasShape(w, h));
  rsAssert(ch2.hasShape(w, h));
  rsAssert(ch3.hasShape(w, h));
  rsAssert(ch4.hasShape(w, h));
  rsFloat32x4 *pIn = img.getPixelPointer(0, 0);
  float *p1 = ch1.getPixelPointer(0, 0);
  float *p2 = ch2.getPixelPointer(0, 0);
  float *p3 = ch3.getPixelPointer(0, 0);
  float *p4 = ch4.getPixelPointer(0, 0);
  for(int i = 0; i < w*h; i++)
  {
    rsFloat32x4 c = pIn[i];  // color
    p1[i] = c[0];            // 1st channel
    p2[i] = c[1];            // 2nd channel
    p3[i] = c[2];            // 3rd channel
    p4[i] = c[3];            // 4th channel
  }
  int dummy = 0;
};
// move into Drawing.h/cpp. if this goes into the rapt library, it should not operate directly on
// the float/rsFloat32x4 types but rather on some type T and rsSimdVector<T, 4>

// todo: combineChannels or merge channels


void renderNewtonFractal()
{
  // User parameters:
  //double xMin   = -2.0;
  //double xMax   = +2.0;
  //double yMin   = -2.0;
  //double yMax   = +2.0;
  double xMin   = -1.9;
  double xMax   = -0.8;
  double yMin   = +0.8;
  double yMax   = +1.9;
  int    w      = 1600;      // image width in pixels
  int    h      = 1600;      // image height
  int    maxIts =  100;      // maximum number of iterations
  double tol    =  1.e-14;   // tolerance in convergence test
  int    smooth =  25;       // number of smoothing passes via gradientify
  // It's actually ok to use a larger number of smoothing passes because the number of iterations 
  // per pass tends to drop a lot. The first pass tends to take around 100 iterations, then second 
  // only around 30, etc. Around pass 20, it's only 3 or 2 iterations per pass.


  //w = 800, h = 800;

  //w = 240; h = 135;   // FHD/8
  //w = 480; h = 270;   // FHD/4
  //w = 960; h = 540; // FHD/2
  // Uncomment for high-quality rendering (takes long - use release build!)
  //w = 1920; h = 1080;   // FHD
  w = 3840; h = 2160; // UHD = 2 * FHD, 4K
  //w = 7860; h = 4320;   // 2 * UHD = 4 * FHD, 8K
  //w *= 3;    h *= 3;    // ..with 3x oversampling

  using Complex = std::complex<double>;
  using Vec2D   = RAPT::rsVector2D<double>;
  using Color   = rsFloat32x4;

  // Define iteration function. The function results from applying the Newton iteration rule 
  // zNew = z - f(z)/f'(z) to the function f(z) = z^4 - 1:
  auto iterFunc = [](Vec2D v, Vec2D p)
  {
    Complex z(v.x, v.y);
    Complex z2 = z*z;

    // Literal Newton iteration formula (uses 1 div, 3 mul, 2 sub):
    Complex f  = z2*z2 - 1.0;             // f(z)  = z^4 - 1
    Complex fp = 4.0 * z2*z;              // f'(z) = 4 * z^3

    //// test with 8 roots - comment out, if you want 4 roots
    //Complex z4 = z2*z2;
    //f  = z4*z4 - 1.0;                     // f(z)  = z^8 - 1
    //fp = 8.0 * z4*z2*z;                   // f'(z) = 8 * z^7
    // no - this doesn't work because the roots array will not match the function. todo: write
    // a function that produces the iteration function from a roots array, then create the roots
    // array programmatically.

    Complex zn = z - f/fp;                // zNew  = z - f(x) / f'(z)
    return Vec2D(zn.real(), zn.imag());

    // Algebraically equivalent (uses 1 div, 4 mul, 1 add):
    //Complex w = (3.0*z2*z2+1.0) / (4.0*z2*z); // (3 z^4 - 1) / (4 z^3)
    //return Vec2D(w.real(), w.imag());

    // I think, for f(z) = z^n - 1, we would get zNew = ((n-1) z^n + 1) / (n z^(n-1)) 
    // -> verify! If correct, implement a function that creates such iterFuncs
  };

  // Define stopping criterion:
  auto stopCriterion = [&](const std::vector<Vec2D>& t)
  {
    size_t N = t.size();
    if(N < 2)
      return false;
    return rsNorm(t[N-1] - t[N-2]) <= tol; // todo: use squared norm (and squared tolerance)
  };



  // Helper function for the coloring algorithms. Returns true, iff a is strictly closer to the
  // reference value r than n:
  auto closer = [](Vec2D a, Vec2D b, Vec2D r)
  {
    double da = rsNorm(r-a);  // todo: use squared norms to get rid of the sqrt
    double db = rsNorm(r-b);
    return da < db;
  };
  std::vector<Vec2D> roots( { Vec2D(1,0), Vec2D(0,1), Vec2D(-1,0), Vec2D(0,-1) });

  // Define coloring function. The parameter t is the trajectory. Here, we use a simple grayscale
  // coloring, target root determines gray value. The coloring is determined by a pair of function,
  // one that processes that tarjectory to produce a 4-float value at each pixel that is called 
  // back for each pixel after convegence with the tarjectory as argument and one post-processing
  // function that takes an image of these previously produced 4-float values and translates them
  // into actual RGBA values. The split is necessary to allow a post-processor to access global
  // features such as the total maximum (over all pixels) of iterations taken in order to normalize
  // things, etc. So, a colorFunc must always be paired with appropriate post-processing function 
  // that knows, how to intepret the values produced by the colorFunc in terms of RGBA values:
  auto colorFunc1 = [&](const std::vector<Vec2D>& t)
  {
    Vec2D vL = rsLast(t);
    int k = findBestMatch(&roots[0], (int) roots.size(), vL, closer);
    return  float(k) / (roots.size()-1);
  };
  auto postProcess1 = [&](rsImage<Color>& img)
  {
    // Do nothing. The colorFunc1 itself already produces the final RGBA values.
  };

  // Another coloring strategy using the target root for the hue and the number of iterations for
  // the lightness, saturation is fixed:
  auto colorFunc2 = [&](const std::vector<Vec2D>& t)
  {
    Vec2D vL = rsLast(t);
    int k = findBestMatch(&roots[0], (int) roots.size(), vL, closer);
    float rootIdx = float(k);           // index of root to which it converged
    float numIts  = float(t.size());    // number of iterations taken
    return rsFloat32x4(rootIdx, numIts, 0.f, 0.f);
  };
  auto colorFunc2_2 = [&](const std::vector<Vec2D>& t)
  {
    // doesn't work yet

    Vec2D vL = rsLast(t);
    int k = findBestMatch(&roots[0], (int) roots.size(), vL, closer);
    float rootIdx = float(k);           // index of root to which it converged
    float distSum = 0.f;                // sum of the distances
    Vec2D root = roots[k];
    for(int i = 0; i < (int)t.size(); i++)
      distSum += (float) rsNorm(root - t[i]);

    // maybe normalize the distSum by the initial distance and/or t.size()
    distSum /= t.size();  // test

    if(t.size() == maxIts)              // set the outliers with huges distance to zero
      distSum = 0.f;

    return rsFloat32x4(rootIdx, distSum, 0.f, 0.f);
  };
  auto postProcess2 = [&](rsImage<Color>& img)
  {
    using AT = RAPT::rsArrayTools;
    using IP = RAPT::rsImageProcessor<float>;
    int w = img.getWidth();
    int h = img.getHeight();

    // Extract root index and number of iterations into a,b (c,d ar unused dummies):
    rsImage<float> a(w,h), b(w,h), c(w,h), d(w,h); 
    splitChannels(img, a, b, c, d);
    // a,b,c,d are genric "color" or just data channels up to our interpretation. We assume that
    // colorFunc2 was used, so a is the root index and b the number of iterations taken.

    // Process numIterations, to be used for lightness L in HSL:
    float *pb = b.getPixelPointer(0, 0);  // pb: pointer to the lightness channel
    AT::normalize(pb, w*h);
    for(int i = 0; i < w*h; i++) {        // ad hoc to remove ugly white diagonal line, replace
      if(pb[i] == 1.f)                    // it by less obstrusive black line
        pb[i] = 0.f; }
    //rsPlotArray(pb, w*h);  // for development
    if(smooth != 0) {
      d.copyPixelDataFrom(b);
      gradientifyFlatRegions(d, b, smooth); }
    IP::gammaCorrection(b, 0.6f);

    // Process root index, to be used for hue H in HSL:
    float hueScaler = 1.f / roots.size();
    AT::scale(a.getPixelPointer(0, 0), w*h, hueScaler);

    // Now the b-image represents our desired L (lightness) values and the a-image contains our
    // desired H (hue) values. Now, we compute the corresponding RGB colors and write them back 
    // into the img as rsFloat32x4:
    float S  = 1.f;                             // saturation is a fixed value
    float A  = 1.f;                             // opacity
    float *H = a.getPixelPointer(0, 0);         // pointer to hues
    float *L = b.getPixelPointer(0, 0);         // pointer to lightnesses
    rsFloat32x4* p = img.getPixelPointer(0, 0);
    for(int i = 0; i < w*h; i++)
    {
      float R, G, B;
      rsColor<float>::hsl2rgb(H[i], S, L[i], &R, &G, &B); 
      p[i] = rsFloat32x4(R, G, B, A);
    }

    // todo: 
    // -we get w completely white diagonal line
    //  -maybe when it didn't converge, we should color the pixel black
    // -maybe form a histogram of the lighness values and from that determine a nonlinear 
    //  mapping function to be applied to the lightness values, maybe it shou duse a smoothed
    //  histogram
    // -increasing contrast and gamma in irfanview seems to help a bit 
    // -is there a way to get rid of the steps by somehow determining a "fractional part" of the 
    //  number of iterations - maybe looking at the last two iterates?
    // -To avoid the stairstep like nature of the lightness, try to find a continuous function that
    //  is similar to the number of iterations. Maybe the sum of distances to the approached root, 
    //  maybe divided by the initial distance and/or with initial distance subtracted off. We want 
    //  smooth color gradients!
    // -could it make sense to use a smooth hue function, too? how could that look like?
    // -or maybe we should desaturate the colors near the boundaries, rationale: no hard steps
    //  at boudaries...but maybe we want hard steps
    // -try a higher order polynomial - we want more hues to get a nice rainbow look
    // -general algo for de-flattening regions with flat color in any image:
    //  -identify flat color regions by finding pixels whose color is equal to the color of all its
    //   neighbors
    //  -for each pixel that belogns to such a flat-color region, assign a new color based on 
    //   linearly interpolating between colors that are adjacent to to the flat region. Let's say
    //   the pixel has coordinates (200, 300) and the flat region around the pixel extends 20 
    //   pixels to the left, 30 to the right, 40 downward and 50 upward, use a bilinear mix between
    //   the colors of the pixels at these 4 positions

    // -Smoothing seems to affect only the outer flat regions, at least with a low-res rendering
    //  of 100x100 - the inner flat regions get misclassified as "rest"
    // -Test rendering at 960x540 and look at the produced AfterStep3.ppm generated by
    //  gradientifyFlatRegions. There are still more or less flat regions. Bumping up our smooth
    //  parameter here doesn't seem to help. I think the grandientify method does not yet work
    //  as intended -> investigate the method experimentally more thoroughly

    // Ideas:
    // 
  
    return;
  };



  // Set up the renderer:
  rsFractalImageRenderer renderer;
  renderer.setIterationFunction(iterFunc);

  //renderer.setColoringFunction(colorFunc1);
  //renderer.setPostProcessingFunction(postProcess1);

  renderer.setColoringFunction(colorFunc2);
  //renderer.setColoringFunction(colorFunc2_2);
  renderer.setPostProcessingFunction(postProcess2);

  renderer.setStoppingCriterion(stopCriterion);
  renderer.setCoordinateRange(xMin, xMax, yMin, yMax);
  renderer.setMaxNumIterations(maxIts);
  renderer.setImageSize(w, h);

  // Render the fractal and post-process the returned raw image by translating the colors to 
  // 24-Bit RGB, etc. and write the result into a .ppm file:
  rsImage<Color> imgRaw = renderer.render();

  // ToDo:
  // -implement oversampling
  // -apply some post-processing, maybe this should be also callback based, i.e. we define a
  //  postProcess function that takes a raw-image as input and produces a final RGB image
  // -maybe symmetrize the image via averaging or crossfading with flipped version
  // -implement better color-spaces and use them in the coloring algo

  rsImage<rsPixelRGB> img = rsConvertImage(imgRaw, true);
  writeImageToFilePPM(img, "NewtonFractalDeg4.ppm");
  rsPrintLine("Newton fractal done");


  // Observations:
  // -When rendering at 1600x1600 and looking at the bottom-left of the PixelClasses.ppm (zoomed 
  //  in), we see some "holes" in the boundary pixel class, i.e. pixels are classified into the 
  //  black ("rest") class rather than in the lightgray ("boundary") class. That's weird! Figure
  //  out, what's going on!
  // -When rendering at really high resolution like UHD with 3x oversampling, the "gradientify"
  //  algorithm seems to become less effective. Although the sharp boundaries between the flat 
  //  regions are smoothed out, some amount of flatness remains. Maybe the absolute size of the 
  //  flat regions (in pixels) determines, how effective gradientification is? Increasing the 
  //  "smooth" variable here and the "maxIts" variable in gradientifyFlatRegions doesn't seem
  //  to help against the flatness either, so it doesn't seem to be an issue of stopping before
  //  convergence - the algo seems to actually converge to that final state.

  // Ideas:
  // -use the number of iteration for L and the root for H, see the frcatals here:
  //    https://www.youtube.com/watch?v=blOARV4lnIM
  //  the white regions are regions of slow convergence
  // -make a composited image taken from different regions of the same fractal, maybe use the same
  //  center but slightly different zoom levels, or shift it a bit left/right/up/down
  //  ...maybe this can be done as post-processing
  // -maybe render raw image with some margins to facilitate post-processing without problems at 
  //  the edges
  // -try it with a higher order polynomial and the permute the roots-array in some way

  // ToDo:
  // -when rendering the fractal at very high resolution, the AfterStep3.ppm and also the final 
  //  result again show flat-looking regions. Check, how many iterations are taken by the 
  //  gradientify algorithm and whether we may stop iterating too early, i.e. before everything 
  //  looks nice.
  // -The AfterSetp3b is useful: reduce green by 80, increase blue by 20 in IrfanViews color 
  //  corrections -> looks good - maybe use as wallpaper

}

void fractal()
{
  //renderNewtonFractal();
  //renderMandelbrot(500, 500);
  //renderMandelbrot(2000, 2000);
  renderMandelbrot(5000, 5000);
}

void parametricCurve2D()
{
  

  // todo: turn into unit test

  // some shorthands:
  using VecN     = std::vector<double>;          // nD vectors
  using Vec2     = rsVector2D<double>;           // 2D vectors
  using Func_1_2 = std::function<Vec2(double)>;  // functions from 1D scalars to 2D vectors
  using Curve2D  = rsParametricCurve2D<double>;  // class to represent 2D curves

  // create and set up the curve object:
  Curve2D c;
  Func_1_2 r = [&](double t) { return Vec2(cos(2*PI*t), sin(2*PI*t)); }; 
  c.setPositionFunction(r);
    // The curve r(t) = (cos(2*pi*t), sin(2*pi*t)) defines a unit circle. This is not a natural 
    // parametrization (a.k.a. parametrization by arc length): the circle is not traversed at unit 
    // speed but at speed 2*pi, so it makes a full revolution, when the parameter t traverses an 
    // interval of unit length, for example from 0 to 1.

  // compute velocities at t = 0, 0.25, 0.5, 0.75:
  double h;   // approximation stepsize for numeric derivatives
  Vec2 v;     // velocity
  h = 1.e-4;  // not yet verified, if this value is best
  v = c.getVelocity(0.0,  h);  // ( 0,     2*pi)
  v = c.getVelocity(0.25, h);  // (-2*pi,  0   )
  v = c.getVelocity(0.5,  h);  // ( 0,    -2*pi)
  v = c.getVelocity(0.75, h);  // ( 2*pi,  0   )

  // compute accelerations:
  Vec2 a;      // acceleration
  h = 1.e-4;
  a = c.getAcceleration(0.0,  h);  // (-4*pi^2,  0     )
  a = c.getAcceleration(0.25, h);  // ( 0,      -4*pi^2)
  a = c.getAcceleration(0.5,  h);  // ( 4*pi^2,  0     )
  a = c.getAcceleration(0.75, h);  // ( 0,       4*pi^2)
  // 4*p^2 = 39.4784176043574
  // for the acceleration, we obtain the best numerical approximation of 7 correct decimal digits
  // with = 1.e-4:
  //   h = 1.e-N, N:   3  4  5  6  7
  //   correct digits: 5  7  7  5  4
  // with 1.e-5, we also get seven correct digits, but the error is a bit worse
  // todo: figure out, if 1.e-4 is also the optimal value for computing the velocities - and how
  // the optimal choice depends on the particular curve

  // sage code to compute velocities and accelerations analytically:
  //  var("t")
  //  r = vector((cos(2*pi*t), sin(2*pi*t))) # position as function of t
  //  v = r.diff(t)                          # velocity
  //  a = v.diff(t)                          # acceleration
  //  v(0),v(1/4),v(1/2),v(3/4),\
  //  a(0),a(1/4),a(1/2),a(3/4)              # print

  double k;  // curvature "kappa"
  k = c.getCurvature(0.0,  h);
  k = c.getCurvature(0.25, h);
  k = c.getCurvature(0.5,  h);
  k = c.getCurvature(0.75, h);
  // todo: figure out optimal value of h, such that the numerical error is minimal


  // re-use v for a generic vector - now we compute the center of the osculating circle
  v = c.getOsculatingCircleCenter(0.0,  h);
  v = c.getOsculatingCircleCenter(0.25, h);
  v = c.getOsculatingCircleCenter(0.5,  h);
  v = c.getOsculatingCircleCenter(0.75, h);
  // should always be 0 - seems to work

  // test arc-length computation:

  VecN t = rsRangeLinear(0.0, 1.0, 65);
  VecN s = c.getArcLengthFunction(t);
  // ok - works in principle, but we should improve accuracy by using trapezoidal integration - 
  // do this in class rsNumericIntegrator - generally move the free functions for numeric calculus
  // into these classes - maybe also drag in the ODE solver from the GNUPlotCPP repo

  VecN t2 = c.getArcLengthParametrization(0.0, 1.0, (int) t.size());

  //rsPlotVectorsXY(t, t2); // should be the identity - and indeed, it is! :-)
  // try something where the curve is traversed by changing speed such as
  // cos(tau*t^2)

  // as examples for 3D curves, draw helix, trefoil knot, 3D Lissaous
  // plot some curves and their evolutes

  int dummy = 0;
}

bool testRotation3D()
{
  using GT = rsGeometricTransforms<double>;

  double A[3][3];  // transformation matrix
  rsVector3D<double> u, v, w;
  bool test = true;

  // test to rotate all 3 unit vectors into all others:
  u.set(1,0,0); v.set(0,1,0); GT::rotationMatrixFromTo(u, v, A); w = A*u; test &= w == v;
  u.set(1,0,0); v.set(0,0,1); GT::rotationMatrixFromTo(u, v, A); w = A*u; test &= w == v;
  u.set(0,1,0); v.set(1,0,0); GT::rotationMatrixFromTo(u, v, A); w = A*u; test &= w == v;
  u.set(0,1,0); v.set(0,0,1); GT::rotationMatrixFromTo(u, v, A); w = A*u; test &= w == v;
  u.set(0,0,1); v.set(1,0,0); GT::rotationMatrixFromTo(u, v, A); w = A*u; test &= w == v;
  u.set(0,0,1); v.set(0,1,0); GT::rotationMatrixFromTo(u, v, A); w = A*u; test &= w == v;

  rsAssert(test);
  return test;
}
// todo: plot a 3D curve - we somehow need to apply a transformation to the 3D vector such that 
// after the transformation, we can just discard the z-coordinate


/*

// y = A*x
template<class T>
void applyMatrix4D(const T A[4][4], const T x[4], T y[4])
{
  rsAssert(x != y, "Cannot be used in place" );
  y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2] + A[0][3]*x[3];
  y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2] + A[1][3]*x[3];
  y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2] + A[2][3]*x[3];
  y[3] = A[3][0]*x[0] + A[3][1]*x[1] + A[3][2]*x[2] + A[3][3]*x[3];
}

// C = A*B
template<class T>
void multiplyMatrices4D(const T A[4][4], const T B[4][4], T C[4][4])
{
  rsAssert(A != C && B != C, "Cannot be used in place" );
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      C[i][j] = T(0);
      for(int k = 0; k < 4; k++)
        C[i][j] += A[i][k] * B[k][j]; }}
}

template<class T>
rsVector2D<T> project(const rsVector3D<T> v, const T A[4][4])
{
  T t1[4], t2[4];  // temporary homogeneous 4D vectors

  // copy 3D input vector into homogeneous 4D vector:
  t1[0] = v.x;   // x
  t1[1] = v.y;   // y
  t1[2] = v.z;   // z
  t1[3] = T(1);  // w = 1

  // apply homogeneous 4x4 matrix:
  applyMatrix4D(A, t1, t2);

  // extract x,y components and put into rsVector2D:
  return rsVector2D<T>(t2[0], t2[1]);
}
// this is not optimized code!

*/

template<class TPix, class TVal>
void plotParametricCurve3D(rsImage<TPix>& img, TPix color, 
  const rsParametricCurve3D<TVal>& crv, const TVal A[4][4], TVal t0, TVal t1, int N)
{
  rsImagePainter<TPix, TVal, TVal> painter(&img);
  painter.setDeTwist(true);

  TVal xMaxPixel = TVal(img.getWidth()  - 1);   // maximum x-coordinate in pixel coordinates
  TVal yMaxPixel = TVal(img.getHeight() - 1);   // same for y-coordinate

  std::vector<TVal> t = rsRangeLinear(t0, t1, N);
    // todo: optionally re-map t to s to get an arc-length parametrization
  for(int n = 0; n < N; n++) {
    rsVector3D<TVal> v3 = crv.getPosition(t[n]);
    rsVector2D<TVal> v2 = project(v3, A);

    // todo: transform to pixel coordinates:
    TVal px = rsLinToLin(v2.x, -1.0, +1.0, TVal(0), xMaxPixel);
    TVal py = rsLinToLin(v2.y, -1.0, +1.0, yMaxPixel, TVal(0));
    // -1..+1 comes from the OpenGL convention of "normalized device coordinates"
    // https://learnopengl.com/Getting-started/Coordinate-Systems
    // https://computergraphics.stackexchange.com/questions/1769/world-coordinates-normalised-device-coordinates-and-device-coordinates

    // maybe use z-coordinate to sclae thickness and/or brightness to give some sense of 3D depth


    painter.paintDot(px, py, color); }


  /*
  // this code would have to be used when we use a line-drawer:
  // rename this function to plot..ViaDots, have a corresponding ...ViaLines function
  rsVector3D<TVal> v3 = crv.getPosition(t[0]);
  rsVector2D<TVal> start = project(v3, A), end;
  for(int n = 1; n < N; n++)
  {
    v3    = crv.getPosition(t[n]);
    end   = project(v3, A);
    start = end;
    drawer.drawLine(start, end);
  }
  */


}
// -maybe t0,t1 should be members of the rsParametricCurve class
// -maybe have an option for re-parametrization by arc-length
// -maybe optionally draw velocity and acceleration vectors...but not - that should be a separate
//  function - and these vectors should be drawn only at some selected values of t

#undef near // some silly header defines these as macros
#undef far
void plotCurve3D()
{
  // tests 3D curve plotting

  // some shortcuts (maybe get rid, if they are used only once):
  using Vec3     = rsVector3D<double>;           // 3D vectors
  using Func_1_3 = std::function<Vec3(double)>;  // functions from 1D scalars to 3D vectors
  using Curve3D  = rsParametricCurve3D<double>;  // class to represent 3D curves


  // image parameters:
  int w = 400;  // width
  int h = 400;  // height

  // drawing parameters:
  float color = 1.f;  // brightness
  int N = 1600;        // number of dots

  // perspective parameters:
  double left   = -2.0; 
  double right  = +2.0; 
  double bottom = -2.0; 
  double top    = +2.0; 
  double near   = -2.0; 
  double far    = +2.0;
  // how are they interpreted? can we find a more convenient parametrization? maybe via 3 vectors
  // (eye, center, up) such as in OpenGL vmath::lookat, vmath::ortho (see OpenGL Prog. Guide, 
  // pg 220)
  // it seems, they do not contain any information
  // see https://bartipan.net/vmath/doc/class_matrix4.html#a0b8035f3d1144444d6835cd60642009d
  // https://bartipan.net/vmath/

  // maybe also look at this:
  // http://glm.g-truc.net/0.9.5/index.html
  // because here:
  // https://stackoverflow.com/questions/17250146/what-and-where-is-vmath-h
  // someone says, vmath.h is buggy

  // alternative perspective parameters (maybe make a boolean switch, which are used):
  Vec3 eye(   1, 1, 0.75);
  Vec3 center(0, 0, 0);
  Vec3 up(    0, 0, 1);

  // default perspective (i think) - look down on the xy plane from z=1 with y-axis being the up
  // direction:
  //Vec3 eye(   0, 0, 1);
  //Vec3 center(0, 0, 0);
  //Vec3 up(    0, 1, 0);
  // todo: figure out, if the function lookAtMatrix4D does what it's supposed to do - darw a simple
  // object for which it is obvious, how it should look like form varius points and see, if it 
  // looks as expected

  double zoom = 0.375;

  // range of curve parameter t:
  double t0 = -2.0;
  double t1 = +2.0;



  // create and set up the curve object:
  Curve3D crv;
  Func_1_3 r;


  // Helix:
  zoom = 0.375; t0 = -2, t1 = 2; N = 2400;
  r = [&](double t) { return Vec3(cos(2*PI*t), sin(2*PI*t), t); }; 

  // Lissajous:
  //zoom = 0.6; t0 = 0; t1 = 1; N = 3200;
  //r = [&](double t) { return Vec3(cos(2*2*PI*t), sin(3*2*PI*t), cos(5*2*PI*t)); };

  // Trefoil knot:
  //zoom = 0.25; t0 = 0, t1 = 2*PI; N = 2400;
  //r = [&](double t) { return Vec3(sin(t)+2*sin(2*t), cos(t)-2*cos(2*t), -sin(3*t)); }; 
  // need other perspective

  // todo: define a Hilbert curve as rsCurve2D object
  // rsVector2D<T> hilbertPolygon(T t, int order); 
  // t goes from 0 to 1, function should return point on Hilbert polygon of given order for given t



  // i think, we should really use the saturating mode

  crv.setPositionFunction(r);

  // create image and transformation matrix:
  rsImageF img(w, h);
  double A[4][4];
  //rsGeometricTransforms<double>::orthographicProjection(A, left, right, bottom, top, near, far);
  lookAtMatrix4D(A, eye, center, up, zoom); // maybe have an additional "zoom" parameter (scalar)
  
  // plot the curve and save result to image file:
  plotParametricCurve3D(img, color, crv, A, t0, t1, N);
  writeImageToFilePPM(img, "Curve3D.ppm");
}

void differentialGeometry()
{
  //testRotation3D();
  //parametricCurve2D();
  plotCurve3D();
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


// color spaces
// https://en.wikipedia.org/wiki/HSL_and_HSV
// https://en.wikipedia.org/wiki/HSL_and_HSV#Disadvantages

// https://en.wikipedia.org/wiki/CIELAB_color_space
// https://en.wikipedia.org/wiki/SRGB
// https://stackoverflow.com/questions/7880264/convert-lab-color-to-rgb
// http://ai.stanford.edu/~ruzon/software/rgblab.html

// how about: 
// saturation =  (max(r,g,b) - min(r,g,b)) / mean(r,g,b)
//            or (max(r,g,b) - min(r,g,b)) / mid( r,g,b)
//            or (max-min) / max
// the last has the property, that saturation is 1 whenever min=0 - this seems to make sense.
// brightness = wr*r + wg*g + wb*b where wr + wg + wb = 1
// hue        = s = min(r,g,b), rs = r-s, gs = g-s, bs= b-s - shifted rgb values
//              if(    bs == 0): between red and green
//              elseif(gs == 0): between red and blue
//              else(->rs == 0): betwen green and blue
// ...is this invertible?
// if hue between red,green -> blue=min etc. - the hue is used to reconstruct which of the r,g,b 
// values are min, max