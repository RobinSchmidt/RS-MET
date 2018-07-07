#include "GraphicsExperiments.h"
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
  static const int N = 80;  // number of data points
  float a = 2.f;
  float b = 7.f;
  float scale = 0.9f;

  // create and set up rsPhaseScopeBuffer object:
  rsPhaseScopeBufferFFD psb;
  psb.setAntiAlias(true);
  psb.setBrightness(100.f);
  psb.setLineDensity(1.f);
  psb.setPixelSpread(0.3f);
  psb.setUseColorGradient(true);
  psb.setSize(400, 400);
  psb.reset();

  // create image:
  float x[N], y[N];
  float s = float(2*PI) / (N-1);
  for(int n = 0; n < N; n++)
  {
    x[n] = scale*sin(s*a*n);
    y[n] = scale*sin(s*b*n);
    psb.processSampleFrame(x[n], y[n]);
  }

  // retrieve and save image:
  psb.getImage();
  writeImageToFilePPM(*psb.getImage(), "PhaseScopeLissajous.ppm");

  //// plot (for reference):
  //GNUPlotter plt;
  //plt.addDataArrays(N, x, y);
  //plt.plot();
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
  //pDrawTriangle = &drawTriangle;
  pDrawTriangle = &drawTriangleAntiAliased;




  /*

  //void drawTriangle(rsImageDrawerFFF& drw, 
  //  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color);

  // parallelogram made from a flat-top and flat-bottom triangle:
  Vec2 
    p1(  2.f,  1.f), 
    p2(  7.f,  1.f), 
    p3(  4.f,  5.f), 
    p4(  9.f,  5.f);
  pDrawTriangle(drw, p1, p2, p3, c);  // flat-top
  pDrawTriangle(drw, p2, p4, p3, c);  // flat-bottom
  //drawTriangleFlatTop(   drw, p1, p2, p3, c);  // p1,p2 above p3, clockwise
  //drawTriangleFlatBottom(drw, p2, p3, p4, c);  // p2 above p3,p4, counterclockwise

  // polygon from 3 general triangles:
  Vec2
    p5(  6.f,  7.f),
    p6( 15.f,  8.f),
    p7(  9.f, 10.f),
    p8(  3.f, 11.f),
    p9( 14.f, 13.f);
  pDrawTriangle(drw, p5, p6, p7, c);
  pDrawTriangle(drw, p6, p9, p7, c);
  pDrawTriangle(drw, p5, p7, p8, c);

  // a general triangle:
  Vec2
    p10(25.f,  1.f),
    p11(31.f, 10.f),
    p12(17.f,  5.f);
  pDrawTriangle(drw, p10, p11, p12, c);
  */

  //pDrawTriangle(drw, V(,), V(,), V(,), c);

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

  //pDrawTriangle(drw, V(,), V(,), V(,), c);  


  // save to file:
  writeImageToFilePPM(img, "PolygonsViaTriangles.ppm");
}

// todo: try anti-aliasing by computing overlap-polygons by clipping the triangle to pixels
// see Foley, page 124ff - iterate over all pixels in the bounding box of the triangle, for each
// perform sutherlnad/hodgeman clipping, compute polygon area and fill the pixel - ultra-slow
// but good for a prototype