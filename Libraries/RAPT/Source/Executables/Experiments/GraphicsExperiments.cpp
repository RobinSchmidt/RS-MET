#include "GraphicsExperiments.h"

// Sources:
// http://hinjang.com/articles/01.html
// https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm
// https://github.com/tobsa/Xiaolin-Wus-Line-Algorithm/blob/master/Xiaolin%20Wu's%20Line%20Algorithm/Source/XiaolinWusLineAlgorithm.cpp

// This is the Wu line drawing algorithm, directly translated from the pseudocode here:
// https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm
// todo: move to a Prototypes.h/cpp pair of files in the Shared folder
inline int   ipart(float x)           { return (int) x;              }
inline int   roundToInt(float x)      { return ipart(x + 0.5f);      }
inline float fpart(float x)           { return x - ipart(x);         }
inline float rfpart(float x)          { return 1 - fpart(x);         }
inline void  swap(float& x, float& y) { float t = x; x = y; y = t;   }
//inline void  plot(ImageF& im, int x, int y, float c){ im(x, y) += c; }
inline float min(float x, float y)    { return x < y ? x : y; }
inline void  plot(ImageF& im, int x, int y, float c){ im(x,y) = min(1.f, im(x,y)+c); }
void drawLineWu(ImageF& img, float x0, float y0, float x1, float y1, float color)
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

// Bresenham line drawing algorithm:
void drawLineBresenham(ImageF& img, int x1, int y1, int x2, int y2, float color)
{
  bool steep = abs(y2 - y1) > abs(x2 - x1);

  if(steep){
    swap(x1, y1);
    swap(x2, y2); }
  if(x1 > x2){
    swap(x1, x2);
    swap(y1, y2); }

  int deltaX = x2 - x1;
  int deltaY = abs(y2 - y1);
  int error  = deltaX / 2;
  int ystep;
  int y = y1;

  if(y1 < y2)
    ystep = 1;
  else
    ystep = -1;

  for(int x = x1; x <= x2; x++){
    if(steep)
      plot(img, y, x, color);
    else
      plot(img, x, y, color);
    error -= deltaY;
    if(error < 0){
      y += ystep;
      error += deltaX; }}
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
  ImageF image(imageWidth, imageHeight);
  ImagePainterFFF painter(&image, nullptr);
  float x1 = margin;
  float y1 = margin;
  int i;

  // create arrays for line endpoints:
  vector<float> x2, y2;
  for(i = 0; i < numLines; i++){    // flat, horizontal'ish
    x2.push_back(imageWidth - margin);
    y2.push_back(margin + i * (imageHeight - margin) / numLines);
  }
  x2.push_back(imageWidth -margin); // 45° diagonal
  y2.push_back(imageHeight-margin);
  for(i = 0; i < numLines; i++){    // steep, vertical'ish
    x2.push_back(margin + i * (imageWidth - margin) / numLines);
    y2.push_back(imageHeight - margin);
  }

  // dotted algorithm:
  for(i = 0; i < x2.size(); i++)
    painter.drawDottedLine(x1, y1, x2[i], y2[i], brightness);
  writeImageToFilePPM(image, "LinesDotted.ppm");

  // Wu algorithm:
  image.clear();
  for(i = 0; i < x2.size(); i++)
    drawLineWu(image, x1, y1, x2[i], y2[i], brightness);
  writeImageToFilePPM(image, "LinesWu.ppm");

  // Bresenham algorithm:
  image.clear();
  for(i = 0; i < x2.size(); i++)
    drawLineBresenham(image, roundToInt(x1), roundToInt(y1), 
      roundToInt(x2[i]), roundToInt(y2[i]), brightness);
  writeImageToFilePPM(image, "LinesBresenham.ppm");
}

