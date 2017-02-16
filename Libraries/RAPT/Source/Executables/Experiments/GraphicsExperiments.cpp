#include "GraphicsExperiments.h"

void lineDrawing()
{
  // Compares different line drawing algorithms. We draw lines of different directions.

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
  float x2, y2;
  int i;

  // draw lines using the algorithm that creates them from dots:
  for(i = 0; i < numLines; i++){ // horizontal'ish
    x2 = imageWidth - margin;
    y2 = margin + i * (imageHeight - margin) / numLines;
    painter.drawDottedLine(x1, y1, x2, y2, brightness);
  }
  for(i = 0; i < numLines; i++){ // vertical'ish
    y2 = imageHeight - margin;
    x2 = margin + i * (imageWidth - margin) / numLines;
    painter.drawDottedLine(x1, y1, x2, y2, brightness);
  }
  painter.drawDottedLine(x1, y1, imageWidth-margin, imageHeight-margin, brightness);
  writeImageToFilePPM(image, "DottedLines.ppm");


  int dummy = 0; 
}

