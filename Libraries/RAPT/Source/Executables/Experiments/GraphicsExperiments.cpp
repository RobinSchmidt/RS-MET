#include "GraphicsExperiments.h"

void drawLineWu(ImageF& image, float x1, float y1, float x2, float y2, float color)
{
  // add Wu algorithm from wikipedia here
}

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
  //float x2, y2;
  int i;

  // create arrays for line endpoints:
  vector<float> xE, yE;
  for(i = 0; i < numLines; i++){    // flat, horizontal'ish
    xE.push_back(imageWidth - margin);
    yE.push_back(margin + i * (imageHeight - margin) / numLines);
  }
  xE.push_back(imageWidth -margin); // 45° diagonal
  yE.push_back(imageHeight-margin);
  for(i = 0; i < numLines; i++){    // steep, vertical'ish
    xE.push_back(margin + i * (imageWidth - margin) / numLines);
    yE.push_back(imageHeight - margin);
  }

  // draw lines using the algorithm that creates them from dots:
  for(i = 0; i < xE.size(); i++)
    painter.drawDottedLine(x1, y1, xE[i], yE[i], brightness);
  writeImageToFilePPM(image, "DottedLines.ppm");

  // draw lines with Wu algorithm:
  image.clear();
  for(i = 0; i < xE.size(); i++)
    drawLineWu(image, x1, y1, xE[i], yE[i], brightness);
  writeImageToFilePPM(image, "WuLines.ppm");
}

