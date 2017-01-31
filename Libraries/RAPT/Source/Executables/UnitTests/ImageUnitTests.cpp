#include "ImageUnitTests.h"

void fillWithCheckerBoardPattern(Image<float>& image)
{
  image.clear();
  for(int x = 0; x < image.getWidth(); x++) {
    for(int y = 0; y < image.getHeight(); y++) {
      if( (x+y) % 2 == 0) // x+y even -> fill pixel
        image(x, y) = 1;
    }
  }
}

void fillWithCross(Image<float>& image)
{
  image.clear();
  int cx = image.getWidth()  / 2;
  int cy = image.getHeight() / 2;
  for(int x = 0; x < image.getWidth(); x++)
    image(x, cy) = 1;
  for(int y = 0; y < image.getHeight(); y++)
    image(cx, y) = 1;
}



bool imagePainterUnitTest()
{
  bool r = true;

  Image<float> image;
  AlphaMask<float> mask;  // maybe use a regular image as mask
  ImagePainter<float, float, float> painter(&image, &mask);

  int imageWidth  = 30;
  int imageHeight = 30; // 50x50 image with 3x3 mask gives an access violation
  int maskSize    = 7;

  // maybe, we should 1st use the simplest case: 1x1 mask

  image.setSize(imageWidth, imageHeight);
  image.clear();

  mask.setSize(maskSize);
  r &= mask.getWidth()  == maskSize;
  r &= mask.getHeight() == maskSize;
  //fillWithCheckerBoardPattern(mask);
  //fillWithCross(mask);
  mask.fillAll(1.f);  // full white

  //painter.paintDotViaMask(0.25f, 0.75f, 1);
  //painter.paintDotViaMask(2.25f, 3.75f, 1);
  //painter.paintDotViaMask(3.25f, 3.75f, 1);

  // draw in center and at all 4 edges:
  float dx = 0.5;
  float dy = 0.5;
  float w  = imageWidth;
  float w2 = w/2;
  float h  = imageHeight;
  float h2 = h/2;
  float b  = 1.0f; // brightness
  painter.paintDotViaMask(    dx,     dy, b);   // top-left
  painter.paintDotViaMask(w2 +dx,     dy, b);   // top-center
  painter.paintDotViaMask(w-1+dx,     dy, b);   // top-right
  painter.paintDotViaMask(    dx, h2 +dy, b);   // center-left
  painter.paintDotViaMask(w-1+dx, h2 +dy, b);   // center-right
  painter.paintDotViaMask(    dx, h-1+dy, b);   // bottom-left
  painter.paintDotViaMask(w2 +dx, h-1+dy, b);   // bottom-center
  painter.paintDotViaMask(w-1+dx, h-1+dy, b);   // bottom-right
  painter.paintDotViaMask(w2 +dx, h2 +dy, b);   // center

  // it seems, things drawn on the left border leak into the right border


  //painter.paintDotViaMask(10.2f, 10.6f, 1);
  // we use the 1000 here, because the painter uses this strange saturating function - maybe, we 
  // should introduce a blend-mode: mix, add, add-and-clip, add-and-saturate, multiply, ...
  // for testing here, we should use either alpha-blend or add-and-clip (should give same results)


  // it seems to work, as long as we are not too close to the edge of the image

  writeImageToFilePPM(mask,  "PaintTestMask.ppm");
  writeImageToFilePPM(image, "PaintTestImage.ppm");

  // ...

  return r;
}