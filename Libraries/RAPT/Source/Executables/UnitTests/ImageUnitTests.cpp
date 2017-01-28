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

bool imagePainterUnitTest()
{
  bool r = true;

  Image<float> image;
  AlphaMask<float> mask;
  ImagePainter<float, float, float> painter(&image, &mask);

  image.setSize(10, 10);
  image.clear();

  mask.setSize(5);
  r &= mask.getWidth()  == 5;
  r &= mask.getHeight() == 5;
  fillWithCheckerBoardPattern(mask);

  painter.paintDotViaMask(4.2f, 4.6f, 1);

  //writeImageFile(image, "PaintTest.ppm");

  // ...

  return r;
}