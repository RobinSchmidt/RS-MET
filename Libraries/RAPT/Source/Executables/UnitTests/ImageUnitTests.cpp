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

  int imageWidth  = 50;
  int imageHeight = 50;
  //int maskWidth  = 13;
  //int maskHeight = 13;

  image.setSize(imageWidth, imageHeight);
  image.clear();

  mask.setSize(7);
  r &= mask.getWidth()  == 7;
  r &= mask.getHeight() == 7;
  //fillWithCheckerBoardPattern(mask);
  fillWithCross(mask);

  //painter.paintDotViaMask(10,    11,    1);
  //painter.paintDotViaMask(10.5f, 11.5f, 1);
  //painter.paintDotViaMask(10.25f, 11.75f, 1);
  painter.paintDotViaMask(10.75f, 11.25f, 1);
  painter.paintDotViaMask(2.75f, 2.25f, 1);                        // top-left - is wrong
  painter.paintDotViaMask(imageWidth-2.75f, imageHeight-2.25f, 1); // bottom-right - looks ok

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