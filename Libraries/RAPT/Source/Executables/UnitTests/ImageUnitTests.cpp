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

  image.setSize(50, 50);
  image.clear();

  mask.setSize(5);
  r &= mask.getWidth()  == 5;
  r &= mask.getHeight() == 5;
  fillWithCheckerBoardPattern(mask);


  //painter.paintDotViaMask(10,    11,    1000);
  painter.paintDotViaMask(10.2f, 10.6f, 1000);
  // we use the 1000 here, because the painter uses this strange saturating function - maybe, we 
  // should introduce a blend-mode: mix, add, add-and-clip, add-and-saturate, multiply, ...
  // for testing here, we should use either alpha-blend or add-and-clip (should give same results)



  writeImageToFilePPM(mask,  "PaintTestMask.ppm");
  writeImageToFilePPM(image, "PaintTestImage.ppm");

  // ...

  return r;
}