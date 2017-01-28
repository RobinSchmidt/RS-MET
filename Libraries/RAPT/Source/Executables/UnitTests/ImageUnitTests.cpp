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

// move this function to the shared code, make a version that supports color:
void writeImageToFilePPM(const Image<float>& img, const char* path)
{
  int w = img.getWidth();
  int h = img.getHeight();
  unsigned char* buf = new unsigned char[w*h*3];

  for(int y = 0; y < h; y++) {
    for(int x = 0; x < w; x++) {
      int i = y*w*3 + x*3;
      unsigned char gray = (unsigned char) (255 * img.getPixelColor(x, y));
      buf[i+0] = gray;
      buf[i+1] = gray;
      buf[i+2] = gray;
    }
  }

  FILE* fd = fopen(path, "wb");  // "wb": write binary
  fprintf(fd, "P6\n%d %d\n255\n", w, h);
  fwrite(buf, 1, w*h*3, fd);

  fclose(fd);
  delete[] buf;
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