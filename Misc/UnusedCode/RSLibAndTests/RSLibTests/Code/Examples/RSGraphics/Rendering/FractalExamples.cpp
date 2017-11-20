#include "FractalExamples.h"

void writeImageToFilePPM(const rsImageRGBA &img, const rsFile &file)
{
  int w = img.getWidth();
  int h = img.getHeight();

  unsigned char* buf = new unsigned char[w*h*3];

  for(int y = 0; y < h; y++)
  {
    for(int x = 0; x < w; x++)
    {
      int offset = y*w*3 + x*3;
      buf[offset + 0] = img.getPixelColor(x, y).r;
      buf[offset + 1] = img.getPixelColor(x, y).g;
      buf[offset + 2] = img.getPixelColor(x, y).b;
    }
  }

  FILE* fd = fopen("test.ppm", "wb");  // \todo use the passed filename
  fprintf(fd, "P6\n%d %d\n255\n", w, h);
  fwrite(buf, 1, w*h*3, fd);
  fclose(fd);
      
  delete [] buf;
}

void createDefaultFractal()
{
  FractalArtGenerator g;
  rsImageRGBA img = g.getOutputImage();

  rsString path = rsGetCurrentApplicationDirectory();
  path += "\\MandelbrotDefault.ppm";
  rsFile file(path);
  writeImageToFilePPM(img, file);
}

void createMandelbrotSet1()
{
  FractalArtGenerator g;
  g.setImageResolution(400, 400);
  g.setFractalIterator(new MandelbrotIterator);
  g.setMaxIterations(20);
  g.setInsideIntensityAlgorithm( new FractalIntensityLastZ);
  g.setOutsideIntensityAlgorithm(new ConstantFractalIntensity(0.0f));
  //g.setInsideIntensityMapper( new IntensityMapperDownUp);
  //g.setOutsideIntensityMapper(new IntensityMapperIdentity);

  rsImageRGBA img = g.getOutputImage();

  rsString path = rsGetCurrentApplicationDirectory();
  path += "\\Mandelbrot1.ppm";
  rsFile file(path);
  writeImageToFilePPM(img, file);
}

void createBuddhaBrot1()
{
  FractalArtGenerator g;
  g.setImageResolution(800, 800);
  g.setFractalIterator(new MandelbrotIterator);
  g.setMaxIterations(20);
  //g.setInsideIntensityAlgorithm( new ConstantFractalIntensity(0.0f));
  g.setInsideIntensityAlgorithm( new FractalIntensityBuddhaBrot);
  g.setOutsideIntensityAlgorithm(new FractalIntensityBuddhaBrot);

  rsImageRGBA img = g.getOutputImage();

  rsString path = rsGetCurrentApplicationDirectory();
  path += "\\Mandelbrot1.ppm";
  rsFile file(path);
  writeImageToFilePPM(img, file);
}