#include "FractalArtGenerator.h"


// we use this subclass for testing the BuddhaBrot rendering technique - for actual production 
// code, this should somehow be done differently
// something is still amiss - there's a sharp circular area which is brighter then the rest
class BuddhaBrotColorizer : public FractalColorizer
{

public:

  BuddhaBrotColorizer()
  {
    colorInsidePixels  = false;  // false, for the normal Buddhabrot, true for Anti-Buddhabrot
    colorOutsidePixels = true;  // vice versa

    // the anti-buddhabrot do not render correctly in release build, in a debug build, it seems
    // to work - there must be some bug
  }

protected:

  // override:
  virtual void computeRawPixelIntensities(rsImage<float> &intensities, rsImage<bool> &insideFlags,
    FractalIterator *iterator, int maxIterations, 
    double xMin, double xMax, double yMin, double yMax)
  {
    insideFlags.fillAll(true);  // pixels initially considered inside... 
    intensities.fillAll(0.f);   // ...with zero intensity

    int seed            = 0;
    int numSamplePoints = 1000000;   // draft
    //int numSamplePoints = 10000000;   // low-quality
    //int numSamplePoints = 100000000;   // medium-quality
    //int numSamplePoints = 1000000000;   // high-quality
    int w = intensities.getWidth();
    int h = intensities.getHeight();
    double sx = w / (xMax-xMin);  // scaler for computing horizontal pixel-index from z.re
    double sy = h / (yMax-yMin); 
    double cx, cy;                // x- and y-value for the complex parameter c
    rsRandomUniform(xMin, xMax, seed);
    for(int k = 1; k <= numSamplePoints; k++)
    {
      cx = rsRandomUniform(xMin, xMax);
      cy = rsRandomUniform(yMin, yMax);

      iterator->setParameter(rsComplexDbl(cx, cy));
      iterator->iterate();

      //insideFlags.setPixelColor(ix, iy, false);

      if(  ( iterator->isDivergenceCertain() && colorOutsidePixels) 
        || (!iterator->isDivergenceCertain() && colorInsidePixels ) )
      {


        // used recorded trajectory to fill pixels:
        for(int i = 1; i <= iterator->iterationCounter; i++)
        {
          double zx = iterator->zTrajectory[i].re;
          double zy = iterator->zTrajectory[i].im;

          // find pixel that gets an increment and incremet it (later: de-interpolate the increment 
          // into 4 pixels):
          double px = sx * (zx - xMin) - 0.5;
          double py = sy * (zy - yMin) - 0.5;

          int    ix = rsRoundToInt(px);
          int    iy = rsRoundToInt(py);

          //int    ix = (int) floor(px+0.5f);  // test: floor(x+0.5) should give a rounded value
          //int    iy = (int) floor(py+0.5f);

          if( intensities.arePixelCoordinatesValid(ix, iy) )
            *intensities.getPointerToPixel(ix, iy) += 1.f;
        }
      }
    }
  }


  bool colorInsidePixels, colorOutsidePixels;

};


//-------------------------------------------------------------------------------------------------
// class FractalArtGenerator

// Construction/destruction:

FractalArtGenerator::FractalArtGenerator()
{
  xMin = -1.5;
  xMax = +0.5;
  yMin = -1.0;
  yMax = +1.0;

  maxIterations = 20;

  w = 800;   
  h = 800;

  iterator  = new MandelbrotIterator;
  colorizer = new FractalColorizer();

 

  //colorizer = new BuddhaBrotColorizer();  // test
}

FractalArtGenerator::~FractalArtGenerator()
{
  delete iterator;
  delete colorizer;
}

// Setup:

void FractalArtGenerator::setImageResolution(int newWidth, int newHeight)
{
  w = newWidth;
  h = newHeight;
}
  
void FractalArtGenerator::setFractalIterator(FractalIterator *newIterator)
{
  if( iterator != newIterator )
  {
    delete iterator;
    iterator = newIterator;
  }
}

void FractalArtGenerator::setMaxIterations(int newMaxIterations)
{
  maxIterations = newMaxIterations;
}

void FractalArtGenerator::setInsideIntensityAlgorithm(FractalIntensityAlgorithm *newAlgorithm)
{
  colorizer->setInsideAlgorithm(newAlgorithm);
}

void FractalArtGenerator::setOutsideIntensityAlgorithm(FractalIntensityAlgorithm *newAlgorithm)
{
  colorizer->setOutsideAlgorithm(newAlgorithm);
}

// Rendering:

rsImageRGBA FractalArtGenerator::getOutputImage()
{
  rsImageRGBA img(w, h);
  colorizer->computeColors(img, iterator, maxIterations, xMin, xMax, yMin, yMax);

   // \todo maybe apply some post-processing here....

  return img;
}

//-------------------------------------------------------------------------------------------------
// subclasses of FractalIterator
/*
// spoils gradients in inside areas
bool MandelbrotIterator::isNonDivergenceCertain() 
{
  if( z != rsComplexDbl(0.0, 0.0) ) // i think, the criteria below are only valid for z0 = 0
    return false; 

  // check, if c is inside the main cardioid (which leads to a convergent sequence):
  double tmp = c.re - 0.25;
  double q   = tmp*tmp + c.im*c.im;
  if( q*(q+tmp) < 0.25*c.im*c.im )
    return true; 

  // check, if c is inside the circle left to the main cardioid (which leads to a 2-periodic
  // sequence):
  tmp = c.re + 1;
  if( tmp*tmp + c.im*c.im < 1.0/16.0 )
    return true;

  return false;
}
*/

//-------------------------------------------------------------------------------------------------
// class FractalColorizer:

FractalColorizer::FractalColorizer()
{
  insideAlgorithm  = new ConstantFractalIntensity(0.0f);
  outsideAlgorithm = new ConstantFractalIntensity(1.0f);

  insideIntensityMapper  = new IntensityMapperIdentity;
  outsideIntensityMapper = new IntensityMapperIdentity;

  // insideColorMapper...

}

FractalColorizer::~FractalColorizer()
{
  delete insideAlgorithm;
  delete outsideAlgorithm;

  delete insideIntensityMapper;
  delete outsideIntensityMapper;
}
  
void FractalColorizer::setInsideAlgorithm(FractalIntensityAlgorithm *newAlgorithm)
{
  if( insideAlgorithm != newAlgorithm )
  {
    delete insideAlgorithm;
    insideAlgorithm = newAlgorithm;
  }
}

void FractalColorizer::setOutsideAlgorithm(FractalIntensityAlgorithm *newAlgorithm)
{
  if( outsideAlgorithm != newAlgorithm )
  {
    delete outsideAlgorithm;
    outsideAlgorithm = newAlgorithm;
  }
}

void FractalColorizer::computeColors(rsImageRGBA &img, FractalIterator *iterator, 
  int maxIterations, double xMin, double xMax, double yMin, double yMax)
{
  int w = img.getWidth();
  int h = img.getHeight();
  rsImage<float> intensities(w, h);
  rsImage<bool>  insideFlags(w, h);

  computePixelIntensities(intensities, insideFlags, iterator, maxIterations, 
    xMin, xMax, yMin, yMax);

  // map intensities to actual colors (for inside and outside points separately):
  for(int iy = 0; iy < h; iy++)
  {
    for(int ix = 0; ix < w; ix++)
    {
      bool        isInside = insideFlags.getPixelColor( ix, iy);
      float       *pFloat  = intensities.getPointerToPixel(ix, iy);
      rsColorRGBA *pColor  = img.getPointerToPixel(     ix, iy);

      // apply intensity and/or color-mapping:
      double intensity;
      if( isInside )
      {
        intensity = insideIntensityMapper->mapValue(*pFloat);
        *pColor   = insideColorMap.mapValue(intensity);
      }
      else
      {
        intensity = outsideIntensityMapper->mapValue(*pFloat);
        *pColor   = outsideColorMap.mapValue(intensity);
      }
      // code-duplication - try to get rid of it
    }
  }
}

void FractalColorizer::computePixelIntensities(
  rsImage<float> &intensities, rsImage<bool> &insideFlags, 
  FractalIterator *iterator, int maxIterations, 
  double xMin, double xMax, double yMin, double yMax)
{
  computeRawPixelIntensities(intensities, insideFlags, iterator, maxIterations, 
    xMin, xMax, yMin, yMax);

  // test:
  rsNormalize(intensities.getPointerToPixel(0, 0), intensities.getNumPixels(), 1.f);

  //normalizePixelIntensities(intensities, insideFlags);
}

void FractalColorizer::normalizePixelIntensities(
  rsImage<float> &intensities, rsImage<bool> &insideFlags)
{
  // retrieve total length and pointers in order to avoid nested loops later:
  int   length        = intensities.getNumPixels();
  float *pIntensities = intensities.getPointerToPixel(0, 0);
  bool  *pInsideFlags = insideFlags.getPointerToPixel(0, 0);

  // find minimum and maximum values for inside and outside points separately:
  float minInside  = RS_INF(float);
  float minOutside = RS_INF(float);
  float maxInside  = 0.f;
  float maxOutside = 0.f;
  for(int k = 0; k < length; k++)
  {
    if( pInsideFlags[k] == true )
    {
      minInside = rsMin(minInside,  pIntensities[k]);
      maxInside = rsMax(maxInside,  pIntensities[k]);
    }
    else
    {
      minOutside = rsMin(minOutside, pIntensities[k]);
      maxOutside = rsMax(maxOutside, pIntensities[k]);
    }
  }

  // map values into 0...1 for inside and outside separately (unless all values are the same in 
  // which case they should remain as is):
  float scalerInside; 
  if( maxInside > minInside )
    scalerInside = 1.f / (maxInside - minInside);
  else
  {
    scalerInside = 1.f;
    minInside    = 0.f;
  }
  float scalerOutside;
  if( maxOutside > minOutside )
    scalerOutside = 1.f / (maxOutside - minOutside);
  else
  {
    scalerOutside = 1.f;
    minOutside    = 0.f;
  }
  for(int k = 0; k < length; k++)
  {
    if( pInsideFlags[k] == true )
      pIntensities[k] = scalerInside  * (pIntensities[k] - minInside);
    else
      pIntensities[k] = scalerOutside * (pIntensities[k] - minOutside);
  }

  // it's all a bit ugly due to the code duplication, but i have no idea how to improve it 
  // meaningfully
}

void FractalColorizer::computeRawPixelIntensities(
  rsImage<float> &intensities, rsImage<bool> &insideFlags,
  FractalIterator *iterator, int maxIterations, 
  double xMin, double xMax, double yMin, double yMax)
{
  int    w  = intensities.getWidth();
  int    h  = intensities.getHeight();
  double dx = (xMax-xMin) / w;
  double dy = (yMax-yMin) / h;
  double x, y;
  insideFlags.fillAll(true);  // pixels initially considered inside... 
  intensities.fillAll(0.f);   // ...with zero intensity

  for(int iy = 0; iy < h; iy++)
  {  
    for(int ix = 0; ix < w; ix++)
    {
      y = yMin + (iy+0.5)*dy;
      x = xMin + (ix+0.5)*dx;
      iterator->setParameter(rsComplexDbl(x, y));
      iterator->iterate();
      if( iterator->isDivergenceCertain() )
      {            
        outsideAlgorithm->updatePixelIntensities(
          intensities, ix, iy, iterator, xMin, xMax, yMin, yMax);
        insideFlags.setPixelColor(ix, iy, false);
      }
      else
      {
        insideAlgorithm->updatePixelIntensities(
          intensities, ix, iy, iterator, xMin, xMax, yMin, yMax);
      }
    }
  }
}
