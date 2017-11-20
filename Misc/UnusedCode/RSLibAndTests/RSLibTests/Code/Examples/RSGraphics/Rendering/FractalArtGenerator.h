#ifndef RS_FRACTALARTGENERATOR_H
#define RS_FRACTALARTGENERATOR_H

#include "../../../Common/TestUtilities.h"

//=================================================================================================

/** An abstract baseclass (and its subclasses) that is used as a strategy object inside the 
FractalArtGenerator class to map grayscale values from the range 0...1 onto itself (but possibly 
nonlinear and/or nonmonotonic). 

\todo get rid of the baseclass/subclasses implementation - just have one class that defines a 
flexible mapping function via breakpoints or via a nonlinear function defined by an expression
that will be evaluated ...hmhm maybe we should provide both options, so we might still need the
baseclass/subclasses approach (one subclas for breakpoints, one for expressions)

*/

class IntensityMapper
{
public:
  virtual float mapValue(float in) = 0;
};
class IntensityMapperIdentity : public IntensityMapper
{
public:
  virtual float mapValue(float in) 
  {
    return in;
  }
};
class IntensityMapperDownUp : public IntensityMapper
{
public:
  virtual float mapValue(float in) 
  {
    if( in <= 0.5f )
      return +1.f - 2.f*in;
    else
      return -1.f + 2.f*in;
  }
};

/** A class for defining a mapping from real numbers to colors. */

class ColorMapRGBA
{
public:
  rsColorRGBA mapValue(double index)
  {
    return rsColorRGBA((float)index); // preliminary
      // \todo: allow the user to define an arbitrary color-gradient via breakpoints. when the 
      // index is 0.0, the leftmost color is returned, when it's 1.0 the rightmost
  }
};

//=================================================================================================

/** Abstract baseclass for iteration recurrences that may lead to either divergent or non-divergent
behavior, based on the value of some complex parameter c. The subclasses will actually implement 
these iteration formulas. Typically, they will be such that the boundary between parameter values 
for which the iteration will and won't diverge, has a self-similar (aka fractal) structure. */

class FractalIterator
{

public:

  /** Constructor. Initializes z, z0 and c to zero. */
  FractalIterator() 
  { 
    z = z0 = c = rsComplexDbl(0.0, 0.0); 
    maxNumIterations = 20;
    zTrajectory = new rsComplexDbl[maxNumIterations+1];
  }

  /** Destructor */
  ~FractalIterator() 
  { 
    delete[] zTrajectory;
  }



  /** Sets the complex parameter which is to interpreted as a point in the complex plane that is
  either inside or outside the set. */
  void setParameter(rsComplexDbl newC) { c = newC; }

  /** Sets the maximum number of iterations ot be taken. */
  void setMaxNumIterations(int newMax) 
  { 
    if( newMax != maxNumIterations )
    {
      delete[] zTrajectory;
      maxNumIterations = newMax; 
      zTrajectory = new rsComplexDbl[maxNumIterations+1];
    }
  } 



  /** Should be overriden by subclasses to return, whether or not divergence can be guaranteed 
  based on the current value of z. For example, in the iteration that leads to the Mandelbrot set, 
  divergence can be guaranteed as soon as the absolute value of z exceeds 2. */
  virtual bool isDivergenceCertain() = 0;

  /** Can be overriden by subclasses to indicate that for the current values of c and z, it can be
  assured that the iteration does not diverge. Sometimes, this can be checked by closed from 
  criteria which can be used to bypass the iteration and thereby optimize the computation. The 
  baseclass implementation just returns false which means the the iteration is entered. */
  virtual bool isNonDivergenceCertain() { return false; }



  /** Resets the iteration variable z to its initial value z0. */
  virtual void reset() 
  { 
    z = z0; 
    iterationCounter = 0;
  }

  /** Performs a number of iterations of the formula up to a maximum and thereby records the 
  trajectory of the z-values into our member array zTrajectory. */
  void iterate()
  {
    reset();
    zTrajectory[0] = z;
    for(iterationCounter = 1; iterationCounter <= maxNumIterations; iterationCounter++)
    {
      advanceZ();
      zTrajectory[iterationCounter] = z;
      if( isDivergenceCertain() )
        break;  // commented out for Buddhabrot - hmm - still wrong but different and interesting
    }
  }


//protected:

  /** Should be overriden to advance the iteration variable z to the next iteration. */
  virtual void advanceZ() = 0;

  rsComplexDbl c;             // the complex parameter
  rsComplexDbl z0;            // initial value for z
  rsComplexDbl z;             // the variable for the iteration
  rsComplexDbl *zTrajectory;  // trajectory of the z-values

  int  maxNumIterations;
  int  iterationCounter;

};

//=================================================================================================

/** Abstract baseclass for algorithms that assign an intensity to each pixel in fractal images. 
Such "intensities" will later be converted to colormap index. */

class FractalIntensityAlgorithm
{

public:

  virtual void updatePixelIntensities(
    rsImage<float> &img, 
    int ix, int iy, // possibly redundant - could be computed from iterator->c
    FractalIterator *iterator,    
    double xMin, double xMax, double yMin, double yMax) = 0;
      // are these parameters actually used? if not, get rid of them

};

//=================================================================================================

/** Wraps two FractalColorizingAlgorithm objects (one for the inside and one for the outside) into
one object and also handles mapping between color-indices to actual colors */

class FractalColorizer
{

public:


  /** \name Construction/Destruction */

  /** Default constructor. Uses solid black for inside points and white for outside points. */
  FractalColorizer::FractalColorizer();

  /** Destructor. */
  FractalColorizer::~FractalColorizer();


  /** \name Setup */

  /** Sets the algorithm that will be used to compute the intensities for inside points. This 
  object takes over ownership over the passed object and will evetually delete it. */
  void setInsideAlgorithm(FractalIntensityAlgorithm *newAlgorithm);

  /** Sets the algorithm that will be used to compute the intensities for outside points. This 
  object takes over ownership over the passed object and will evetually delete it. */
  void setOutsideAlgorithm(FractalIntensityAlgorithm *newAlgorithm);


  /** \name Color computation */

  /** Computes and sets a color for every pixel in the passed image. */
  void computeColors(rsImageRGBA &img, 
    FractalIterator *iterator, int maxIterations, 
    double xMin, double xMax, double yMin, double yMax);

  /** Computes an "intensity" (to be used as index into a colormap) for every pixel in the passed 
  image. The value of this index is normalized to the range 0...1 and the boolean flag in the 
  bool-image determines whether a point is considred as inside or outside the set. */
  virtual void computePixelIntensities(rsImage<float> &intensities, rsImage<bool> &insideFlags,
    FractalIterator *iterator, int maxIterations, 
    double xMin, double xMax, double yMin, double yMax);

  /** Computes a number for every pixel in the passed image that will later be converted to an 
  index in a colormap. This value here is raw, which means, it's not necessarily normalized to
  the range 0...1 and no transfer function is applied. */
  virtual void computeRawPixelIntensities(rsImage<float> &intensities, rsImage<bool> &insideFlags,
    FractalIterator *iterator, int maxIterations, 
    double xMin, double xMax, double yMin, double yMax);

protected:

  /** Normalizes raw pixel intensity values into the range 0...1 (or 0...-1 for outside pixels) by
  subtracting the minium value and then scaling by 1/(max-min). */
  //void normalizePixelIntensities(float *buffer, int length);
  void normalizePixelIntensities(rsImage<float> &intensities, rsImage<bool> &insideFlags);


  FractalIntensityAlgorithm *insideAlgorithm,       *outsideAlgorithm;
  IntensityMapper           *insideIntensityMapper, *outsideIntensityMapper;
  ColorMapRGBA               insideColorMap,         outsideColorMap;

};

//=================================================================================================

/** A class for rendering fractal art images. This class is the high-level facade that is supposed 
to be used by client code. It wraps an object for the actual fractal iteration formula and an 
object for the coloring rules. */

class FractalArtGenerator
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  FractalArtGenerator();

  /** Destructor. */
  virtual ~FractalArtGenerator();


  /** \name Setup */

  /** Sets the resolution of the final output image in pixels. */
  void setImageResolution(int newWidth, int newHeight);

  /** Sets the iteration object to be used. The FractalArtGenerator takes over ownership of the 
  object and will take care for its deletion. */
  void setFractalIterator(FractalIterator *newIterator);

  /** Sets the maximum number of iterations that will be taken until a point is considered 
  non-divergent (i.e. inside the set). Higher values tend to classify more points as being outside 
  the set, so higher values lead to a kind of erosion compared to images rendered with lower 
  values. Also, higher values lead to more accurate classification with respect to the true set. 
  Put differently, high values reduce the chance of misclassification, and a misclassification 
  will always be a false positive - i.e. a point is classified to be in the set when in fact it is 
  not. */
  void setMaxIterations(int newMaxIterations);

  /** Sets the intensity-computation algorithm to be used for inside points. The 
  FractalArtGenerator takes over ownership of the object and will take care for its deletion. */
  void setInsideIntensityAlgorithm(FractalIntensityAlgorithm *newAlgorithm);

  /** Sets the intensity-computation algorithm to be used for outside points. The 
  FractalArtGenerator takes over ownership of the object and will take care for its deletion. */
  void setOutsideIntensityAlgorithm(FractalIntensityAlgorithm *newAlgorithm);


  /** \name Rendering: */

  /** Renders the output image and returns it (rendering may take a long time). */
  rsImageRGBA getOutputImage();

protected:

  double xMin, xMax, yMin, yMax; // determines the visible area of the xy-plane

  int w, h;  // image width and height -> maybe use an rsImageRGBA member which manages these

  int maxIterations;

  FractalIterator  *iterator;
  FractalColorizer *colorizer;

};






//=================================================================================================
// concrete subclasses for fractal iterators:

class MandelbrotIterator : public FractalIterator
{

public:

  MandelbrotIterator() 
  { 
    bailOutRadiusSquared = 4.0; 
    //bailOutRadiusSquared = 8.0;  // test
    //bailOutRadiusSquared = 1000.0;  // test
  }

  /** Sets the radius of z, at which the sequence is considered to be certainly divergent. 2.0 is
  the smallest radius after which the seuqence is certain to diverge, but you can choose a larger 
  value for greater flexibility in the visual output. */
  void setBailOutRadius(double newRadius) { bailOutRadiusSquared = newRadius*newRadius; }

  /** Applies the iteration formula z <- z^2 + c to our iteration variable z. */
  virtual void advanceZ() { z = z*z + c; }

  /** Checks the divergence criterion which is abs(z) > R for some radius R >= 2. */
  virtual bool isDivergenceCertain() { return z.re*z.re + z.im*z.im > bailOutRadiusSquared; }

  /** Performs some closed-form formula criteria, which, when true, indicate that the point is
  inside the Mandelbrot set without having to resort to iteration. */
  //virtual bool isNonDivergenceCertain();

protected:

  double bailOutRadiusSquared;

};

//=================================================================================================
// concrete subclasses for fractal intensity calculation algorithms:

class ConstantFractalIntensity : public FractalIntensityAlgorithm
{

public:

  ConstantFractalIntensity::ConstantFractalIntensity(float valueToUse)
  {
    value = valueToUse;
  }

  void updatePixelIntensities(rsImage<float> &img, 
    int ix, int iy, 
    FractalIterator *iterator,    
    double xMin, double xMax, double yMin, double yMax)
  {
    img.setPixelColor(ix, iy, value);
  }

  float value;

};

// for algorithms based on the last value of z, the maximum number of iterations has a great impact
// on the intensity value - compute inside intesities with 20 and 21 as maximum number of 
// iterations and compare. that could be used to use the rendition with 20 for luminance and the
// rendition with 21 for hue (and maybe a third one with 22 for saturation)
class FractalIntensityLastZ : public FractalIntensityAlgorithm
{
public:
  void updatePixelIntensities(rsImage<float> &img, 
    int ix, int iy, 
    FractalIterator *iterator,    
    double xMin, double xMax, double yMin, double yMax)
  {
    rsComplexDbl z = iterator->z;
    double       a = 0.25;  // to avoid division by zero, can be made a user-parameter later

    float value = (float) (1.0 / (a + z.re*z.re + z.im*z.im)); 
       // later allow different formulas here, maybe include an expression evaluator with a,b,c,d
       // parameters

    img.setPixelColor(ix, iy, value);
  }
};


// the results are not yet as expected - probabl, we indeed need to choose values of c 
// probabilistically -
class FractalIntensityBuddhaBrot : public FractalIntensityAlgorithm
{
public:

  void updatePixelIntensities(rsImage<float> &img, 
    int ix, int iy, 
    FractalIterator *iterator,    
    double xMin, double xMax, double yMin, double yMax)
  {
    if( !iterator->isDivergenceCertain() )
      return;

    int w = img.getWidth();
    int h = img.getHeight();
    double sx = w / (xMax-xMin);  // scaler for computing horizontal pixel-index from z.re
    double sy = h / (yMax-yMin);  // 

    for(int i = 0; i <= iterator->iterationCounter; i++)
    {
      double zx = iterator->zTrajectory[i].re;
      double zy = iterator->zTrajectory[i].im;

      // find pixel that gets an increment and incremet it (later: de-interpolate the increment 
      // into 4 pixels):
      double px = sx * (zx - xMin) - 0.5;
      double py = sy * (zy - yMin) - 0.5;
      ix = rsRoundToInt(px);
      iy = rsRoundToInt(py);
      if( img.arePixelCoordinatesValid(ix, iy) )
        *img.getPointerToPixel(ix, iy) += 1.f;
    }
  }

};





#endif
