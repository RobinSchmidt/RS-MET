#ifndef RAPT_PHASESCOPEBUFFER_H_INCLUDED
#define RAPT_PHASESCOPEBUFFER_H_INCLUDED

/** Implements the buffering for a phasescope analyzer. */

template<class TSig, class TPix, class TPar> // signal, pixel, parameter types
class PhaseScopeBuffer
{

public:

  /** Constructor. */
  PhaseScopeBuffer();

  /** Sets the sample rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the frame rate. */
  void setFrameRate(TPar newFrameRate);
  // Replace the frame rate with an integer number decayGranularity that gives the number of 
  // samples after which a multiplication of the whole buffer with the decayFactor is 
  // automatically triggered. Reasonable values should be around 1000. At 44.1kHz sample rate, the 
  // decay update rate would be 44.1 Hz which is equal to a reasonable frame rate. With lower 
  // values (around 100), we would apply the decay more often, so we would see less artifacts like
  // lines having a constant gray value for some length, then switch, etc. - the color gradients on
  // the display will be smoother, especially with low glow times

  /** Sets the overall brightness. This parameter, together with the sample rate, determines the 
  weight by which new dot are added in. */
  void setBrightness(TPix newBrightness);

  /** Sets the time it takes for "color" to decay away. */
  void setDecayTime(TPar newDecayTime);
    // rename to setPixelDecayTime

  /** Sets the size of the matrix into which we buffer the incoming samples. This should correspond 
  to the pixel size of the display. */
  void setSize(int newWidth, int newHeight);

  /** Sets the maximum size of the image into which we draw. The class will pre-allocate an 
  appropriate amount of memory for the pixels and subsequent call to setSize will not reallocate 
  memory but just interpret the preallocated memory differently 8and possibly leaving part of it 
  unused. If you call setSize with a size larger than the pre-allocated maximum size, memory 
  re-allocation will occur (and the new maximum size will be increased). When the user resizes the 
  display from the GUI thread while we write samples into the buffer in the audio thread, we don't 
  want to let memory reallocation take place... */
  void setMaxSizeWithoutReAllocation(int newMaxWidth, int newMaxHeight);

  /** Switches anti-aliasing on/off. */
  void setAntiAlias(bool shouldAntiAlias);

  /** Sets up the density of the lines the connect our actual datapoints. It actually determines the 
  number of artificial datapoints that are inserted (by linear interpolation) between our actual 
  incoming datapoints. If set to zero, it will just draw the datapoints as dots. */
  void setLineDensity(TPar newDensity);

  /** Sets a limit of the number of artificial datapoints that may be inserted per drawing 
  operation. This is important to keep the cpu usage under control. */
  void setLineDensityLimit(int newMaxNumDotsPerLine);

  /** Sets up a weight by which each pixel is not only accumulated into its actual place but also 
  into the neighbouring pixels to add more weight or thickness. It should be a value between 0 
  and 1. */
  void setPixelSpread(TPar newSpread);
  // maybe rename to setDotSpread

  /** Converts the raw left- and right signal amplitude values to the matrix indices, where the 
  data should be written. This is the xy-pixel coordinates (kept still as real numbers), where the 
  display is to be illuminated in response to the given amplitude values. */
  void toPixelCoordinates(TSig &x, TSig &y);

  /** Accepts one input sample frame for buffering. */
  void processSampleFrame(TSig left, TSig right);

  /** Applies our pixel decay-factor to the matrix of buffered values. This assumed to be called at 
  the frame rate. */
  virtual void applyPixelDecay();

  /** Resets the internal buffer to all zeros. */
  void reset();

  /** Returns a pointer to our image that we use as buffer. */
  ImageResizable<TPix> *getImage() { return &image; }

  /** Returns the sample rate. */
  inline TPar getSampleRate() { return sampleRate; }

  /** Returns the frame rate. */
  inline TPar getFrameRate() { return frameRate; }

  /** Returns the width in pixels. */
  inline int getWidth() { return image.getWidth(); }

  /** Returns the height in pixels. */
  inline int getHeight() { return image.getHeight(); }

protected:

  /** Adds a line to the given x,y coordinates (in pixel coordinates). The starting point of the 
  line are the old pixel coordinates xOld, yOld (member variables). It takes into account our line 
  density - when it's set to zero, it will just draw a dot at the new given position. */
  virtual void addLineTo(TSig x, TSig y);

  /** Draws a line by inserting a number of dots along the line. The number is proportional to the 
  given density parameter and to the Euclidean distance between the two endpoints (i.e. the length 
  of the line). The color will be scaled inversely proportional to the length, such that the total
  amount of color added to the picture is independent of the length. The maxNumDots parameter
  is for restricting the number of dots that are used which might be important in realtime 
  situations. scaleByNumDots ...
  \todo: maybe make this color scaling optional  */
  void drawDottedLine(TSig x1, TSig y1, TSig x2, TSig y2, TPix color, TPar density = 1, 
    int maxNumDots = 0, bool scaleByNumDots = false, TPar minDotDistance = 1);
  // rename to drawLineDotted

  /** Updates the pixel decay factor according to the settings of frame rate and desired decay 
  time. */
  virtual void updateDecayFactor();

  /** Updates the pixel insertion factor (i.e. the weight by which new dots are multiplied when 
  they get added in according to the settings of sample rate and brightness parameter. */
  void updateInsertFactor();


  TPar sampleRate;
  TPar frameRate;
  TPar decayTime;      // pixel illumination time
  TPar decayFactor;    // factor by which pixels decay (applied at frameRate)
  TPar lineDensity;    // density of the artificial points between actual datapoints
  int  maxDotsPerLine;
  TPar thickness;      // line (or dot) thickness from 0 to 1. 0: one pixel, 1: 3 pixels
                       // maybe rename to spread or weight or something

  TPix brightness;     // determines weight by which dots are added in
  TPix insertFactor;   // factor by which are pixels "inserted" (applied at sampleRate)

  TSig xOld, yOld;     // pixel coordinates of old datapoint (one sample ago)

  // members for actual painting on an image:
  //Image<TPix> image;
  ImageResizable<TPix> image;
  //ImagePainter<TPix, float, TSig> painter; // float: weight-type for alpha mask
  //ImagePainter<TPix, TPar, TSig> painter;  // using TPar for the 2nd TWgt template parameter
  //                                         // might not be ideal
  ImagePainter<TPix, TSig, TSig> painter;  // using TPar for the 2nd TWgt template parameter
                                           // might not be ideal

};

//=================================================================================================

/** Extends the basic PhaseScopeBuffer class with some more artistic features such as an alpha-mask 
rendered dot, blurring between frames, etc. These features are factored out into a subclass to keep 
the baseclass lean. */

template<class TSig, class TPix, class TPar> // signal, pixel, parameter types
class PhaseScopeBuffer2 : public PhaseScopeBuffer<TSig, TPix, TPar>
{

public:

  /** Constructor. */
  PhaseScopeBuffer2();

  /** Switches between simple and alpha-mask based dot drawing. */
  void setUseAlphaMask(bool shouldUseMask);

  /** Sets up a dependency of the pixel decay time on the current brightness value of the pixel. 
  When it's zero, the pixel decay time is independent of the brightness such that the decay behaves
  the same way as in the baseclass. When it's greater than zero, bright pixels will decay faster 
  than dark pixels, when it's less than zero, bright pixels will decay more slowly than dark 
  pixels. It can be used to make the decay initially fast and later slow of vice versa. */
  void setPixelDecayByValue(TPar newDecayByValue);

  /** Sets a dependency of the pixel decay time on the global average brightness of the whole 
  screen. This may help to avoid an overcrowding of the canvas by making it decay faster when its
  fuller. */
  void setPixelDecayByAverage(TPar newDecayByAverage);

  /** Sets the brightness for the additional lines drawn by the line-drawer. */
  void setLineBrightness(TPar newBrightness);

  /** Sets the width of the lines in pixels. */
  void setLineWidth(TPar newWidth);

  /** Selects one of the line-profile functions from LineDrawer::lineProfiles. */
  void setLineProfile(int newProfile);

  /** Switches dot drawing on/off */
  void setDrawDots(bool shouldDrawDots) { drawDots = shouldDrawDots; }

  /** Switches line drawing on/off */
  void setDrawLines(bool shouldDrawLines) { drawLines = shouldDrawLines; }

  /** Overriden to take into account our new decay-by-value and decay-by average features. */
  virtual void applyPixelDecay() override;

  /** Alpha mask used for drawing a "dot", a public member, such that we don't need to implement
  lots of delegations here for setting it up (maybe move to protected and provide access 
  functions). */
  AlphaMask<TSig> dotMask; 

protected:

  virtual void updateDecayFactor() override;
  virtual void addLineTo(TSig x, TSig y) override;

  TPar decayByValue = 0;     // dependency of pixel decay on current pixel value
  TPar decayByAverage = 0;   // dependency of pixel decay on global average brightness
  TPar decayFactorAt1;
  TPix lineBrightness = 0;
  bool drawLines = false;
  bool drawDots  = true;

  LineDrawer<TPix, TSig, TSig> lineDrawer;
};

#endif