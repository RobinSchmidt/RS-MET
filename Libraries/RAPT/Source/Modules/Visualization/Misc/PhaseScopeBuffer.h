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

  /** Sets the size of the matrix into which we buffer the incoming samples. This should correspond 
  to the pixel size of the display. */
  void setSize(int newWidth, int newHeight);

  /** Switches anti-aliasing on/off. */
  void setAntiAlias(bool shouldAntiAlias);

  /** Sets up the density of the lines the connect our actual datapoints. It actually determines the 
  number of artificial datapoints that are inserted (by linear interpolation) between our actual 
  incoming datapoints. If set to zero, it will just draw the datapoints as dots. */
  void setLineDensity(TPar newDensity);

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
  void bufferSampleFrame(TSig left, TSig right);

  /** Applies our pixel decay-factor to the matrix of buffered values. This assumed to be called at 
  the frame rate. */
  void applyPixelDecay();

  /** Resets the internal buffer to all zeros. */
  void reset();

  /** Returns a pointer to our image that we use as buffer. */
  Image<TPix> *getImage() { return &image; }

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
  void addLineTo(TSig x, TSig y);

  /** Updates the pixel decay factor according to the settings of frame rate and desired decay 
  time. */
  void updateDecayFactor();

  /** Updates the pixel insertion factor (i.e. the weight by which new dots are multiplied when 
  they get added in according to the settings of sample rate and brightness parameter. */
  void updateInsertFactor();


  TPar sampleRate;
  TPar frameRate;
  TPar decayTime;      // pixel illumination time
  TPar decayFactor;    // factor by which pixels decay (applied at frameRate)
  TPar lineDensity;    // density of the artificial points between actual datapoints
  TPar thickness;      // line (or dot) thickness from 0 to 1. 0: one pixel, 1: 3 pixels
                       // maybe rename to spread or weight or something

  TPix brightness;     // determines weight by which dots are added in
  TPix insertFactor;   // factor by which are pixels "inserted" (applied at sampleRate)

  TSig xOld, yOld;     // pixel coordinates of old datapoint (one sample ago)

  // members for actual painting on an image:
  Image<TPix> image;
  ImagePainter<TPix, TPar, TSig> painter;

};

//=================================================================================================

/** Extends the basic PhaseScopeBuffer class with some more artistic features such as an alpha-mask 
rendered dot, blurring bewteen frames, etc. These features are factored out into a subclass to keep 
the baseclass lean. */

template<class TSig, class TPix, class TPar> // signal, pixel, parameter types
class PhaseScopeBuffer2 : public PhaseScopeBuffer<TSig, TPix, TPar>
{

public:

  /** Constructor. */
  PhaseScopeBuffer2();

  /** Sets the dot size in pixels. */
  void setDotSize(TPar newSize);

  /** Sets the amount of blur for the dot from 0 (sharp) to 1 (maximally blurry). */
  void setDotBlur(TPar newBlur);

  /** Switches between simple and alpha-mask based dot drawing. */
  void setUseAlphaMask(bool shouldUseMask);


protected:


  AlphaMask<TPar> dotMask; // alpha mask used for drawing a "dot"

};

#endif