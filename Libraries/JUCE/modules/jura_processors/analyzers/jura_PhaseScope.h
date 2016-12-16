#ifndef jura_PhaseScope_h
#define jura_PhaseScope_h
  

/** Implements the buffering for a phasescope analyzer (maybe move this class to RAPT). */

class JUCE_API PhaseScopeBuffer
{

public:

  /** Constructor. */
  PhaseScopeBuffer();

  /** Destructor. */
  virtual ~PhaseScopeBuffer();

  /** Sets the sample rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the frame rate. */
  void setFrameRate(double newFrameRate);

  /** Sets the time it takes for "color" to decay away. */
  void setDecayTime(double newDecayTime);

  /** Sets the size of the matrix into which we buffer the incoming samples. This should correspond 
  to the pixel size of the display. */
  void setSize(int newWidth, int newHeight);

  /** Switches anti-aliasing on/off. */
  void setAntiAlias(bool shouldAntiAlias);

  /** Converts the raw left- and right signal amplitude values to the matrix indices, where the 
  data should be written. This is the xy-pixel coordinates (kept still as real numbers), where the 
  display is to be illuminated in response to the given amplitude values. */
  void convertAmplitudesToMatrixIndices(double &x, double &y);

  /** Accepts one input sample frame for buffering. */
  void bufferSampleFrame(double left, double right);

  /** Applies our pixel decay-factor to the matrix of buffered values. This assumed to be called at 
  the frame rate. */
  void applyPixelDecay();

  /** Resets the internal buffer to all zeros. */
  void reset();

  /** Returns the buffered value at the give xy pixel position. No bounds checking is done - you must
  make sure that the indices are valid. */
  inline float getValueAt(int x, int y) { return buffer[x][y]; }

  /** Returns a pointer to our internally stored data matrix. */
  float** getDataMatrix() { return buffer; }

protected:

  /** Adds a dot into our data matrix at the given position (given in matrix-index (i.e. pixel-) 
  coordinates using bilinear deinterpolation for anti-aliiasing. */
  void addDot(float x, float y);

  /** Like addDot but without anti-aliasing (we just round the coordinates to the nearest 
  integer). */
  void addDotFast(float x, float y);

  /** Allocates the memory for our data matrix buffer. */
  void allocateBuffer();

  /** Frees the memory for our data matrix buffer. */
  void freeBuffer();

  /** Updates the pixel decay factor according to the settings of frame rate and desired decay 
  time. */
  void updateDecayFactor();

  double sampleRate;
  double frameRate;
  double decayTime;
  float  decayFactor;    // factor by which pixels decay (applied at frameRate)
  float  insertFactor;   // factor by which are pixels "inserted" (applied at sampleRate)
  int    width, height;
  bool   antiAlias;      // flag to switch anti-aliasing on/off

  // the actual matrix-shaped buffer (maybe use a kind of MatrixView class that wraps a std::vector 
  // later):
  float *bufferFlat;
  float **buffer;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeBuffer)
};

//=================================================================================================

/** Implements a phasescope analyzer */

class JUCE_API PhaseScope : public jura::AudioModule
{

public:

  PhaseScope(CriticalSection *lockToUse);

  /** Sets the desired pixel size. */
  void setPixelSize(int width, int height);

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

  inline Colour getColourAt(int x, int y) 
  { 
    uint8 c = 255;
    const Colour baseColor(c, c, c, c);  // make user selectable member later
    return baseColor.withAlpha(phaseScopeBuffer.getValueAt(x, y));
  }

protected:

  PhaseScopeBuffer phaseScopeBuffer;

  friend class PhaseScopeDisplay;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope)
};

//=================================================================================================

/** Implements the GUI display for the phase scope. 

\todo Maybe this class should be derived from the CoordinateSystem baseclass and use the angular 
and radial grids from there. */

class JUCE_API PhaseScopeDisplay : public AudioModuleEditor, public Timer
{

public:

  PhaseScopeDisplay(jura::PhaseScope *newPhaseScopeToEdit);

  virtual void resized() override;
  virtual void paint(Graphics &g)	override;
  virtual void timerCallback() override;

protected:

  PhaseScope *phaseScope;

  Image image;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeDisplay)
};



#endif 