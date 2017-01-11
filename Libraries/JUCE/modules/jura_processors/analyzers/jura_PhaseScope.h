#ifndef jura_PhaseScope_h
#define jura_PhaseScope_h
  
//=================================================================================================

/** Implements a phasescope analyzer. 

\todo
-maybe implement a linear decay in addition to the exponential decay: 
 new = expDecay * old - linDecay (the result must be clipped at zero)
-line density should be either 0 or 1, intermediate values are not useful
-independently adjustable line-brightness and dot-brightness  
 -the same thing also for the thicknesses
-make line-thickness dependent on dot-distance via parameter
-implement different drawing modes (the current one, and some based on juce's vector drawing engine
 with lines, circles, etc.)
-implement rainbow mode (hue-rotation, we may choose a basic hue and rotation speed which may be 
 set to zero)
*/

class JUCE_API PhaseScope : public jura::AudioModule, public jura::ImageUpdater
{

public:

  PhaseScope(CriticalSection *lockToUse);

  /** Creates the parameters to control the drawing. 
  \todo: maybe make this function a virtual member function of the AudioModule baseclass to be 
  overriden - the call it from the constructor there. */
  virtual void createParameters();

  /** Sets the desired pixel size for the display. This function will update the internal pixel 
  size of the buffer according to the passed values and the desired rescaling factor. " */
  void setDisplayPixelSize(int width, int height);

  // parameter setup functions (to be used for the callbacks from the parameters):
  void setBrightness(double newBrightness);
  void setAfterGlow(double newGlow);
  void setLineDensity(double newDensity);
  void setPixelSpread(double newSpread);
  void setPixelScale(double newFactor);
  void setAntiAlias(bool shouldAntiAlias);
  void setFrameRate(double newRate);
  //void setDrawingMode(int newMode);
  void setRainbowMode(bool shouldUseRainbowColors); // maybe provide more modes and a function 
    // setColorMode(int newMode) - can have different settings: fixed color, hue rotation, 
    // alternating colors, colormapped values, etc.

  //inline void triggerPixelDecay() { needsPixelDecay = true; }  // obsolete

  // inquiry functions:
  inline double getFrameRate() { return phaseScopeBuffer.getFrameRate(); }
  inline int getInternalPixelWidth() { return phaseScopeBuffer.getWidth(); }
  inline int getInternalPixelHeight() { return phaseScopeBuffer.getHeight(); }

  // overriden from AudioModule baseclass:
  AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

  /** Returns the color that should be used for this frame. We assume here that this function is 
  called at frame rate (it will update the internal color-period counter, so it msut be called at 
  the correct rate - otherwise the color period will be wrong). */
  Colour getAndUpdateColor();

protected:

  /** Updates the size of the internal buffer according to the settings of displayWidth, 
  displayHeight and pixelScale. */
  void updateBufferSize();

  /** Updates the image for the scope picture, i.e. writes/converts the content of the 
  phaseScopeBuffer member into the image member. */
  void updateScopeImage();

  /** Updates our repaint triggering interval - needs to be called whenever frame rate or sample
  rate changes. It's the number of samples between two frames. */
  void updateRepaintInterval();

  //bool needsPixelDecay;   // flag to indicate that a pixel decay should be triggered
  bool rainbow;           // indicates usage of rainbow colors (i.e. hue rotation)
  double colorPeriod;     // period for one complete color change cycle (in seconds)
  double colorCounter;

  double pixelScale;      // scale factor between internal and external pixel sizes of the image
  int displayWidth;       // display width in pixels
  int displayHeight;      // display height in pixels


  int repaintIntervalInSamples;
  int repaintCounter;

  juce::Image image;

  // this object is reponsible for drawing the incoming data onto a virtual screen:
  RAPT::PhaseScopeBuffer<double, float, double> phaseScopeBuffer;

  friend class PhaseScopeDisplay;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope)
};

//=================================================================================================

/** Implements the GUI display for the phase scope. 

\todo Maybe this class should be derived from the CoordinateSystem baseclass and use the angular 
and radial grids from there. 

\todo BUG: there's a flickering of the lines...could this be related to threading issues?
..it's more obvious with faster decay times - maybe it's because the decay is applied in the GUI 
thread whereas accumulation is done in the audio-thread? instead of calling applyPixelDecay in the 
GUI thread we could set a flag in the phaseScope audio module and the apply the decay there
-> done - this seems to help indeed but also seems to introduce tearing artifacts */

class JUCE_API PhaseScopeDisplay : public Component, public ImageUpdateListener,
  public ChangeListener, public ChangeBroadcaster
  /*public Timer*/
{

public:

  PhaseScopeDisplay(jura::PhaseScope *newPhaseScopeToEdit);
  virtual ~PhaseScopeDisplay();

  virtual void resized() override;
  virtual void paint(Graphics &g)	override;
  //virtual void timerCallback() override;
  virtual void imageWasUpdated(juce::Image* image) override;
  virtual void changeListenerCallback(ChangeBroadcaster *source) override;

protected:

  PhaseScope *phaseScope;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeDisplay)
};

//=================================================================================================

/** Implements a GUI editor for the XY phase scope. */

class JUCE_API PhaseScopeEditor : public AudioModuleEditor
{

public:

  PhaseScopeEditor(jura::PhaseScope *newPhaseScopeToEdit);

  virtual void createWidgets();

  //void setDisplayPixelSize(int newWidth, int newHeight);

  virtual void resized() override;

protected:

  PhaseScope *scope;
  PhaseScopeDisplay display;
  int widgetMargin;

  // Widgets:
  RSlider *sliderBrightness, *sliderAfterglow, *sliderPixelSpread, *sliderPixelScale, 
    *sliderLineDensity, *sliderFrameRate;
  RButton *buttonAntiAlias, *buttonRainbow;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeEditor)
};

#endif 