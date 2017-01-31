#ifndef jura_PhaseScope_h
#define jura_PhaseScope_h
  
//=================================================================================================

/** Implements a phasescope analyzer. */

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
  void setDotLimit(double newLimit);
  void setPixelSpread(double newSpread);
  void setPixelScale(double newFactor);
  void setAntiAlias(bool shouldAntiAlias);
  void setFrameRate(double newRate);

  // inquiry functions:
  inline double getFrameRate() { return phaseScopeBuffer.getFrameRate(); }
  inline int getInternalPixelWidth() { return phaseScopeBuffer.getWidth(); }
  inline int getInternalPixelHeight() { return phaseScopeBuffer.getHeight(); }

  // overriden from AudioModule baseclass:
  virtual AudioModuleEditor *createEditor() override;
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

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

  double pixelScale;      // scale factor between internal and external pixel sizes of the image
  int displayWidth;       // display width in pixels
  int displayHeight;      // display height in pixels

  int repaintIntervalInSamples;
  int repaintCounter;

  bool bypassPixelDecay;

  // this object is reponsible for drawing the incoming data onto a virtual screen:
  //RAPT::PhaseScopeBuffer<double, float, double> phaseScopeBuffer;
  RAPT::PhaseScopeBuffer2<double, float, double> phaseScopeBuffer;

  juce::Image image;       // image for the display
  jura::ColorMap colorMap; // the color map to translate the buffered data matrix to colors

  friend class PhaseScopeDisplay;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope)
};

//=================================================================================================

/** Implements the GUI display for the phase scope. 

\todo BUG: there's a flickering of the lines...could this be related to threading issues?
..it's more obvious with faster decay times - maybe it's because the decay is applied in the GUI 
thread whereas accumulation is done in the audio-thread? instead of calling applyPixelDecay in the 
GUI thread we could set a flag in the phaseScope audio module and the apply the decay there
-> done - this seems to help indeed but also seems to introduce tearing artifacts */

class JUCE_API PhaseScopeDisplay : public Component, public ImageUpdateListener,
  public ChangeListener, public ChangeBroadcaster
{

public:

  PhaseScopeDisplay(jura::PhaseScope *newPhaseScopeToEdit);
  virtual ~PhaseScopeDisplay();

  virtual void resized() override;
  virtual void paint(Graphics &g)	override;
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
  virtual void resized() override;

protected:

  PhaseScope *scope;
  PhaseScopeDisplay display;
  int widgetMargin;

  // Widgets:
  RSlider *sliderBrightness, *sliderAfterglow, *sliderPixelSpread, *sliderPixelScale, 
    *sliderLineDensity, *sliderDotLimit, *sliderFrameRate;
  RButton *buttonAntiAlias;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeEditor)
};



//=================================================================================================

/** Extends the basic PhaseScope by some more artistic features such as a customizable dot, 
blurring, etc. */

class JUCE_API PhaseScope2 : public jura::PhaseScope
{

public:

  PhaseScope2(CriticalSection *lockToUse);

  // additional setup functions:
  void setDotSize(double newSize);
  void setDotBlur(double newBlur);

  // overriden from PhaseScope baseclass:
  virtual void createParameters() override;
  virtual AudioModuleEditor *createEditor() override;

protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope2)
};


/** GUI editor for the extended PhaseScope */

class JUCE_API PhaseScopeEditor2 : public PhaseScopeEditor
{

public:

  PhaseScopeEditor2(jura::PhaseScope2 *newPhaseScopeToEdit);


  virtual void createWidgets() override;
  virtual void resized() override;

protected:


  // additional widgets:
  RSlider *sliderDotSize, *sliderDotBlur;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeEditor2)
};

#endif 