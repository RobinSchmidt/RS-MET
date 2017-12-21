#ifndef jura_PhaseScope_h
#define jura_PhaseScope_h
  
//=================================================================================================

/** Implements a phasescope analyzer. */

class JUCE_API PhaseScope : public jura::AudioModule, public jura::ImageUpdater
{

public:

  PhaseScope(CriticalSection *lockToUse);

  virtual ~PhaseScope();

  /** Creates the parameters to control the drawing. 
  \todo: maybe make this function a virtual member function of the AudioModule baseclass to be 
  overriden - the call it from the constructor there. */
  virtual void createParameters();

  /** Sets the desired pixel size for the display. This function will update the internal pixel 
  size of the buffer according to the passed values and the desired rescaling factor. " */
  void setDisplayPixelSize(int width, int height);

  // parameter setup functions (to be used for the callbacks from the parameters):
  void setBrightness(double newBrightness);
  void setAfterGlow(double newGlow);  // rename to setPixelDecayTime
  void setLineDensity(double newDensity);
  void setDotLimit(double newLimit);
  void setPixelSpread(double newSpread);
  void setPixelScale(double newFactor);
  void setAntiAlias(bool shouldAntiAlias);
  void setFrameRate(double newRate);
  void setScaleX(double newScale);
  void setScaleY(double newScale);
  //void setShearX(double newShear);
  //void setShearY(double newShear);
  void setRotation(double degrees);
  //void setShiftX(double newShift);
  //void setShiftY(double newShift);
  void setOneDimensionalMode(bool shouldBe1D);

  void setScanningFrequency(double newFrequency);
  void setNumCyclesShown(int newNumCycles);
  //void setZoom(double newZoom);
  void setSyncMode(bool shouldSync);
    // todo: maybe have the function bodies here in the header file - they are all trivial 
    // delegations
   

  // inquiry functions:
  inline double getFrameRate() { return phaseScopeBuffer->getFrameRate(); }
  inline int getInternalPixelWidth() { return phaseScopeBuffer->getWidth(); }
  inline int getInternalPixelHeight() { return phaseScopeBuffer->getHeight(); }
  const ColorMap& getColorMap() { return colorMap; }
  LoadableColorMap* getColorMapPointer() { return &colorMap; }


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

  bool bypassPixelDecay;        // not needed anymore - Elan wanted this at some point

  // this object is reponsible for drawing the incoming data onto a virtual screen:
  //RAPT::PhaseScopeBuffer<double, float, double> phaseScopeBuffer;
  //RAPT::PhaseScopeBuffer2<double, float, double> phaseScopeBuffer;
  RAPT::rsPhaseScopeBuffer<double, float, double> *phaseScopeBuffer;
    // maybe declare a pointer to a PhaseScopeBuffer here and in subclass PhaseScope2 declare a 
    // pointer to PahseScopeBuffer2 and in the constructor re-assign the inherited pointer - so we
    // don't have the overhead of PhaseScopeBuffer2 here.

  juce::Image image;       // image for the display
  //jura::ColorMap colorMap; // the color map to translate the buffered data matrix to colors
  jura::LoadableColorMap colorMap; // the color map to translate the buffered data matrix to colors

  friend class PhaseScopeDisplay;
  //friend class PhaseScopeEditor2;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope)
};

//=================================================================================================

/** Implements the GUI display for the phase scope. 

\todo BUG: there's a flickering of the lines...could this be related to threading issues?
..it's more obvious with faster decay times - maybe it's because the decay is applied in the GUI 
thread whereas accumulation is done in the audio-thread? instead of calling applyPixelDecay in the 
GUI thread we could set a flag in the phaseScope audio module and the apply the decay there
-> done - this seems to help indeed but also seems to introduce tearing artifacts 

maybe we should have a version of it which uses OpenGL - here's a tutorial:
http://duriansoftware.com/joe/An-intro-to-modern-OpenGL.-Table-of-Contents.html
*/

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

class JUCE_API PhaseScopeEditor : public AudioModuleEditor, public ParameterObserver
{

public:

  PhaseScopeEditor(jura::PhaseScope *newPhaseScopeToEdit);
  virtual ~PhaseScopeEditor();

  virtual void createWidgets();
  virtual void resized() override;
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

protected:

  PhaseScope *scope;
  PhaseScopeDisplay display;
  int widgetMargin;

  // Widgets:
  RSlider *sliderBrightness, *sliderAfterglow, *sliderPixelSpread, *sliderPixelScale, 
    *sliderLineDensity, *sliderDotLimit, *sliderFrameRate;

  RButton *buttonAntiAlias, *button1D, *buttonSync;

  AutomatableSlider *sliderScaleX, *sliderScaleY, *sliderShearX, *sliderShearY,
    *sliderRotation, *sliderShiftX, *sliderShiftY,
    *sliderScanFreq, *sliderNumCycles; // , *sliderZoom;
  // use regular (non-automatable) sliders



  ColorMapLoader* colorMapLoader;

  //RButtonPainter3D buttonPainter; // only temporary, experimental
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeEditor)
};



////=================================================================================================
//
///** Extends the basic PhaseScope by some more artistic features such as a customizable dot, 
//blurring, etc. 
//todo: rename to PrettyScope
//
//Give me these controls
//1. min-pixel-distance between dots
//2. maximum draw calls per frame (one control that limits both dots and lines)
//3. Let's start creating controls that are resolution-dependent. OH Actually if all size controls scale with resolution I don't have to limit the scope size, because then CPU usage is tied to resolution which is very good. So, give me percentages of resolution for dot size and line size rather than pixel size.
//4. Give me a master brightness control (so I can control dot/line brightness simultaneously)
//
//*/
//
//class JUCE_API PhaseScope2 : public jura::PhaseScope
//{
//
//public:
//
//  PhaseScope2(CriticalSection *lockToUse);
//
//  // additional setup functions:
//  void setPixelDecayByValue(double newDecayByValue);
//  void setPixelDecayByAverage(double newDecayByAverage);
//  // todo: add smear/blur-functions setLeft/Right/Up/DownSmear...
//
//  void setDrawDots(bool shouldDraw);
//  void setUseBigDot(bool shouldUseBigDot);
//  void setDotSize(double newSize);
//  void setDotBlur(double newBlur);
//  void setDotInnerSlope(double newSlope);
//  void setDotOuterSlope(double newSlope);
//
//  void setDrawLines(bool shouldDraw);
//  void setLineBrightness(double newBrightness);
//  void setLineWidth(double newWidth);
//  void setLineProfile(int newProfile);
//
//
//  // overriden from PhaseScope baseclass:
//  virtual void createParameters() override;
//  virtual AudioModuleEditor *createEditor() override;
//
//protected:
//
//  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScope2)
//};
//
//
///** GUI editor for the extended PhaseScope */
//
//class JUCE_API PhaseScopeEditor2 : public PhaseScopeEditor, public RSliderListener
//{
//
//public:
//
//  PhaseScopeEditor2(jura::PhaseScope2 *newPhaseScopeToEdit);
//
//
//  virtual void createWidgets() override;
//  virtual void resized() override;
//  virtual void paint(Graphics& g) override;
//
//  virtual void rSliderValueChanged(RSlider* slider) override;
//
//protected:
//
//  /** Updates the image for previewing the dot. */
//  void updatePreviewDot();
//
//  // additional widgets:
//  RSlider *sliderDecayByValue, *sliderDecayByAverage;
//  RButton *buttonBigDot, *buttonDrawDots, *buttonDrawLines;
//  RSlider *sliderDotSize, *sliderDotBlur, *sliderDotInnerSlope, *sliderDotOuterSlope,
//    *sliderLineBrightness, *sliderLineWidth;
//  RComboBox *boxLineProfile;
//
//  // image for previewing the dot:
//  RAPT::AlphaMask<float> dotPreviewMask;          
//  juce::Image dotPreviewImage; 
//    // add a preview for the line-profile, too
//
//
//  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseScopeEditor2)
//};

#endif 