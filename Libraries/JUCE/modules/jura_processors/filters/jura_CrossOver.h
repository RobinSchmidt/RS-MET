#ifndef jura_CrossOver_h
#define jura_CrossOver_h

/** This class wraps rosic::CrossOver into a rosof::AudioModule to facilitate its use as plugIn. */

class CrossOverAudioModule : public AudioModule
{

  friend class CrossOverPlotEditor;
  friend class CrossOverModuleEditor;

public:

  CrossOverAudioModule(CriticalSection *newPlugInLock, rosic::CrossOver4Way *crossOverToWrap = nullptr);

  virtual ~CrossOverAudioModule();

  AudioModuleEditor* createEditor() override;


  //virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setSampleRate(double newSampleRate)
  {
    ScopedLock scopedLock(*lock);
    wrappedCrossOver->setSampleRate(newSampleRate);
  }

  /** Overriden to deal with the multichannel stuff. */
  virtual void processBlock(AudioSampleBuffer& buffer, int startSample, int length)
  {
    ScopedLock scopedLock(*lock);

    if(buffer.getNumChannels() != 8)
      return;
    else
    {
      //wrappedCrossOver->processBuffer(buffer.getSampleData(0, startSample), length);
        // fast but works only for buffers arranged continuously in memory...

      /*
      int c, n;
      double sampleFrame[8];

      float *pointers[8];
      for(int c=0; c<8; c++)
        pointers[c] = buffer.getSampleData(c, 0);

      for(n=0; n<length; n++)
      {
        sampleFrame[0] = pointers[0][n];
        sampleFrame[1] = pointers[1][n];
        wrappedCrossOver->processSampleFrame(sampleFrame);
        for(c=0; c<8; c++)
          pointers[c][n] = (float) sampleFrame[c];
      }
      */
      float *pointers[8];
      for(int c=0; c<8; c++)
        pointers[c] = buffer.getWritePointer(c, startSample);
      wrappedCrossOver->processBuffer(pointers, length);
    }
  }

  virtual void reset()
  {
    ScopedLock scopedLock(*lock);
    wrappedCrossOver->resetBuffers();
  }

protected:

  void createStaticParameters();

  rosic::CrossOver4Way *wrappedCrossOver;
  bool wrappedCrossOverIsOwned = false;

  // adapter functions for the callbacks (boilerplate):
  void setBandActive_0_0(bool shouldBeActive)        { wrappedCrossOver->setBandActive(shouldBeActive, 0, 0); }
  void setCrossoverFrequency_0_0(double newFrequency) { wrappedCrossOver->setCrossoverFrequency(newFrequency, 0, 0); }
  void setSlope_0_0(int newSlope)                     { wrappedCrossOver->setSlope(newSlope, 0, 0); }

  void setBandActive_1_0(bool shouldBeActive)        { wrappedCrossOver->setBandActive(shouldBeActive, 1, 0); }
  void setCrossoverFrequency_1_0(double newFrequency) { wrappedCrossOver->setCrossoverFrequency(newFrequency, 1, 0); }
  void setSlope_1_0(int newSlope)                     { wrappedCrossOver->setSlope(newSlope, 1, 0); }

  void setBandActive_1_1(bool shouldBeActive)        { wrappedCrossOver->setBandActive(shouldBeActive, 1, 1); }
  void setCrossoverFrequency_1_1(double newFrequency) { wrappedCrossOver->setCrossoverFrequency(newFrequency, 1, 1); }
  void setSlope_1_1(int newSlope)                     { wrappedCrossOver->setSlope(newSlope, 1, 1); }


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/**

This class plots the frequency responses of a rosic::CrossOver object and allows for
editing parameters like the cutoff frequencies.

*/

class CrossOverPlotEditor	: virtual public SpectrumDisplayOld, public ParameterObserver, public ChangeBroadcaster
{

  /** Enumeration of the handles that can be grabbed and dragged by the mouse.  */
  enum dragHandles
  {
    NONE = 0,
    FREQUENCY_1_1,
    FREQUENCY_2_1,
    FREQUENCY_2_2
  };

public:

  //-------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  //CrossOverPlotEditor(const juce::String& name = juce::String(T("CrossOverPlotEditor")));
  CrossOverPlotEditor(CriticalSection *newPlugInLock, CrossOverAudioModule* newCrossOverModuleToEdit);

  /** Destructor. */
  virtual ~CrossOverPlotEditor();

  //-------------------------------------------------------------------------------------------------------------------------------------
  // parameter handling: \todo: should we ever use this class when dynamically assigning CrossOverAudioModule objects to this editor,
  // the functions below will need more sophisticated implementations - we must take care of always checking against NULL and stuff - see
  // EqualizerPlotEditor for reference

  /** Passes a pointer the the actual rosic::CrossOver object which is to be edited. */
  //virtual void setCrossOverToEdit(rosic::CrossOver4Way* newCrossOverToEdit);

  /** Assigns a rojue::Parameter object to one of the on/off switches for observation and manipulation. */
  virtual void assignParameterOnOff(int treeLevel, int indexInLevel, Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to one of the crossover frequencies for observation and manipulation. */
  virtual void assignParameterFreq(int treeLevel, int indexInLevel, Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to one of the slopes for observation and manipulation. */
  virtual void assignParameterSlope(int treeLevel, int indexInLevel, Parameter* parameterToAssign);

  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted) {}

  //-------------------------------------------------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the currently selected crossover via it's level inside the tree and its index inside the tree-level. Each level (start
  counting at 0) has 2^level indices, so the 1st level has only 1 index, the 2nd has 2, the 3rd has 4, etc. (remark: currently only 2
  levels are present but it makes sense to already use this structure to support later extensions). It will assign the
  reference-parameters to -1 if none is selected. */
  virtual void getSelectedTreeLevelAndIndex(int &treeLevel, int &indexInLevel);

  //-------------------------------------------------------------------------------------------------------------------------------------
  // callbacks:

  /** This method is called when one of the assigned rosic::AutomatableParameters has been changed - we override it here in the subclass
  to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendMessage = false);

  /** Overrides the changeListetnerCcallback in order to receive messages which this object sends to itself. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Overrides mouseMove in order to update the cursor according to what is under the mouse. */
  virtual void mouseMove(const MouseEvent &e);

  /** Overrides mouseDown for assigning currentDragHandle. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseDrag for adjusting one of the crossover frequencies. */
  virtual void mouseDrag(const MouseEvent& e);

  /** Overrides mouseUp to reset the currentDragHandle to NONE. */
  virtual void mouseUp(const MouseEvent& e);

  /** Overrides the resized-method. */
  virtual void resized();

  /** Updates the frequency response plot. */
  virtual void updatePlot();



protected:

  /** Returns the handle for mouse grab/drag under the specified position (in pixels) as one of the values in enum dragHandles. */
  virtual int getDragHandleAt(int x, int y);

  /** Does the setup of the filter according to some new mouse position) */
  virtual void setupFilterAccordingToMousePosition(double mouseX, double mouseY);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handles. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL, XmlElement *targetSVG = NULL);

  /** Colourizes the background with colours that represent the center frequencies of each band. */
  virtual void colourizeBackground(Graphics &g, juce::Image *targetImage);

  /** Draws the indicator that marks one of the crossovers as selected. */
  virtual void drawSelectionIndicator(Graphics &g, juce::Image *targetImage);

  /** Draws a vertical line (which serves as handle to adjust one of the crossover frequencies) at the given frequency. If the 'active'
  flag is false, the line will be drawn dashed to indicate that the respective crossover is turned off. */
  virtual void drawVerticalLineAtFrequency(Graphics &g, juce::Image* targetImage, double frequency, int treeLevel, bool active);

  /** Draws a triangle (which serves as handle to switch one of the crossovers on/off) at the given frequency. If the 'active' flag is
  false, the triangle will be drawn non-filled to indicate that the respective crossover is turned off. */
  virtual void drawTriangleSwitchAtFrequency(Graphics &g, juce::Image* targetImage, double frequency, int treeLevel, bool active);

  /** Restricts the x-coordinate (in pixels) such that drag-handle with given index can't be pulled through other drag-handles. */
  virtual double restrictDragHandleX(double x, int dragHandleIndex);

  /** Returns the x-coordinate (in pixels) of one of the drag-handles. @see dragHandles */
  virtual double getDragHandleX(int dragHandleIndex);

  /** Returns the colour that should be used for the low band's background. */
  virtual Colour getLowBandColour();
  virtual Colour getLowMidBandColour();    ///< @see getLowBandColour
  virtual Colour getMidBandColour();       ///< @see getLowBandColour
  virtual Colour getHighMidBandColour();   ///< @see getLowBandColour
  virtual Colour getHighBandColour();      ///< @see getLowBandColour


  CriticalSection      *plugInLock;              // mutex to access the edited AudioModule object
  CrossOverAudioModule *crossOverModuleToEdit;   // the edited AudioModule object

  // the parameters which wil cause re-plotting and therefore must be listened to:
  Parameter *onOff21Parameter, *onOff22Parameter, *freq11Parameter, *freq21Parameter, *freq22Parameter,
    *slope11Parameter, *slope21Parameter, *slope22Parameter;

  // magnitude response display stuff:
  int    numBins;
  double *frequencies, *magnitudes1, *magnitudes2, *magnitudes3, *magnitudes4;
  double **allMagnitudes;

  int selectedLevel, selectedIndex;
  int currentlyDraggedHandle;


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class CrossOverModuleEditor : public AudioModuleEditor //, public RComboBoxObserver
{

public:

  CrossOverModuleEditor(CriticalSection *newPlugInLock, CrossOverAudioModule* newCrossOverAudioModule);

  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  /** Makes currently required widgets visible and currently not required widgets invisible. */
  virtual void updateWidgetVisibility();

  CrossOverAudioModule *crossOverModuleToEdit;
  CrossOverPlotEditor  *plotEditor;

  // 1st index: level in tree, 2nd index: index within level. 4 bands -> 3 crossover-frequencies, will later become 7 (for 8 bands)
  // ...we may want to use a Tree class then.
  RSlider *frequency11Slider, *frequency21Slider, *frequency22Slider;
  RSlider *slope11Slider, *slope21Slider, *slope22Slider;
  RButton *monoButton;


  juce_UseDebuggingNewOperator;
};

#endif
