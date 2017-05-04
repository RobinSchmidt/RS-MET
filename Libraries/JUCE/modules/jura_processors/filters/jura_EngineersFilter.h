#ifndef jura_EngineersFilter_h
#define jura_EngineersFilter_h


/** This class wraps rosic::EngineersFilter into a rosof::AudioModule to facilitate its use as 
plugIn. */

class EngineersFilterAudioModule : public AudioModule
{

  friend class EngineersFilterModuleEditor;

public:

  /** Constructor to use when you want to wrap an existing rosic::EngineersFilter object without 
  this AudioModule taking ownership (suitable, if the object already exists a smember of some 
  higher level dsp object). */
  EngineersFilterAudioModule(CriticalSection *newPlugInLock, rosic::EngineersFilter *sciFilterToWrap);

  /** Constructor to use ehen there's no existing rosic::EngineersFilter object to be wrapped. in 
  this case, we'll create one here and take over ownership  (i.e. will also delete it in our 
  destructor). */
  EngineersFilterAudioModule(CriticalSection *newPlugInLock);

  virtual ~EngineersFilterAudioModule();






  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedEngineersFilter->setSampleRate(newSampleRate);
  }

  // obsolete?
  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    for(int n = 0; n < numSamples; n++)
    {
      double dL = left[n];
      double dR = right[n];
      wrappedEngineersFilter->getSampleFrameDirect1(&dL, &dR);
      left[n]  = (float)dL;
      right[n] = (float)dR;
    }
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedEngineersFilter->getSampleFrameDirect1(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

  virtual void reset()
  {
    wrappedEngineersFilter->reset();
  }


protected:

  void initializeAutomatableParameters();

  rosic::EngineersFilter *wrappedEngineersFilter;

  bool wrappedEngineersFilterIsOwned = false;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class plots the frequency responses of a rosic::EngineersFilter object and allows for 
editing parameters like the cutoff frequencies. */

class EngineersFilterPlotEditor	: virtual public SpectrumDisplayOld, public ChangeBroadcaster 
  //, public ParameterObserver, 
{

  /** Enumeration of the handles that can be grabbed and dragged by the mouse.  */
  enum dragHandles
  {
    NONE = 0,
    LOW_MID,
    MID_HIGH
  };

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  EngineersFilterPlotEditor(const juce::String& name = "EngineersFilterPlotEditor");   

  /** Destructor. */
  virtual ~EngineersFilterPlotEditor(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::EngineersFilter object which is to be edited. */
  virtual void setEngineersFilterToEdit(rosic::EngineersFilter* newEngineersFilterToEdit);

  /** Assigns a rojue::Parameter object to the cutoff frequency of the lowpass for observation and 
  manipulation. */
  //virtual void assignParameterLowFreq(Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to the slope of the lowpass for observation and 
  manipulation. */
  //virtual void assignParameterLowSlope(Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to the cutoff frequency of the highpass for observation and 
  manipulation. */
  //virtual void assignParameterHighFreq(Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to the slope of the highpass for observation and 
  manipulation. */
  //virtual void assignParameterHighSlope(Parameter* parameterToAssign);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  //virtual void parameterChanged(Parameter* parameterThatHasChanged);
  //virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted);

  /** This method is called when one of the assigned rosic::AutomatableParameters has been changed. 
  We override it here in the subclass to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendMessage = false);

  /** Overrides the changeListetnerCcallback in order to receive messages which this object sends
  to itself. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Overrides mouseMove in order to update the cursor according to what is under the mouse. */
  virtual void mouseMove(const MouseEvent &e);

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop up 
  when the right button is clicked for 
  MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e);

  /** Overrides mouseUp to reset the currentDragHandle to NONE. */
  virtual void mouseUp(const MouseEvent& e);

  /** Overrides the resized-method. */
  virtual void resized();

  /** Updates the frequency response plot. */
  virtual void updatePlot();


protected:

  /** Returns the handle for mouse grab/drag under the specified position (in pixels) as one of 
  the values in enum dragHandles. */
  virtual int getDragHandleAt(int x, int y);

  /** Does the setup of the filter according to some new mouse position) */
  virtual void setupFilterAccordingToMousePosition(double mouseX, double mouseY);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handles. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL, 
    XmlElement *targetSVG = NULL);

  /** Pointer to the actual rosic::EngineersFilter object which is being edited. */
  rosic::EngineersFilter* sciFilterToEdit;

  // the parameters which wil cause re-plotting and therefore must be listened to:
  //Parameter *lowFreqParameter, *lowSlopeParameter, *highFreqParameter, *highSlopeParameter;

  // magnitude response display stuff:
  int    numBins;
  double *frequencies, *magnitudes;
  //double *lowpassMagnitudes, *bandpassMagnitudes, *highpassMagnitudes;
  //double **allMagnitudes;

  int currentlyDraggedHandle;


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Editor for EngineersFilterAudioModule */

class EngineersFilterModuleEditor : public AudioModuleEditor, public RComboBoxObserver, 
  public RSliderListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  EngineersFilterModuleEditor(CriticalSection *newPlugInLock, 
    EngineersFilterAudioModule* newEngineersFilterAudioModule);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();


protected:

  /** Makes currently required widgets visible and currently not required widgets invisible. */
  virtual void updateWidgetVisibility();

  EngineersFilterAudioModule *sciFilterModuleToEdit;

  EngineersFilterPlotEditor  *plotEditor;

  // use RNamedComboBox:
  //RLabel    *modeLabel, *methodLabel;
  //RComboBox *modeComboBox, *methodComboBox;

  RNamedComboBox *modeComboBox, *methodComboBox;


  RSlider   *frequencySlider, *orderSlider, *bandwidthSlider, *gainSlider, *rippleSlider,
    *rejectionSlider, *inRippleSlider, *outRippleSlider;

  // *transitionWidthSlider, *passbandRippleSlider, *stopBandAttenuationSlider, *stopBandRippleSlider;
  /*
  RLabel    *bandsLabel, *lowLabel, *highLabel;
  RComboBox *numBandsComboBox;
  RButton   *monoButton;
  RSlider   *lowFrequencySlider, *lowSlopeSlider, *highFrequencySlider, *highSlopeSlider;
  */

  juce_UseDebuggingNewOperator;
};

#endif 
