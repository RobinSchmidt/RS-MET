#ifndef jura_EngineersFilter_h
#define jura_EngineersFilter_h

// to quickly switch between old and new implementation during debugging - get rid later
typedef rosic::rsEngineersFilterStereo EFLT; // elliptic bandpass (among others) doesn't work
//typedef rosic::rsEngineersFilterOld EFLT;

/** This class wraps rosic::EngineersFilter into a rosof::AudioModule to facilitate its use as 
plugIn. */

class EngineersFilterAudioModule : public AudioModule
{

  friend class EngineersFilterModuleEditor;


public:

  /** Constructor to use when you want to wrap an existing rosic::EngineersFilter object without 
  this AudioModule taking ownership (suitable, if the object already exists a smember of some 
  higher level dsp object). */
  EngineersFilterAudioModule(CriticalSection *newPlugInLock, EFLT *sciFilterToWrap);

  /** Constructor to use ehen there's no existing rosic::EngineersFilter object to be wrapped. in 
  this case, we'll create one here and take over ownership  (i.e. will also delete it in our 
  destructor). */
  EngineersFilterAudioModule(CriticalSection *newPlugInLock);

  void init();

  virtual ~EngineersFilterAudioModule();

  AudioModuleEditor* createEditor(int type) override;

  void setSampleRate(double newSampleRate) override
  {
    wrappedEngineersFilter->setSampleRate(newSampleRate);
  }

  void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedEngineersFilter->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
      //wrappedEngineersFilter->getSampleFrameDirect1(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

  void processStereoFrame(double* left, double* right) override
  {
    wrappedEngineersFilter->getSampleFrameStereo(left, right);
  }

  void reset() override
  {
    wrappedEngineersFilter->reset();
  }


protected:

  void createParameters();

  EFLT *wrappedEngineersFilter;
  bool wrappedEngineersFilterIsOwned = false;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(EngineersFilterAudioModule)
};

//=================================================================================================

/** This class plots the frequency responses of a rosic::EngineersFilter object and allows for 
editing parameters like the cutoff frequencies. */

class EngineersFilterPlotEditor	: virtual public rsSpectrumPlot, public ChangeBroadcaster 
  , public ParameterObserver
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
  virtual void setEngineersFilterToEdit(EFLT* newEngineersFilterToEdit);

  // functions to assign the parameters, we observe (in order to update the plot) and possibly 
  // later also to manipulated them:
  void assign(Parameter*& target, Parameter* p) { target = p; p->registerParameterObserver(this);} 
  void assignParameterMode(     Parameter* p) { assign(modeParam,      p); }
  void assignParameterMethod(   Parameter* p) { assign(methodParam,    p); }
  void assignParameterOrder(    Parameter* p) { assign(orderParam,     p); }
  void assignParameterFrequency(Parameter* p) { assign(freqParam,      p); }
  void assignParameterBandwidth(Parameter* p) { assign(bandwidthParam, p); }
  void assignParameterGain(     Parameter* p) { assign(gainParam,      p); }
  void assignParameterRipple(   Parameter* p) { assign(rippleParam,    p); }
  void assignParameterRejection(Parameter* p) { assign(rejectionParam, p); }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void parameterChanged(Parameter* param) override;
  //virtual void parameterChanged(Parameter* parameterThatHasChanged);
  //virtual void parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted);

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
  //rosic::rsEngineersFilterStereo* sciFilterToEdit;
  EFLT* sciFilterToEdit;

  // the parameters which wil cause re-plotting and therefore must be listened to:
  //Parameter *lowFreqParameter, *lowSlopeParameter, *highFreqParameter, *highSlopeParameter;
  Parameter *modeParam, *methodParam, *orderParam, *freqParam, *bandwidthParam, *gainParam, 
    *rippleParam, *rejectionParam;

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

class EngineersFilterModuleEditor : public AudioModuleEditor, public ParameterObserver
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  EngineersFilterModuleEditor(CriticalSection *newPlugInLock, 
    EngineersFilterAudioModule* newEngineersFilterAudioModule);
  // todo: get rid of newPlugInLock parameter

  virtual ~EngineersFilterModuleEditor();

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void resized() override;
  virtual void parameterChanged(Parameter* param) override;
  virtual void updateWidgetsAccordingToState() override;


protected:

  virtual void createWidgets();

  /** Makes currently required widgets visible and currently not required widgets invisible. */
  virtual void updateWidgetVisibility();

  EngineersFilterAudioModule *sciFilterModuleToEdit;
  EngineersFilterPlotEditor  *plotEditor;

  RNamedComboBox *modeComboBox, *methodComboBox;

  RSlider   *frequencySlider, *orderSlider, *bandwidthSlider, *gainSlider, *rippleSlider,
    *rejectionSlider;

  juce_UseDebuggingNewOperator;
};

#endif 
