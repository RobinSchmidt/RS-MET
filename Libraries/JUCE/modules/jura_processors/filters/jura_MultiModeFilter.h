#ifndef jura_MultiModeFilter_h
#define jura_MultiModeFilter_h

class MultiModeFilterAudioModule : public AudioModule
{

  friend class MultiModeFilterModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultiModeFilterAudioModule(CriticalSection *newPlugInLock, 
    rosic::MultiModeFilter *newMultiModeFilterToWrap);

  //---------------------------------------------------------------------------------------------
  // overrides:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedMultiModeFilter->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    jassertfalse; // not yet implemented
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::MultiModeFilter *wrappedMultiModeFilter;


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class plots the frequency response of a rosic::MultiModeFilter object and allowas for 
editing parameters like th cutoff frequency and resonance by dragging around some node-point.

\todo: implement new locking strategy and parameter handling */

class MultiModeFreqResponseEditor	: virtual public SpectrumDisplayOld, public ParameterObserver, 
  public ChangeBroadcaster
{

public:

  //-------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultiModeFreqResponseEditor(const juce::String& name = juce::String("MultiModeFreqResponseEditor"));   

  /** Destructor. */
  virtual ~MultiModeFreqResponseEditor(); 

  //-------------------------------------------------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::MoogyFilter object which is to be edited. Make 
  sure to call this function again with a NULL-pointer when the object get deleted for some 
  reason. */
  virtual void setFilterToEdit(rosic::MultiModeFilter* newFilterToEdit);

  /** Assigns a Parameter object to the frequency (horizontal axis) for observation and manipulation. */
  virtual void assignParameterFreq(Parameter* parameterToAssign);

  /** Assigns a Parameter object to the resonance (vertical axis) for observation and manipulation. */
  virtual void assignParameterReso(Parameter* parameterToAssign);

  /** Assigns a Parameter object to the resonance (vertical axis) for observation and manipulation. */
  virtual void assignParameterQ(Parameter* parameterToAssign);

  /** Assigns a Parameter object for the filter gain for observation. */
  virtual void assignParameterGain(Parameter* parameterToAssign);

  /** Assigns a Parameter object for the filter morph for observation. */
  virtual void assignParameterMorph(Parameter* parameterToAssign);

  /** Un-Assigns a previously Parameter object to the horizontal axis. */
  virtual void unAssignParameterFreq();

  /** Un-Assigns a previously Parameter object to the verical axis. */
  virtual void unAssignParameterReso();

  /** Un-Assigns a previously Parameter object to the verical axis. */
  virtual void unAssignParameterQ();

  /** Un-Assigns a previously Parameter object for the filter gain. */
  virtual void unAssignParameterGain();

  /** Un-Assigns a previously Parameter object for the filter morph. */
  virtual void unAssignParameterMorph();

  //-------------------------------------------------------------------------------------------------------------------------------------
  // callbacks:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted);

  /** This method is called when one of the assigned rosic::AutomatableParameters has been changed - we override it here in the subclass 
  to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendMessage = false);

  /** Overrides the changeListetnerCcallback in order to receive messages which this object sends to itself. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop up when the right button is clicked for 
  MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e);

  /** Overrides the resized-method. */
  virtual void resized();

  /** Updates the frequency response plot. */
  virtual void updatePlot();



protected:

  /** Does the setup of the filter according to some new mouse position) */
  virtual void setupFilterAccordingToMousePosition(double mouseX, double mouseY);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handle. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL, XmlElement *targetSVG = NULL);

  /** Converts a resonance value to an y-coordinate in components/image coordinates. */
  double resoToY(double reso, juce::Image *targetImage = NULL);

  /** Converts an y-coordinate in components/image coordinates to a resonance value. */
  double yToReso(double y, juce::Image *targetImage = NULL);

  /** Converts a Q-value to an y-coordinate in components/image coordinates. */
  double qToY(double q, juce::Image *targetImage = NULL);

  /** Converts an y-coordinate in components/image coordinates to a Q-value. */
  double yToQ(double y, juce::Image *targetImage = NULL);

  /** Radius of the dot-handle to be drawn. */
  float dotRadius;

  /** Pointer to the actual rosic::MultiModeFilter object which is being edited. */
  rosic::MultiModeFilter* filterToEdit;

  // the parameters which wil cause re-plotting and therefore must be listened to:
  Parameter* freqParameter;
  Parameter* resoParameter;
  Parameter* qParameter;
  Parameter* gainParameter;
  Parameter* morphParameter;

  // magnitude response display stuff:
  int    numBins;
  double *frequencies, *magnitudes;

  
  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class MultiModeFilterModuleEditor : public AudioModuleEditor, public RComboBoxObserver, 
  public RSliderListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MultiModeFilterModuleEditor(CriticalSection *newPlugInLock, 
    MultiModeFilterAudioModule* newMultiModeFilterAudioModule);

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::MoogyFilter object which is to be edited. Make 
  sure to call this function again with a NULL-pointer when the object get deleted for some 
  reason. */
  //virtual void setFilterToEdit(rosic::MultiModeFilter* newFilterToEdit);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  //virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

  //---------------------------------------------------------------------------------------------
  // public data members:  \todo: move to protected

  // the widgets:
  RSlider *freqSlider, *freqByKeySlider, *freqByVelSlider, *resoSlider, *qSlider, 
    *driveSlider, *orderSlider, *gainSlider, *morphSlider, *transitionSlider, 
    *preAllpassSlider, *makeUpSlider;
  //*freq2ScaleSlider, *freq2OffsetSlider, *q2ScaleSlider, *gain2ScaleSlider;

  //RLabel*    modeLabel;
  //RComboBox* modeComboBox;
  RNamedComboBox *modeComboBox;
  RButton        *twoStagesButton;

  MultiModeFreqResponseEditor *frequencyResponseDisplay;

protected:

  /** Updates the arrangement of the widgets according to the chosen filter mode. */
  virtual void updateWidgetArrangement();

  /** Arranges the widgets that are common for all filter types (i.e. the mode-combo-box, 
  the frequency slider and it's key- and vel-slider. */
  //virtual void arrangeCommonWidgets();

  /** Arranges the widgets for Moogish Lowpass mode. */
  virtual void arrangeWidgetsForMoogishLowpassMode();

  /** Arranges the widgets for the first order filter modes which don't have a gain parameter, 
  these are: Lowpass 6 dB/oct, Highpass 6 dB/oct, Allpass 1st order. */
  virtual void arrangeWidgetsForFirstOrderWithoutGain();

  /** Arranges the widgets for the second order filter modes which don't have a gain parameter, 
  these are: Lowpass 12 dB/oct, Highpass 12 dB/oct, Bandpass 2*6 dB/oct, Bandstop 2*6 dB/oct, 
  Allpass 2nd order. */
  virtual void arrangeWidgetsForSecondOrderWithoutGain();

  /** Arranges the widgets for the first order filter modes which have a gain parameter, 
  these are: Low Shelv 1st order, High Shelv 1st order. */
  virtual void arrangeWidgetsForFirstOrderWithGain();

  /** Arranges the widgets for the second order filter modes which have a gain parameter, 
  these are: Low Shelv 2nd order, High Shelv 2nd order, Peak/Dip. */
  virtual void arrangeWidgetsForSecondOrderWithGain();

  /** Arranges the widgets for morphable modes mode. */
  virtual void arrangeWidgetsForMorphableMode();

  /** Pointer to the actual rosic::MultiModeFilter object which is being edited. */
  rosic::MultiModeFilter* filterToEdit;

  juce_UseDebuggingNewOperator;
};

#endif 
