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

#endif 
