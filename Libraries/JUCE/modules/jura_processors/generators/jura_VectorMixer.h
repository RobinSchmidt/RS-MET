#ifndef jura_VectorMixer_h
#define jura_VectorMixer_h

class VectorMixerAudioModule : public AudioModule
{

  friend class VectorMixerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  VectorMixerAudioModule(CriticalSection *newPlugInLock, 
    rosic::VectorMixer *newVectorMixerToWrap);

  //---------------------------------------------------------------------------------------------
  // automation:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  /*
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);
  */

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutL = *inOutR = 0.0;
  }


protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::VectorMixer *wrappedVectorMixer;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/**

This class implements a user interface XY-pad for a rosic::VectorMixer object.

\todo turn this into a general purpose widget with callback and obesrver
callback: xyPadChanged(double newX, double newY, bool xWasChanged, bool yWasChanged);

*/

class VectorMixerPad	: virtual public CoordinateSystemOld, public ParameterObserver, 
  public ChangeBroadcaster
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  VectorMixerPad(rosic::VectorMixer* newVectorMixerToEdit, 
    const juce::String& name = juce::String("XYPad")  );   

  /** Destructor. */
  virtual ~VectorMixerPad(); 

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::MoogyFilter object which is to be edited. Make 
  sure to call this function again with a NULL-pointer when the object get deleted for some 
  reason. */
  virtual void setVectorMixerToEdit(rosic::VectorMixer* newVectorMixerToEdit);

  /** Assigns a rojue::Parameter object to the horizontal axis for observation and 
  manipulation. */
  virtual void assignParameterX(Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to the vertical axis for observation and 
  manipulation. */
  virtual void assignParameterY(Parameter* parameterToAssign);

  /** Un-Assigns a previously rojue::Parameter object for the horizontal axis. */
  virtual void unAssignParameterX();

  /** Un-Assigns a previously rojue::Parameter object for the vertical axis. */
  virtual void unAssignParameterY();

  /** This method is called when one of the assigned rojue::Parameters has been changed
  - we override it here in the subclass to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendMessage = false);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted) {}

  /** Overrides the changeListetnerCcallback in order to receive messages which this object sends
  to itself. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop 
  up when the right button is clicked for MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e);

  /** Overrides the resized-method. */
  virtual void resized();

  /** Updates the vector plot. */
  //virtual void updatePlot();

  /** Overrides the inherited method from the CoordinateSystem base-class. */
  virtual void drawCoordinateSystem(Graphics &g, juce::Image* targetImage = NULL, 
    XmlElement* targetSVG = NULL);

protected:

  /** Does the setup of the VectorMixer according to some new mouse position. */
  virtual void setupVectorMixerAccordingToMousePosition(double mouseX, double mouseY);

  /** Pointer to the actual rosic::VectorMixer object which is being edited. */
  rosic::VectorMixer* vectorMixerToEdit;

  /** Radius of the dot-handle to be drawn. */
  float dotRadius = 8.f;

  Colour dotColour;

  // the parameters which will cause re-plotting and therefore must be listened to:
  Parameter* xParameter;
  Parameter* yParameter;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

//class VectorMixerModuleEditor : virtual public AudioModuleEditor, virtual public VectorMixerPad
class VectorMixerModuleEditor : public AudioModuleEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. You must pass a valid (non NULL) pointer to an VectorMixerAudioModule 
  object which will be accessed by this editor. */
  VectorMixerModuleEditor(CriticalSection *newPlugInLock, 
    VectorMixerAudioModule* newVectorMixerAudioModule);

  //---------------------------------------------------------------------------------------------
  // overrides to avoid ambigous access errors and dominance warnings:
  /*
  virtual void mouseDown(       const MouseEvent &e);
  virtual void mouseDrag(       const MouseEvent &e);
  virtual void mouseMove(       const MouseEvent &e);
  virtual void mouseEnter(      const MouseEvent &e);
  virtual void mouseExit(       const MouseEvent &e);
  virtual void mouseUp(         const MouseEvent &e);
  virtual void mouseDoubleClick(const MouseEvent &e);
  virtual void mouseWheelMove(  const MouseEvent &e, 
  float wheelIncrementX, float wheelIncrementY);
  virtual void paint(Graphics &g);
  */
  virtual void resized();

  //virtual void updateWidgetsAccordingToState();
  //virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& name);
  //virtual XmlElement* getStateAsXml(const juce::String& stateName = juce::String(T("VectorPadState"))) const;

protected:

  VectorMixerPad *vectorPad;

  juce_UseDebuggingNewOperator;
};

#endif 
