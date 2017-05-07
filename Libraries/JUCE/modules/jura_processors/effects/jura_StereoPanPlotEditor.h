#ifndef jura_StereoPanPlotEditor_h
#define jura_StereoPanPlotEditor_h

/** This class plots the curves for the pan law and lets the user adjust the pan-position. */

class StereoPanPlotEditor : virtual public CurveFamilyPlotOld, public ParameterObserver, 
  public ChangeBroadcaster
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  StereoPanPlotEditor(const juce::String& name = juce::String("StereoPanPlotEditor"));

  /** Destructor. */
  virtual ~StereoPanPlotEditor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::StereoPan object which is to be edited. */
  virtual void setStereoPanToEdit(rosic::StereoPan* newStereoPanToEdit);

  /** Assigns a rojue::Parameter object to the pan-parameter for observation and manipulation. */
  virtual void assignParameterPan(Parameter* parameterToAssign);

  /** Assigns a rojue::Parameter object to the pan-law parameter for observation. */
  virtual void assignParameterPanLaw(Parameter* parameterToAssign);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted);

  /** This method is called when one of the assigned rosic::Parameter has been changed
  - we override it here in the subclass to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendMessage = false);

  /** Overrides the changeListetnerCallback in order to receive messages which this object sends
  to itself. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Overrides the mouseEnter callback in order to show the description in the dedicated field
  when the mouse enters the widget. */
  //virtual void mouseEnter(const MouseEvent &e);

  /** Overrides mouseMove in order to update the cursor according to what is under the mouse. */
  //virtual void mouseMove(const MouseEvent &e);

  /** Overrides the mouseExit callback in order to make the description disappear when the mouse
  leaves the widget. */
  //virtual void mouseExit(const MouseEvent &e);

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop
  up when the right button is clicked for MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e);

  /** Overrides mouseUp to reset the currentDragHandle to NONE. */
  //virtual void mouseUp(const MouseEvent& e);

  /** Overrides the resized-method. */
  virtual void resized();

  /** Updates the frequency response plot. */
  virtual void updatePlot();


protected:

  /** Does the setup of the filter according to some new mouse position) */
  virtual void setupPanAccordingToMousePosition(double mouseX);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handles. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL,
    XmlElement *targetSVG = NULL);

  /** Pointer to the actual rosic::StereoPan object which is being edited. */
  rosic::StereoPan* stereoPanToEdit;

  // the parameters which wil cause re-plotting and therefore must be listened to:
  Parameter *panParameter, *panLawParameter;

  // plot stuff:
  int    numValues;
  double *p, *gLL, *gRL, *gLR, *gRR;
  double **allGains;

  juce_UseDebuggingNewOperator;
};

#endif  
