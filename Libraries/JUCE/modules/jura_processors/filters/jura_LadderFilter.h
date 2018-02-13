#ifndef jura_LadderFilter_h
#define jura_LadderFilter_h
  
//=================================================================================================

/** Wraps two RAPT LadderFilter instances into a single object, for stereo processing. In order to
make it more interesting than just instantiating the RAPT::LadderFilter template for a kind of 
stereo signal type (which would be also possible), we allow the two cutoff frequencies of both 
channels to be different by introducing a stereo spread parameter. 

\todo: 
-rename the class to LadderStereo
-make the GUI editor work for the basic ladder filter (not only for this stereo version)
 (maybe we have to make this a subclass of the underlying RAPT::LadderFilter)
-let the user freely adjust the mixing coefficients for the filters stages (via the DraggableNumber
 widget) - the mode-combobox then selects between various predefined settings
*/

class JUCE_API Ladder : public jura::ModulatableAudioModule
  // : public jura::AudioModuleWithMidiIn // why WithMidiIn?
{

public:

  Ladder(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr);
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createStaticParameters();
    // maybe rename to createParameters

  /** Creates the GUI editor (returns an object of an appropriate subclass of AudioModuleEditor) */
  AudioModuleEditor *createEditor() override;

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

  // target functions for callbacks that are called on parameter changes:
  void setCutoff(double newCutoff);
  void setResonance(double newResonance);
  void setMode(int newMode);
  void setStereoSpread(double newSpreadInSemitones);
  void setMidSideMode(bool shouldBeInMidSideMode);

  // for the magnitude plot/editor later...
  double getCutoff()    { return cutoff; }
  double getResonance() { return 0.0;    }  // preliminary

  double getMagnitudeAt(double frequency); // returns the magnitude at the given frequency

  void getMagnitudeResponse(const double *frequencies, double *magnitudes, int numBins, 
    bool inDecibels = true);

protected:

  // some members to store the parameter values:
  double cutoff      = 1000.0;
  double spread      = 0.0;
  double freqFactorL = 1.0;
  double freqFactorR = 1.0;
  bool   midSideMode = false;

  // embedded core DSP objects from the RAPT library:
  //rsLadderDD ladderL, ladderR;
  rosic::rsLadderFilterStereo wrappedLadder;
    // use the one with different cutoffs for the channels later (maybe)

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Ladder)
};

//=================================================================================================
// the magnitude response plot/editor:
// this is the old, deprecated version - now we use the new rsLadderPlotEditor - when transition is 
// finished, this old class can be deleted

class JUCE_API LadderSpectrumEditor : public rsSpectrumPlot, public ParameterObserver, 
  public ChangeBroadcaster // why is this a changeBroadcaster? may this be obsolete? all widgets
    // (i.e. the sliders and ths plot-editor) actually observe the Parameter object, so we should 
    // not need any other mechanism to sync the widgets
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  LadderSpectrumEditor(const juce::String& name = juce::String("LadderSpectrumEditor"));   

  /** Destructor. */
  virtual ~LadderSpectrumEditor(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::MoogyFilter object which is to be edited. Make 
  sure to call this function again with a NULL-pointer when the object get deleted for some 
  reason. */
  virtual void setFilterToEdit(jura::Ladder* newFilterToEdit);

  /** Assigns a Parameter object to the frequency (horizontal axis) for observation and 
  manipulation. */
  virtual void assignParameterFreq(Parameter* parameterToAssign);

  /** Assigns a Parameter object to the resonance (vertical axis) for observation and 
  manipulation. */
  virtual void assignParameterReso(Parameter* parameterToAssign);

  /** Un-Assigns a previously Parameter object to the horizontal axis. */
  virtual void unAssignParameterFreq();

  /** Un-Assigns a previously Parameter object to the verical axis. */
  virtual void unAssignParameterReso();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted) override;

  /** This method is called when one of the assigned rosic::AutomatableParameters has been changed.
  We override it here in the subclass to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendMessage = false);

  /** Overrides the changeListenerCallback in order to receive messages which this object sends 
  to itself. */
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop up 
  when the right button is clicked for MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e) override;

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e) override;

  /** Overrides the resized-method. */
  virtual void resized() override;

  /** Updates the frequency response plot. */
  virtual void updatePlot();

protected:

  /** Does the setup of the filter according to some new mouse position) */
  virtual void setupFilterAccordingToMousePosition(double mouseX, double mouseY);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handle. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL, 
    XmlElement *targetSVG = NULL) override;

  /** Converts a resonance value to an y-coordinate in components/image coordinates. */
  double resoToY(double reso);

  /** Converts an y-coordinate in components/image coordinates to a resonance value. */
  double yToReso(double y);

  /** Radius of the dot-handle to be drawn. */
  float dotRadius;

  /** Pointer to the actual jura::Ladder object which is being edited. */
  jura::Ladder* filterToEdit;

  // the parameters which wil cause re-plotting and therefore must be listened to:
  Parameter* freqParameter;
  Parameter* resoParameter;

  // magnitude response display stuff:
  int    numBins;
  double *frequencies, *magnitudes;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LadderSpectrumEditor)
};


//=================================================================================================

class JUCE_API rsLadderPlotEditor : public rsVectorPad
{

public:

  rsLadderPlotEditor(jura::Ladder* ladderModuleToEdit);
  virtual ~rsLadderPlotEditor() {}

  virtual void parameterChanged(Parameter* p) override;
  //virtual void paintOverChildren (Graphics& g) override;
  virtual void resized() override;

protected:

  jura::Ladder* ladderToEdit; 

  rsFunctionPlot* freqRespPlot; // frequency response plot (maybe factor out a pointer to rsPlot)

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsLadderPlotEditor)
};


//=================================================================================================

/** This is the GUI editor class for the jura::Ladder audio processor. */

class JUCE_API LadderEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  LadderEditor(jura::Ladder *newLadderToEdit);
  virtual void resized() override;
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;

protected:

  Ladder *ladderToEdit;

  LadderSpectrumEditor *frequencyResponseDisplay; // old - delete when transition is done

  rsLadderPlotEditor* plotEditor;   // the new plot editor

  //AutomatableSlider *cutoffSlider, *resonanceSlider, *spreadSlider;
  ModulatableSlider *cutoffSlider, *resonanceSlider, *spreadSlider;
  AutomatableComboBox *modeComboBox;
  //RComboBox *modeComboBox;
  //RButton *invertButton;  

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LadderEditor)
};

#endif 