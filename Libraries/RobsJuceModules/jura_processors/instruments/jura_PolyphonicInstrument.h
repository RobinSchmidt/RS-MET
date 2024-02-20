#ifndef jura_PolyphonicInstrument_h
#define jura_PolyphonicInstrument_h

/** This class wraps rosic::PolyphonicInstrument into a rosof::AudioModule to facilitate its 
use as plugIn. */

class PolyphonicInstrumentAudioModule : public AudioModuleWithMidiIn
{

  friend class PolyphonicInstrumentModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  PolyphonicInstrumentAudioModule(CriticalSection *newPlugInLock, 
    rosic::PolyphonicInstrument *instrumentToWrap);
  // maybe deprecate this...

  /** Constructor that doesn't need a pointer to a valid rosic::PolyphonicInstrument object to
  be passed. Instead, use setInstrumentToWrap some time soon after the constructor. */
  PolyphonicInstrumentAudioModule(CriticalSection *newPlugInLock);

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  /** Sometimes, it may not be possible to pass a valid pointer to the underlying instrument 
  to the constructor because the object will have to be created in the constructor of the 
  subclass. In this case, use the other constructor that doesn't need a pointer to be passed,
  create the object in you subclass constructor and then call this function. */
  void setInstrumentToWrap(rosic::PolyphonicInstrument *instrumentToWrap);

  /** Overrides the parameterChanged() method of the indirect AutomationListener base class in
  order to respond to interesting automation events. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
    // get rid of this

  /** Returns the state (i.e. the settings of all relevant parameters) in form of an
  XmlElement. */
  //virtual XmlElement* getStateAsXml(XmlElement* xmlElementToStartFrom = NULL);

  /** Recalls a state (i.e. the settings of all relevant parameters) from an XmlElement. */
  //virtual void setStateFromXml(const XmlElement& xmlState);

  /** Overriden to call allNotesOff before restoring the state. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate) override
  {
    if(underlyingRosicInstrument != NULL)
      underlyingRosicInstrument->setSampleRate(newSampleRate);
  }

  virtual void setBeatsPerMinute(double newBpm) override
  {
    if(underlyingRosicInstrument != NULL)
      underlyingRosicInstrument->setBeatsPerMinute(newBpm);
  }

  virtual void noteOn(int noteNumber, int velocity) override
  {
    if(underlyingRosicInstrument != NULL)
      underlyingRosicInstrument->noteOn(noteNumber, velocity, 0);
  }

  virtual void noteOff(int noteNumber) override
  {
    if(underlyingRosicInstrument != NULL)
      underlyingRosicInstrument->noteOn(noteNumber, 0, 0);
  }

  virtual void setPitchBend(int pitchBendValue) override
  {
    if( underlyingRosicInstrument != NULL )
    {
      double wheelValueMapped = (double) (pitchBendValue-8192) / 8192.0; // check this
      underlyingRosicInstrument->setPitchBend(wheelValueMapped);
    }
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR)
  {
    if(underlyingRosicInstrument != NULL)
      underlyingRosicInstrument->getSampleFrameStereo(outL, outR);
  }

protected:

  /** Fills the array of automatable parameters. */
  void createParameters();

  rosic::PolyphonicInstrument *underlyingRosicInstrument = nullptr;

  friend class PolyphonicInstrumentEditor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyphonicInstrumentAudioModule)
};

//=================================================================================================

class PolyphonicInstrumentEditor : public AudioModuleEditor, public ChangeBroadcaster, 
  public RSliderListener, public TuningFileManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Constructor. */
  PolyphonicInstrumentEditor(CriticalSection *newPlugInLock,   // get rid of the lock parameter
    PolyphonicInstrumentAudioModule* newInstrumentToEdit);

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the colours for the preset section widgets. */
  virtual void setPresetSectionColourScheme(const WidgetColourScheme& newColourScheme);

  /** Sets the colours for the tuning section widgets. */
  virtual void setTuningSectionColourScheme(const WidgetColourScheme& newColourScheme);

  /** Sets the text colour for the info field. */
  virtual void setInfoFieldTextColour(const Colour newColour);

  /** Attaches this editor to the actual plugin which is to be edited. */
  virtual void setInstrumentToEdit(rosic::PolyphonicInstrument* newInstrumentToEdit);

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  //virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  //virtual void  rLabelTextChanged(RLabel *rLabelThatHasChanged);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void updateWidgetsAccordingToState();
  virtual void resized();

protected:

  //-----------------------------------------------------------------------------------------------
  // \name Misc





  rosic::PolyphonicInstrument* instrumentEngine; 
  // The underlying rosic::PolyphonicInstrument object. We acquire the plugInLock whenever
  // this object or the inherited moduleToEdit is accessed

  // The widgets for the global instrument parameters:
  RSlider *levelSlider, *levelByKeySlider, *levelByVelSlider, *midSideRatioSlider, 
    *numVoicesSlider, *compSlider, *masterTuneSlider, *wheelRangeSlider, *glideTimeSlider;
  RButton *tuningMinusButton, *tuningPlusButton, *tuningLoadButton, *glideButton;
  RTextField  *tuningLabel, *tuningFileNameLabel;  // ToDo: use a FileSelectionBox instead

  // ToDo:
  // -Add a slider for the pass-through level. Maybe it should be scaled on percent and go from
  //  -100...+100 or -200...+200. Maybe it should use a sinh-function for scaling


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PolyphonicInstrumentEditor)
};


#endif
