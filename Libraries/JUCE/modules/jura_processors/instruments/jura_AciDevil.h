#ifndef jura_AciDevil_h
#define jura_AciDevil_h


/** This class wraps rosic::AciDevil into a rosof::AudioModule to facilitate its use as plugIn. */

class AciDevilAudioModule : public AudioModuleWithMidiIn
{

  friend class AciDevilModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor for wrapping an existing rosic::AciDevil object. */
  AciDevilAudioModule(CriticalSection *newPlugInLock, rosic::AciDevil *aciDevilToWrap);

  /** Constructor that creates our own owned rosic::AciDevil object. */
  AciDevilAudioModule(CriticalSection *newPlugInLock);

  /** Init code called form both constructors. */
  void init();

  virtual ~AciDevilAudioModule();

  AudioModuleEditor* createEditor() override;


  //---------------------------------------------------------------------------------------------
  // parameter settings:


  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedAciDevil->setSampleRate(newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  //virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  //  bool markAsClean);

  //virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  ///** Calculates a stereo-ouput frame. */
  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  //{
  //  *inOutL = *inOutR = wrappedAciDevil->getSample();
  //}
  //virtual void processBlockStereo(float *left, float *right, int numSamples)
  //{
  //  for(int n=0; n<numSamples; n++)
  //    left[n] = right[n] = (float)wrappedAciDevil->getSample();
  //}
  //// maybe the two functions above are obsolote now


  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      inOutBuffer[0][n] = inOutBuffer[1][n] = wrappedAciDevil->getSample();
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    *left = *right = wrappedAciDevil->getSample();
  }

  //---------------------------------------------------------------------------------------------
  // others:

  virtual void noteOn(int noteNumber, int velocity) override
  {
    wrappedAciDevil->noteOn(noteNumber, velocity, 0.0);
  }

  virtual void noteOff(int noteNumber) override
  {
    wrappedAciDevil->noteOn(noteNumber, 0, 0.0);
  }

  virtual void setPitchBend(int pitchBendValue) override
  {
    double wheelValueMapped = (double) (pitchBendValue-8192) / 8192.0; // check this
    wrappedAciDevil->setPitchBend(wheelValueMapped);
  }

//virtual void allNotesOff();



/** Overrides setMidiController which is inherited from both base-classes - and we simply we
pass through the function call to both of them here. */
//virtual void setMidiController(int controllerNumber, int controllerValue); 


protected:

  void createParameters();

  AcidSequencerAudioModule *sequencerModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::AciDevil *wrappedAciDevil;
  bool wrappedAciDevilIsOwned = false;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class AciDevilModuleEditor : public AudioModuleEditor 
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AciDevilModuleEditor(CriticalSection *newPlugInLock, AciDevilAudioModule* newAciDevilAudioModule);


  //---------------------------------------------------------------------------------------------
  // others:

  /** Overrides resized(). */    
  virtual void paint(Graphics &g);

  /** Overrides resized(). */    
  virtual void resized();

protected:

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  AciDevilAudioModule *aciDevilModuleToEdit;

  juce::Rectangle<int>  globalRectangle, oscRectangle, filterRectangle, filterEnvRectangle, 
    ampRectangle, sequencerRectangle;

  RTextField *globalLabel, *oscLabel, *filterLabel, *filterEnvLabel, *ampLabel, *sequencerLabel, 
    *normalLabel, *accentLabel, *subOscLabel;

  RSlider *masterLevelSlider, *accentSlider, *slideTimeSlider, *waveformSlider, *subOscLevelSlider, 
    *subOscWaveformSlider, *cutoffSlider, *resonanceSlider, *envModSlider, *normalDecaySlider, 
    *normalAttackSlider, *accentDecaySlider, *accentAttackSlider, *upwardFractionSlider, 
    *ampDecaySlider, *ampSustainSlider, *ampReleaseSlider, *distortionDriveSlider;

  RTextField *filterModeLabel;

  RComboBox *filterModeBox;
  RButton   *sequencerButton; 

  AcidSequencerModuleEditor *sequencerEditor;

  juce_UseDebuggingNewOperator;
};

#endif
