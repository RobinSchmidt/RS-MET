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

  AudioModuleEditor* createEditor(int type) override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedAciDevil->setSampleRate(newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

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

protected:

  void createParameters();

  AcidSequencerAudioModule *sequencerModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::AciDevil *wrappedAciDevil;
  bool wrappedAciDevilIsOwned = false;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AciDevilAudioModule)
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

  virtual void createWidgets();

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  AciDevilAudioModule *aciDevilModuleToEdit;

  juce::Rectangle<int>  globalRectangle, oscRectangle, filterRectangle, filterEnvRectangle, 
    ampRectangle, distRectangle, sequencerRectangle;

  RTextField *globalLabel, *oscLabel, *filterLabel, *filterEnvLabel, *ampLabel, *distLabel, 
    *sequencerLabel, *normalLabel, *accentLabel, *subOscLabel;

  rsAutomatableSlider *masterLevelSlider, *accentSlider, *slideTimeSlider, *waveformSlider, 
    *pulseWidthSlider, *subOscLevelSlider, *subOscWaveformSlider, *cutoffSlider, *resonanceSlider, 
    *envModSlider, *normalDecaySlider, *normalAttackSlider, *accentDecaySlider, 
    *accentAttackSlider, *upwardFractionSlider, *ampDecaySlider, *ampSustainSlider, 
    *ampReleaseSlider, *distortionDriveSlider;
    // maybe we can use an array to reduce the number of member variables? ah - no, we need all the
    // separate pointers in resized

  RTextField *filterModeLabel;

  rsAutomatableComboBox *filterModeBox;
  //RButton   *sequencerButton; 

  AcidSequencerModuleEditor *sequencerEditor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AciDevilModuleEditor)
};

#endif
