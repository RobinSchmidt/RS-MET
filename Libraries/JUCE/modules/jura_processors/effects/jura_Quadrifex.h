#ifndef jura_Quadrifex_h
#define jura_Quadrifex_h

class QuadrifexModuleEditor;

class QuadrifexAudioModule : public AudioModule
{

  friend class QuadrifexModuleEditor;

public:

  QuadrifexAudioModule(CriticalSection *newPlugInLock, rosic::Quadrifex *quadrifexToWrap);

  virtual ~QuadrifexAudioModule();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the editor which edits this object - call this in the constructor of the editor
  with 'this' and in the destructor with 'NULL' as parameter. We must maintain a pointer to the
  editor here in order to inform it about changes of the assigned effect-algorithms in the slots
  to make the sub-editors consistent with this AudioModule. */
  void setEditor(QuadrifexModuleEditor *newEditor);

  /** Sets up the given slot to the given algorithm. */
  void setEffectAlgorithm(int slotIndex, int newAlgorithmIndex);

  //---------------------------------------------------------------------------------------------
  // others:

  /*
  virtual void acquireLock()
  {
    mutex.enter();
    wrappedQuadrifex->acquireLock();
  }
  virtual void releaseLock()
  {
    wrappedQuadrifex->releaseLock();
    mutex.exit();
  }
  */

  /** Be careful to not de-reference the pointer returned pointer after making it invalid via a
  call to setEffectAlgorithm. */
  virtual AudioModule* getEffectAudioModule(int slotIndex)
  {
    return effectModules[slotIndex];
  }

//---------------------------------------------------------------------------------------------
// overrides:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  virtual void setSampleRate(double newSampleRate)
  {
    ScopedLock scopedLock(*lock);
    wrappedQuadrifex->setSampleRate(newSampleRate);
  }
  virtual void setBeatsPerMinute(double newBpm)
  {
    ScopedLock scopedLock(*lock);
    wrappedQuadrifex->setTempoInBPM(newBpm);
  }
  virtual void trigger()
  {
    ScopedLock scopedLock(*lock);
    wrappedQuadrifex->trigger();
  }
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    ScopedLock scopedLock(*lock);
    *inOutL = *inOutR = 0.0;
    //wrappedQuadrifex->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    ScopedLock scopedLock(*lock);
    wrappedQuadrifex->processBlock(left, right, numSamples);
  }

  /*
  virtual void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
  {
    ScopedLock scopedLock(*plugInLock);
    acquireLock();
    if( wrappedQuadrifex == NULL || buffer.getNumChannels() < 1 )
    {
      jassertfalse;
      releaseLock();
      return;
    }
    float *left, *right;
    left  = buffer.getSampleData(0, 0);
    if( buffer.getNumChannels() < 2 )
      right = buffer.getSampleData(0, 0);
    else
      right = buffer.getSampleData(1, 0);
    wrappedQuadrifex->processBlock(left, right, buffer.getNumSamples());
    releaseLock();
  }
  */

  virtual void reset()
  {
    ScopedLock scopedLock(*lock);
    wrappedQuadrifex->reset();
  }

protected:

  virtual void initializeAutomatableParameters();

  /** Converts an index for a routing into a string. */
  virtual juce::String slotRoutingIndexToString(int index);

  /** Converts a string into an index for a slot routing. */
  virtual int stringToSlotRoutingIndex(const juce::String& routingString);

  /** Converts an index for an effect-algorithm into a string. */
  virtual juce::String effectAlgorithmIndexToString(int index);

  /** Converts a string into an index for an effect-algorithm. */
  virtual int stringToEffectAlgorithmIndex(const juce::String& algoString);

  rosic::Quadrifex      *wrappedQuadrifex;
  QuadrifexModuleEditor *editor;

  RoutingMatrixAudioModule *matrixModule;
  AudioModule              *effectModules[rosic::Quadrifex::numEffectSlots];

  //CriticalSection mutex; // we must get rid of this mutex - it causes now deadlocks

  // for storing and recalling the states when an effect is replaced:
  XmlElement *bitCrusherStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *chorusStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *combBankStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *combResonatorStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *combStereoizerStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *compressorStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *dualTwoPoleFilterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *equalizerStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *expanderStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *flangerStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *formantShifterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *fourPoleFilterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *frequencyShifterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *harmonicsStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *ladderFilterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *limiterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *modulatedAllpassStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *noiseGateStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *noisifierStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *phaserStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *phaseStereoizerStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *pingPongEchoStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *pitchShifterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *reverbStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *ringModulatorStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *simpleDelayStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *sineOscillatorStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *singleSidebandModulatorStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *slewRateLimiterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *slopeFilterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *stereoPanStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *stereoWidthStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *tremoloStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *twoPoleFilterStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *vibratoStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *wahWahStates[rosic::Quadrifex::numEffectSlots];
  XmlElement *waveShaperStates[rosic::Quadrifex::numEffectSlots];

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

#endif 
