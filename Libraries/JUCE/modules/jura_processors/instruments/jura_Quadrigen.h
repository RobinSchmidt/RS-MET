#ifndef jura_Quadrigen_h
#define jura_Quadrigen_h

class QuadrigenModuleEditor;

class QuadrigenAudioModule : public AudioModule
{

  friend class QuadrigenModuleEditor;

public:

  QuadrigenAudioModule(CriticalSection *newPlugInLock, rosic::Quadrigen *quadrigenToWrap);

  virtual ~QuadrigenAudioModule();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the editor which edits this object - call this in the constructor of the editor
  with 'this' and in the destructor with 'NULL' as parameter. We must maintain a pointer to the
  editor here in order to inform it about changes of the assigned generator-algorithms in the slots
  to make the sub-editors consistent with this AudioModule. */
  void setEditor(QuadrigenModuleEditor *newEditor);

  /** Sets up the given slot to the given algorithm. */
  void setGeneratorAlgorithm(int slotIndex, int newAlgorithmIndex);

  //---------------------------------------------------------------------------------------------
  // others:

  virtual void acquireLock()
  {
    mutex.enter();
    wrappedQuadrigen->acquireLock();
  }

  virtual void releaseLock()
  {
    wrappedQuadrigen->releaseLock();
    mutex.exit();
  }

  /** Be careful to not de-reference the pointer returned pointer after making it invalid via a
  call to setGeneratorAlgorithm. */
  virtual AudioModule* getGeneratorAudioModule(int slotIndex)
  {
    return generatorModules[slotIndex];
  }

  //---------------------------------------------------------------------------------------------
  // overrides:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

  virtual void setSampleRate(double newSampleRate)
  {
    acquireLock();
    wrappedQuadrigen->setSampleRate(newSampleRate);
    releaseLock();
  }
  virtual void setBeatsPerMinute(double newBpm)
  {
    acquireLock();
    wrappedQuadrigen->setTempoInBPM(newBpm);
    releaseLock();
  }
  virtual void trigger()
  {
    acquireLock();
    wrappedQuadrigen->trigger();
    releaseLock();
  }
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    acquireLock();
    *inOutL = *inOutR = 0.0;
    //if( wrappedQuadrigen != NULL )
    //  wrappedQuadrigen->getSampleFrameStereo(inOutL, inOutR);
    releaseLock();
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    acquireLock();
    wrappedQuadrigen->processBlock(left, right, numSamples);
    releaseLock();
  }

  virtual void reset()
  {
    acquireLock();
    wrappedQuadrigen->reset();
    releaseLock();
  }

protected:

  virtual void initializeAutomatableParameters();

  /** Converts an index for an generator-algorithm into a string. */
  virtual juce::String generatorAlgorithmIndexToString(int index);

  /** Converts a string into an index for an generator-algorithm. */
  virtual int stringToGeneratorAlgorithmIndex(const juce::String& algoString);

  rosic::Quadrigen      *wrappedQuadrigen;
  QuadrigenModuleEditor *editor;

  RoutingMatrixAudioModule *matrixModule;
  AudioModule              *generatorModules[rosic::Quadrigen::numGeneratorSlots];

  CriticalSection mutex; // for what is this? check this... can't we use the inherited lock?

  // for storing and recalling the states when an generator is replaced:
  XmlElement *oscillatorStereoStates[rosic::Quadrigen::numGeneratorSlots];

  juce_UseDebuggingNewOperator;
};

//=============================================================================================

#endif 
