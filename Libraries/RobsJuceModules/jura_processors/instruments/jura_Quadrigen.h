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

class QuadrigenModuleEditor : public AudioModuleEditor
{

  friend class QuadrigenAudioModule;

public:

  QuadrigenModuleEditor(CriticalSection *newPlugInLock, QuadrigenAudioModule* newQuadrigenAudioModule);

  virtual ~QuadrigenModuleEditor();

  virtual void initializeColourScheme();

  virtual void acquireLock()
  {
    mutex.enter();
    if( quadrigenModuleToEdit != NULL )
      quadrigenModuleToEdit->acquireLock();
  }
  virtual void releaseLock()
  {
    if( quadrigenModuleToEdit != NULL )
      quadrigenModuleToEdit->releaseLock();
    mutex.exit();
  }

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton  *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox  *rComboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void mouseDown(const MouseEvent& e);
  virtual void paint(Graphics &g);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  /** Returns the index of the generator slot that corresponds to the drawn slot-box at the given 
  position if any, -1 if none. */
  virtual int getSlotIndexAtPixelPosition(int x, int y);

  /** Selects a new generator algorithm for one of the slots. */
  virtual void setGeneratorAlgorithm(int slotIndex, int newAlgorithmIndex);

  /** Calls createEditorForSlotIfNeeded() for all slots. */
  //virtual void createEditorsIfNeeded();

  /** Checks, whether the current editor for the given slot is of the approriate kind and if not,
  it deletes it and creates the approriate one. It will also update the widgets of the 
  (new or old) editor. */
  //virtual void createEditorForSlotIfNeeded(int slotIndex);

  /** Removes all the child editors and nulls the pointers. */
  virtual void removeChildEditors();

  /** Removes the child editor in the given slot and nulls the pointer. */
  virtual void removeChildEditorInSlot(int slotIndex);

  /** Creates an editor for the given slot (slotIndex) of the given kind (algorithmIndex). It 
  will also try to associate the editor with the corresponding AudiModule and update the widgets 
  according to the state of the AudioModule. */
  virtual void createEditorForSlot(int slotIndex, int algorithmIndex);

  QuadrigenAudioModule *quadrigenModuleToEdit;

  //juce::Rectangle leftRectangle1, leftRectangle2, leftRectangle3, eqBandParamRectangle;
  juce::Rectangle<int> globalRectangle;
  juce::Rectangle<int> slotRectangles[rosic::Quadrigen::numGeneratorSlots];
  //RComboBox *generatorSelectComboBoxes[4];

  RoutingMatrixModuleEditor *matrixEditor;
  AudioModuleEditor         *moduleEditors[rosic::Quadrigen::numGeneratorSlots];

  CriticalSection mutex;

  int oldAlgorithmIndices[rosic::Quadrigen::numGeneratorSlots];

  juce_UseDebuggingNewOperator;
};

#endif 
