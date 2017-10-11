#ifndef jura_Quadrifex_h
#define jura_Quadrifex_h

class QuadrifexModuleEditor;

class QuadrifexAudioModule : public AudioModule
{

  friend class QuadrifexModuleEditor;

public:

  QuadrifexAudioModule(CriticalSection *newPlugInLock, 
    rosic::Quadrifex *quadrifexToWrap = nullptr);

  virtual ~QuadrifexAudioModule();

  AudioModuleEditor* createEditor() override;

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the editor which edits this object - call this in the constructor of the editor
  with 'this' and in the destructor with 'NULL' as parameter. We must maintain a pointer to the
  editor here in order to inform it about changes of the assigned effect-algorithms in the slots
  to make the sub-editors consistent with this AudioModule. 
  \todo get rid of this pointer - use observer pattern, i.e. let editor register itself as 
  observer
  */
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

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    ScopedLock scopedLock(*lock);
    if(numChannels != 2)
      return;
    wrappedQuadrifex->processBlock(inOutBuffer[0], inOutBuffer[1], numSamples);
  }

  //virtual void processStereoFrame(double *left, double *right) override;
  // needs to be overriden, when we want to use the mod-system

  /*
  // older (now obsolete) audio callbacks:
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

  rosic::Quadrifex *wrappedQuadrifex;
  bool wrappedQuadrifexIsOwned = false;



  QuadrifexModuleEditor *editor;
  //....uuuhhhh....this is dangerous. the juce doc says that we should not keep a pointer to
  // the editor - we need to get rid of this


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

/** This class ...  */

class QuadrifexRoutingDiagram : /*public RectangleComponent,*/ public ChangeBroadcaster,
  public CoordinateSystemOld
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  QuadrifexRoutingDiagram(CriticalSection *newPlugInLock);

  /** Destructor. */
  virtual ~QuadrifexRoutingDiagram();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Use this function to inform the object about a new selected algorithm. */
  virtual void setAlgorithmIndex(int slotIndex, int newAlgorithmIndex);

  /** Sets up a new routing and causes a repaint. */
  virtual void setSlotRouting(int newRouting);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the internal variable that represents the selected algorithm. Should be used to
  retrieve the selected algorithm when receiving change-messages form this object (which are send
  out when the algorithm in one of the slots should be changed from this object. */
  virtual int getAlgorithmIndex(int slotIndex) const { return algorithmIndices[slotIndex]; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void mouseEnter(const MouseEvent &e);
  virtual void mouseExit(const MouseEvent &e);
  virtual void mouseDown(const MouseEvent& e);
  virtual void mouseDrag(const MouseEvent& e);
  virtual void paint(Graphics &g);
  //virtual void resized();

protected:

  /** Returns the index of the effect slot that corresponds to the drawn slot-box at the given
  position if any, -1 if none. */
  virtual int getSlotIndexAtPixelPosition(int x, int y);


  virtual void drawRoutingDiagram(Graphics &g);

  virtual void drawSlotBox(Graphics &g, float x, float y, float w, float h, int slotIndex,
    int algoIndex);



  juce::Rectangle<int> slotBoxes[rosic::Quadrifex::numEffectSlots];

  int slotRouting;

  // get rid of them - consolidate in parent QuadrifexModuleEditor:
  int algorithmIndices[rosic::Quadrifex::numEffectSlots];
  int oldAlgorithmIndices[rosic::Quadrifex::numEffectSlots];

  CriticalSection *plugInLock; // mutex to access the edited AudioModule object 
                               // maybe rename this to "lock"

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class QuadrifexModuleEditor : public AudioModuleEditor, public RPopUpMenuObserver,
  public RComboBoxObserver // obsolete?
{

  friend class QuadrifexAudioModule;
  friend class QuadrifexRoutingDiagram;

public:

  // construction/destruction:
  QuadrifexModuleEditor(CriticalSection *newPlugInLock, 
    QuadrifexAudioModule* newQuadrifexAudioModule);
  virtual ~QuadrifexModuleEditor();

  // setup:
  virtual void initializeColourScheme();

  // callbacks:
  //virtual void rButtonClicked(RButton  *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox  *rComboBoxThatHasChanged);
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged);
  virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void mouseDown(const MouseEvent& e);
  virtual void paint(Graphics &g);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();


protected:

  /** Opens the menu to select an effect algorithm for the given effect-slot at the given
  position (with respect to this Component). */
  virtual void openEffectSelectionMenuForSlot(int slotIndex, juce::Point<int> menuPosition);

  /** Returns the index of the effect slot that corresponds to the drawn slot-box at the given
  position if any, -1 if none. */
  virtual int getSlotIndexAtPixelPosition(int x, int y);

  /** Selects a new effect algorithm for one of the slots. */
  virtual void setEffectAlgorithm(int slotIndex, int newAlgorithmIndex);

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

  /** Sets up the position and size of the popup-editor(s) that may (or may not) be present in
  the editor for the given slot index. */
  virtual void setupPopupEditors(int slotIndex);

  QuadrifexAudioModule *quadrifexModuleToEdit;

  RTextField *routingLabel, *permutationLabel;
  RComboBox  *routingComboBox, *permutationComboBox;
  RSlider    *dryWetSlider, *wetLevelSlider, *feedbackSlider;

  //Rectangle leftRectangle1, leftRectangle2, leftRectangle3, eqBandParamRectangle;
  juce::Rectangle<int> globalRectangle;
  juce::Rectangle<int> slotRectangles[rosic::Quadrifex::numEffectSlots];
  //RComboBox *effectSelectComboBoxes[4];

  QuadrifexRoutingDiagram   *routingDiagram;
  RoutingMatrixModuleEditor *matrixEditor;
  AudioModuleEditor         *moduleEditors[rosic::Quadrifex::numEffectSlots];
  EffectSelectionPopup      *effectSelectionPopup;

  // CriticalSection mutex;

  int oldAlgorithmIndices[rosic::Quadrifex::numEffectSlots];

  int slotForWhichMenuIsOpen;

  juce_UseDebuggingNewOperator;
};

#endif 
