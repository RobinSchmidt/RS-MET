#ifndef jura_WaveTable_h
#define jura_WaveTable_h

//===============================================================================================
// class StandardWaveformRendererAudioModule:

class StandardWaveformRendererAudioModule : public AudioModule
{
public:
  StandardWaveformRendererAudioModule(CriticalSection *newPlugInLock,
    rosic::StandardWaveformRenderer *newStandardWaveformRendererToWrap);
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
protected:
  virtual void initializeAutomatableParameters();
  rosic::StandardWaveformRenderer *wrappedStandardWaveformRenderer;
  juce_UseDebuggingNewOperator;
};

//===============================================================================================
// class WaveformBufferAudioModule:

class WaveformBufferAudioModule : public AudioModule, public AudioFileManager
{
public:
  WaveformBufferAudioModule(CriticalSection *newPlugInLock, 
    rosic::WaveformBuffer *newWaveformBufferToWrap);
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;
  virtual bool loadFile(const juce::File& fileToLoad) override;
  virtual bool saveToFile(const juce::File& fileToSaveTo) override;
  virtual bool setAudioData(AudioSampleBuffer* newBuffer, const juce::File& underlyingFile,
    bool markAsClean) override;
  virtual void setWaveformFromFile(const juce::String &fileToLoadFrom);
protected:
  rosic::WaveformBuffer *wrappedWaveformBuffer;
  juce_UseDebuggingNewOperator;
};


//===============================================================================================
// class WaveformRendererAudioModule:

class WaveformRendererAudioModule : public AudioModule
{
  friend class WaveformRendererEditor;
public:
  WaveformRendererAudioModule(CriticalSection *newPlugInLock, 
    rosic::WaveformRenderer *newWaveformRendererToWrap);
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
protected:
  virtual void initializeAutomatableParameters();
  rosic::WaveformRenderer *wrappedWaveformRenderer;
  StandardWaveformRendererAudioModule *standardRendererModule;
  WaveformBufferAudioModule           *waveformBufferModule;
  juce_UseDebuggingNewOperator;
};

//===============================================================================================
// class WaveTableAudioModule:

/**

This class wraps rosic::WaveTable into a rosof::AudioModule to facilitate its use as
plugIn or sub-module inside a plugin.

*/

class WaveTableAudioModule : public AudioModule //, public ChangeListener
{

  friend class WaveTableModuleEditorPopUp;
  friend class WaveTableModuleEditorCompact;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WaveTableAudioModule(CriticalSection *newPlugInLock, rosic::WaveTable *newWaveTableToWrap);

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  /** We override this function inherited from StateManager to trigger an update of the waveform
  buffers whenever some parameter changes that should affect the wavefrom buffer (this does not
  happen automatically because changes in the WaveformRenderer do not go through its embedding
  WaveTable object). */
  virtual void markStateAsDirty() override;

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  // maintain them only temporarily for the transition of the Quadrifex-presets:
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  virtual void setWaveformFromFile(const juce::String &fileToLoadFrom);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override {}

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  // wrapped rosic object:
  rosic::WaveTable *wrappedWaveTable;

  // child modules:
  WaveformRendererAudioModule *rendererModule;

  //FileManager audioFileManager;  // keep only temporarily

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class StandardWaveformEditor:

class StandardWaveformEditor : public AudioModuleEditor, public ChangeBroadcaster, 
  public RComboBoxObserver
{
public:
  StandardWaveformEditor(CriticalSection *newPlugInLock, 
    StandardWaveformRendererAudioModule* newRendererModuleToEdit); 
  virtual void rComboBoxChanged(RComboBox *rComboBoxChanged);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  RNamedComboBox *shapeComboBox;
};


//=================================================================================================
// class WaveformBufferEditor:

class WaveformBufferEditor : public AudioModuleEditor, public ChangeBroadcaster
{
public:
  WaveformBufferEditor(CriticalSection *newPlugInLock, 
    WaveformBufferAudioModule* newWaveformBufferModuleToEdit); 
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  FileSelectionBox          *fileSelectionBox;
  WaveformBufferAudioModule *waveformBufferModuleToEdit;
};


//=================================================================================================
// class WaveformRendererEditor:

class WaveformRendererEditor : public AudioModuleEditor, public ChangeBroadcaster, 
  public RComboBoxObserver
{
public:
  WaveformRendererEditor(CriticalSection *newPlugInLock, 
    WaveformRendererAudioModule* newWaveformRendererModuleToEdit); 
  virtual void rComboBoxChanged(RComboBox *rComboBoxChanged);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();
  virtual void updateWidgetVisibility();
  juce_UseDebuggingNewOperator;
protected:
  RComboBox               *modeComboBox;
  StandardWaveformEditor  *standardEditor;
  WaveformBufferEditor    *bufferEditor;
  rosic::WaveformRenderer *renderer;
};


//=================================================================================================
// class WaveTableModuleEditor:

class WaveTableModuleEditorPopUp : public AudioModuleEditor, public ChangeBroadcaster, 
  public RComboBoxObserver
{

  friend class WaveTableModuleEditorCompact;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  WaveTableModuleEditorPopUp(CriticalSection *newPlugInLock, 
    WaveTableAudioModule* newWaveTableModuleToEdit); 
  virtual ~WaveTableModuleEditorPopUp();

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rComboBoxChanged(RComboBox *rComboBoxChanged);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();
  virtual void updatePlot();

protected:

  // pointers to the edited objects (wrapped and non-wrapped):
  rosic::WaveTable     *waveTableToEdit;
  WaveTableAudioModule *waveTableModuleToEdit;

  // sub-editor for the waveform-renderer:
  WaveformRendererEditor *rendererEditor;

  // plot and related stuff:
  rsDataPlot  *waveformDisplay;
  rsPlot *emptyDisplay;
  double *xValues, *yValuesL, *yValuesR;
  int    numSamplesInPlot;

  // widgets:
  // all the manipulation widgets come here.....or we make an extra class for them

  // other ideas: stereo-shift, invert left and right separately

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class WaveTableModuleEditorCompact

class WaveTableModuleEditorCompact : public AudioModuleEditor
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  WaveTableModuleEditorCompact(CriticalSection *newPlugInLock, 
    WaveTableAudioModule* newWaveTableModuleToEdit);  

  virtual ~WaveTableModuleEditorCompact();  

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the bounds of the popup editor relative to the top-left position of the 
  edit-button. */
  virtual void setPopUpEditorBounds(int x, int y, int w, int h);

  /** Sets the headline text for both, the compact editor and the embedded popup editor. */
  virtual void setHeadlineText(const juce::String& newHeadlineText);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the bounds of the edit-button. */
  //virtual Rectangle getEditButtonBounds() const { return editButton->getBounds(); }

  //---------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();
  virtual void updatePlot();

protected:

  // pointers to the edited objects (wrapped and non-wrapped):
  rosic::WaveTable     *waveTableToEdit;
  WaveTableAudioModule *waveTableModuleToEdit;

  // the popup editor that opens when clicking on the 'Edit' button:
  WaveTableModuleEditorPopUp *popUpEditor;

  // widgets:
  RTextField *waveformLabel;
  RButton    *editButton;

  // plot and related stuff:
  rsDataPlot  *waveformDisplay;
  rsPlot *emptyDisplay;
  double *xValues, *yValuesL, *yValuesR;
  int    numSamplesInPlot;

  // bounds of the big editor relative to the top-left position of the edit-button:
  int popUpEditorX, popUpEditorY, popUpEditorW, popUpEditorH; 

  juce_UseDebuggingNewOperator;
};

#endif 
