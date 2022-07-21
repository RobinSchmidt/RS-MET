#pragma once


/** An SFZ engine that keeps track of which .sfz file is currently loaded and allows for skipping
through a directory of .sfz files. This functionality is inherited from jura::FileManager whereas 
the actual sfz playback functionality is inherited from rosic::Sampler::rsSamplerEngine2. */

class JUCE_API SfzPlayer : public jura::FileManager, public rosic::Sampler::rsSamplerEngine2
{

public:

  /** Constructor. Configures the .sfz and .wav directories which are to be used. */
  SfzPlayer();

  /** Overriden from jura::FileManager to try to load the given file which is assumed to be an .sfz
  file. If all works well, on return, we will have loaded the given sfz file and our lastValidSfz 
  member will have stored the content of the file and the function will return true. It it goes 
  wrong, we will revert to our lastValidSfz and return false. */
  bool loadFile(const juce::File& fileToLoad) override;

  /** Opens a file saving dialog by which the user can save the current content of our lastValidSfz
  variable into a file. */
  bool saveToFile(const juce::File& fileToSaveTo) override;

  /** Tries to load an .sfz file whose relative (to our sfzRootDir) path is given by the passed 
  string and reports whether or not this was successful.  */
  bool loadFile(const juce::String& relativePath);


  bool setupFromSfzString(const juce::String& newSfz);
  // this may be called from the editor or some widget near it to try to update the instrument 
  // according to the given string...


protected:

  using Engine     = rosic::Sampler::rsSamplerEngine2;  
  using ReturnCode = rosic::Sampler::rsReturnCode;

  /** Sets up the member variables that define where the app expects sfz-files, samples, etc. 
  Called in constructor. */
  virtual void setupDirectories();


  /** The root directory where we expect all the sfz files to be. */
  juce::String sfzRootDir;
  // I think, sfzRootDir is redundant with the std::string in the rsSamplerEngine2 baseclass, so
  // maybe try to get rid.

  //juce::String sampleRootDir;
  // This is (or should be) also stored in the Engine baseclass object (as std::string)

  /** When trying to update according to a new sfz string or file, we may fail because the sfz may
  be malformed, files may be missing, etc. In such a case, we'll revert to the last known valid
  sfz string which we keep in this variable. */
  juce::String lastValidSfz;

  // Maybe we need a boolean flag to indicate, whether the current state was saved and/or our
  // lastValidSfz is in sync with the most recently loaded or saved sfz-file. The current state
  // of the engine should always be in sync with our lastValidSfz member...we should actually make
  // that a class invariant, I think.

};

//=================================================================================================

/** A sampler with functionality roughly based on the sfz specification. It has jura::FileManager
as baseclass to keep track of the currently loaded .sfz file. The editor has also FileManager as
baseclass and uses it, as usual, for keeping track of the loaded .xml file - .xml load/save is 
managed by the GUI but .sfz load/save is managed by the AudioModule. That's a bit like with the 
wavefiles on oscillator modules (right?). */

class JUCE_API SamplerModule : public jura::AudioModuleWithMidiIn /*, public jura::FileManager*/
{

  // Oh no - we should not derive from FileManager. We are already an indirect subclass of it due
  // to AudioModule which is a StateFileManager...what can we do?
  // Maybe we need to make a class SfzEngine that derives from jura::FileManager and 
  // rosic::Sampler::SamplerEngine2 and this call is then responsible for keeping track of the sfz
  // file

public:

  SamplerModule(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr);



  // overriden from AudioModule baseclass:
  AudioModuleEditor* createEditor(int type) override;

  void setSampleRate(double newSampleRate) override;
  void setGain(double newGain);

  bool setupFromSfzString(const juce::String& newSfz) { return engine.setupFromSfzString(newSfz); }



  int getNumActiveLayers() const { return engine.getNumActiveLayers(); }


  void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;


  // Midi Handling:
  void noteOn(int key, int vel) override;
  void noteOff(int key) override;
  //void handleMidiMessage(MidiMessage message) override;

  // Audio Processing:
  void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  void processStereoFrame(double *left, double *right) override;

  void reset() override;


protected:

  /** Creates the parameters of the sampler that sit directly in the xml file. They have nothing to
  do with the opcodes defined in the sfz files and they also do not alter the sfz settings that the
  engine works with. Anything that is controlled by these parameters is either post-processing step
  of the engine's output or controls a global behavior that is not specified in the sfz 
  specification such as the resampling quality, the selection whether opcodes should work 
  accumulatively or overridingly, polyphony, a global gain, etc.. */
  virtual void createParameters();





  void setBusMode(bool shouldAccumulate);

  // Shorthands for convenience:
  using ReturnCode = rosic::Sampler::rsReturnCode;
  using Event      = rosic::Sampler::rsMusicalEvent<float>;

  jura::SfzPlayer engine;  // maybe rename to sfzPlayer


  friend class SamplerEditor;  // maybe try to get rid

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerModule)
};


//=================================================================================================

/** Editor for SamplerAudioModule */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor, public jura::FileManager, 
  public juce::Timer, public juce::CodeDocument::Listener 
{

public:

  SamplerEditor(SamplerModule* samplerToEdit);

  virtual ~SamplerEditor();

  // Overrides                                   // overriden from...
  bool loadFile(const juce::File& f) override;   //   FileManager
  bool saveToFile(const juce::File& f) override; //   FileManager
  void timerCallback() override;                 //   Timer
  void resized() override;                       //   AudioModuleEditor
  void codeDocumentTextInserted(const String &newText, int insertIndex) override; // CodeDocument::Listener 
  void codeDocumentTextDeleted(int startIndex, int endIndex) override;            // CodeDocument::Listener 
  void rButtonClicked(RButton *buttonThatWasClicked) override;



protected:

  virtual void createWidgets();

  virtual void setCodeIsParsed(bool isParsed);

  virtual void setCodeIsSaved(bool isSaved);

  virtual void setCodeIsDirty() { setCodeIsParsed(false); setCodeIsSaved(false); }

  void parseCurrentEditorContent();  // maybe should return a bool?
  void saveCurrentEditorContent();   // dito?


  SamplerModule* samplerModule = nullptr;

  jura::RTextField *instrumentLabel;
  jura::MeteringDisplayWithText *layersMeter;

  // SFZ text editor and adjacent widgets:
  jura::FileSelectionBox *sfzFileLoader;
  jura::RClickButton     *parseButton;
  juce::CodeDocument sfzDoc;           // Declare doc before the editor because the editor holds a 
  juce::CodeEditorComponent sfzEditor; // reference to it (-> order of construction/destruction)
  //juce::CodeTokeniser sfzTokenizer;
  // ToDo: implement this, see:
  //   https://docs.juce.com/master/classCodeTokeniser.html
  // The CodeTokeniser class is abstract. We need to make a subclass rsSfzTokenizer. Maybe someone
  // else already did that? Check open-source SFZ sampler projects. When we have that, we need to 
  // pass a pointer to our tokenizer to the constructor of the editor

  // Maybe have an RClickButton for Reparse/Parse/Update ...hwoever we want to call it. maybe it 
  // should get automatically highlighted, as soon as text was edited such that the engine is out
  // of date. The FileLoader should also show a "Dirty" star next to the filenam, when the current
  // state is not save inot a file



  /*

  RTextField *numLayersLabel, *numLayersOfLabel, ;
  RDraggableNumber *maxNumLayersSlider;

  RTextField *cpuLoadLabel, *cpuLoadField, *ramLoadLabel, *ramLoadField;  
  // todo: diskLoad - maybe ram should also show the total occupation (by all apps) and the 
  // remaining available
  */

  // Flags to indicate whether the current content of our sfz code is parsed and/or saved to disk:
  bool codeIsParsed = false;
  bool codeIsSaved  = false;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerEditor)
};
// maybe it should derive from juce::Timer to periodically update the levelMeter and the 
// numLayersField