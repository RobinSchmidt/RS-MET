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

  /** Tries to set up the sfz engine according to the given string which is supposed to be an sfz
  instrument definition and reports whether or not thsi was successful. It may fail due to 
  malformed sfz code, missing sample files, etc. In case of failure, the old instrument will be 
  restored. the function needs to know whether the string comes from an .sfz file on disk or from
  some in memory buffer in order decide if the inherited state of the FileManager's activeFile 
  should be considered dirty or not. If we update from an internal memory buffer, the sfz 
  file-widget will show a asterisk next to the filename to mark it as "dirty". If we update from a
  file, no asterisk will be shown. */
  bool setupFromSfzString(const juce::String& newSfz, bool stringComesFromFile);
  // this may be called from the editor or some widget near it to try to update the instrument 
  // according to the given string...

  const juce::String& getCurrentSfz() const { return lastValidSfz; }


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

/** A subclass of jura::Parameter that can be connected to an SFZ opcode ...tbc...

ToDo: 
-Maybe later derive from AutomatableParameter or some other Parameter subclass to enable more 
 features. I think ModulatableParameter would not be adequate though because modulation is managed
 within the sampler engine itself. But how about smoothing? */

/*
// Code currently commented - I'm not yet sure, if that's the right way to handle it. Maybe it can
// be deleted at some stage.
class JUCE_API SfzOpcodeParameter : public jura::Parameter 
{

public:

  using jura::Parameter::Parameter;        // Inherit constructors
  using Opcode = rosic::Sampler::Opcode;   // Shorthand for convenience

  // Inquiry:
  bool isGlobalSetting()  const { return groupIndex == -1 && regionIndex == -1; }
  bool isGroupSetting()   const { return groupIndex != -1 && regionIndex == -1; }
  bool isRegionSetting()  const { return groupIndex != -1 && regionIndex != -1; }


protected:

  int    groupIndex  = -1;
  int    regionIndex = -1;
  Opcode opcode      = Opcode::Unknown;

  // Maybe add later:
  // int location = -1; // position in the sfz code as character index of the first character of 
     // the opcode. Maybe it could be convenient to store line and column, too

  // Maybe instead of storing group- and region indices, we should store pointers to the actual
  // object. Rationale: when the suer edits the sfz file, the indices may become out of date 
  // whereas a pointer may still remain valid. But actually, re-parsing will rebuild the SFZ 
  // datastructure anyway, so actually no, the pointers would be invalidated anyway. Hmm...Maybe
  // we need to figure out the new indices based on the new and old location of the opcode in the 
  // code whenever the code changes?
};
*/

//=================================================================================================

/** A sampler with functionality roughly based on the sfz specification. It has jura::FileManager
as baseclass to keep track of the currently loaded .sfz file. The editor has also FileManager as
baseclass and uses it, as usual, for keeping track of the loaded .xml file - .xml load/save is 
managed by the GUI but .sfz load/save is managed by the AudioModule. That's a bit like with the 
wavefiles on oscillator modules (right?). */

class JUCE_API SamplerModule : public jura::AudioModuleWithMidiIn, public jura::FileManagerListener
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

  bool setupFromSfzString(const juce::String& newSfz, bool stringComesFromFile) 
  { return sfzPlayer.setupFromSfzString(newSfz, stringComesFromFile); }



  /** Returns the number of layers that are currently playing. Used for displaying that information
  on the GUI.  */
  int getNumActiveLayers() const { return sfzPlayer.getNumActiveLayers(); }

  /** Returns the number of SfzOpcodeParameter objects that we allocate on construction. The 
  objects can be connected to SFZ opcodes at runtime in order to enable the user to control the 
  value of a given opcode via a slider in the GUI. */
  int getNumOpcodeParameters() const { return numOpcodeParams; }


  void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  void activeFileChanged(FileManager *fileMan) override; 

  /** Overriden from AudioModuleWithMidiIn to handle the SfzOpcodeParameters in a different way 
  than usual */
  void parameterChanged(Parameter* param) override;


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

  jura::SfzPlayer sfzPlayer;

  static const int numOpcodeParams = 8;  // Preliminary. Maybe have more later


  friend class SamplerEditor;  // Maybe try to get rid

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerModule)
};

//=================================================================================================

/** A class for represneting the tree nodes in the tree-view for showing the sfz structure. Nodes
can represent either structuring elements like the <group> or <region> tags in the sfz spec or 
opcodes. Opcodes are leaf nodes, structural elements typically not (unless they are empty). */

class SfzTreeViewNode : public jura::RTreeViewNode
{

public:


  using jura::RTreeViewNode::RTreeViewNode;


  /** Type for the user-data that can be stored at the tree-nodes. It contains the information that
  is required to find the corresponding setting in the SfzInstrument datastructure. 
  ToDo: Maybe also store information to find it in the sfz-code such as line/column/location. Maybe
  drag the class out of SfzTreeViewNode and get rid of SfzTreeViewNode. */
  struct Data
  {
    Data(){}

    /** The type of the data stored at the nodes depends on the type of the node. In order to be 
    able to tell, which type it is, we define an enum. */
    enum class Type
    {
      group,
      region,
      playbackSetting,    // e.g. tune, volume, ...
      modulationRouting,  // e.g. lfo3_cutoff2, adsr2_volume1
      unknown             // maybe get rid - type should always be known
    };

    /** We define a sum-type of the passible datatypes that can be stored at the nodes here. 
    ToDo: maybe use std::variant instead (but that requires C++17). */
    union Variant  // find better name
    {
      Variant(){}

      rosic::Sampler::PlaybackSetting   playbackSetting;
      rosic::Sampler::ModulationRouting modRouting;
    };

    Type type = Type::unknown;
    Variant data;
    int groupIndex  = -1;
    int regionIndex = -1;
  };


  Data data;

protected:




  //Type type = Type::unknown;
  // Maybe store data like group/region index (if applicable) etc. Data, such that we can find the
  // corresponding settings in an rosic::Sampler::SfzInstrument. We may need to retrieve pointers
  // to regions and settings...we'll see...



  //jura::RWidget* widget = nullptr;


};
// Maybe we don't need ot make a subclass of RTreeViewNode. Instead, attach the additional data
// in the data pointer that RTreeViewNode defines for exactly this purpose. Maybe the data stored
// should be Type tag, and an PlaybackSetting object which we may leave empty when it's no 
// applicable

//=================================================================================================

/** A subclass of RTreeView for displaying and manipulating the structure of an sfz instrument 
definition. Because sfz instruments have 3 hierarchy levels (global, group, region), the tree has
3 levels, too. Leaf nodes are the individual opcodes. 

ToDo:
When the user clicks on such an opcode node, we'll show some appropriate widgets in some other area 
of the GUI to manipulate the data of the node. In particular, we'll show some textbox that 
describes the type of opcode in more detail along with its range, unit, default value, etc and, 
importantly, a slider to manipulate the numeric value or a combobox to switch between the available
options - whatever is appropriate to the given opcode at hand...tbc... 

When the user clicks on a group node, we may show the number of regions in this group and perhaps 
some other relevant information. In case of a region node, we may show key- and velocity range,
number of opcodes (well...maybe that's not very useful...maybe there's more interesting data to 
show - but we need to fill the space with something). */

class SfzTreeView : public jura::RTreeView
{

public:

  SfzTreeView();

  //using SfzInstrument = rosic::Sampler::SfzInstrument;

  /** Builds (or updates) the internal tree from the given sfz instrument datastructure. This needs
  to be called from the editor, whenever the instrument was changed by some other widget, e.g. the
  code editor to update the contents of the TreeView. */
  void buildTreeFromSfz(const rosic::Sampler::SfzInstrument& sfz);

  void clearTree();


protected:

  SfzTreeViewNode rootNode;  // Manages the lifetimes of all its child-nodes

  //jura::SfzPlayer* player = nullptr; 
  // hmm..nahhh...Let's try to avoid the coupling to this class as long as possible
};


//=================================================================================================

/** Editor for SamplerAudioModule. It features a code editor to edit the sfz file some metering
widgets showing the current system load, etc...tbc...  */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor, 
  public juce::Timer, public juce::CodeDocument::Listener, public jura::FileManagerListener,
  public jura::RTreeViewObserver
{

public:

  SamplerEditor(SamplerModule* samplerToEdit);

  virtual ~SamplerEditor();

  // Overrides                                   // overriden from...
  void timerCallback() override;                 //   Timer
  void resized() override;                       //   AudioModuleEditor
  void treeNodeClicked(RTreeView *treeView, RTreeViewNode *node, const MouseEvent &mouseEvent, 
    int clickPosition) override;
  void treeNodeChanged(RTreeView *treeView, RTreeViewNode *node) override;
  void codeDocumentTextInserted(const String &newText, int insertIndex) override; // CodeDocument::Listener 
  void codeDocumentTextDeleted(int startIndex, int endIndex) override;            // CodeDocument::Listener 

  void rButtonClicked(RButton *buttonThatWasClicked) override;
  void activeFileChanged(FileManager *fileMan) override;                    // FileManagerListener



protected:

  virtual void createWidgets();

  virtual void setCodeIsParsed(bool isParsed);

  virtual void setCodeIsSaved(bool isSaved);

  void setCodeIsDirty() { setCodeIsParsed(false); setCodeIsSaved(false); }
  void setCodeIsClean() { setCodeIsParsed(true);  setCodeIsSaved(true);  }


  /** Tries to parse the current content of the code editor and set up the engine accordingly. If
  parsing fails, the engine will revert into its previous settings. ...tbc... ToDo: explain 
  rationale of this behavior. */
  void parseCodeEditorContent();  // maybe should return a bool?
  //void saveCodeEditorContent();   // dito?

  /** Updates the TreeView with the sfz structure accoring to the instrument definition that is
  currently loaded in the engine. */
  void updateTreeView();

  SamplerModule* samplerModule = nullptr;

  jura::RTextField *instrumentLabel;
  jura::MeteringDisplayWithText *layersMeter;


  // SFZ tree view and adjacent widgets:
  SfzTreeView* sfzTree;
  //SfzTreeNodeWidgetSet* nodeWidgets;


  // SFZ text editor and adjacent widgets:
  jura::FileSelectionBox *sfzFileLoader;
  jura::RLabeledTextEntryField *sfzStatusField;
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