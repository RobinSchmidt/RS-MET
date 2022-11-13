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
  // rename to getCurrentSfzString

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

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzPlayer)
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
// may be oboslete. we are doing it differently...
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


  using Opcode     = rosic::Sampler::Opcode;
  using ReturnCode = rosic::Sampler::rsReturnCode;

  // ToDo: return a ReturnCode instead of an int:

  int setRegionSetting(int groupIndex, int regionIndex, Opcode type, float value, int index)
  { return sfzPlayer.setRegionSetting(groupIndex, regionIndex, type, value, index); }

  int setGroupSetting(int groupIndex, Opcode type, float value, int index)
  { return sfzPlayer.setGroupSetting(groupIndex, type, value, index); }

  int setInstrumentSetting(Opcode type, float value, int index)
  { return sfzPlayer.setInstrumentSetting(type, value, index); }





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
  //using ReturnCode = rosic::Sampler::rsReturnCode;
  using Event      = rosic::Sampler::rsMusicalEvent<float>;

  jura::SfzPlayer sfzPlayer;

  static const int numOpcodeParams = 8;  // Preliminary. Maybe have more later


  friend class SamplerEditor;  // Maybe try to get rid

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerModule)
};

//=================================================================================================
// The Mediator pattern infrastructure for coordinating the interaction of different GUI elements:

class SamplerInterfaceMediator;

/** Baseclass for the GUI components of the sampler. These components coordinate with one another 
through the use of a mediator. */

class SamplerInterfaceComponent : public jura::MediatedColleague
{

public:

  SamplerInterfaceComponent() {}
  virtual ~SamplerInterfaceComponent() {}

  /** Enumeration to be used in the patchChanged callback. */
  enum class PatchChangeType
  {
    opcodeValueChanged,
    //opcodeAdded,
    //opcodeRemoved,
    //opcodeReplaced,
    //regionAdded,
    //regionRemoved,
    //groupAdded,
    //groupRemoved

    unknown
  };


  class PatchChangeInfo : public jura::rsMessageData
  {
  public:
    PatchChangeInfo() {}


    PatchChangeType type = PatchChangeType::unknown;
    int groupIndex  = -1;
    int regionIndex = -1;



    rosic::Sampler::PlaybackSetting oldSetting;
    // maybe we need it to be a sum-type of PlaybackSetting and ModulationRouting ...just like
    // we have done in SfzTreeViewNode. Maybe factor that "Variant" internal class out and use it
    // here too

    float newValue = 0.f;  // Value that we want to set the oldSetting to



    void reset()
    {
      type = PatchChangeType::unknown;
      groupIndex  = -1;
      regionIndex = -1;
      oldSetting.reset();
      newValue = 0.f;
    }
    // move to .cpp

  };


  /** Subclasses must override this to deal with changes of the patch structure. */
  virtual void handlePatchUpdate(const PatchChangeInfo& info) = 0;

  /** Overrides the inherited callback from jura::MediatedColleague to try to dynamically cast the
  passed messageData into a PatchChangeInfo and, if successful, calls the more specific 
  handlePatchUpdate with it (which subclasses must override). If the cast fails, it will trigger an 
  assertion. This should not happen and indicates a bug somewhere else. */
  void handleMediatorNotification(MediatedColleague* originator,
    int messageCode = 0, rsMessageData* messageData = nullptr) override;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerInterfaceComponent)
};

/** The mediator class that coordinates the interaction between the various components of the
sampler-engine's graphical user interface. */

class SamplerInterfaceMediator : public jura::Mediator
{

  
public:

  //using jura::Mediator::Mediator;  // inherit constructor...doesn't seem to work

  SamplerInterfaceMediator() {}      // ...so let's define an empty constructor ourselves


protected:

  std::vector<SamplerInterfaceComponent*> colleagues;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerInterfaceMediator)
};

//=================================================================================================

/** Class for representing the data at a node in the SFZ patch structure. You can think of the node
as a tree-node where <group> and <region> nodes are intermediate levels and the opcodes are the 
leaf nodes. This class stores the information that is required to facilitate the implementation of
certain GUI elements such as the overlay widgets in the TreeView. The information stored here can 
also be used to find the corresponding setting in the SfzInstrument datastructure which is needed 
in order to update the data inside the engine such that parameter changes are audible.
ToDo:
-Maybe also store information to find it in the sfz-code such as line/column/location. 
-Maybe drag the class out of SfzTreeViewNode and get rid of SfzTreeViewNode. But when we build the
 tree, this information is actually not available...hmm... */

class SfzNodeData
{

public:



  /** The type of the data stored at the nodes depends on the type of the node. In order to be 
  able to tell, which type it is, we define an enum. */
  enum class Type
  {
    group,
    region,
    playbackSetting,    // e.g. tune, volume, ...
    modulationRouting,  // e.g. lfo3_cutoff2, adsr2_volume1
    unknown             // maybe get rid, maybe type should always be known
  };


  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  SfzNodeData();
  //SfzNodeData(){}  
  // Try to enforce usage of the factory methods by making the constructor private to ensure that 
  // only valid nodes can be created 

  SfzNodeData(const SfzNodeData& other);


  static SfzNodeData createEmptyNode();

  static SfzNodeData createGroupNode(int groupIndex);

  static SfzNodeData createRegionNode(int groupIndex, int regionIndex);

  static SfzNodeData createPlaybackSettingNode(int groupIndex, int regionIndex, 
    const rosic::Sampler::PlaybackSetting& playbackSetting);

  static SfzNodeData createModulationRoutingNode(int groupIndex, int regionIndex, 
    const rosic::Sampler::ModulationRouting& playbackSetting);


  //-----------------------------------------------------------------------------------------------
  // \name Setup



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns true, if this is a node for an opcode, i.e. a leaf-node. Among other things, this 
  determines whether or not it makes sense to show an editing widget for the node to the user. */
  bool isOpcodeNode() const;

  using OpcodeFormat = rosic::Sampler::OpcodeFormat;

  /** Returns the format of the data that is stored at this node. The return type is a type 
  used in the SfzCodeBook in rosic. Typical values are: Boolean, Natural (unsigned int), Integer,
  Float, String. This information is used to decide, what kind of widget should be displayed to 
  edit the data (button, slider, combobox, text-field, etc.).  */
  OpcodeFormat getOpcodeFormat() const;



  rosic::Sampler::PlaybackSetting getPlaybackSetting() const
  {
    RAPT::rsAssert(type == Type::playbackSetting);
    // Retrieving a meaningless PlaybackSetting is probably a bug

    return data.playbackSetting;
  }

  int getGroupIndex() const { return groupIndex; }

  int getRegionIndex() const { return regionIndex; }



//protected:  // todo: make protected

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


private:

  //SfzNodeData() {} 


  //JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzNodeData)
  // Doesn't compile when this is not commented out. VS-compiler says:
  //   SfzNodeData::SfzNodeData(const jura::SfzNodeData &)  member function already defined or 
  //   declared
  // Figure out why and fix it!
};

//=================================================================================================

/** A class for representing the tree nodes in the tree-view for showing the sfz structure. Nodes
can represent either structuring elements like the <group> or <region> tags in the sfz spec or 
opcodes. Opcodes are leaf nodes, structural elements typically not (unless they are empty). */

class SfzTreeViewNode : public jura::RTreeViewNode
{

public:


  //using jura::RTreeViewNode::RTreeViewNode;

  SfzTreeViewNode() {}  // maybe make protected

  SfzTreeViewNode(const juce::String& nodeText, const SfzNodeData& nodeData) 
    : RTreeViewNode(nodeText), data(nodeData) {}


  static SfzTreeViewNode* createGroupNode(const juce::String& nodeText, int groupIndex)
  {
    return new SfzTreeViewNode(nodeText, SfzNodeData::createGroupNode(groupIndex));
  }

  static SfzTreeViewNode* createRegionNode(const juce::String& nodeText, int groupIndex, 
    int regionIndex)
  {
    return new SfzTreeViewNode(nodeText, SfzNodeData::createRegionNode(groupIndex, regionIndex));
  }

  static SfzTreeViewNode* createPlaybackSettingNode(const juce::String& nodeText, int groupIndex, 
    int regionIndex, const rosic::Sampler::PlaybackSetting& playbackSetting)
  {
    return new SfzTreeViewNode(nodeText, 
      SfzNodeData::createPlaybackSettingNode(groupIndex, regionIndex, playbackSetting));
  }

  static SfzTreeViewNode* createModulationRoutingNode(const juce::String& nodeText, int groupIndex, 
    int regionIndex, const rosic::Sampler::ModulationRouting& playbackSetting)
  {
    return new SfzTreeViewNode(nodeText, 
      SfzNodeData::createModulationRoutingNode(groupIndex, regionIndex, playbackSetting));
  }



  /** Returns true, if this is a node for an opcode, i.e. a leaf-node. Among other things, this 
  determines whether or not it makes sense to show an editing widget for the node to the user. */
  bool isOpcodeNode() { return data.isOpcodeNode(); }

  using OpcodeFormat = rosic::Sampler::OpcodeFormat;

  /** Returns the format of the data that is stored at this node. The return type is a type 
  used in the SfzCodeBook in rosic. Typical values are: Boolean, Natural (unsigned int), Integer,
  Float, String. This information is used to decide, what kind of widget should be displayed to 
  edit the data (button, slider, combobox, text-field, etc.).  */
  OpcodeFormat getOpcodeFormat() { return data.getOpcodeFormat(); }

  const SfzNodeData& getNodeData() const { return data; }


  SfzNodeData data; // make protected!

protected:




  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzTreeViewNode)
};
// Maybe we don't need ot make a subclass of RTreeViewNode. Instead, attach the additional data
// in the data pointer that RTreeViewNode defines for exactly this purpose. Maybe the data stored
// should be Type tag, and an PlaybackSetting object which we may leave empty when it's no 
// applicable

//=================================================================================================

/** A class that holds different types of widgets for the different kinds of SFZ opcodes. It makes
one of those widgets visible at a time - whatever is appropriate for the opcode in question. 
For continuous parameter opcodes (e.g. cutoff=1000) it will show a slider, for choice opcodes (e.g. 
fil_type=lpf_2p) it will show a combo box, for freeform text parameters (e.g. sample="Guitar.wav") a
text entry field. */

class SfzOpcodeWidgetSet : public jura::WidgetSet, public SamplerInterfaceComponent
  ,public jura::RSliderListener, public jura::RComboBoxObserver
  , public jura::RTextEntryFieldObserver
{

public:

  SfzOpcodeWidgetSet();

  virtual~SfzOpcodeWidgetSet() {}

  /** Client cod calls this to set up the kind of opcode/setting that this widget set edits. This 
  will determine whther we will show a slider, combobox, etc. and the range of the slider, the 
  available options in the combobox, etc. */
  //void setSettingToEdit(int groupIndex, int regionIndex, 
  //  const rosic::Sampler::PlaybackSetting& setting);

  /** Client code calls this to set up the kind of node that this widget set edits. This will 
  determine whther we will show a slider, combobox, etc. and the range of the slider, the available
  options in the combobox, etc. */
  void setSfzNodeToEdit(const SfzNodeData& nodeData);




  void handlePatchUpdate(const PatchChangeInfo& info) override;

  // Overriden callbacks for the widgets:
  void rSliderValueChanged(RSlider* s) override;
  void rComboBoxChanged(RComboBox* cb) override;
  void textChanged(RTextEntryField *tf) override;


  enum class WidgetMode { slider, button, chooser, text, none };

protected:

  bool wantsExponentialSlider(rosic::Sampler::Opcode op) const;

  // Overriden juce callbacks:
  void resized() override;


  void setWidgetMode(WidgetMode newMode);

  void createWidgets();

  void updateVisibilities();


  jura::RSlider*         slider;     // For continuous parameters, e.g. cutoff=1000
  jura::RComboBox*       comboBox;   // For choice parameters, e.g. fil_type=lpf_2p
  jura::RTextEntryField* textField;  // For freeform string parameters, e.g. sample="Guitar.wav"

  WidgetMode mode = WidgetMode::none;

  PatchChangeInfo patchChangeInfo;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzOpcodeWidgetSet)
};


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

class SfzTreeView : public jura::RTreeView, public jura::SamplerInterfaceComponent
{

public:

  SfzTreeView();

  //using SfzInstrument = rosic::Sampler::SfzInstrument;
  using Opcode = rosic::Sampler::Opcode;

  /** Builds (or updates) the internal tree from the given sfz instrument datastructure. This needs
  to be called from the editor, whenever the instrument was changed by some other widget, e.g. the
  code editor to update the contents of the TreeView. */
  void buildTreeFromSfz(const rosic::Sampler::SfzInstrument& sfz);

  void clearTree();

  void handlePatchUpdate(const PatchChangeInfo& info) override;


  SfzTreeViewNode* findNode(const PatchChangeInfo& info);


  void mouseMove(const MouseEvent& e) override;
  void setMediator(Mediator *newMediator) override;





protected:

  jura::RTreeViewNode* getGroupNode(jura::RTreeViewNode* parent, int groupIndex);
  jura::RTreeViewNode* getRegionNode(jura::RTreeViewNode* parent, int regionIndex);
  jura::RTreeViewNode* getOpcodeNode(jura::RTreeViewNode* parent, const PatchChangeInfo& info);


  SfzTreeViewNode rootNode;  // Manages the lifetimes of all its child-nodes

  void createWidgets();
  //void updateVisibilities();

  SfzOpcodeWidgetSet *overlayWidgets = nullptr;

  /** Hides all the overlay widgets, i.e. calls juce::Component::setVisible(false) on them. */
  void hideOverlayWidgets();

  /** Shows the overlay widget that is suitable for the given tree node at the given y-coordinate, 
  i.e. a slider for a continuous parameter, a box for a choice parameter, etc. */
  void showOverlayWidget(SfzTreeViewNode* node, int y);


  //jura::SfzPlayer* player = nullptr; 
  // hmm..nahhh...Let's try to avoid the coupling to this class as long as possible

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzTreeView)
};

//=================================================================================================

/** A subeditor for displaying and manipulating an sfz opcode. It will show the name of the 
opcode and a little help text that tells the user about what the opcode does, what its nominal 
range is, which unit it has, etc. It will also show a widget that is appropriate to manipulate the
actual parameter value. */

class SfzOpcodeEditor : public jura::Editor, public jura::SamplerInterfaceComponent, 

  // these baseclasses may be obsolete after refactoring:
  public jura::RSliderListener, public jura::RButtonListener, public jura::RComboBoxObserver,
  public jura::RTextEntryFieldObserver
{



public:

  //using namespace rosic::Sampler;

  SfzOpcodeEditor();
  virtual ~SfzOpcodeEditor() {}



  void setSettingToEdit(int groupIndex, int regionIndex, 
    const rosic::Sampler::PlaybackSetting& setting);



  void handlePatchUpdate(const PatchChangeInfo& info) override;

  // Overriden callbacks for the widgets:
  void rSliderValueChanged(RSlider* s) override;
  void rButtonClicked(RButton* b) override;
  void rComboBoxChanged(RComboBox* cb) override;
  void textChanged(RTextEntryField *tf) override;
  // also obsolete after refactoring

  void setMediator(Mediator *newMediator) override;


  // Overriden juce callbacks:
  void resized() override;




protected:

  //bool wantsExponentialSlider(const rosic::Sampler::PlaybackSetting& setting);




  void createWidgets();

  void updateVisibilities();

  jura::RTextField *opcodeField;     // Shows name/syntax of active opcode
  jura::RTextField *helpField;       // Shows a description text for active opcode

  // This stuff shall be factored out into SfzOpcodeWidgetSet and then we will use an object of 
  // this class here:
  enum class WidgetMode { slider, button, chooser, text, none };
  void setWidgetMode(WidgetMode newMode);
  jura::RSlider* slider;             // Sets continuous parameters
  jura::RButton* button;             // Sets boolean parameters
  jura::RComboBox* comboBox;         // Sets choice parameters
  jura::RTextEntryField* textField;  // Sets freeform string parameters
  WidgetMode mode = WidgetMode::none;
  PatchChangeInfo patchChangeInfo;
  // In our widget callbacks, we (re)assign the value field in the embedded 
  // rosic::Sampler::PlaybackSetting member and then send out a change message to the mediator 
  // which will, in turn, will inform all colleagues about the desired change ...maybe it should
  // actually perform the change before broadcasting the change message? Or shall some other object
  // be responsible for actually performing the change? Maybe the editor itself? Not sure yet...


  //SfzPlayer* sfzPlayer = nullptr; 
  // Pointer to the underlying SFZ player. Should be assigned after construction and remain valid
  // for the whole lifetime of the SfzOpcodeEditor object. 
  // ...not yet sure, if we need this - avoiding it would promote a looser coupling


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SfzOpcodeEditor)
};
// Later, this should communicate with an SfzEditorMediator


//=================================================================================================

/** A subclass of juce::CodeEditor specifically for editing sfz files and for interacting with 
other sfz manipulation components  */

class SfzCodeEditor : public juce::CodeEditorComponent, public jura::SamplerInterfaceComponent
{

public:

  using juce::CodeEditorComponent::CodeEditorComponent;

  void handlePatchUpdate(const PatchChangeInfo& info) override;

  /** For the given PatchChangeInfo, this function figures out the starting position and length of
  the code segment that is responsible for this setting and stores them in the output parameters
  "position" and "length". If no fitting segment is found, it will assign both outputs to -1. */
  void findCodeSegment(const PatchChangeInfo& info, int* startPos, int* endPos);

protected:


  //juce::CodeTokeniser sfzTokenizer;
  // ToDo: implement this, see:
  //   https://docs.juce.com/master/classCodeTokeniser.html
  // The CodeTokeniser class is abstract. We need to make a subclass rsSfzTokenizer. Maybe someone
  // else already did that? Check open-source SFZ sampler projects. When we have that, we need to 
  // pass a pointer to our tokenizer to the constructor of CodeEditorComponent

};


//=================================================================================================

/** Editor for SamplerAudioModule. It features a code editor to edit the sfz file some metering
widgets showing the current system load, etc...tbc...  */

class JUCE_API SamplerEditor : public jura::AudioModuleEditor, 
  public juce::Timer, public juce::CodeDocument::Listener, public jura::FileManagerListener,
  public jura::RTreeViewObserver, public SamplerInterfaceComponent
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

  void handlePatchUpdate(const PatchChangeInfo& info) override;



protected:

  virtual void createWidgets();

  virtual void connectGuiElementsToMediator();



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

  void updateVisibilities();

  void makePlayWidgetsVisible(bool shouldBeVisible);
  void makeEditWidgetsVisible(bool shouldBeVisible);


  SamplerModule* samplerModule = nullptr;

  // Buttons to toggle between the different pages:
  jura::RRadioButton *playButton, *editButton;
  jura::RRadioButtonGroup guiPageButtonGroup;


  // Realtime metering widgets for the "Play" page:
  jura::MeteringDisplayWithText *layersMeter;
  // ToDo: meters for CPU and RAM usage (maybe later disk traffic, when DFD is implemented), 
  // output level, MIDI activity indicator, maybe RAM meter should also show the total occupation
  // (by all apps) and the remaining available, if that is possible


  // SFZ tree view and adjacent widgets:
  jura::RTextField* structureField;  // Says: "Patch Structure", placed above the TreeView
  SfzTreeView* sfzTree;


  // SFZ text editor and adjacent widgets:
  jura::FileSelectionBox *sfzFileLoader;
  jura::RClickButton *parseButton;
  juce::CodeDocument sfzDoc;           // Declare doc before the editor because the editor holds a
  //juce::CodeEditorComponent sfzEditor; // reference to it (-> order of construction/destruction)
  jura::SfzCodeEditor sfzEditor; 


  SfzOpcodeEditor* opcodeEditor;

  // The mediator object that coordinates the interactions between the different parts of the GUI:
  SamplerInterfaceMediator guiMediator;



  // Flags to indicate whether the current content of our sfz code is parsed and/or saved to disk:
  bool codeIsParsed = false;
  bool codeIsSaved  = false;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SamplerEditor)
};
// maybe it should derive from juce::Timer to periodically update the levelMeter and the 
// numLayersField
