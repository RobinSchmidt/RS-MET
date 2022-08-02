#ifndef jura_Liberty_h
#define jura_Liberty_h

//=================================================================================================

/** This is a class for keeping track of the state of the user interface of the modular synth 
(i.e., which panel is open, what is the scroll/zoom-setting, etc.).  */

class LibertyInterfaceState
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  LibertyInterfaceState();

  //-----------------------------------------------------------------------------------------------
  // persistence:

  //  ...getAsXml, setFromXml

  enum panels
  {
    STRUCTURE_PANEL = 0,
    INTERFACE_PANEL,
    HYBRID_PANEL
  };

  int activePanel;

};

//=================================================================================================

/** This class wraps romos::Liberty into a jura::AudioModule to facilitate its use as plugIn. */

class LibertyAudioModule : public PolyphonicInstrumentAudioModule, public ActionBroadcaster
//class LibertyAudioModule : public AudioModule
{

  friend class LibertyEditor;
  friend class LibertyInterfaceMediator;
  friend class ModularBlockDiagramPanel;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  //LibertyAudioModule(CriticalSection *newPlugInLock, romos::Liberty *modularSynthToWrap);

  LibertyAudioModule(CriticalSection *newPlugInLock);

  void init(); // called from constructors, encapsulates their common code

  virtual ~LibertyAudioModule();

  AudioModuleEditor* createEditor(int type) override;

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedLiberty->setSampleRate(newSampleRate);
  }

  //-----------------------------------------------------------------------------------------------
  // persistence:

  /** Adds module-type specific state data to an existing xml-element, if any. For example, a 
  constant-module may have its value stored here, or a multi-stage equalizer may store the number 
  of biquad stages and settings per stage, etc. */
  static void writeModuleTypeSpecificStateDataToXml(romos::Module *module, XmlElement* xmlState);

  /** Restores module-type specific state data from a xml-element, if any. 
  @see writeModuleTypeSpecificStateDataToXml */
  static void restoreModuleTypeSpecificStateDataFromXml(romos::Module *module, 
    const XmlElement& xmlState);

  /** Returns the state of the passed module as a (pointer to) an XmlElement. The caller is 
  responsible for eventually deleting the XmlElement or adding it as child-element to some other 
  XmlElement in which case it will be deleted together the parent-element. The second argument 
  specifies whether or not the external connections of the module should be stored - when the user 
  stores a container from the diagram-panel, we don't want this - otherwise, normally yes. */
  static XmlElement* getModuleStateAsXml(romos::Module *module, bool withExternalConnections);

  /** Sets up the state of the passed module from the passed XmlElement. First, it calls 
  createAndSetupEmbeddedModulesFromXml, then it calls createConnectionsFromXml making the 
  restoration of the state a two-pass process through the XmlElement. */
  static void setModuleStateFromXml(const XmlElement& xmlState, romos::Module *module);

  /** Creates all the child modules (recursively) for the passed module and sets theri state from 
  the passed XmlElement. This function  does not wire up the connections. */
  static void createAndSetupEmbeddedModulesFromXml(const XmlElement& xmlState, 
    romos::Module *module);

  /** Assuming that the child modules have been already properly created, it creates all the 
  connections in the passed module (and recursively inside its children, if any. */
  static void createConnectionsFromXml(const XmlElement& xmlState, romos::Module *module);

  /** Returns the state of the whole instrument as XmlElement. */
  virtual XmlElement* getStateAsXml(const juce::String &stateName, bool markAsClean) override;

  /** Restores the state of the whole instrument from the passed XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;

  /** Calling this functions inside setStateFromXml is a quick and dirty ad-hoc solution to restore 
  the names and positions of the top-level I/O modules. */
  virtual void restoreTopLevelInOutStates(const XmlElement& xmlState);

  //  \todo override the xml get/set methods to include the GUI state

  /** Returns a pointer to the interface state object that is maintained here. This pointer can be 
  used by the GUI to set and retrieve interface settings from the LibertyAudioModule object. We 
  maintain the interface state here for storing it in the presets. ...nah- we don't need that 
  function - the GUI can access the member directly */
  //LibertyInterfaceState* getInterfaceState() { return &interfaceState; }

  //-----------------------------------------------------------------------------------------------
  // event-handling

  /** Triggers a note-on event. */
  virtual void noteOn(int noteNumber, int velocity) override;

  /** Triggers a note-off event. */
  virtual void noteOff(int noteNumber) override;


  //-----------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    if(wrappedLiberty->isSilent())
      return;

    double tmpL, tmpR;
    for(int n = 0; n < numSamples; n++)
    {
      wrappedLiberty->getSampleFrameStereo(&tmpL, &tmpR);
      inOutBuffer[0][n] += tmpL; 
      inOutBuffer[1][n] += tmpR;
    }
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    double tmpL, tmpR;
    wrappedLiberty->getSampleFrameStereo(&tmpL, &tmpR);
    *left  += tmpL;
    *right += tmpR;
  }

  /*
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedLiberty->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    if(wrappedLiberty->isSilent())
    {
      rosic::fillWithZeros(left, numSamples);
      rosic::fillWithZeros(right, numSamples);
    }
    else
    {
      wrappedLiberty->getBlockOfSampleFramesStereo(left, right, numSamples);
      //for(int n = 0; n < numSamples; n++)
      //{
      //  double dL = (double) left[n];
      //  double dR = (double) right[n];
      //  wrappedLiberty->getSampleFrameStereo(&dL, &dR);
      //  left[n]  = (float) dL;
      //  right[n] = (float) dR;
      //}
    }
  }
  */

  //-----------------------------------------------------------------------------------------------
  // others:

  virtual void reset() override
  {
    wrappedLiberty->reset();
  }

protected:

  /** Pointer to the underlying object which is wrapped. */
  romos::Liberty *wrappedLiberty;
  bool wrappedLibertyIsOwned = false;

  LibertyInterfaceState interfaceState; // maintains info about open panels, scroll-positions, etc.

  juce::File macroDirectory;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** A mix-in class to mix into the basic jura widget classes. 

ToDo:
-Docmuent intention - why do we need this? Maybe it can be useful in other contexts as well. If so,
 rename it with a more general name, e.g. NamedParameterWidget or soemthing  */

class LibertyModuleWidget
{

public:

  LibertyModuleWidget(const juce::String &name) 
  { 
    widgetParameterName = name; 
  }

  juce::String getWidgetParameterName() const   
  { 
    return widgetParameterName; 
  }

  juce_UseDebuggingNewOperator;

protected:

  juce::String widgetParameterName;

};

/** RSlider, augmented by the LibertyModuleWidget mix-in class. */
class LibertySlider : public RSlider, public LibertyModuleWidget
{
public:
  LibertySlider(const juce::String &name) : LibertyModuleWidget(name), RSlider(name) { }
  juce_UseDebuggingNewOperator;
};

/** RComboBox, augmented by the LibertyModuleWidget mix-in class. */
class LibertyComboBox : public RComboBox, public LibertyModuleWidget
{
public:
  LibertyComboBox(const juce::String &name) : LibertyModuleWidget(name), RComboBox(name) { }
  juce_UseDebuggingNewOperator;
};

/** RNamedComboBox, augmented by the LibertyModuleWidget mix-in class. */
class LibertyNamedComboBox : public RNamedComboBox, public LibertyModuleWidget
{
public:
  LibertyNamedComboBox(const juce::String &name) : LibertyModuleWidget(name), 
    RNamedComboBox(name, name) { }
  juce_UseDebuggingNewOperator;
};

/** RTextEntryField, augmented by the LibertyModuleWidget mix-in class. */
class LibertyTextEntryField : public RTextEntryField, public LibertyModuleWidget
{
public:
  LibertyTextEntryField(const juce::String &name) 
    : LibertyModuleWidget(name), RTextEntryField(name) { }
  juce_UseDebuggingNewOperator;
};


/** RLabeledTextEntryField, augmented by the LibertyModuleWidget mix-in class. */
class LibertyLabeledTextEntryField : public RLabeledTextEntryField, public LibertyModuleWidget
{
public:
  LibertyLabeledTextEntryField(const juce::String &name) 
    : LibertyModuleWidget(name), RLabeledTextEntryField(name) 
  { 
    removeChildWidget(entryField, true, true);
    entryField = new LibertyTextEntryField(name);
    addChildWidget(entryField);
  }

  juce_UseDebuggingNewOperator;
};

// todo: make classes LibertyButton, etc.
// we also need to override setValueFromString in RButton, RComboBox, RTextEntryField, etc. then

//=================================================================================================
// class ModulePropertiesEditor:

/** This class is used to edit some properties of modules in the modular synthesizer such as the 
name, the polyphony setting, etc. */

class ModulePropertiesEditor : public Editor, public RSliderListener, public RComboBoxObserver, 
  public RTextEntryFieldObserver, public RButtonListener, public ActionListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor.  */  
  ModulePropertiesEditor(LibertyAudioModule *newLiberty, romos::Module* newModuleToEdit);
  /*ModulePropertiesEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit);*/


  /** Destructor. */
  virtual ~ModulePropertiesEditor();


  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rSliderValueChanged(RSlider         *rSlider);
  virtual void rComboBoxChanged(   RComboBox       *comboBoxThatHasChanged);
  virtual void textChanged(        RTextEntryField *rTextEntryFieldThatHasChanged);
  virtual void rButtonClicked(     RButton         *buttonThatWasClicked);
  virtual void resized();

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Updates the widgets from the moduleToEdit member. */
  virtual void updateWidgetsFromModuleState();


protected:


  /** Called internally by the callbacks for the specific widgets rSliderValueChanged, 
  rComboBoxChanged, etc.. The function tries to cast the widget to a LibertyWidget, reads out the 
  parameter-name and sets the parameter with given in the moduleToEdit (provided, that this 
  moduleToEdit has parameters, which it should when this function is called). */
  virtual void widgetChanged(RWidget *widgetThatHasChanged);


  virtual void actionListenerCallback(const String& message) override;


  CriticalSection* plugInLock;
  LibertyAudioModule* libertyModule;
  romos::Module* moduleToEdit;

  RTextField *moduleTypeLabel, *moduleTypeField;
  RButton    *polyButton;


  // preset load/save section

  // non-editable fields (just for info):
  // String moduleKind;
  // int x, y;
  // numAudioInputs, numChildModules, numAudioOutputs, numEventInputs, numEventOutputs

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// the main GUI classes

class LibertyInterfaceMediator;

/** Baseclass for the GUI components of the modular system. These components coordinate with one 
another through the use of a mediator. */

class LibertyInterfaceComponent : public MediatedColleague
{

  //friend class LibertyEditor;

public:

  /** Enumeration to be used in the moduleChanged callback (? is that comment out of date? where is
  this callback -> figure out and update comment). */
  enum moduleChangeTypeIndices
  {
    MODULE_NAME,               // the module's name in the diagram was changed
    //POSITION,                // the module's position in the diagram was changed
    NUM_CHILDREN,              // the number of child modules changed due to adding/removing
    NUM_CONNECTIONS,
    POLYPHONY,                 // polyphony of one or more of the child-modules changed
    //NUM_AUDIO_INPUTS,        // the number of audio inputs was changed
    //NUM_AUDIO_OUTPUTS,       // the number of audio outputs was changed
    //NUM_EVENT_INPUTS,        // the number of event inputs was changed
    //NUM_EVENT_OUTPUTS,       // the number of event outputs was changed
    //CONNECTIVITY,            // the connectivity between any two child modules was changed due to adding/removing wires
    //SELECTION,               // the set of selected modules has changed

    //EDITED_MODULE_CONTAINER,   // module container shown in the editor changed

    CONTAINER_SHOWN_IN_DIAGRAM,
    MODULE_TO_SHOW_EDITOR_FOR


    //MODULE_FOCUS             // the focused moduled was changed
  };

  /** Constructor.  */  
  LibertyInterfaceComponent(LibertyInterfaceMediator *interfaceMediatorToUse);

  LibertyInterfaceMediator* getInterfaceMediator() const;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class coordinates the different components in the user interface for the modular 
system. */

class LibertyInterfaceMediator : public Mediator
{

  friend class ModularStructureTreeView;
  friend class ModulePropertiesEditorHolder;
  friend class ModularBlockDiagramPanel;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor.  */  
  LibertyInterfaceMediator(CriticalSection *newPlugInLock, 
    LibertyAudioModule* newLibertyModuleToEdit);

  /** Destructor. */
  virtual ~LibertyInterfaceMediator();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Chooses a new ModuleContainer object to be shown and edited in the block diagram. */
  void setContainerToShowInDiagram(romos::ContainerModule* containerToShow);
  // \todo pass the pointer to the colleague that requests the setEd....in order to pass it through 
  // to the other colleagues

  /** Sets the module for which the editor will be shown. */
  void setModuleToShowEditorFor(romos::Module *module);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the toplevel module which represents the outermost module that can be edited. */
  romos::TopLevelModule*  getTopLevelModule() const { return topLevelModule; }

  /** Returns the container that is currently being edited/shown in the diagram editor. */
  romos::ContainerModule* getContainerShownInDiagram() const { return containerShownInDiagram; }

  /** Returns the module that is currently marked as selected (might be NULL in case nothing is 
  selected or also when multiple modules are selected) */
  romos::Module* getModuleToShowEditorFor() const { return moduleToShowEditorFor; }


protected:

  CriticalSection         *plugInLock; 
  LibertyAudioModule      *modularSynthModuleToEdit;
  romos::TopLevelModule   *topLevelModule;
  romos::ContainerModule  *containerShownInDiagram;
  romos::Module           *moduleToShowEditorFor;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class is used to visually represent a modular patch as tree and provides facilities for 
selection of modules. */

class ModularStructureTreeView : public LibertyInterfaceComponent, public RTreeView
{

public:

  /** Constructor.  */  
  ModularStructureTreeView(LibertyInterfaceMediator *interfaceMediatorToUse);

  /** Destructor.  */  
  virtual ~ModularStructureTreeView();

  /** Callback that is called from the mediator when another GUI component has changed something 
  that might be interest here. */
  void mediatorHasSentNotification(MediatedColleague *originatingColleague, 
    int messageCode = 0, void* messageData = nullptr) override;

  /** Overriden from RTreeView to .... */
  void nodeClicked(RTreeViewNode *nodeThatWasClicked, const MouseEvent &mouseEvent, 
    int clickPosition) override;

protected:

  /** Returns a pointer to the module that is at the passed x-position or NULL if none. */
  //virtual romos::ModuleProperties* getModuleAt(int x) const; // maybe later

  /** Re-builds the whole tree of nodes. */
  virtual void rebuildTree();

  /** Creates a subtree for the passed moduleToCreateSubTreeFor and hangs it into the passed 
  parentNodeToUse. */
  virtual void createAndHangInSubTree(RTreeViewNode *parentNodeToUse, 
    romos::Module *moduleToCreateSubTreeFor);

  /** Highlights the tree node that is associated with the focused module (and un-highlights all 
  others). */
  virtual void updateNodeHighlighting();

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class is used to hold the properties editor for the selected module. We need this wrapper 
because the content component must be created and deleted dynamically according to the selected 
module. */

class ModulePropertiesEditorHolder : public ColourSchemeComponent, public LibertyInterfaceComponent
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor.  */  
  ModulePropertiesEditorHolder(LibertyInterfaceMediator *interfaceMediatorToUse);

  /** Destructor.  */  
  virtual ~ModulePropertiesEditorHolder();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Callback that is called from the mediator when another GUI component has changed something 
  that might be interest here. */
  void mediatorHasSentNotification(MediatedColleague *originatingColleague, 
    int messageCode = 0, void* messageData = nullptr) override;

  void paint(Graphics &g) override;
  void resized() override;


protected:

  /** Creates the properties editor for the currently selected module (or a default top-level 
  editor, if none is selected). */
  virtual void createPropertiesEditorForSelectedModule();

  ModulePropertiesEditor *currentEditor;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class ModularBlockDiagramPanel:

/** This class is used to visually represent a block-diagram for an object of class 
romos::ModuleContainer.

\todo factor out a class that is applicable to handle plugin-routing, too (like the JUCE plugin 
host) - decouple it from the romos stuff ...maybe call it FlowchartEditor and move it over to 
jura framework....maybe also have a class BlockDiagramElement that wraps the romos::Module class

\todo when the focused module gets deleted (from outside, namely in the tree), set the new focus 
to its parent 

maybe rename to ModularPatcher

*/

class ModularBlockDiagramPanel : public LibertyInterfaceComponent, public ColourSchemeComponent, 
  public RTreeViewObserver, public RPopUpMenuObserver, public RTextEntryFieldObserver 
{

public:

  /** Different styles to draw the grid. */
  enum gridStyles
  {
    NO_GRID,
    GRID_LINES,
    DOTTED_GRID
  };

  /** Enumeration of the different actions that are possible on tne selection. */
  enum actionsOnSelection
  {
    EDIT_NAME,
    SAVE_CONTAINER,
    EXPORT_TO_CODE,
    DELETE_SELECTION,
    SET_POLYPHONIC,
    SET_MONOPHONIC,
    CONTAINERIZE,
    UNCONTAINERIZE
    // duplicate, copy, minimize numInputs, ...
  };

  enum actionsOnFreeArea
  {
    LOAD_CONTAINER = 0 
    //LOAD_CONTAINER = romos::ModuleTypeRegistry::NUM_MODULE_TYPES + 1  
    // preliminary, to avoid clashes between module type-indices and additional possible actions 
    // which are shown in the same TreeView
  };


  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor.  */  
  //ModularBlockDiagramPanel(CriticalSection *newPlugInLock, LibertyAudioModule* newLibertyModuleToEdit);
  ModularBlockDiagramPanel(LibertyInterfaceMediator *interfaceMediatorToUse);

  /** Destructor.  */  
  virtual ~ModularBlockDiagramPanel();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Allows an outlying class to inform this panel about the available screen estate, so we can 
  adjust the size here to at least fill this available space. */
  virtual void setAvailabeSizeForCanvas(int w, int h);

  /** Updates the bounds of this component. The passed values are interpreted as minimum required 
  bounds. The function itself will possibly extend these bounds according to the space needed to 
  accomodate for all the modules that we must draw. */
  virtual void updateCanvasBounds(int x, int  y, int w, int h); 

  //-----------------------------------------------------------------------------------------------
  // inqiury:

  /** Informs whether the passed module is currently selected or not */
  bool isModuleSelected(romos::Module *module) const 
  { 
    return rosic::containsElement(selectedModules, module);
    //return selectedModules.hasElement(module); 
  }

  virtual int getMinXInPixels() const;
  virtual int getMinYInPixels() const;
  virtual int getMaxXInPixels() const;
  virtual int getMaxYInPixels() const;
  //virtual int getMinXInPinDistances() const;
  //virtual int getMinYInPinDistances() const;
  //virtual int getMaxXInPinDistances() const;
  //virtual int getMaxYInPinDistances() const;

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void mouseExit(       const MouseEvent &e);
  virtual void mouseMove(       const MouseEvent &e);
  virtual void mouseDown(       const MouseEvent &e);
  //virtual void mouseDoubleClick(const MouseEvent &e);
  virtual void mouseDrag(       const MouseEvent &e);
  virtual void mouseUp(         const MouseEvent &e);

  virtual void paint(Graphics &g);
  virtual void paintOverChildren(Graphics &g);

  //virtual void resized();
  //virtual void updateDiagram();

  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged);

  //virtual void rDialogBoxChanged(RDialogBox* dialogBoxThatHasChanged);
  //virtual void rDialogBoxOKClicked(RDialogBox* dialogBoxThatWantsToAcceptAndLeave);
  //virtual void rDialogBoxCancelClicked(RDialogBox* dialogBoxThatWantsToBeCanceled);

  virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged);

  virtual void treeNodeClicked(RTreeView *treeView, RTreeViewNode *nodeThatWasClicked, 
    const MouseEvent &mouseEvent, int clickPosition);
  virtual void treeNodeChanged(RTreeView *treeView, RTreeViewNode *nodeThatHasChanged);
  //virtual void drawCoordinateSystem(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

  /** Callback that is called from the mediator when another GUI component has changed something 
  that might be interest here. */
  void mediatorHasSentNotification(MediatedColleague *originatingColleague, int messageCode = 0, 
    void* messageData = nullptr) override;

  //virtual void modulePropertiesChanged(romos::Module* module, int whichProperty);

protected:

  /** Open a popup menu to let the user insert modules. */
  virtual void openModuleInsertionMenu(int x, int y);

  /** Opens a popup menu to let the user choose an action to perform on the set of currently 
  selected modules such as deleting, duplicating, copying, containerizing, etc. */
  virtual void openActOnSelectionMenu(int x, int y);

  /** Opens the dialog for showing and editing the selected module's properties - this should 
  be called only when a single module is currently selected, otherwise it would be confusing. 
  However, in any case, it shows the properties of the first module in our array of selected 
  modules. */
  //virtual void openPropertiesDialog();

  /** Assuming that a single module is selected, this function opens a text entry field which 
  allows the user to change the name of the selected module. The field will appear directly on 
  top of the module. */
  virtual void openModuleNameEntryField();

  /** Opens the dialog for loading (and iserting) containers that are stored on disk. */
  virtual void openContainerLoadDialog();

  /** Assuming that a single module is selected and this module is a container, this function 
  opens a dialog to allow the user to save the container to disk. */
  virtual void openContainerSaveDialog();

  /** Opens a dialog to export the selected module to a code fragment that can be used to 
  programmatically create the module. */
  virtual void openExportToCodeDialog();

  /** Inserts a module of the kind given by the identifier (@see romos::moduleIdentifiers) at the 
  given coordinates. */
  //virtual void insertModule(int moduleIdentifer, int xInPinDistances, int yInPinDistances);
    // old

  virtual void insertModule(const juce::String& typeName, int xInPinDistances, int yInPinDistances);
    // new


  /** Adds the passed module to our array of selected modules and takes care to also select all 
  connections to and from that module. */
  //virtual void addModuleToSelection(romos::Module *moduleToAdd);

  /** Removes the passed module from our array of selected modules (if present) and takes care to 
  also de-select all connections to and from that module unless these connections also involve 
  other modules that remain selected. */
  //virtual void removeModuleFromSelection(romos::Module *moduleToRemove);

  // make these functions generic so as to work also for the lasso-array...

  /** First, it de-selects all selected modules and connections and then it selects the passed 
  module (with its connection) only. It also notifies the mediator. */
  virtual void selectSingleModule(romos::Module *moduleToSelect);

  /** Adds the passed module to the passed array and takes care to also add all connections that 
  belong to the module to the passed array of connections. */
  virtual void addModuleWithConnectionsToArray(romos::Module *moduleToAdd, 
    std::vector<romos::Module*> &moduleArray, 
    std::vector<romos::AudioConnection> &connectionArray);

  /** Removes the passed module from the passed array (if present) and takes care to also remove 
  all connections belonging to the module from the passed connection-array unless these connections 
  also involve other modules that remain in the array. */
  virtual void removeModuleWithConnectionsFromArray(romos::Module *moduleToRemove, 
    std::vector<romos::Module*> &moduleArray, std::vector<romos::AudioConnection> &connectionArray);

  /** Sets the polyphony for the selected modules. */
  virtual void setPolyphonyForSelection(bool shouldBePolyphonic);

  /** Deletes all modules and connections that are currently selected. */
  virtual void deleteSelection();

  /** Puts all modules that are currently selected into a container module. */
  virtual void containerizeSelection();

  /** Extracts all container modules that are currently selected and puts their content modules as 
  direct child-modules of the focused one. */
  virtual void unContainerizeSelection();

  /** Minimizes the number of input pins of all selected modules by investigating which pins are 
  superfluos, reconfiguring the connections accordingly and deleting the now obsolete pins. */
  virtual void minimizeNumberOfInputs();

  /** Fills the availableModulesTreeView with all the available modules. */
  virtual void fillAvailableModulesTreeView();

  /** Updates our array of audio connections - we maintain them here in an array in order to 
  deteremine when a mouse is over a connection. */
  virtual void updateAudioConnectionArray();

  /** Draws a grid to show the alignment positions. */
  virtual void drawGrid(Graphics &g);

  /** Draws the actual flowchart diagram. */
  virtual void drawDiagram(Graphics &g);

  /** Draws the graphical representation of the passed module at its proper position on the 
  Graphics canvas. */
  virtual void drawModule(Graphics &g, romos::Module *moduleToDraw);

  /** Draws the input pins including their names of the given module. You should also pass the 
  rectangle into which the module is drawn. */
  virtual void drawInputPins(Graphics &g, romos::Module *module, 
    juce::Rectangle<int> moduleRectangle);

  /** Draws the output pins including their names of the given module. */
  virtual void drawOutputPins(Graphics &g, romos::Module *module, 
    juce::Rectangle<int> moduleRectangle);

  /** Draws the incoming connections for the passed module. */
  //virtual void drawIncomingConnectionsForModule(Graphics &g, romos::Module *module, Rectangle<int> moduleRectangle);
  virtual void drawIncomingConnectionsForModule(Graphics &g, romos::Module *module);

  /** Assigns the passed reference variables to the x, y, coordinates and width and height that 
  are desired to draw the rectangle for the module on the block diagram. The includingPins flag 
  indicates whether or not the extended width of the module including the pins should be returned 
  (the pins stick out a bit). */
  virtual void getRectangleForModuleInPixels(romos::Module *module, int &x, int &y, int &w, int &h, 
    bool includingPins) const;

  /** Returns the rectangle for the module on the block diagram. The includingPins flag indicates 
  whether or not the extended width of the module including the pins should be returned (the pins 
  stick out a bit). */
  virtual juce::Rectangle<int> getRectangleForModuleInPixels(romos::Module *module, 
    bool includingPins) const;

  /** Returns the height of the rectangle in which the title of the module appears - this may vary, 
  depending on the module - some have no title at all. */
  virtual int getModuleTitleHeightInPixels(romos::Module *module) const; // { return 2*t + 2*m + hn; }

  /** Returns true if the rendering of the passed module is inside the passed rectangle, false 
  otherwise. */
  virtual bool isModuleInsideRectangle(romos::Module *module, juce::Rectangle<int> rectangle, 
    bool includingPins) const;

  /** Returns the offset by which the module must be shifted upward such that the midpoint of the 
  first pin exactly aligns with a grid-line. \todo make this depend on the particular module - 
  some have no titles, etc. */
  virtual int getOffsetY(romos::Module *module) const 
  { 
    return getModuleTitleHeightInPixels(module) + m + (int) floor(0.5*smallFontHeight); 
  }

  /** Returns the width that is required to properly draw the module. The includingPins flag 
  indicates whether or not the extended width of the module including the pins should be returned
  (the pins stick out a bit). */
  virtual int getRequiredWidthForModuleInPixels(romos::Module *module, bool includingPins) const;

  /** Returns the height that is required to properly draw the module (including title) with all 
  its input and output pins. */
  virtual int getRequiredHeightForModuleInPixels(romos::Module *module) const;

  /** Returns the bounds (in pixels) of a pin of some Module. You need to specify the kind of the 
  pin (AUDIO/EVENT), its direction (INCOMING/OUTGOING), and its index. You also need to pass the 
  module in question and the rectangle in which it is drawn.  */
  virtual juce::Rectangle<int> getPinBounds(int kindOfPin, int direction, int pinIndex, 
    romos::Module *module, juce::Rectangle<int> moduleRectangle) const;

  /** Returns the center coordinates (in pixels) of a pin of some Module. You need to specify the 
  kind of the pin (AUDIO/EVENT), its direction (INCOMING/OUTGOING), and its index. You also need 
  to pass the module in question and the rectangle in which it is drawn. The reference parameters 
  xPin and yPin will then be assigned to the center coordinates of the pin. */
  virtual void getPinCenterCoordinates(int kindOfPin, int direction, int pinIndex, 
    romos::Module *module, juce::Rectangle<int> moduleRectangle, int &xPin, int &yPin) const;

  /** Returns a Line object that represents the line to be drawn (in pixels) for the passed 
  conncetion. */
  virtual juce::Line<float> getLineForConnection(romos::AudioConnection connection) const;

  /** Returns an array with all modules inside the rectangle (given in in pixels). */
  virtual std::vector<romos::Module*> getModulesInRectangle(juce::Rectangle<int> rectangle) const;

  /** Returns a pointer to the module that is at the passed position or NULL if none. The 
  considerPins flag indicates whether or not the extended width of the module including the pins 
  should be considered (the pins stick out a bit). */
  virtual romos::Module* getModuleAtPixels(int x, int y, bool considerPins) const;

  /** Returns a pointer to the connection that is at or near the passed position or NULL if none. 
  \todo factor out baseclass Connection to treat event connections the same way. */
  virtual romos::AudioConnection getConnectionAtPixels(int x, int y) const;

  /** Assigns the 3 reference variables to the kind, direction (@see romos::connectionKinds and 
  romos::connectionDirections ) and index identifiers of the pin that is found at the pixel 
  location given by x, y. The function assumes that the module that sits at the given location and 
  its rectangle have already been found and are passed via the "module" and moduleRectangle 
  parameters. If there is no pin to be found at the given pixel location, the reference parameters 
  will be assigned to -1 and the function will return false, true otherwise. */
  virtual bool getPinPropertiesAtPixels(int x, int y, romos::Module* module, 
    juce::Rectangle<int> moduleRectangle, int &kindOfPin, int &directionOfPin, 
    int &indexOfPin) const;

  /** Converts a coordinate given in pin-distances to the corresponding coordinate in pixels. */
  virtual int inPixels(int pinDistances) const { return pinDistance * pinDistances; }

  /** Converts a coordinate given in pixels to the corresponding coordinate in pin-distances. */
  virtual int inPinDistances(int pixels) const 
  { 
    return roundToInt(pixels / (double) pinDistance); 
  }

  /** Returns a pixel value (x- or y-coordinate) that is at the nearest grid position with respect 
  to the passed coordinate. */
  virtual int snapPixelPositionToGrid(int pixelPosition) const 
  { 
    return inPixels(inPinDistances(pixelPosition)); 
  }

  //-----------------------------------------------------------------------------------------------
  // data members:

  RTreeView               *availableModulesTreeView;
  RTreeViewNode           *treeRootNode;
  RPopUpMenu              *actOnSelectionMenu;
  RTextEntryField         *nameEntryField;

  romos::AudioConnection  tmpAudioConnection;
  //romos::AudioConnection  *tmpAudioConnection;
  //romos::EventConnection  *tmpEventConnection;

  std::vector<romos::Module*> selectedModules;
  std::vector<romos::Module*> modulesInLasso;
  std::vector<romos::AudioConnection> selectedAudioConnections;
  std::vector<romos::AudioConnection> audioConnectionsInLasso;
  std::vector<romos::AudioConnection> allAudioConnections;

  int  mouseDownX, mouseDownY;
  int  selectionOffsetX, selectionOffsetY; // to draw the selected modules at an offsetted position during drag (in pixels)
  int  availableWidth, availableHeight;    // we store these to let the canvas always fill this size
  int  gridStyle;

  int m, t, s, bigFontHeight, normalFontHeight, smallFontHeight; 

  // drawing parameters: margin, thickness, pin stickout, normal height, small height
  // make them static const - or maybe it's good to have them tweakable at runtime...

  int pinDistance;
  int arrowLength, arrowHeadLength;  // for the arrows in the I/O modules

  Colour wireColour;
  Colour highlightColour;  // preliminarily used for selections, etc.

  jura::RectangleComponent *lassoComponent;
  jura::RectangleComponent *pinHighlighter; 
  // later use a custom component that can assume other shapes (ellipses, etc.) - maybe define 
  // ShapeComponent

  // pointers for coveniently accessing the bitmap fonts:
  const BitmapFont* bigFont    = &BitmapFontRoundedBoldA10D0::instance;
  const BitmapFont* normalFont = &BitmapFontRoundedBoldA9D0::instance;
  const BitmapFont* smallFont  = &BitmapFontRoundedA7D0::instance;

  juce_UseDebuggingNewOperator;
};

// we need ComponentScrollContainer - then uncomment...

//=================================================================================================

/** This is the editor for the whole modular synthesizer. */

class LibertyEditor : public AudioModuleEditor //, public LibertyInterfaceComponent  // public PolyphonicInstrumentEditor
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  LibertyEditor(CriticalSection *newPlugInLock, LibertyAudioModule* newLibertyAudioModule);

  /** Denstructor. */
  virtual ~LibertyEditor();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  //virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  /** Creates the proprties editor for the currently selected module. */
  //virtual void createPropertiesEditorForSelectedModule();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  LibertyAudioModule *modularSynthAudioModule;

  LibertyInterfaceMediator *interfaceMediator;

  ModularStructureTreeView     *structureTreeView;   
  ModularBlockDiagramPanel     *blockDiagramPanel;
  ComponentScrollContainer     *diagramScrollContainer;
  ModulePropertiesEditorHolder *moduleEditorHolder;

  //ModulePropertiesEditor   *modulePropertiesEditor;

  //RRadioButton  *structureButton, *interfaceButton; // *hybridButton;
  //RRadioButtonGroup radioButtonGroup;

  juce_UseDebuggingNewOperator;
};

#endif
