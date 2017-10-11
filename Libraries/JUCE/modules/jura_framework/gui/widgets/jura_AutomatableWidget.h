#ifndef jura_AutomatableWidget_h
#define jura_AutomatableWidget_h

class rsModulationConnectionWidget;

/** A component for setting up the modulations of some ModulationTarget. */

class JUCE_API rsModulationSetup : public ColourSchemeComponent, public RButtonListener, 
  public rsDeletionRequester, public rsGarbageCollector, public RPopUpMenuObserver, 
  public RTextEntryFieldObserver
{

public:

  /** Constructor. You must pass the widget that you want to assign modulations to (actually, they
  will be assigned to its underlying Parameter assumed to be a ModulatableParameter) and you must 
  also pass a pointer to the MetaParameterManager object that should be used for attaching a
  meta-parameter to the modulation-depths of the connections to be made. */
  rsModulationSetup(AutomatableWidget* widgetToModulate, MetaParameterManager* metaManager);

  /** Destructor. */
  virtual ~rsModulationSetup();

  // callbacks:
  //virtual void paint(Graphics& g) override;
  virtual void resized() override;
  virtual void rButtonClicked(RButton *button) override;
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;
  virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged) override;

protected:

  /** Adds a modulation connection between the ModulationTarget (underlying our widget member) and 
  the source with given index among the available and not yet connected sources. */
  void addConnection(int index);

  /** Removes a modulation connection between the ModulationTarget (underlying our widget member) 
  and the source with given index among the connected sources. */
  void removeConnection(int index);

  /** Creates or removes widget-sets (depth-sliders, etc.) for the modulation the existing 
  moduldation connections. */
  void updateConnectionWidgetsArray();

  /** Shows the popup menu with the available (and not yet connected) ModulationSources. */
  void showConnectableSourcesPopUp();

  /** Shows the popup menu with the connected ModulationSources (for removal). */
  void showRemovableSourcesPopUp();

  /** Returns true, if we have a slider associated with the passed parameter. */
  bool hasSlider(MetaControlledParameter* p);

  /** Adds a slider for the passed parameter to our amountSliders array. You must also pass a 
  pointer to the ModulationConnection whose depth this parameter controls so the slider's popup
  menu can access it in order to set up some other aspects of the connection as well. */
  void addWidgetsForConnection(ModulationConnection* c); 

  /** Removes the widgets for the connection with given index. */
  void removeWidgetsForConnection(int index); 

  /** Clears the array of connections widgets. */
  void clearConnectionWidgets();

  /** Updates the size in order to provide space for all required widgets. */
  void updateSize();

  // functions for getting/setting the clipping limits:
  void setClipMin(double newMin);
  void setClipMax(double newMax);
  double getClipMin();
  double getClipMax();

  AutomatableWidget* widget;         // our owner widget
  MetaParameterManager* metaManager; // used for meta-controlling modulation amounts
  //ModulationManager*    modManager;  // needed for debugging

  // owned widgets:
  RTextField* modulationsLabel;
  std::vector<rsModulationConnectionWidget*> connectionWidgets;

  RButton *addButton, *removeButton;
  RClickButtonNotifyOnMouseUp* closeButton;
  RLabeledTextEntryField *clipMinField, *clipMaxField;

  RPopUpMenu *connectableSourcesPopUp = nullptr; // created when needed the first time
  RPopUpMenu *removableSourcesPopUp   = nullptr; // ditto

  static const int sliderHeight = 16, sliderDistance = 4;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsModulationSetup)
};


//=================================================================================================

/** RWidget subclass that adds automation facilities either via MIDI or the host automation system 
using setParameter(). If you want to make a widget automatable, derive it from some widget class 
and also from this class, for example, like: 
class JUCE_API AutomatableSlider : public RSlider, public AutomatableWidget
\todo: maybe move into jura_processors module (it's relevant only for audio plugins). 

maybe rename to something that fits better (we now have also lumped in the modulation stuff)

*/

class JUCE_API AutomatableWidget : virtual public RPopUpMenuObserver, public rsDeletor
{

public:

  AutomatableWidget(RWidget *widgetToWrap);

  virtual ~AutomatableWidget();

  /** Returns true if the popup menu is currently open, false otherwise. */
  bool isPopUpOpen();

  /** Tries to cast the Parameter that is underlying the wrapped widget into an 
  AutomatableParameter and returns the pointer to it. Note that this may return a nullptr, when 
  the Parameter is not of type AutomatableParameter. */
  AutomatableParameter* getAutomatableParameter();

  /** Similar to @see getAutomatbleParameter. */
  MetaControlledParameter* getMetaControlledParameter();

  /** Similar to @see getAutomatbleParameter. */
  ModulatableParameter* getModulatableParameter();

  /** Returns a pointer to the MetaParameterManager that is used by the underlying parameter, if 
  any - a nullptr otherwise. */
  MetaParameterManager* getMetaParameterManager();

  // overrides:
  void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;
  void rPopUpMenuDismissed(RPopUpMenu* menuThatwasDismissed) override;


protected:

  /** Enumeration of the identifiers to used as return-values for the right-click popup menu. */
  enum popUpIds
  {
    ENTER_VALUE = 1,  // nope - these two should be available in the slider baseclass already
    DEFAULT_VALUE,    

    MIDI_ASSIGN,
    MIDI_LEARN,
    MIDI_MIN,
    MIDI_MAX,
    MIDI_REVERT,

    META_ATTACH,
    META_DETACH,

    MODULATION_SETUP
  };

  /** Clears the popup-menu and then calls createPopUpMenuItems() */
  virtual void updatePopUpMenu();

  /** Populates the right-click popup menu with items, according to the settings of this RSlider. */
  virtual void addPopUpMenuItems();
  // rename to populateRightClickPopupMenu...or something

  /** Adds the MIDI related items to the popoup menu (if applicable). */
  virtual void addPopUpMidiItems();

  /** Adds the MetaParameter related items to the popup menu (if applicable). */
  virtual void addPopUpMetaItems();

  /** Adds the Modulation related items to the popup menu (if applicable). */
  virtual void addPopUpModulationItems();

  /** Opens the PopupMenu that appears on right clicks. */
  virtual void openRightClickPopupMenu();

  /** Closes the popup menu. */
  virtual void closePopUp();


  virtual void showModulationSetup();

  /** This is called from the modulation setup popup window when the user clicks on its 
  close-button. We don't actually delete the object here, we just make it invisible, so we can 
  reuse the same object when the modulation setup is opened again. It's deleted in our 
  destructor. */
  virtual void deleteObject(rsDeletionRequester* objectToDelete) override;

  RWidget *wrappedWidget;                 // widget that is being made automatable
  bool popUpIsOpen = false;               // maybe we can get rid of this?
  RPopUpMenu *rightClickPopUp = nullptr;  // object created when it's needed for the 1st time
  rsModulationSetup* modSetup = nullptr;  // ditto for modulation setup

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableWidget)
};

//=================================================================================================

class JUCE_API AutomatableSlider : public RSlider, public AutomatableWidget
{

public:

  AutomatableSlider();
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;
  virtual void mouseDown(const MouseEvent& e) override;

protected:

  virtual void addPopUpMenuItems() override;
  virtual void addPopUpEnterValueItem();
  virtual void addPopUpDefaultValueItems();

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableSlider)
};


class JUCE_API AutomatableComboBox : public RComboBox, public AutomatableWidget
{

public:

  AutomatableComboBox();
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void parameterChanged(Parameter* p) override;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableComboBox)
};

class JUCE_API AutomatableButton : public RButton, public AutomatableWidget
{

public:

  AutomatableButton(const juce::String& buttonText = juce::String::empty);
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void parameterChanged(Parameter* p) override;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableButton)
};

//=================================================================================================

/** A special slider subclass specifically for being used for the depth-sliders in a modulation
setup. They need some special additional options in their popup menu such as facilities to set 
min/max values and the modulation mode. */

class JUCE_API rsModulationDepthSlider : public AutomatableSlider
{

public:

  /** Constructor. You must pass a pointer to the ModulationConnection object on which this slider
  operates such that we can use our popup menu to set up the modulation mode (absolute vs 
  relative), too. */
  rsModulationDepthSlider(ModulationConnection* connection) : modConnection(connection) {}

  virtual ~rsModulationDepthSlider() {}

  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;

  virtual void addPopUpMenuItems() override;
  virtual void addPopUpMinMaxAndModeItems();

protected:

  // for getting and setting the min/max values of the mod-depth slider and for the clipping of the
  // modulated value and getting/setting relative mode (maybe these should be inlined):
  inline double getModDepthMin() { return assignedParameter->getMinValue(); }
  inline double getModDepthMax() { return assignedParameter->getMaxValue(); }
  inline int    getModMode()     { return modConnection->getMode(); }

  inline void setModDepthMin(double newMin)  { assignedParameter->setMinValue(newMin); }
  inline void setModDepthMax(double newMax)  { assignedParameter->setMaxValue(newMax); }
  inline void setModMode(    int    mode)    { modConnection->setMode(mode); }

  /** Additional item ids for this subclass. */
  enum popUpIds2
  {
    MOD_DEPTH_MIN = popUpIds::MODULATION_SETUP+1,
    MOD_DEPTH_MAX,
    MOD_MODE_ABSOLUTE,
    MOD_MODE_RELATIVE,
    MOD_MODE_EXPONENTIAL
  };

  ModulationConnection* modConnection = nullptr; // needs to be assigned in constructor

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsModulationDepthSlider)
};

//=================================================================================================

/** A class that encapsulates the depth-slider and remove-button for a ModulationConnection. 

\todo:
-add an on/off button, maybe also solo - i think, this needs a new member in ModulationConnection
 isActive ...hmm...but this will cause overhead...we'll see
*/

class JUCE_API rsModulationConnectionWidget : public RWidget, public rsDeletionRequester
{

  friend class rsModulationSetup;

public:

  rsModulationConnectionWidget(ModulationConnection* connection, rsGarbageCollector* deletor);

  virtual ~rsModulationConnectionWidget() {}

  void resized() override;
  void paint(Graphics& g) override {}

protected:

  rsModulationDepthSlider* depthSlider;
  RClickButtonNotifyOnMouseUp* removeButton;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsModulationConnectionWidget)
};

//=================================================================================================

/** A slider subclass specifically for sliders that have modulatable parameters assigned. It may
paint itself differently, depending on the modulations that are applied to its underlying 
parameter. */

class JUCE_API ModulatableSlider : public AutomatableSlider, public ModulationTargetObserver
{

public:

  ModulatableSlider() {}
  ~ModulatableSlider();
  void modulationsChanged() override;
  void assignParameter(Parameter* parameterToAssign) override;
  void paint(Graphics& g) override;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulatableSlider)
};


#endif   
