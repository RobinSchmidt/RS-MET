#ifndef jura_AutomatableWidget_h
#define jura_AutomatableWidget_h


class AutomatableWidget;
class AutomatableSlider;

/** A component for setting up the modulations of some ModulationTarget. */

class JUCE_API rsModulationSetup : public ColourSchemeComponent, public RButtonListener, 
  public rsDeletionRequester, public RPopUpMenuObserver
{

public:

  /** Constructor. */
  rsModulationSetup(AutomatableWidget* widgetToModulate);

  /** Destructor. */
  virtual ~rsModulationSetup();

  // callbacks:
  virtual void paint(Graphics& g) override;
  virtual void resized() override;
  virtual void rButtonClicked(RButton *button) override;
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;

protected:

  /** Adds a modulation connection between the ModulationTarget (underlying our widget member) and 
  the source with given index among the available and not yet connected sources. */
  void addConnection(int index);

  /** Creates or removes sliders for the modulation amounts according to the existing 
  connections. */
  void updateAmountSliderArray();

  /** Shows the popup menu with the available (and not yet connected) ModulationSources. */
  void showConnectableSourcesPopUp();

  /** Returns true, if we have a slider associated with the passed parameter. */
  bool hasSlider(MetaControlledParameter* p);

  /** Adds a slider for the passed parameter to our amountSliders array. */
  void addSliderFor(MetaControlledParameter* p);

  /** Updates the size in order to provide space for all required widgets. */
  void updateSize();

  // pointer to our owner:
  AutomatableWidget* widget;

  // owned widgets:
  RTextField* modulationsLabel;
  std::vector<AutomatableSlider*> amountSliders;
  RButton* addButton;
  RClickButton* closeButton;
  RPopUpMenu *sourcesPopUp = nullptr; // created when needed the first time

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

  // overrides:
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;


protected:

  /** Enumeration of the identifiers to used as return-values for the right-click popup menu. */
  enum rightClickPopUpItemIdentifiers
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


    //MODULATOR_CONNECT // maybe factor out modulation related stuff into ModulatableSlider class
                      // ...hmm but maybe not, because we may wnat to modulate other kinds of 
                      // widgets like DraggableNumber
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
  bool popUpIsOpen = false;
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


#endif   
