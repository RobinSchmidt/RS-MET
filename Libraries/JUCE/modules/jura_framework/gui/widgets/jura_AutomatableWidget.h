#ifndef jura_AutomatableWidget_h
#define jura_AutomatableWidget_h

/** RWidget subclass that adds automation facilities. 
Maybe move into jura_processors module
*/

class JUCE_API AutomatableWidget : public RWidget, public RPopUpMenuObserver
{

public:

  AutomatableWidget();
  ~AutomatableWidget();

  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;


protected:

  /** Enumeration of the identifiers to used as return-values for the right-click popup menu. */
  enum rightClickPopUpItemIdentifiers
  {
    ENTER_VALUE = 1,
    DEFAULT_VALUE,
    MIDI_ASSIGN,
    MIDI_LEARN,
    MIDI_MIN,
    MIDI_MAX,
    MIDI_REVERT
  };

  /** Clears the popup-menu and then calls createPopUpMenuItems() */
  virtual void updatePopUpMenu();

  /** Populates the right-click popup menu with items, according to the settings of this RSlider. */
  virtual void addPopUpMenuItems();
  // rename to populateRightClickPopupMenu...or something

  // called from createPopUpMenuItems:
  virtual void addPopUpMidiItems();

  /** Opens the PopupMenu that appears on right clicks. */
  virtual void openRightClickPopupMenu();

  /** Opens a modal field for manually entering a value and returns the value entered. */
  virtual double openModalNumberEntryField();

  RPopUpMenu *rightClickPopUp = nullptr; // object created when it's needed for the 1st time

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableWidget)
};

//=================================================================================================



#endif   
