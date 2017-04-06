#ifndef jura_AutomatableWidget_h
#define jura_AutomatableWidget_h

/** RWidget subclass that adds automation facilities. 
Maybe move into jura_processors module  */

class JUCE_API AutomatableWidget : virtual public RPopUpMenuObserver
{

public:

  AutomatableWidget(RWidget *widgetToWrap);

  virtual ~AutomatableWidget();

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

  /** Returns a pointer to the AutomatableParameter object that is assigned to the wrapped 
  widget.  */
  AutomatableParameter* getParameter();


  RPopUpMenu *rightClickPopUp = nullptr; // object created when it's needed for the 1st time
  RWidget *wrappedWidget;                // widget that is being made automatable

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableWidget)
};

//=================================================================================================

/** A slider that can be automated via MIDI or the host automation system (using setParameter()) */

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

};


#endif   
