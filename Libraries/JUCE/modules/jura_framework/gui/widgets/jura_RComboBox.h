#ifndef jura_RComboBox_h
#define jura_RComboBox_h

// todo: drag out the Item class 

class RComboBox;

/** Observer class for RComboBoxes */

class JUCE_API RComboBoxObserver
{

public:

  virtual ~RComboBoxObserver() {}

  /** Callback that gets called when the user has selected a new item in a combobox. */
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) = 0;

};

//=================================================================================================

/**

A component that lets the user choose from a drop-down list of choices. The combo-box has a list of 
text strings, each with an associated id number, that will be shown in the drop-down list when the 
user clicks on the component. The currently selected choice is displayed in the combo-box, and this 
can either be read-only text, or editable. To find out when the user selects a different item or 
edits the text, you can register a RComboBoxObserver to receive callbacks. @see RComboBoxObserver

*/

class JUCE_API RComboBox : public RTextField, public RPopUpMenuObserver, public RPopUpOwner
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates a combo-box. On construction, the text field will be empty, so you should call the 
  setSelectedId() or setText() method to choose the initial value before displaying it. */
  RComboBox (const juce::String& componentName = String::empty);

  /** Destructor. */
  virtual ~RComboBox();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Adds an item to be shown in the drop-down list. */
  virtual void addItem(int itemResultId, const juce::String& itemText, bool isEnabled = true, 
    bool isTicked = false);

  /** This allows items in the drop-down list to be selectively disabled. When you add an item, 
  it's enabled by default, but you can call this method to change its status. If you disable an 
  item which is already selected, this won't change the current selection - it just stops the user 
  choosing that item from the list. */
  void setItemEnabled(int index, bool shouldBeEnabled);
    // \todo: implement this

  /** Changes the text for an existing item. */
  void setItemText(int index, const juce::String& newText);

  /** Removes all the items from the drop-down list. If this call causes the content to be 
  cleared, then a change-message will be broadcast unless dontSendChangeMessage is true. 
  @see addItem, removeItem, getNumItems */
  void clear(const bool dontSendChangeMessage = false);

  /** Selects one of the items and optionally sends out a notification to our observers.  */
  void selectItemByIndex(int indexToSelect, bool sendNotification);

  /** Assuming that the passed text matches one of the item-texts, this function sets the contents 
  of the combo-box's text field to the respective item. If there is no matching item, a debug-break 
  will be issued. */
  void selectItemFromText(const juce::String& textToSelect, bool sendNotification);

  /** Overriden from RWidget - calls selectItemFromText. */
  void setStateFromString(const juce::String &stateString, bool sendChangeMessage);

  /** Registers an observer that will be called when the box's content changes. */
  virtual void registerComboBoxObserver(RComboBoxObserver* const observerToRegister);

  /** Deregisters a previously-registered observer. */
  virtual void deRegisterComboBoxObserver(RComboBoxObserver* const observerToDeRegister);

  /** Overrides a RWidget::assignParameter in order to set up the available options to those that 
  are defined in the Parameter. */
  virtual void assignParameter(Parameter* parameterToAssign);

  /** Sets the maximum width and height for the popup menu that appears on clicking. */
  virtual void setMaxPopUpSize(int newMaxWidth, int newMaxHeight)
  {
    maxPopUpWidth = newMaxWidth; maxPopUpHeight = newMaxHeight;
  }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of items that have been added to the list. Note that this doesn't include 
  headers or separators. */
  //virtual int getNumItems() const { return popUpMenu->getNumItems(); }
  // too unspecific, deprecated in favor of:

  /** Returns the number of items that are in the 1st level of the tree. For a flat array of 
  choices, this equals the number of possible choices. */
  virtual int getNumTopLevelItems() const { return popUpMenu->getNumTopLevelItems(); }
   // this returns the same number that was formerly returned by getNumItems

  /** Returns the number of selectable items, that is the number of leaf nodes in the tree. */
  virtual int getNumSelectableItems() const { return popUpMenu->getNumSelectableItems(); }


  /** Returns the text for one of the items in the list. Note that this doesn't include headers or 
  separators. @param index  the item's
  index from 0 to (getNumItems() - 1) */
  virtual const juce::String getItemText(const int index) const;

  /** Returns the index of the item that is currently active, -1 if none. */
  //virtual int getSelectedItemIndex() const { return currentIndex; }

  /** Returns a (pointer to) the selected item if any, otherwise a NULL pointer will be returned. */
  virtual RTreeViewNode* getSelectedItem() const { return popUpMenu->getSelectedItem(); }

  /** Returns the identifier of the selected item if any, otherwise -1. */
  virtual int getSelectedItemIdentifier() const;

  /** Returns the text of the selected item if any, otherwise "nothing selected" will be 
  returned . */
  virtual juce::String getSelectedItemText() const;

  /** Overriden from RWidget - returns the state (defined as the currently selected text) as 
  string. */
  virtual juce::String getStateAsString() const;

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overriden ffrom RPopUpOwner. */
  virtual void rPopUpDismissedByClickOnOwner(ROwnedPopUpComponent *popUp);

  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged);

  virtual void updateWidgetFromAssignedParameter(bool sendNotification = false);
  virtual void mouseDown(const MouseEvent&);
  virtual void paint(Graphics&);

protected:

  /** Sends out the notification to our observers that the selected item was changed. */
  virtual void sendComboBoxChangeNotifications();

  /** Opens the popup-menu with the items, allowing the user to select one. */
  virtual void openPopUp();

  /** Opens the popup-menu with the items, if currently open. */
  virtual void closePopUp();

  RPopUpMenu *popUpMenu;
  juce::Array<RComboBoxObserver*> comboBoxObservers;

  bool dontOpenPopUpOnNextMouseClick;

  int maxPopUpWidth = INT_MAX, maxPopUpHeight = INT_MAX;

private:

  RComboBox (const RComboBox&);
  const RComboBox& operator= (const RComboBox&);
  virtual void setText(const juce::String &newText) { RTextField::setText(newText); }

  juce_UseDebuggingNewOperator
};

//=================================================================================================
// class RNamedComboBox:

/** A ComboBox that has a name field (realized as RLabel) in front of it or above.

\todo: override setDescriptionField, colour-scheme etc. to pass the values through to the child 
widgets...or maybe don't let this class have child-components (i.e. the RTextField) - instead, 
handle drawing differently in paint

*/

class JUCE_API RNamedComboBox : public RComboBox
{

public:

  enum nameLabelPositions
  {
    INVISIBLE,
    LEFT_TO_BOX,
    ABOVE_BOX
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  RNamedComboBox(const juce::String& componentName, const juce::String& comboBoxName);
  virtual ~RNamedComboBox();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the name of the combobox - this is the text that appears in front of (or above) the 
  actual combobox. */
  virtual void setComboBoxName(const juce::String& newName) { nameLabel->setText(newName); }

  /** Sets the width (in pixels) of the name-label which appears in front of or above the box. */
  virtual void setNameLabelWidth(int newWidth)
  {
    nameLabelWidth = newWidth;
    resized();
    //repaint(); 
  }

  /** Sets the of the name-label relative to the box @see: nameLabePositions. */
  virtual void setNameLabelPosition(int newPosition)
  {
    nameLabelPosition = newPosition;
    repaint();
  }

  /** Overrides RWidget's setColourScheme in order to pass the new scheme through to the label. */
  virtual void setColourScheme(const WidgetColourScheme& newColourScheme);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the width (in pixels) of the name-label which appears in front of or above the box. */
  virtual int getNameLabelWidth() const
  {
    return nameLabelWidth;
  }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void paint(Graphics &g);
  virtual void resized();

protected:

  /** Overriden to alter positioning. */
  virtual void openPopUp();

  RTextField* nameLabel;
  int nameLabelWidth, nameLabelPosition;

private:

  /** Moved to private area, because it should not be used for this subclass. Instead, the 
  two-parameteric form of this function should be used. */
  //virtual void setColourScheme(const WidgetColourScheme& newColourScheme) {}

  juce_UseDebuggingNewOperator
};

#endif  
