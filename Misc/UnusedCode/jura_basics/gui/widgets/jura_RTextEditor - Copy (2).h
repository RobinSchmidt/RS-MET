#ifndef jura_RTextEditor_h
#define jura_RTextEditor_h

const int textEditorLineSpacing = 2;

class RTextEditor;
class RTextHolderComponent;

//=================================================================================================

/** Receives callbacks from a RTextEditor component when it changes. */

class JUCE_API RTextEditorListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Destructor. */
  virtual ~RTextEditorListener()  {}

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Called when the user changes the text in some way. */
  virtual void rTextEditorTextChanged(RTextEditor& editor) = 0;

  /** Called when the user presses the return key. */
  virtual void rTextEditorReturnKeyPressed(RTextEditor& editor) = 0;

  /** Called when the user presses the escape key. */
  virtual void rTextEditorEscapeKeyPressed(RTextEditor& editor) = 0;

  /** Called when the text editor loses focus. */
  virtual void rTextEditorFocusLost(RTextEditor& editor) = 0;

};

//=================================================================================================

/** A component containing text that can be edited. A RTextEditor can either be in single- or 
multi-line mode, and supports mixed fonts and colours. @see RTextEditorListener, Label

\todo remove unneeded aspects: password - done
\todo factor multiline/scrolling functionality into subclass
\todo factor undo/redo into subclass
\todo make all text uniform ...mmm ...maybe not
\todo perhaps get rid of the permissibleCharacters restirction (allow for all characters)

*/

class JUCE_API RTextEditor : public RWidget //, public SettableTooltipClient
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates a new, empty text editor. */
  RTextEditor(const juce::String& componentName = juce::String::empty);

  /** Destructor. */
  virtual ~RTextEditor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Puts the editor into either multi- or single-line mode. By default, the editor will be in 
  single-line mode, so use this if you need a multi-line editor. See also the 
  setReturnKeyStartsNewLine() method, which will also need to be turned on if you want a
  multi-line editor with line-breaks. @see isMultiLine, setReturnKeyStartsNewLine */
  void setMultiLine (const bool shouldBeMultiLine, const bool shouldWordWrap = true);

  /** Changes the behaviour of the return key. If set to true, the return key will insert a 
  new-line into the text; if false it will trigger a call to the 
  RTextEditorListener::textEditorReturnKeyPressed() method. By default this is set to false, and 
  when true it will only insert new-lines when in multi-line mode (see setMultiLine()). */
  void setReturnKeyStartsNewLine (const bool shouldStartNewLine);

  /** Indicates whether the tab key should be accepted and used to input a tab character, or 
  whether it gets ignored. By default the tab key is ignored, so that it can be used to switch 
  keyboard focus between components. */
  void setTabKeyUsedAsCharacter (const bool shouldTabKeyBeUsed) throw();

  /** Changes the editor to read-only mode. By default, the text editor is not read-only. If you're 
  making it read-only, you might also want to call setCaretVisible (false) to get rid of the caret. 
  The text can still be highlighted and copied when in read-only mode. 
  @see isReadOnly, setCaretVisible */
  void setReadOnly (const bool shouldBeReadOnly);

  /** Makes the caret visible or invisible. By default the caret is visible. @see setCaretColour, 
  setCaretPosition */
  void setCaretVisible (const bool shouldBeVisible) throw();

  /** Enables/disables a vertical scrollbar. (This only applies when in multi-line mode). When the 
  text gets too long to fit in the component, a scrollbar can appear to allow it to be scrolled. 
  Even when this is enabled, the scrollbar will be hidden unless it's needed. By default the 
  scrollbar is enabled. */
  void setScrollbarsShown (bool shouldBeEnabled) throw();

  /** Allows a right-click menu to appear for the editor. (This defaults to being enabled). If 
  enabled, right-clicking (or command-clicking on the Mac) will pop up a menu of options such as 
  cut/copy/paste, undo/redo, etc. */
  void setPopupMenuEnabled (const bool menuEnabled) throw();

  /** Sets the font to use for newly added text. This will change the font that will be used next 
  time any text is added or entered into the editor. It won't change the font of any existing 
  text - to do that, use applyFontToAllText() instead. @see applyFontToAllText */
  void setFont(const BitmapFont* newFont) throw();

  /** Applies a font to all the text in the editor. This will also set the current font to use for 
  any new text that's added. */
  void applyFontToAllText(const BitmapFont* newFont);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns true if the editor is in multi-line mode. */
  bool isMultiLine() const throw();

  /** Returns the value set by setReturnKeyStartsNewLine(). See setReturnKeyStartsNewLine() for 
  more info. */
  bool getReturnKeyStartsNewLine() const throw() { return returnKeyStartsNewLine; }

  /** Returns true if the tab key is being used for input. @see setTabKeyUsedAsCharacter */
  bool isTabKeyUsedAsCharacter() const throw() { return tabKeyUsed; }

  /** Returns true if the editor is in read-only mode. */
  bool isReadOnly() const throw();

  /** Returns true if the caret is enabled. @see setCaretVisible */
  bool isCaretVisible() const throw() { return caretVisible; }

  /** Returns true if scrollbars are enabled. @see setScrollbarsShown */
  bool areScrollbarsShown() const throw() { return scrollbarVisible; }

  /** Returns true if the right-click menu is enabled. @see setPopupMenuEnabled */
  bool isPopupMenuEnabled() const throw() { return popupMenuEnabled; }

  /** Returns true if a popup-menu is currently being displayed. */
  bool isPopupMenuCurrentlyActive() const throw() { return menuActive; }

  /** Returns the font that's currently being used for new text. @see setFont */
  const BitmapFont* getFont() const throw();

  /** If set to true, focusing on the editor will highlight all its text. (Set to false by 
  default). This is useful for boxes where you expect the user to re-enter all the text when they 
  focus on the component, rather than editing what's already there. */
  void setSelectAllWhenFocused (const bool b) throw();

  /** Sets limits on the characters that can be entered. */
  void setInputRestrictions(const int maxTextLength, 
    const juce::String& allowedCharacters = juce::String::empty) throw();

  /** When the text editor is empty, it can be set to display a message. This is handy for things 
  like telling the user what to type in the box - the string is only displayed, it's not taken to 
  actually be the contents of the editor. */
  void setTextToShowWhenEmpty(const juce::String& text, const Colour& colourToUse) throw();

  /** Changes the size of the scrollbars that are used. Handy if you need smaller scrollbars for a 
  small text box. */
  void setScrollBarThickness (const int newThicknessPixels);

  /** Shows or hides the buttons on any scrollbars that are used. */
  //void setScrollBarButtonVisibility (const bool buttonsVisible);

  /** Registers a listener to be told when things happen to the text. */
  void addListener (RTextEditorListener* const newListener) throw();

  /** Deregisters a listener. */
  void removeListener (RTextEditorListener* const listenerToRemove) throw();

  /** Returns the entire contents of the editor. */
  const juce::String getText() const throw();

  /** Returns a section of the contents of the editor. */
  const juce::String getTextSubstring(const int startCharacter, 
    const int endCharacter) const throw();

  /** Returns true if there are no characters in the editor. This is more efficient than calling 
  getText().isEmpty(). */
  bool isEmpty() const throw();

  /** Sets the entire content of the editor. */
  void setText(const juce::String& newText, const bool sendTextChangeMessage = true);

  /** Inserts some text at the current cursor position. If a section of the text is highlighted, it 
  will be replaced by this string, otherwise it will be inserted. To delete a section of text, you 
  can use setHighlightedRegion() to highlight it, and call insertTextAtCursor (String::empty). */
  void insertTextAtCursor (juce::String textToInsert);

  /** Deletes all the text from the editor. */
  void clear();

  /** Deletes the currently selected region, and puts it on the clipboard. */
  void cut();

  /** Copies any currently selected region to the clipboard. */
  void copy();

  /** Pastes the contents of the clipboard into the editor at the cursor position. */
  void paste();

  /** Moves the caret to be in front of a given character. */
  void setCaretPosition (const int newIndex) throw();

  /** Returns the current index of the caret. */
  int getCaretPosition() const throw();

  /** Selects a section of the text. */
  void setHighlightedRegion (int startIndex, int numberOfCharactersToHighlight) throw();

  /** Returns the first character that is selected. If nothing is selected, this will still return 
  a character index, but getHighlightedRegionLength() will return 0. @see setHighlightedRegion, 
  getHighlightedRegionLength */
  int getHighlightedRegionStart() const throw()  { return selectionStart; }

  /** Returns the number of characters that are selected. */
  int getHighlightedRegionLength() const throw() 
  { 
    return jmax (0, selectionEnd - selectionStart); 
  }

  /** Returns the section of text that is currently selected. */
  const juce::String getHighlightedText() const throw();

  /** Finds the index of the character at a given position. The co-ordinates are relative to the 
  component's top-left. */
  int getTextIndexAt (const int x, const int y) throw();

  /** Returns the total width of the text, as it is currently laid-out. This may be larger than the 
  size of the RTextEditor, and can change when the RTextEditor is resized or the text changes. */
  int getTextWidth() const throw();

  /** Returns the maximum height of the text, as it is currently laid-out. This may be larger than 
  the size of the RTextEditor, and can change when the RTextEditor is resized or the text 
  changes. */
  int getTextHeight() const throw();

  /** Changes the size of the gap at the top and left-edge of the editor. By default there's a gap 
  of 4 pixels. */
  void setIndents (const int newLeftIndent, const int newTopIndent) throw();

  /** Changes the size of border left around the edge of the component. */
  void setBorder(const BorderSize<int>& border) throw();

  /** Returns the size of border around the edge of the component. */
  const BorderSize<int> getBorder() const throw();

  /** Used to disable the auto-scrolling which keeps the cursor visible. If true (the default), 
  the editor will scroll when the cursor moves offscreen. If set to false, it won't. */
  void setScrollToShowCursor(const bool shouldScrollToShowCursor) throw();

  // used internally:
  virtual void paint (Graphics& g);
  virtual void paintOverChildren (Graphics& g);
  virtual void mouseDown (const MouseEvent& e);
  virtual void mouseUp (const MouseEvent& e);
  virtual void mouseDrag (const MouseEvent& e);
  virtual void mouseDoubleClick (const MouseEvent& e);
  //virtual void mouseWheelMove(const MouseEvent& e, float wheelIncrementX, float wheelIncrementY);
  virtual void mouseWheelMove(const MouseEvent &e, const MouseWheelDetails &wheel);
  bool keyPressed(const KeyPress& key);
  void focusGained(FocusChangeType cause);
  void focusLost(FocusChangeType cause);
  void resized();
  void enablementChanged();
  void colourChanged();

  //  maybe they should be all made virtual, so subclasses can override them...

protected:

  /** This adds the items to the popup menu. By default it adds the cut/copy/paste items, but you 
  can override this if you need to replace these with your own items. If you want to add your own 
  items to the existing ones, you can override this, call the base class's addPopupMenuItems() 
  method, then append your own items. When the menu has been shown, performPopupMenuAction() will 
  be called to perform the item that the user has chosen. The default menu items will be added 
  using item IDs in the range 0x7fff0000 - 0x7fff1000, so you should avoid those values for your 
  own menu IDs. If this was triggered by a mouse-click, the mouseClickEvent parameter will be a 
  pointer to the info about it, or may be null if the menu is being triggered by some other means.
  @see performPopupMenuAction, setPopupMenuEnabled, isPopupMenuEnabled */
  virtual void addPopupMenuItems(PopupMenu& menuToAddTo, const MouseEvent* mouseClickEvent);

  /** This is called to perform one of the items that was shown on the popup menu. If you've 
  overridden addPopupMenuItems(), you should also override this to perform the actions that you've 
  added. If you've overridden addPopupMenuItems() but have still left the default items on the 
  menu, remember to call the superclass's performPopupMenuAction() so that it can perform the 
  default actions if that's what the user clicked on.
  @see addPopupMenuItems, setPopupMenuEnabled, isPopupMenuEnabled */
  virtual void performPopupMenuAction (const int menuItemID);

  /** Scrolls the minimum distance needed to get the caret into view. */
  void scrollToMakeSureCursorIsVisible() throw();

  /** Used internally to dispatch a text-change message. */
  void textChanged() throw();

  /** Counts the number of characters in the text. This is quicker than getting the text as a 
  string if you just need to know the length. */
  int getTotalNumChars() throw();

  /** Begins a new transaction in the UndoManager. */
  void newTransaction() throw();

  /** Used internally to trigger an undo or redo. */
  void doUndoRedo (const bool isRedo);

  /** Can be overridden to intercept return key presses directly */
  virtual void returnPressed();

  /** Can be overridden to intercept escape key presses directly */
  virtual void escapePressed();

  // internal functions:
  void moveCaret(int newCaretPos) throw();
  void moveCursorTo(const int newPosition, const bool isSelecting) throw();
  void handleCommandMessage(int commandId);

  //===============================================================================================

private:

  // internal functions:
  void coalesceSimilarSections() throw();
  void splitSection(const int sectionIndex, const int charToSplitAt) throw();
  void clearInternal(UndoManager* const um) throw();
  void insert(const juce::String& text, const int insertIndex, BitmapFont const* font, 
    const Colour& colour, UndoManager* const um, const int caretPositionToMoveTo) throw();
  void reinsert(const int insertIndex, const juce::Array<void*>& sections) throw();
  void remove(const int startIndex, int endIndex, UndoManager* const um, 
    const int caretPositionToMoveTo) throw();
  void getCharPosition(const int index, int& x, int& y, int& lineHeight) const throw();
  int indexAtPosition(const int x, const int y) throw();
  int findWordBreakAfter(const int position) const throw();
  int findWordBreakBefore(const int position) const throw();
  void drawContent(Graphics& g);
  void updateTextHolderSize() throw();
  //float getWordWrapWidth() const throw();
  int getWordWrapWidth() const throw();
  void timerCallbackInt();
  void repaintCaret();
  void repaintText(int textStartIndex, int textEndIndex);

  // data members:
  Viewport             *viewport;
  RTextHolderComponent *textHolder;
  BorderSize<int> borderSize;
  bool readOnly                   : 1;
  bool multiline                  : 1;
  bool wordWrap                   : 1;
  bool returnKeyStartsNewLine     : 1;
  bool caretVisible               : 1;
  bool popupMenuEnabled           : 1;
  bool selectAllTextWhenFocused   : 1;
  bool scrollbarVisible           : 1;
  bool wasFocused                 : 1;
  bool caretFlashState            : 1;
  bool keepCursorOnScreen         : 1;
  bool tabKeyUsed                 : 1;
  bool menuActive                 : 1;

  UndoManager undoManager;
  int cursorX, cursorY, cursorHeight;
  int maxTextLength;
  int selectionStart, selectionEnd;
  int leftIndent, topIndent;
  unsigned int lastTransactionTime;
  BitmapFont const *currentFont;
  int totalNumChars, caretPosition;
  juce::Array<void*> sections;
  juce::String textToShowWhenEmpty;
  Colour       colourForTextWhenEmpty;

  enum
  {
    notDragging,
    draggingSelectionStart,
    draggingSelectionEnd
  } dragType;

  juce::String allowedCharacters;  // factor this into a subclass - allow all by default...
  juce::Array<RTextEditorListener*> listeners;


  // misc stuff:
  friend class RTextEditorInsertAction;
  friend class RTextEditorRemoveAction;
  friend class RTextHolderComponent;
  friend class RTextEditorViewport;

  RTextEditor (const RTextEditor&);
  const RTextEditor& operator= (const RTextEditor&);

  juce_UseDebuggingNewOperator
};

#endif  
