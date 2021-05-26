#ifndef jura_RTextField_h
#define jura_RTextField_h

//=================================================================================================
// class RTextField:

/** This class implements a simple textfield that can be shown on the GUI. */

class JUCE_API RTextField : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RTextField(const juce::String& initialText = juce::String());

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the text to be shown. */
  virtual void setText(const juce::String& newText);

  /** Sets up the text to be shown. */
  virtual void setTextFromStdString(const std::string& newText) { setText(juce::String(newText)); }

  /** Sets the justification for positioning the text. */
  virtual void setJustification(const Justification& newJustification);

  /** Overriden from RWidget - sets the text to the passed stateString. */
  virtual void setStateFromString(const juce::String &stateString, bool sendChangeMessage = true);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the text that is shown. */
  virtual juce::String getText() const { return text; }

  /** Returns the x-coordinate at which the left border of the leftmost character is drawn. */
  virtual int getTextPixelPositionX() const;

  /** Returns the y-coordinate at which the top border of the highest characters is drawn. */
  virtual int getTextPixelPositionY() const;

  /** Retruns the width in pixels such that the text fits ideally into the textfield. */
  virtual int getWidthToFitText() const;

  /** Overriden from RWidget - returns the state (defined here as the text). */
  virtual juce::String getStateAsString() const;

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides assignParameter() inherited from RWidget in order to update the text. */
  virtual void assignParameter(Parameter* parameterToAssign) override;

  /** Overrides parameterChanged() inherited from RWidget in order to update the text. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  /** Overrides updateWidgetFromAssignedParameter() inherited from RWidget. */
  virtual void updateWidgetFromAssignedParameter(bool sendChangeMessage = false) override;

  /** Paints the text field. */
  virtual void paint(Graphics& g) override;

protected:

  /** Margin for the text from the left (or right, depending on justification) border of the 
  component (in pixels). */
  static const int horizontalMargin = 4;

  // data members:
  juce::String  text;
  Justification justification;

  //static const BitmapFont* font; // moved doen to RWidget

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RTextEntryField (and associated observer class):

class RTextEntryField;

class JUCE_API RTextEntryFieldObserver
{

public:

  virtual ~RTextEntryFieldObserver() {}

  /** This is called whenever the user has entered new text and finished his typing action via
  enter, escape or losing focus. */
  virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged) = 0;

  /** This is called whenever the user types in a new character, deletes a character or portion of
  the text or pastes a portion of text from the clipboard but did not yet finish typing by
  pressing enter, escape or losing focus. You may override this when you want to be informed about
  text-changes during typing - the default implementation does nothing. */
  virtual void somethingWasTypedIn(RTextEntryField *rTextEntryFieldThatHasChanged) { }

  // \todo: may rename these methods to avoid clashes with TextEditor if necessary
};


/** This class implements a textfield, the content of which may be edited by the user. Moreover,
the class sends out notifications when the text was changed to any registered
RTextEntryFieldObserver.

\todo: implement auto scrolling (such that caret is always visible) */

class JUCE_API RTextEntryField : public RTextField, public Timer
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RTextEntryField(const juce::String& initialText = juce::String());

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the permitted characters and optionally deletes the non-permitted ones from the current
  text. You can pass an empty string to indicate that all characters are permitted (which is the
  default). */
  virtual void setPermittedCharacters(const juce::String& newCharacters,
    bool deleteNonPermittedCharsFromText = true);

  /** Registers an observer that will get notified about events via the respective callbacks. */
  virtual void registerTextEntryFieldObserver(RTextEntryFieldObserver *observerToRegister);

  /** De-registers a previously registered observer. */
  virtual void deRegisterTextEntryFieldObserver(RTextEntryFieldObserver *observerToDeRegister);

  /** Selects all text. */
  virtual void selectAll();

  /** De-selects selected text (if any). */
  virtual void deSelect();

  /** You can used this function to tell the text-entry field that the current text is invalid.
  This will cause the widget to be displayed with a reddish background colour such that user knows
  that something is wrong. */
  virtual void markTextAsInvalid(bool shouldBeMarkedInvalid);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the pixel-position of the top-left pixel of the caret. */
  virtual Point<int> getCaretPixelPosition() const;

  /** Returns the pixel x-coordinate of the character with given index in our displayed string. */
  virtual int characterIndexToPixelPosition(int index) const;

  /** Returns the index of the character in our displayed string that is at or closest to the given
  pixel x-coordinate. */
  virtual int pixelPositionToCharacterIndex(int pixelX) const;

  /** Returns the text before the start of the selected area. */
  virtual juce::String getTextBeforeSelection() const;

  /** Returns the text inside the selected area. */
  virtual juce::String getTextInsideSelection() const;

  /** Returns the text after the end of the selected area. */
  virtual juce::String getTextAfterSelection() const;

  /** Returns the substring before the caret (or before the selection) and the substring after the
  caret (or after the selection) separately which are called head and tail here respectively. This
  is useful for inserting soemthing in between and thereby possibly replacing the selected section,
  for example. */
  virtual void getHeadAndTailString(juce::String &head, juce::String &tail) const;

  /** Returns true if the selection is empty (i.e. nothing is selected), false otherwise. */
  virtual bool isSelectionEmpty() const { return selectionEnd <= selectionStart; }

  /** Overriden from RTextField in order to optionally send out a change message. */
  virtual void setStateFromString(const juce::String &stateString, bool sendChangeMessage = true);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overriden to loose focus when the entry filed is modal and the user clicks somewhere else. */
  virtual void inputAttemptWhenModal();

  /** Overriden to edit the text (insert/delete characters, copy/paste, etc.). */
  virtual bool keyPressed(const KeyPress &key);

  /** Overriden to start the blinking of the caret. */
  virtual void focusGained(FocusChangeType cause);

  /** Overriden to stop the blinking of the caret. */
  virtual void focusLost(FocusChangeType cause);

  /** Implements the timerCallback in order to make the caret periodically visible/invisible (i.e.
  blinking). */
  virtual void timerCallback();

  /** Overriden to position the caret and de-select any selected portion of the text. */
  virtual void mouseDown(const MouseEvent &e);

  /** Overriden to select the whole text. */
  virtual void mouseDoubleClick(const MouseEvent &e);

  /** Overriden to select a substring. */
  virtual void mouseDrag(const MouseEvent &e);

  /** Paints the text-field taking into account selection, caret-visibility, etc. */
  virtual void paint(Graphics& g);

protected:

  /** Called whenever the user has finished typing to notify our observers about the text-change. */
  virtual void sendTextChangeNotificationIfChanged();

  /** Called whenever the user has typed in a new character or deletes or pastes portions of text. */
  virtual void sendTypingNotification();


  // data members:

  static const int blinkInterval = 400;  // blink interval for the caret in milliseconds

  int  caretPosition, selectionStart, selectionEnd;
  bool caretVisible, replaceMode, textInvalid;

  juce::String oldText;
  juce::String permittedCharacters;

  juce::Array<RTextEntryFieldObserver*> textEntryFieldObservers;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RLabeledTextEntryField:

/** This class implements a text entry field with an associated field for a label. */

class JUCE_API RLabeledTextEntryField : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RLabeledTextEntryField(const juce::String& labelText = juce::String(),
    const juce::String& entryFieldText = juce::String());

  /** Destructor. */
  virtual ~RLabeledTextEntryField();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the text to be shown in the label. */
  virtual void setLabelText(const juce::String &newText) { labelField->setText(newText); }

  /** Sets up the text to be shown in the text entry box. */
  virtual void setEntryFieldText(const juce::String &newText) { entryField->setText(newText); }

  /** Sets the justification for positioning the text in the label. */
  virtual void setLabelJustification(const Justification& newJustification)
  { labelField->setJustification(newJustification); }

  /** Sets the justification for positioning the text in the text entry box. */
  virtual void setEntryFieldJustification(const Justification& newJustification)
  { entryField->setJustification(newJustification); }

  /** Sets the same description for both, label and entry field. */
  virtual void setDescription(const juce::String &newDescription);

  /** Sets the width of the label. */
  virtual void setLabelWidth(int newWidth);

  /** Overriden from RWidget - sets the text in the entry field to the passed stateString. */
  virtual void setStateFromString(const juce::String &stateString, bool sendChangeMessage = true);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the text that is shown in the label. */
  virtual juce::String getLabelText() const { return labelField->getText(); }

  /** Returns the text that is shown in the text entry field. */
  virtual juce::String getEntryFieldText() const { return entryField->getText(); }

  /** Returns a pointer to the embedded RTextEntryField. This may be needed, for example, to
  register observers for the entry field. */
  virtual RTextEntryField* getTextEntryField() const { return entryField; }

  /** Overriden from RWidget - returns the state (defined here as the text inside the entry
  field). */
  virtual juce::String getStateAsString() const { return entryField->getStateAsString(); }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void resized();

protected:

  RTextField      *labelField;
  RTextEntryField *entryField;
  int labelWidth;

  juce_UseDebuggingNewOperator;
};

// \todo maybe make an RTextFieldWithModalEntryField, add comments and separator-lines in the cpp
// file

#endif
