#ifndef jura_RButton_h
#define jura_RButton_h  // maybe rename to RButtons

class RButton;

//=================================================================================================

/** A class to receive callbacks when a button is clicked. */

class JUCE_API RButtonListener
{

public:

  virtual ~RButtonListener() {}

  /** Called when the button is clicked. */
  virtual void rButtonClicked(RButton* button) = 0;

};

//=================================================================================================

/** Baseclass for custom painting of buttons. In a concrete subclass, you override the
paintButton method to actually paint the button. To assign a custom painter to a button, you would
call myButton->setPainter(myButtonPainter) where myButtonPainter is asumed to be a pointer to an 
object of your RButtonPainter subclass. */

class JUCE_API RButtonPainter
{
public:
  RButtonPainter() {}
  virtual ~RButtonPainter() {}

  /** You need to override this method to actually paint the button using the passed Graphics 
   object. */
  virtual void paint(Graphics& g, RButton* button) = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RButtonPainter)
};

//=================================================================================================

/** This is a class for representing buttons an a GUI. */

class JUCE_API RButton : public RWidget
{

public:

  enum buttonSymbols
  {
    NO_SYMBOL = 0,
    PLUS,
    MINUS,
    ARROW_UP,    // rename to TRIANGLE_UP etc.
    ARROW_DOWN,
    ARROW_LEFT,
    ARROW_RIGHT,
    PLAY,
    SKIP_FORWARD,
    SKIP_BACK,
    MUTE,
    LOOP,
    CLOSE,

    NUM_SYMBOLS
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructs a button with a symbol. */
  RButton(int newSymbolIndex);

  /** Constructs a button with text. */
  RButton(const juce::String& buttonText = juce::String("RButton"));

  /** Destructor. */
  virtual ~RButton();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Chooses a new symbol for this button. */
  virtual void setSymbolIndex(int newSymbolIndex);

  /** Changes the button's text. */
  void setButtonText(const juce::String& newText) throw(); // why throw?

  /** Decides, whether or not this button should change its on/off state when clicked. */
  void setClickingTogglesState(const bool shouldToggle);

  /** A button has an on/off state associated with it, and this changes that. By default buttons
  are 'off' and for simple buttons that you click to perform an action you won't change this. Toggle
  buttons, however will want to change their state when turned on or off. */
  void setToggleState(const bool shouldBeOn, const bool sendChangeNotification);

  /** Overriden from RWidget - sets the toggle-state to "off" when the string-as-integer is 0 and
  sets it to "on" when the
  string-as-integer is != 0 (presumably "1"). */
  virtual void setStateFromString(const juce::String &stateString, bool sendChangeMessage = true);

  /** Sets up a painter object that can be used for custom painting. */
  virtual void setPainter(RButtonPainter* painterToUse) { painter = painterToUse; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the text displayed in the button. */
  const juce::String getButtonText() const throw() { return text; }

  /** Returns true if the button in 'on'. */
  bool getToggleState() const throw() { return isOn; }

  /** Returns "1" when the toggle-state is "on", otherwise "0". */
  virtual juce::String getStateAsString() const;

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Registers a listener to receive events when this button's state changes. */
  void addRButtonListener(RButtonListener* const newListener) throw();

  /** Removes a previously-registered button listener */
  void removeRButtonListener(RButtonListener* const listener) throw();

  /** Overrides the method inherited from RWidget */
  virtual void updateWidgetFromAssignedParameter(bool sendChangeMessage = false);

  /** Paints the button. */
  virtual void paint(Graphics &g);

  /** @internal */
  virtual void mouseDown(const MouseEvent& e);

  /** @internal */
  virtual void enablementChanged();

protected:

  /** Sends out a message to our listeners that this button has been clicked. */
  void sendClickMessage();

  /** Overrides the inherited clicked callback in order to update an assigned Parameter
  (if any). */
  virtual void clicked();

  /** Draws the symbol onto the button (if any). */
  virtual void drawSymbol(Graphics &g) const;


  juce::String text;
  juce::Array<RButtonListener*> buttonListeners;
  RButtonPainter *painter = nullptr;

  bool isOn;
  bool clickTogglesState;
  int  symbolIndex;

private:

  RButton(const RButton&);
  const RButton& operator=(const RButton&);

  juce_UseDebuggingNewOperator;
};


//=================================================================================================
// class RClickButton:

/** This class implements a button that paints as active when mouse is down, and switches back into
inactive state when the mouse is up again. maybe rename to RBangButton */

class JUCE_API RClickButton : public RButton
{

public:

  RClickButton(int newSymbolIndex);

  RClickButton(const juce::String& buttonText = juce::String("RButton"));

  virtual void mouseDown(const MouseEvent& e);

  virtual void mouseUp(const MouseEvent& e);


  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RClickButtonNotifyOnMouseUp:

/** This class is like an RClickButton but it sends out a click message only on mouse-up events.
Moreover, these mouse-up events must occur inside this button after a mouse-down event occured.
This behaviour is desirable for OK/Cancel buttons on dialog boxes, for example. */

class JUCE_API RClickButtonNotifyOnMouseUp : public RClickButton
{

public:

  RClickButtonNotifyOnMouseUp(int newSymbolIndex);

  RClickButtonNotifyOnMouseUp(const juce::String& buttonText = juce::String("RButton"));

  virtual void mouseDown(const MouseEvent& e);

  virtual void mouseUp(const MouseEvent& e);

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RClickButtonWithAutoRepeat

/** This class is like an RClickButton but it sends out click messages repeatedly as long as it is
clicked. */

class JUCE_API RClickButtonWithAutoRepeat : public RClickButton, public Timer
{

public:

  RClickButtonWithAutoRepeat(int newSymbolIndex);

  RClickButtonWithAutoRepeat(const juce::String& buttonText = juce::String("RButton"));

  virtual void mouseDown(const MouseEvent& e);

  virtual void mouseUp(const MouseEvent& e);

  virtual void timerCallback();

protected:

  int initialDelay, timeInterval; // both values in ms \todo maybe provide setters and getters for these

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RRadioButton and RRadioButtonGroup:

class RRadioButtonGroup;

/** This class implements a button that can be used in a group of mutually exclusively pressed
buttons. That means, only one at a time canbe in 'pressed' state (this kind of behavior is also
known as "radio-button"). */

class JUCE_API RRadioButton : public RButton
{

  //friend class RRadioButtonGroup;

public:

  RRadioButton(int newSymbolIndex);

  RRadioButton(const juce::String& buttonText = juce::String("RButton"));

  /** Sets the radio-group of which this button should become a member. */
  virtual void addToRadioButtonGroup(RRadioButtonGroup *newGroupToUse);
  // maybe rename to setRadioButtonGroup

  /** Overriden clicked to make sure that all other buttons in the same radio-group are going to be
  siwtched off. */
  virtual void clicked();

  /** Overriden to make sure that all other buttons in the same radio-group are going to be
  switched off. */
  void setToggleState(const bool shouldBeOn, const bool sendChangeNotification);

protected:

  RRadioButtonGroup *radioGroupToUse;

  juce_UseDebuggingNewOperator;
};

class JUCE_API RRadioButtonGroup
{

public:

  virtual ~RRadioButtonGroup() {}

  virtual void addButtonToRadioGroup(RRadioButton *buttonToAdd);

  virtual void removeButtonFromRadioGroup(RRadioButton *buttonToRemove);

  /** Toggles the passed button (which is assumed to be a member of this group) on and toggles all
  other buttons in the group off. */
  virtual void toggleRadioButtonOn(RRadioButton *buttonToToggleOn, bool sendNotifications);

  /** Returns true if the passed button pointer points to a member of the radio-group. */
  virtual bool isButtonMemberOfGroup(RButton *buttonToCheck);

protected:

  juce::Array<RRadioButton*> radioButtons;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RHyperlinkButton

class JUCE_API RHyperlinkButton : public RButton
{

public:

  /** Creates a RHyperlinkButton with given text and refering to the given URL. */
  RHyperlinkButton(const juce::String& linkText, const URL& linkURL);

  /** Destructor. */
  ~RHyperlinkButton();

  /** Changes the URL that the button will trigger. */
  void setURL(const URL& newURL) throw();

  /** Returns the URL that the button will trigger. */
  const URL& getURL() const throw() { return url; }

protected:

  void clicked();

  void paint(Graphics& g);

  URL url;

private:

  RHyperlinkButton(const RHyperlinkButton&);
  const RHyperlinkButton& operator= (const RHyperlinkButton&);

  juce_UseDebuggingNewOperator
};

//=================================================================================================

/** A button painter class that gives the button a pseudo 3D'ish look.
todo: test it... */

class JUCE_API RButtonPainter3D : public RButtonPainter
{
public:
  RButtonPainter3D(){}
  virtual void paint(Graphics& g, RButton *button) override;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RButtonPainter3D)
};

// todo: make a class RButtonImageSet with members for up to 4 images: up, down, upMouseOver,
// downMouseOver

#endif
