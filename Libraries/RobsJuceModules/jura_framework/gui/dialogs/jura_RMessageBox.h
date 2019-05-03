#ifndef jura_RMessageBox_h
#define jura_RMessageBox_h

/** This class is a simple message box to display messages to the user - it can make itself 
invisible by clicking on an OK button, so it can be used and re-used as a member of some class 
that needs to display messages. */

class RMessageBox : public ColourSchemeComponent, public RButtonListener
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructs a button with a symbol. */
  RMessageBox();

  /** Destructor. */
  virtual ~RMessageBox();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets the text for the headline of this box. */
  virtual void setHeadlineText(const juce::String& newText);

  /** Sets the text for the body of this box. */
  virtual void setBodyText(const juce::String& newText);

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides rButtonClicked() in order to make this box invisible in response to a click on
  the OK button. */
  virtual void rButtonClicked(RButton* button);

  /** Overrides resized() */
  virtual void resized();

  //=============================================================================================
  juce_UseDebuggingNewOperator;

protected:

  RButton     *okButton;
  RTextEditor *bodyTextField;
  RTextField  *headlineTextField;

};

#endif  
