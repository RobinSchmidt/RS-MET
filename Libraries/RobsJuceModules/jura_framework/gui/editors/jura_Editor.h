#ifndef jura_Editor_h
#define jura_Editor_h

/** This class is a component, intended to serve as base-class for all components that represent 
some kind of editor or sub-editor, for example an envelope-editor inside an editor for a 
synthesizer. */

class JUCE_API Editor : public ColourSchemeComponent
{

public:

  enum headlineStyles
  {
    NO_HEADLINE = 0,
    SUB_HEADLINE,
    MAIN_HEADLINE,

    NUM_HEADLINE_STYLES
  };

  enum headlinePositions
  {
    TOP_CENTER = 0,
    TOP_LEFT,
    TOP_RIGHT,
    BOTTOM_LEFT,
    BOTTOM_RIGHT,

    NUM_HEADLINE_POSITIONS
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  Editor(const juce::String& newEditorName = juce::String("Editor"));

  /** Destructor. */
  virtual ~Editor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Lets the outlying class change the headline. */
  virtual void setHeadlineText(const juce::String& newHeadlineText);

  /** Sets the style in which the headline is to be drawn. @see headlineStyles */
  virtual void setHeadlineStyle(int newHeadlineStyle);

  /** Sets the position where the headline is to be drawn. @see headlinePositions */
  virtual void setHeadlinePosition(int newHeadlinePosition);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the headline as String. */
  virtual const juce::String getHeadlineString();

  /** Returns the y-coordinate of the bottom of the of the headline. */
  virtual int getHeadlineBottom();

  /** Returns the x-coordinate of the right edge of the of the headline. */
  virtual int getHeadlineRight();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  //virtual void rButtonClicked(RButtonBase *buttonThatWasClicked) {}
  //virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged) {}
  //virtual void rLabelTextChanged(RLabel* rLabelThatHasChanged) {}
  //virtual void rSliderValueChanged(RSlider *sliderThatHasChanged) {}
  void paint(Graphics &g) override;
  void resized() override;

protected:

  /** Adds a child-editor to this one (and optionally also adds it as child-component to this one). 
  Using this method rather than Component's addChildComponent allows for looping through all the 
  sub-editors widgets to set up their colours, description-fields
  etc. */
  virtual void addChildEditor(Editor* editorToAdd, bool addAsChildComponent = true, 
    bool makeVisible = true);

  /** Removes a child-editor from this one and optionally deletes the object */
  virtual void removeChildEditor(Editor *editorToRemove, bool deleteObject);

  /** Draws the headline for the editor according to the chosen headline-style. */
  virtual void drawHeadline(Graphics& g);

  /** A flag to indicate whether this editor should be drawn with an enclosing rectangle or not. */
  bool drawWithEnclosingRectangle;

  /** An array of all the child-editors that have been added by addChildEditor(). */
  juce::Array<Editor*, CriticalSection> childEditors;

  juce::String headlineText;
  int          headlineStyle, headlinePosition, headlineX, headlineY;

  juce::Array<juce::Rectangle<int> > guiLayoutRectangles;

  // these are not always needed, therefore nullptr by default:
  RButton *closeButton;               

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/**

This class extends the editor class by adding facilities to save/restore a state to/from a file

*/

class EditorWithStateFile : virtual public Editor, public StateFileManager, public ChangeListener
{

public:


  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  EditorWithStateFile(const juce::String& name = juce::String("EditorWithStateFile"));  

  /** Destructor. */
  virtual ~EditorWithStateFile();  

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void resized();
  virtual int  changeListenerCallback(void* objectThatHasChanged);

  /** Override this in your subclass to update your widgets. */
  //virtual void stateChanged() = 0;

  /** Updates the state-widget set. Override this in your subclass to update your custom widgets 
  there as well (but always call the baseclass method). */
  //virtual void updateWidgetsAccordingToState();

protected:

  StateLoadSaveWidgetSet* stateWidgetSet;

  juce_UseDebuggingNewOperator;
};



#endif  
