#ifndef jura_StateLoadSaveWidgetSet_h
#define jura_StateLoadSaveWidgetSet_h

/** This class is a component that puts together the widgets that are needed for loading, saving 
and displaying a currently loded state/preset file. It will send out a change-message whenever a 
preset new preset is loaded. 

\todo: let a popup menu appear, when the user clicks on the text-field which displays the current
preset - in this menu, the available presets (in the current folder) are shown and the user can 
select a new one. */

class JUCE_API StateLoadSaveWidgetSet : public WidgetSet, virtual public RButtonListener, 
  virtual public StateWatcher, virtual public ChangeBroadcaster
{

public:

  enum layouts
  {
    ONE_LINE,
    LABEL_AND_BUTTONS_ABOVE
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  StateLoadSaveWidgetSet(const juce::String& newStateLoadSaveWidgetSetName = 
    juce::String("StateLoadSaveWidgetSet"));

  /** Destructor. */
  virtual ~StateLoadSaveWidgetSet();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Updates the text displayed in the field tat shows the name of the curren t state. */
  virtual void updateStateNameField();

  /** Sets the juce::Label in which the descriptions will appear. */
  virtual void setDescriptionField(RTextField* newDescriptionField);

  /** Chooses one of the layouts for arranging the widgets around the filename filed. */
  virtual void setLayout(int newLayout);

  /** Makes the label with the preset 'headline' visible or invisible. */
  virtual void setPresetLabelVisible(bool shouldBeVisible)
  {
    stateLabel->setVisible(shouldBeVisible);
  }

  /*
  // for debug:
  virtual void setColourScheme(const WidgetColourScheme& newColourScheme) override
  {
    WidgetSet::setColourScheme(newColourScheme);
  }

  virtual void setColourSchemeFromXml(const XmlElement* widgetColours) override
  {
    WidgetSet::setColourSchemeFromXml(widgetColours);
  }
  */
  

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
  void rButtonClicked(RButton *buttonThatWasClicked) override;

  /** Implements the purely virtual method inherited from StateWatcher. */
  void stateDirtyFlagChanged(StateManager* stateManager) override;

  /** Overrides paint in order to possibly draw the headline. */
  //void paint(Graphics &g) override;

  /** Overrides the paint-method of the component base-class in order not to draw the outline
  on top of the child components. */
  //void paintOverChildren(Graphics &g) override;

  /** Overrides the resized()-method of the component base-class. */
  void resized() override;

  //-----------------------------------------------------------------------------------------------
  // public data members:

  RTextField   *stateLabel, *stateFileNameLabel;
  RClickButton *stateMinusButton, *statePlusButton, *stateLoadButton, *stateSaveButton;

protected:

  int layout;

  juce_UseDebuggingNewOperator;
};

#endif  
