#ifndef rojue_StateLoadSaveWidgetSet_h
#define rojue_StateLoadSaveWidgetSet_h

#include "../../file_management/rojue_StateFileManager.h"
#include "rojue_RButton.h"
#include "rojue_RTextField.h"
#include "../editors/rojue_ColourSchemeComponent.h"

namespace rojue
{

  /**

  This class is a component that puts together the widgets that are needed for loading, saving and displaying a currently loded 
  state/preset file. It will send out a change-message whenever a preset new preset is loaded.

  */

  class StateLoadSaveWidgetSet : virtual public ColourSchemeComponent, virtual public RButtonListener, 
    virtual public StateWatcher, virtual public DescribedItem, virtual public ChangeBroadcaster
  {

  public:

    enum layouts
    {
      ONE_LINE,
      LABEL_AND_BUTTONS_ABOVE
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    StateLoadSaveWidgetSet(const juce::String& newStateLoadSaveWidgetSetName = juce::String(T("StateLoadSaveWidgetSet")));  

    /** Destructor. */
    virtual ~StateLoadSaveWidgetSet();  

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Updates the text displayed in the field tat shows the name of the curren t state. */
    virtual void updateStateNameField();

    /** Sets the juce::Label in which the descriptions will appear. */
    virtual void setDescriptionField(RTextField* newDescriptionField);

    /** Chooses one of the layouts for arranging the widgets around the filename filed. */
    virtual void setLayout(int newLayout);

    /** Makes the label with the preset 'headline' visible or invisible. */
    virtual void setPresetLabelVisible(bool shouldBeVisible)
    { stateLabel->setVisible(shouldBeVisible); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // callbacks:

    /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
    virtual void rButtonClicked(RButton *buttonThatWasClicked);

    /** Implements the purely virtual method inherited from StateWatcher. */
    virtual void stateDirtyFlagChanged(StateManager* stateManager);

    /** Overrides piant in order to possibly draw the headline. */
    virtual void paint(Graphics &g);

    /** Overrides the paint-method of the component base-class in order not to draw the outline 
    on top of the child components. */
    virtual void paintOverChildren(Graphics &g);

    /** Overrides the resized()-method of the component base-class. */
    virtual void resized();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // state management:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // public data members:

    RTextField   *stateLabel, *stateFileNameLabel;
    RClickButton *stateMinusButton, *statePlusButton, *stateLoadButton, *stateSaveButton;

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    int layout;

  };

}

#endif  
