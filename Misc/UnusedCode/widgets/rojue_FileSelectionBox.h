#ifndef rojue_FileSelectionBox_h
#define rojue_FileSelectionBox_h

#include "../../file_management/rojue_FileManager.h"
#include "rojue_WidgetSet.h"
#include "rojue_RButton.h"
#include "rojue_RTextField.h"

namespace rojue
{

  /** 
  
  A box for selecting a file with associated widgets. Works together with a FileManager object. Derived from ChangeBroadcaster, it sends 
  out a change-message whenever the active file changes. 
  
  */

  class FileSelectionBox : public WidgetSet, public RButtonListener, public FileManagerListener, public ChangeBroadcaster
  {

  public:

    /** Positions where the label for the box can appear. */
    enum labelPositions
    {
      LABEL_INVISIBLE,
      LABEL_LEFT,
      LABEL_ABOVE
    };

    /** Positions where the load/save/+/- buttons can appear. */
    enum buttonPositions
    {
      BUTTONS_ABOVE,
      BUTTONS_RIGHT,
      BUTTONS_BELOW
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - you must pass a valid (non-NULL) FileManager object. */
    FileSelectionBox(const juce::String& componentName, FileManager *fileManagerToUse);

    /** Destructor. */
    virtual ~FileSelectionBox();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets the position where the label appears relative to the actual box. @see: labelPositions */
    virtual void setLabelPosition(int newPosition) { labelPosition = newPosition; resized(); }

    /** Sets the position where the load/save/+/- buttons appear relative to the actual box. @see: buttonPositions */
    virtual void setButtonsPosition(int newPosition) { buttonsPosition = newPosition; resized(); }

    /** Sets the height of the box that shows the filename (typically, the box should be either 16 or 20 pixels high, default is 16). */
    virtual void setBoxHeight(int newHeight) { boxHeight = newHeight; resized(); }

    /** Selects whether or not the save-button should be visible. */
    virtual void setSaveButtonVisible(bool shouldBeVisible) { saveButton->setVisible(shouldBeVisible); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // callbacks:

    virtual void rButtonClicked(RButton *buttonThatWasClicked);
    virtual void activeFileChanged(FileManager *fileManagerThatHasChanged);
    virtual void resized();

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator

  protected:

    /** Updates the content shown in the box according to the current file in the associated FileManager. */
    virtual void updateBoxContent();

    // data members:
    FileManager *fileManager;
    RTextField  *fileLabel, *fileNameBox;
    RButton     *loadButton, *saveButton, *plusButton, *minusButton;
    int         labelPosition, buttonsPosition, boxHeight;

  };

}

#endif  
