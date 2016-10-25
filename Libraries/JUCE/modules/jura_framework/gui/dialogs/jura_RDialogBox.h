#ifndef jura_RDialogBox_h
#define jura_RDialogBox_h

//#include "../editors/rojue_Editor.h"
//#include "../widgets/rojue_Widgets.h"

class RDialogBox;

/** A class for receiving events from a RDialogBox.  */

class JUCE_API RDialogBoxListener
{

public:

  /** Destructor. */
  virtual ~RDialogBoxListener() {}

  /** Called when the user has changed something relevant in a RDialogBox. This can be useful to 
  (temporarily) apply the settings from the box immediately. Here, the 'temporarily' means that you
  should have some means to revert to the old settings in case the box is left via the 
  cancel-button. */
  virtual void rDialogBoxChanged(RDialogBox* dialogBoxThatHasChanged) = 0;

  /** Called when the user has clicked on the OK-button. The caller should then take care of 
  applying the settings that have been made in the box (for example by retrieving them via getters 
  from some custom RDialogBox subclass) and the deleting the box (or making it invisible). */
  virtual void rDialogBoxOKClicked(RDialogBox* dialogBoxThatWantsToAcceptAndLeave) = 0;

  /** Called when the user has clicked on the cancel-button. The caller should then take care of 
  deleting the box (or making it invisible) without applying the settings that have been made in 
  the box. */
  virtual void rDialogBoxCancelClicked(RDialogBox* dialogBoxThatWantsToBeCanceled) = 0;

};

//=================================================================================================

/** This class serves as baseclass for various dialogs. */

class JUCE_API RDialogBox : virtual public Editor, public RButtonListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RDialogBox();

  /** Destructor. */
  virtual ~RDialogBox();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Registers a listener that will be called when the box's content changes. */
  virtual void addListener(RDialogBoxListener* const listenerToAdd);

  /** Deregisters a previously-registered listener. */
  virtual void removeListener(RDialogBoxListener* const listenerToRemove);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void resized();

protected:

  /** Calls 'dialogBoxChanged' for all our registered listeners with this box as argument. */
  virtual void sendChangeNotification();

  /** Calls 'dialogBoxOKClicked' for all our registered listeners with this box as argument. */
  virtual void sendOKClickedNotification();

  /** Calls 'dialogBoxCancelClicked' for all our registered listeners with this box as argument. */
  virtual void sendCancelClickedNotification();

  juce::Array<RDialogBoxListener*, CriticalSection> listeners;
  RClickButtonNotifyOnMouseUp *cancelButton, *okButton;


  juce_UseDebuggingNewOperator;
};

#endif
