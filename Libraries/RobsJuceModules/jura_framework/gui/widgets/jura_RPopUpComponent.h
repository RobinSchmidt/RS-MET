#ifndef jura_RPopUpComponent_h
#define jura_RPopUpComponent_h


/** This class implements a component that can appear like pop up on the desktop like, for example,
for a popup menu. The component is to be used by passing a content-component
(via setContentComponent) which will contain the actual stuff. After the sontent component has been
passed, the popup can be made to appear modally or non-modally via the respective show...
methods. */

class JUCE_API RPopUpComponent : public RWidget  //, public ComponentMovementWatcher
{

  friend class RComboBox; // preliminray

public:

  enum attachPositions
  {
    BELOW = 0,
    ABOVE,
    RIGHT,
    LEFT
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates an empty popup component. */
  RPopUpComponent();

  /** Destructor. */
  virtual ~RPopUpComponent();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** This function is to be used to pass the actual content-component of the popup component. The
  RPopUpComponent class will take over responsibility for eventually deleting the content
  component. */
  virtual void setContentComponent(Component *newContentComponent);

  /** Alternative to setContentComponent(Component*) for RWidgets to do handle colorschemens
  etc.. */
  virtual void setContentWidget(RWidget *newContentComponent);

  /** Decides whether or not the menu should be dismissed when it looses kayboard focus. */
  virtual void setDismissOnFocusLoss(bool shouldDismiss) { dismissOnFocusLoss = shouldDismiss; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns true when the popup is currently open, false otherwise. */
  virtual bool isOpen() const { return open; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void inputAttemptWhenModal();
  virtual void focusLost(FocusChangeType cause);
  virtual bool canModalEventBeSentToComponent(const Component* targetComponent);

  /*
  virtual void mouseEnter(const MouseEvent &e);  // overriden to make sure that highlighting appears, if appropriate
  virtual void mouseExit( const MouseEvent &e);  // overriden to make sure that highlighting disappears
  */
  virtual void mouseDown(const MouseEvent &e);  // overriden to update the selected item
  /*
  virtual void mouseMove( const MouseEvent &e);  // overriden to update the highlighting
  virtual void paint(Graphics&);
  */

  //-----------------------------------------------------------------------------------------------
  // others:

  // experimental:
  virtual void showAtNonModal(int screenX, int screenY, int width = 0, int height = 0);

  /** Shows the menu modally at a specified location on the screen. The passed coordinates are used
  to position the top-left pixel of the menu. You may also pass a width for the menu - if you leave
  that at the default value (zero), the menu will work out the required width itself (by using
  getRequiredWidth()). The same goes likewise for the height. When the function has returned, the
  user can retrieve the selected item via getSelectedItem(). */
  virtual void showAt(bool showModally, int screenX, int screenY, int width = 0, int height = 0);

  /** Similar to showAt, but uses the current position of the mouse for the top-left corner. */
  virtual void showAtMousePosition(bool showModally, int width = 0, int height = 0);

  /** Removes the component from the desktop and exits the modal loop if necessary. */
  virtual void dismiss();

protected:

  Component *contentComponent;
  bool      open;
  bool      dismissOnFocusLoss;

  juce_UseDebuggingNewOperator;
};


//=================================================================================================
// class RPopUpOwner

class ROwnedPopUpComponent; // forward decalaration

/**
This class serves as basclass for all Component subclasses that need to keep informed about the
question whether a popup was dismissed due to a click of the owner-component. The owner component
may be interested in that because it will receive a mouseDown callback right after the popup
disappeared - normally it would then open the popup, but if it knows that this callback is due to
the click while the popup was still open, it should perhaps not bring it up immediately again.
RPopUp, will call popUpDismissedByClickOnOwner whenever such a situation occurs. The owner may then
implement a way to avoid to respond to the very next mousclick on it. A bit clunky, admittedly.
*/

class JUCE_API RPopUpOwner
{

public:

  virtual ~RPopUpOwner() {}

  /** Callback that gets called whenever the owned popup was dismissed due to a click on its owner.
  Subclasses can override this to treat this case differently from the cases where the popup is
  dismissed due to a click somewhere else. Distinguishing these cases may be necessary because in
  the first case, the owner will receive a mousDown callback immediately after the dismissal which
  probably should be ignored instead of bringing up the popup (again). */
  virtual void rPopUpDismissedByClickOnOwner(ROwnedPopUpComponent *popUp) {}

  juce_UseDebuggingNewOperator
};


//=================================================================================================
// class ROwnedPopUpComponent:

/**

A popup component that can be attached to another component in order to follow its movements.

*/

class JUCE_API ROwnedPopUpComponent : public RPopUpComponent, public ComponentMovementWatcher
{

public:


  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates an empty popup component. */
  ROwnedPopUpComponent(Component *ownerComponent);

  //-----------------------------------------------------------------------------------------------
  // others:


  virtual void focusLost(FocusChangeType cause);

  virtual void componentMovedOrResized(bool wasMoved, bool wasResized);
  virtual void componentPeerChanged();
  virtual void componentVisibilityChanged();

  virtual bool canModalEventBeSentToComponent(const Component* targetComponent);

  /** Shows the popup attached to another component as - for example - in comboboxes. */
  virtual void show(bool showModally, int attachPosition, int width = 0, int height = 0,
    int xOffset = 0, int yOffset = 0);


  juce_UseDebuggingNewOperator;

protected:

  virtual bool isMouseDownOnOwner();

  Component *ownerComponent;

};

#endif
