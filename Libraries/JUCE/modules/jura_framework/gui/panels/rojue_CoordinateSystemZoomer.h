#ifndef rojue_CoordinateSystemZoomer_h
#define rojue_CoordinateSystemZoomer_h

#include "rojue_CoordinateSystem.h"
#include "../widgets/rojue_RButton.h"
#include "../widgets/rojue_RScrollBar.h"
#include "../editors/rojue_ColourSchemeComponent.h"

namespace rojue
{

  /**

  This a class for a CoordinateSystem with zooming and scrolling capabilities. It has a pointer 
  to the actual CoordinateSystem (or some subclass thereof) object for which the zooming and 
  scrolling is required. This pointer needs to be set up via the setCoordinateSystem()-method - 
  otherwise it will be NULL and emptyness will be shown.
  
  \todo: y-zoom on shift/or ctrl-wheel

  */

  class CoordinateSystemZoomer : virtual public ColourSchemeComponent, public RButtonListener, public RScrollBarListener

  {

  public:

    enum verticalMouseWheelModes
    {
      noResponseToVerticalMouseWheel = 0,
      forwardToCoordinateSystem,
      horizontalScrollViaVerticalMouseWheel,
      verticalScrollViaVerticalMouseWheel,
      horizontalZoomViaVerticalMouseWheel,
      verticalZoomViaVerticalMouseWheel,
      zoomViaVerticalMouseWheel
    };

    CoordinateSystemZoomer();  
    ///< Constructor.

    virtual ~CoordinateSystemZoomer();  
    ///< Destructor.

    //-----------------------------------------------------------------------------------------------
    // setup:

    virtual void setCoordinateSystem(CoordinateSystem* newSystemToBeShown);
    /**< This function should be used to pass the CoordinateSystem-object (or an 
    object of some subclass thereof) which should be shown, zoomed and 
    scrolled. */

    virtual void invalidateCoordinateSystemPointer();
    /**< This function should be called when the CoordinateSystem, to which we
    hold a pointer here, was deleted for some reason. */

    virtual void setWidgetDescriptionField(RTextField* newDescriptionField);
    /**< Sets the juce::Label in which the descriptions for the widgets will appear. */

    virtual void setVerticalMouseWheelMode(int newMode);
    /**< Chooses, in which way the zoomer should respond to the vertical 
    mousewheel. For the available modes, see verticalMouseWheelModes. */

    virtual void setRelativeMargins(double newRelativeMarginLeft, double newRelativeMarginRight,
      double newRelativeMarginTop, double newRelativeMarginBottom);
    /**< Sets up the relative margins which should be used in the maximally zoomed out state. A 
    coordinate-system has typically a maximum meaningful range to be displayed, for example, a 
    unit circle can be displayed with a coordinate-system range of x: -1.5...+1.5, y: -1.5...+1.5 
    and it would not be very meaningful to plot it in a system with much larger range.  
    But we might want to have some margins around this range. */

    //-----------------------------------------------------------------------------------------------
    // appearance setup:

    virtual void setZoomerSize(int newZoomerSize) { zoomerSize = newZoomerSize; };

    virtual int getZoomerSize() { return zoomerSize; };
    /**< Returns the size of the zoom- and scroll-widgets. This is the edge-length
    of the buttons, and the thickness of the scrollbars. */

    virtual void hideScrollBarX(bool shoulBeHidden);
    /**< Can be used to hide the ScrollBar for the x-axis (or make it visible again). */

    virtual void hideScrollBarY(bool shoulBeHidden);
    /**< Can be used to hide the ScrollBar for the x-axis (or make it visible again). */

    //-----------------------------------------------------------------------------------------------
    // callbacks:

    virtual void rButtonClicked(RButton *button);
    /**< Implements the rButtonClicked()-method of the ButtonListener base-class. */

    virtual void scrollBarMoved(RScrollBar *scrollBarThatHasMoved, const double newRangeStart);
    /**< Implements the scrollBarMoved()-method of the ScrollbarListener base-class. */

    virtual void mouseWheelMove(const MouseEvent& e,  
      float  wheelIncrementX, float  wheelIncrementY);
    /**< Overrides the mouseWheelMove()-method of the indirect MouseListener base-class to allow 
    zooming by mouswheel. */

    // additionally, we have to override all the other mouse-callbacks in order 
    // to forward them to the CoordinateSystem:
    virtual void mouseMove       (const MouseEvent& e);
    virtual void mouseEnter      (const MouseEvent& e);
    virtual void mouseExit       (const MouseEvent& e);
    virtual void mouseDown       (const MouseEvent& e);
    virtual void mouseDrag       (const MouseEvent& e);
    virtual void mouseUp         (const MouseEvent& e);
    virtual void mouseDoubleClick(const MouseEvent& e);

    /** Overrides paint (inherited from ColourSchemeComponent) in order to suppress drawing of the 
    gradient background. */
    virtual void paint(Graphics &g);

    //-----------------------------------------------------------------------------------------------
    // CoordinateSystem manipulation:

    virtual void shiftX(double shiftAmount);
    /**< Causes the CoordinateSystem to be shifted left or right and updates the scrollBarX. */

    virtual void shiftY(double shiftAmount);
    /**< Causes the CoordinateSystem to be shifted up or down and updates the scrollBarY. */

    virtual void zoomX(double zoomFactor, double relativeCenter = 0.5);
    /**< Zooms the x-coordinate in or out (according to whether the zoomFactor is larger or smaller 
    than 1) and updates the scrollBarX. The relative center determines around which (relative) 
    reference x-coordinate will be zoomed - this is useful for taking the mouse-position into 
    account when zooming via mouse-wheel ("zoom to mouse-coursor" instead of "zoom to center*). 
    */

    virtual void zoomY(double zoomFactor, double relativeCenter = 0.5);
    /**< Zooms the y-coordinate in or out @see: zoomX(). */

    virtual void zoomToAllX();
    /**< Zooms x-axis out maximally.  */

    virtual void zoomToAllY();
    /**< Zooms y-axis out maximally.  */

    virtual void zoomToAllXY();
    /**< Zooms out maximally to show all.  */

    //virtual void setCurrentRangeAndUpdateScrollBarsX(double newMinX, double newMaxX);
    /**< Sets the current range and update the scrollbars, avoiding a double redraw. May be 
    deprecated when the CoordinateSystem only redraws itself on actual range-changes */

    //virtual void setCurrentRangeAndUpdateScrollBarsY(double newMinY, double newMaxY);
    /**< @see setCurrentRangeAndUpdateScrollBarsX */


    //-----------------------------------------------------------------------------------------------
    // mostly for internal use:

    virtual void alignWidgetsToCoordinateSystem();
    /**< Causes the widgets (the zoom-buttons and scrollbars) to align themselves according to the 
    position and size of the CoordinateSystem object, which is being manipulated */

    virtual void updateScrollbars();
    /**< Updates the positions and sizes of the scrollbar thumbs. */

    //===============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    virtual double transformToScrollBarCoordinateX(double x);
    /**< Transforms from an x-value given in the CoordinateSystem's coordinates to normalized 
    coordinates in the range 0...1 to be used by the horizontal scrollbar. */

    virtual double transformFromScrollBarCoordinateX(double x);
    /**< Transforms from a normlized x-value in the range 0...1 to the CoordinateSystem's 
    coordinates. */

    virtual double transformToScrollBarCoordinateY(double y);
    /**< Transforms from an y-value given in the CoordinateSystem's coordinates to normalized 
    coordinates in the range 0...1 to be used by the vertical scrollbar. */

    virtual double transformFromScrollBarCoordinateY(double y);
    /**< Transforms from a normlized y-value in the range 0...1 to the CoordinateSystem's 
    coordinates. */

    double zoomFactorPerClickX, zoomFactorPerClickY;
    /**< The factor by which is zoomed in/out on one click on a plus/minus button for x- and 
    y-axis seperately. */

    double shiftPerClickX, shiftPerClickY;
    /**< The amount by which the CoordinateSystem is shifted left/right or up/down on one click on 
    an arrow-button on the scrollbar. */

    int verticalMouseWheelMode;
    /**< The mode of the(vertical) mousewheel - can scroll or zoom horizontally, vertically, both or 
    none. */

    double relativeMarginLeft, relativeMarginRight, relativeMarginTop, relativeMarginBottom;
    /**< Some margins for the mximum zoomed out level (in percent). */

    bool supressZoomability;
    /**< a flag which allows to suppress all the zooming and scrolling widgets, leaving more space 
    for the actual coordinate system. */

    int zoomerSize;
    /**< thickness of the scrollbars and edge-length of the buttons. */

    bool scrollBarIsHiddenX, scrollBarIsHiddenY;
    /**< Flags to indicate that the scrollbars should be hidden. */

    // the colour-scheme for the scrollbars:
    Colour scrollBarBackgroundColour;
    Colour scrollBarOutlineColour;
    Colour scrollBarThumbColour;

    // embedded components:
    CoordinateSystem* theCoordinateSystem;
    RButton*          zoomInButtonX;
    RButton*          zoomOutButtonX;
    RButton*          zoomToAllButtonX;
    RScrollBar*       scrollBarX;
    RButton*          zoomInButtonY;
    RButton*          zoomOutButtonY;
    RButton*          zoomToAllButtonY;
    RScrollBar*       scrollBarY;
    RButton*          zoomToAllButtonXY;

  };

}

#endif 