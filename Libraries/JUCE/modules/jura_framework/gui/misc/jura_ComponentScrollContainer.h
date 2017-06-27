#ifndef jura_ComponentScrollContainer_h
#define jura_ComponentScrollContainer_h

/** This class is component which can hold another component of possibly larger sizer and  
allows for scrolling around in that other component. */

class ComponentScrollContainer : public ColourSchemeComponent, public RScrollBarListener, 
  public ComponentListener //ComponentMovementWatcher 
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. The passed component will be considered as a child of "this" component and will
  be deleted with "this" if the optional 2nd boolean parameter is true. */
  ComponentScrollContainer(Component *contentComponentToScroll, 
    bool deleteScrolleeInDestructor = true);

  // override the constructor to take a ColourSchemeComponent....

  /** Denstructor. */
  virtual ~ComponentScrollContainer();

  //-----------------------------------------------------------------------------------------------
  // callbacks:
  virtual void scrollBarMoved(RScrollBar* scrollBarThatHasMoved, const double newRangeStart);
  virtual void paint(Graphics &g);
  virtual void paintOverChildren(Graphics &g);
  virtual void resized();
  virtual void componentMovedOrResized(Component &component, bool wasMoved, bool wasResized);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Updates the visibility and bounds of the horizontal and vertical scrollbars depending on 
  whether either of the two or both are actually needed. */
  virtual void updateScrollBarBoundsAndVisibility();

protected:

  Component  *contentComponent;
  RScrollBar *leftRightScrollBar, *upDownScrollBar;
  RWidget    *bottomRightCoverage;  // when  both scrollbars are visible, we typically want to cover the bottom right corner
  int scrollBarThickness;
  int xOffset, yOffset;  // offset for the content component with respect to the (0,0) coordinate of "this" component
  bool autoDeleteScrollee = true;

  //int xView, yView;

  juce_UseDebuggingNewOperator;
};


#endif  
