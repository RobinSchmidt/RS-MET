#ifndef jura_RScrollBar_h
#define jura_RScrollBar_h

class RScrollBar;

/** A class for receiving events from an RScrollBar. You can register an RScrollBarListener with
an RScrollBar using the RScrollBar::addListener() method, and it will be called when the bar's
position changes. */

class JUCE_API RScrollBarListener
{

public:

  /** Destructor. */
  virtual ~RScrollBarListener() {}

  /** Called when an RScrollBar is moved. */
  virtual void scrollBarMoved(RScrollBar* scrollBarThatHasMoved, const double newRangeStart) = 0;

};

//=================================================================================================

/** A scrollbar component. To use a scrollbar, set up its total range using the setRangeLimits() 
method - this sets the range of values it can represent. Then you can use setCurrentRange() to 
change the position and size of the scrollbar's 'thumb'. Registering an RScrollBarListener with the 
scrollbar will allow you to find out when the user moves it, and you can use the 
getCurrentRangeStart() to find out where they moved it to.

\todo: implement auto-repeat for the scroll-in-pages functionality (is this what the inherited 
Timer was for?) 
-we perhaps need an AutoRepeatingButton class 

*/

class JUCE_API RScrollBar : public RWidget, public RButtonListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates an RScrollbar. */
  RScrollBar(const bool isVertical);

  /** Destructor. */
  virtual ~RScrollBar();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the minimum and maximum values that the bar will move between. The bar's thumb will 
  always be constrained so that the top of the thumb will be >= minimum, and the bottom of the 
  thumb <= maximum. @see setCurrentRange */
  virtual void setRangeLimits(const double minimum, const double maximum) throw();

  /** Changes the position of the scrollbar's 'thumb'. This sets both the position and size of the 
  thumb - to just set the position without changing the size, you can use setCurrentRangeStart(). 
  If this method call actually changes the scrollbar's position, it will trigger a call to 
  RScrollBarListener::scrollBarMoved() for all the listeners that are registered. */
  virtual void setCurrentRange(double newStart, double newSize) throw();

  /** Moves the bar's thumb position. This will move the thumb position without changing the thumb 
  size. Note that the maximum thumb start position is 
  getMaximumRangeLimit() - getCurrentRangeSize(). If this method call actually changes the 
  scrollbar's position, it will trigger a call to RScrollBarListener::scrollBarMoved() for all the 
  listeners that are registered. */
  virtual void setCurrentRangeStart(double newStart) { setCurrentRange(newStart, rangeSize); }

  /** Sets the amount by which the up and down buttons will move the bar. The value here is in 
  terms of the total range, and is added or subtracted from the thumb position when the user clicks 
  an up/down (or left/right) button. */
  virtual void setSingleStepSize(const double newSingleStepSize) 
  { 
    singleStepSize = newSingleStepSize; 
  }

  /** Moves the scrollbar by a number of single-steps. This will move the bar by a multiple of its 
  single-step interval (as specified using the setSingleStepSize() method). A positive value here 
  will move the bar down or to the right, a negative value moves it up or to the left. */
  virtual void moveScrollbarInSteps(const int howManySteps)
  {
    setCurrentRangeStart(rangeStart + howManySteps * singleStepSize);
  }

  /** Moves the scroll bar up or down in pages. This will move the bar by a multiple of its current 
  thumb size, effectively doing a page-up or down. A positive value here will move the bar down or 
  to the right, a negative value moves it up or to the left. */
  virtual void moveScrollbarInPages(const int howManyPages) 
  { 
    setCurrentRangeStart (rangeStart + howManyPages * rangeSize); 
  }

  /** Scrolls to the top (or left). */
  virtual void scrollToTop() { setCurrentRangeStart(min); }

  /** Scrolls to the bottom (or right). */
  virtual void scrollToBottom() { setCurrentRangeStart(max - rangeSize); }

  /** Changes the delay before the up and down buttons autorepeat when they are held down. For an 
  explanation of what the parameters are for, see Button::setRepeatSpeed(). 
  @see Button::setRepeatSpeed */
  virtual void setButtonRepeatSpeed (const int initialDelayInMillisecs, 
    const int repeatDelayInMillisecs, const int minimumDelayInMillisecs = -1) throw();

  /** Registers a listener that will be called when the scrollbar is moved. */
  virtual void addListener(RScrollBarListener* const listener) throw();

  /** Deregisters a previously-registered listener. */
  virtual void removeListener(RScrollBarListener* const listener) throw();

  /** Sets up the colour-scheme - overriden in order to pass the colorscheme through to the 
  buttons. */
  virtual void setColourScheme(const WidgetColourScheme& newColourScheme);

  /** Sets up the colour-scheme - overriden in order to pass the description field through to the 
  buttons. */
  virtual void setDescriptionField(RTextField *newDescriptionField);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the lower value that the thumb can be set to. This is the value set by 
  setRangeLimits().  */
  virtual double getMinimumRangeLimit() const throw() { return min; }

  /** Returns the upper value that the thumb can be set to. This is the value set by 
  setRangeLimits(). */
  virtual double getMaximumRangeLimit() const throw() { return max; }

  /** Returns the position of the top of the thumb. @see setCurrentRangeStart */
  virtual double getCurrentRangeStart() const throw() { return rangeStart; }

  /** Returns the current size of the thumb. @see setCurrentRange */
  virtual double getCurrentRangeSize() const throw() { return rangeSize; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton* buttonThatWasClicked);
  virtual bool keyPressed(const KeyPress& key);

  virtual void mouseWheelMove(const MouseEvent &e, const MouseWheelDetails &wheel);
  //virtual void mouseWheelMove(const MouseEvent& e, float wheelIncrementX, float wheelIncrementY);

  virtual void mouseDown(const MouseEvent& e);
  virtual void mouseDrag(const MouseEvent& e);
  virtual void mouseUp(const MouseEvent& e);
  virtual void resized();
  virtual void paint(Graphics& g);

protected:

  /** Notifies our listeneres that the scrollbar has moved or its size has changed. */
  virtual void notifyListeners();

  /** Returns the thumb's bound in pixels. */
  virtual void getThumbBoundsInPixels(int &x, int &y, int &w, int &h);

  /** Applies the indents for drawing the scrollbar's thumb to the passed pixel-bounds - these 
  bounds should represent the thumb without indents. */
  virtual void applyIndents(int &x, int &y, int &w, int &h);

  /** Returns the maximum length for the thumb which is the full length of the scrollbar (width or 
  height) minus the space for the buttons. */
  virtual int getMaxThumbSizeInPixels();

  /** Returns the current length of the thumb. */
  virtual int getThumbSizeInPixels();

  /** Returns the size of the up/down left/right buttons in pixles. */
  virtual int getButtonSize();

  // data members:

  double min, max, rangeStart, rangeSize, singleStepSize, dragStartRange;
  int    dragStartMousePos;
  bool   vertical, isDraggingThumb;

  RClickButtonWithAutoRepeat *forwardButton, *backwardButton;

  juce::Array<RScrollBarListener*, CriticalSection> listeners;


  RScrollBar(const RScrollBar&);
  const RScrollBar& operator= (const RScrollBar&);
  juce_UseDebuggingNewOperator;
};

#endif   
