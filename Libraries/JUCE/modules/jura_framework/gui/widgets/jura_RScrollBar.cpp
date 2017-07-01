//#include "rojue_RScrollBar.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RScrollBar::RScrollBar(const bool isVertical)
{ 
  vertical          = isVertical;
  min               = 0.0;
  max               = 1.0;
  rangeStart        = 1.0;
  rangeSize         = 0.1;
  singleStepSize    = 0.1;
  isDraggingThumb   = false;
  dragStartRange    = 0.0;
  dragStartMousePos = 0;; 

  if( vertical )
  {
    addAndMakeVisible( forwardButton  = new RClickButtonWithAutoRepeat(RButton::ARROW_DOWN) ); 
    addAndMakeVisible( backwardButton = new RClickButtonWithAutoRepeat(RButton::ARROW_UP)   ); 
    setDescription(String("Scroll vertically"));
    forwardButton->setDescription(String("Scroll down"));
    backwardButton->setDescription(String("Scroll up"));
  }
  else
  {
    addAndMakeVisible( forwardButton  = new RClickButtonWithAutoRepeat(RButton::ARROW_RIGHT) ); 
    addAndMakeVisible( backwardButton = new RClickButtonWithAutoRepeat(RButton::ARROW_LEFT)  ); 
    setDescription(String("Scroll horizontally"));
    forwardButton->setDescription(String("Scroll right"));
    backwardButton->setDescription(String("Scroll left"));
  }

  //forwardButton->setClickingTogglesState(false); // re-activate
  forwardButton->addRButtonListener(this);

  //backwardButton->setClickingTogglesState(false); // re-activate
  backwardButton->addRButtonListener(this);

  setButtonRepeatSpeed(100, 50, 10);
}

RScrollBar::~RScrollBar()
{  
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void RScrollBar::setRangeLimits(const double newMinimum, const double newMaximum) throw()
{
  min = newMinimum;
  max = newMaximum;
  jassert(max >= min);        // these can't be the wrong way round!
  setCurrentRangeStart(rangeStart);
  repaint();
}

void RScrollBar::setCurrentRange(double newStart, double newSize) throw()
{
  newSize  = jlimit(0.0, max-min,     newSize);
  newStart = jlimit(min, max-newSize, newStart);
  if(rangeStart != newStart || rangeSize != newSize) 
  {
    rangeStart = newStart;
    rangeSize  = newSize;
    repaint();
    notifyListeners(); // maybe make this conditional
  }
}

void RScrollBar::setButtonRepeatSpeed(const int initialDelayInMillisecs, 
  const int repeatDelayInMillisecs, const int minimumDelayInMillisecs) throw()
{
  //forwardButton->setRepeatSpeed( initialDelayInMillisecs,  repeatDelayInMillisecs,  minimumDelayInMillisecs);
  //backwardButton->setRepeatSpeed(initialDelayInMillisecs,  repeatDelayInMillisecs,  minimumDelayInMillisecs);
  // re-activate these
}

void RScrollBar::addListener(RScrollBarListener* const listener) throw()
{
  listeners.getLock().enter();
  listeners.addIfNotAlreadyThere(listener);
  listeners.getLock().exit();
}

void RScrollBar::removeListener(RScrollBarListener* const listener) throw()
{
  listeners.getLock().enter();
  listeners.removeFirstMatchingValue(listener);
  listeners.getLock().exit();
}

void RScrollBar::setColourScheme(const WidgetColourScheme& newColourScheme)
{
  RWidget::setColourScheme(newColourScheme);
  forwardButton->setColourScheme(newColourScheme);
  backwardButton->setColourScheme(newColourScheme);
}

void RScrollBar::setDescriptionField(RTextField *newDescriptionField)
{
  RWidget::setDescriptionField(newDescriptionField);
  forwardButton->setDescriptionField(newDescriptionField);
  backwardButton->setDescriptionField(newDescriptionField);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RScrollBar::rButtonClicked(RButton* buttonThatWasClicked)
{
  if( buttonThatWasClicked == forwardButton )
    moveScrollbarInSteps(1);
  else if( buttonThatWasClicked == backwardButton )
    moveScrollbarInSteps(-1);
}

void RScrollBar::notifyListeners()
{
  const double value = getCurrentRangeStart();
  listeners.getLock().enter();
  for(int i=0; i<listeners.size(); i++)
    listeners[i]->scrollBarMoved(this, value);
  listeners.getLock().exit();
}

void RScrollBar::resized()
{
  if( vertical )
  {
    int buttonSize = getWidth();
    forwardButton->setBounds( 0, getHeight()-buttonSize, buttonSize, buttonSize);
    backwardButton->setBounds(0, 0,                      buttonSize, buttonSize);
  }
  else
  {
    int buttonSize = getHeight();
    forwardButton->setBounds( getWidth()-buttonSize, 0, buttonSize, buttonSize);
    backwardButton->setBounds(0,                     0, buttonSize, buttonSize);
  }
}

void RScrollBar::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());

  g.setColour(getHandleColour());
  int x, y, w, h;
  getThumbBoundsInPixels(x, y, w, h);
  applyIndents(x, y, w, h);
  g.fillRect(x, y, w, h);

  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), RWidget::outlineThickness);

  //if( !isEnabled() )
  //  g.fillAll(Colours::lightgrey.withAlpha(0.75f));
}

void RScrollBar::mouseDown (const MouseEvent& e)
{
  dragStartRange = rangeStart;
  int x, y, w, h;
  getThumbBoundsInPixels(x, y, w, h);

  isDraggingThumb = false;
  if( vertical )
  {
    dragStartMousePos = e.y;
    if( dragStartMousePos < y )
      moveScrollbarInPages(-1);
    else if( dragStartMousePos >= y+h )
      moveScrollbarInPages(1);
    else
      isDraggingThumb = true;
  }
  else
  {
    dragStartMousePos = e.x;
    if( dragStartMousePos < x )
      moveScrollbarInPages(-1);
    else if( dragStartMousePos >= x+w )
      moveScrollbarInPages(1);
    else
      isDraggingThumb = true;
  }
}

void RScrollBar::mouseDrag(const MouseEvent& e)
{
  int thumbSize    = getThumbSizeInPixels();
  int maxThumbSize = getMaxThumbSizeInPixels();
  if( thumbSize == maxThumbSize )
    return;
  if(isDraggingThumb)
  {
    int    deltaPixels = ((vertical) ? e.y : e.x) - dragStartMousePos;
    double newStart    = dragStartRange + deltaPixels * 
                         ((max-min)-rangeSize) / (maxThumbSize - thumbSize);
    setCurrentRangeStart(newStart);
  }

  /*
  // old (buggy)
  if(isDraggingThumb)
  {
    const int deltaPixels = ((vertical) ? e.y : e.x) - dragStartMousePos;
    setCurrentRangeStart(dragStartRange + deltaPixels * ((max-min)-rangeSize) / (getMaxThumbSizeInPixels() - getThumbSizeInPixels()));
  }
  */
}

void RScrollBar::mouseUp(const MouseEvent&)
{
  isDraggingThumb = false;
  repaint();
}

void RScrollBar::mouseWheelMove(const MouseEvent &e, const MouseWheelDetails &wheel)
{
  //float increment = vertical ? wheelIncrementY : wheelIncrementX;  // old
  float increment = vertical ? wheel.deltaY : wheel.deltaX;
  if(increment < 0)
    increment = jmin(increment * 10.0f, -1.0f);
  else if(increment > 0)
    increment = jmax(increment * 10.0f, 1.0f);
  setCurrentRangeStart(rangeStart - singleStepSize*rangeSize * increment);
}
// old:
//void RScrollBar::mouseWheelMove(const MouseEvent&, float wheelIncrementX, float wheelIncrementY)
//{
//  float increment = vertical ? wheelIncrementY : wheelIncrementX;
//  if(increment < 0)
//    increment = jmin(increment * 10.0f, -1.0f);
//  else if(increment > 0)
//    increment = jmax(increment * 10.0f, 1.0f);
//  setCurrentRangeStart(rangeStart - singleStepSize * increment);
//}

bool RScrollBar::keyPressed (const KeyPress& key)
{
  if( !isVisible() )
    return false;

  if(key.isKeyCode(KeyPress::upKey) || key.isKeyCode(KeyPress::leftKey))
    moveScrollbarInSteps(-1);
  else if(key.isKeyCode(KeyPress::downKey) || key.isKeyCode(KeyPress::rightKey))
    moveScrollbarInSteps(1);
  else if(key.isKeyCode(KeyPress::pageUpKey))
    moveScrollbarInPages(-1);
  else if(key.isKeyCode(KeyPress::pageDownKey))
    moveScrollbarInPages(1);
  else if(key.isKeyCode(KeyPress::homeKey))
    scrollToTop();
  else if(key.isKeyCode(KeyPress::endKey))
    scrollToBottom();
  else
    return false;

  return true;
}

void RScrollBar::getThumbBoundsInPixels(int &x, int &y, int &w, int &h)
{
  double normalizedRange = rangeSize  / (max-min);
  double normalizedStart = rangeStart / (max-min);
  double pixelStart      = getButtonSize() + normalizedStart * getMaxThumbSizeInPixels();  
  double pixelSize       =                   normalizedRange * getMaxThumbSizeInPixels();

  if( pixelSize < getButtonSize() )
  {
    double normalizedMid = normalizedStart + 0.5 * normalizedRange;
    pixelSize            = getButtonSize();
    normalizedRange      = pixelSize / (double) getMaxThumbSizeInPixels();
    normalizedStart      = normalizedMid - 0.5 * normalizedRange;
    pixelStart           = getButtonSize() + normalizedStart * getMaxThumbSizeInPixels();  
    if( pixelStart < getButtonSize() )
      pixelStart = getButtonSize();
    else if( pixelStart > getMaxThumbSizeInPixels() )
      pixelStart = getMaxThumbSizeInPixels();
  }

  if( vertical )
  {
    x = 0;
    y = (int) floor(pixelStart);
    w = getWidth();
    h = (int) ceil(pixelSize); 
  }
  else
  {
    x = (int) floor(pixelStart); 
    y = 0;
    w = (int) ceil(pixelSize); 
    h = getHeight();
  }
}

void RScrollBar::applyIndents(int &x, int &y, int &w, int &h)
{
  int indent = 2;
  if( vertical )
  {
    x += indent + RWidget::outlineThickness;
    y += indent;
    w -= 2*indent + 2*RWidget::outlineThickness;
    h -= 2*indent;
  }
  else 
  {
    x += indent;
    y += indent + RWidget::outlineThickness;
    w -= 2*indent;
    h -= 2*indent + 2*RWidget::outlineThickness;
  }
}

int RScrollBar::getMaxThumbSizeInPixels()
{
  if( vertical )
    return getHeight() - 2*getButtonSize();
  else
    return getWidth() - 2*getButtonSize();
}

int RScrollBar::getThumbSizeInPixels()
{
  int x, y, w, h;
  getThumbBoundsInPixels(x, y, w, h);
  if( vertical )
    return h;
  else
    return w;
}

int RScrollBar::getButtonSize()
{
  if( vertical )
    return getWidth();
  else
    return getHeight();
}
