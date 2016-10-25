#include "rojue_ComponentScrollContainer.h"
using namespace rojue;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

ComponentScrollContainer::ComponentScrollContainer(Component *contentComponentToScroll)
{
  contentComponent = contentComponentToScroll;
  contentComponent->addComponentListener(this);
  addAndMakeVisible(contentComponent);

  scrollBarThickness = 16;
  xOffset            = 0;
  yOffset            = 0;

  leftRightScrollBar = new RScrollBar(false);
  leftRightScrollBar->setRangeLimits( 0.0, jmax(contentComponent->getWidth(), 1));
  leftRightScrollBar->setCurrentRange(0.0, jmax(contentComponent->getWidth(), 1));
  leftRightScrollBar->setSingleStepSize(1.0);
  leftRightScrollBar->addListener(this);
  addWidget(leftRightScrollBar);

  upDownScrollBar = new RScrollBar(true);
  upDownScrollBar->setRangeLimits( 0.0, jmax(contentComponent->getHeight(), 1));
  upDownScrollBar->setCurrentRange(0.0, jmax(contentComponent->getHeight(), 1));
  upDownScrollBar->setSingleStepSize(1.0);
  upDownScrollBar->addListener(this);
  addWidget(upDownScrollBar);

  bottomRightCoverage = new RWidget();
  addWidget(bottomRightCoverage);
}

ComponentScrollContainer::~ComponentScrollContainer()
{
  contentComponent->removeComponentListener(this);
  deleteAllChildren();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void ComponentScrollContainer::scrollBarMoved(RScrollBar* scrollBarThatHasMoved, const double newRangeStart)
{
  if( scrollBarThatHasMoved == upDownScrollBar )
    yOffset = -roundDoubleToInt(newRangeStart);
  else if( scrollBarThatHasMoved == leftRightScrollBar )
    xOffset = -roundDoubleToInt(newRangeStart);
  contentComponent->setTopLeftPosition(xOffset, yOffset); // we'll see if thsi works
}

void ComponentScrollContainer::paint(Graphics &g)
{
  // overriden to suppress the gradient drawing that is inherited from ColourSchemeComponent
}

void ComponentScrollContainer::paintOverChildren(Graphics &g)
{
  // overriden to suppress the outline drawing that is inherited from ColourSchemeComponent
  // maybe make this optional - that is, let the user choose whether or not to draw th outline
}

void ComponentScrollContainer::resized()
{
  Component::resized();
  updateScrollBarBoundsAndVisibility();
}

void ComponentScrollContainer::componentMovedOrResized(Component &component, bool wasMoved, bool wasResized)
{
  // thsi stuff does not yet work...
  /*
  // new:
  //xOffset = contentComponent->getX();
  //yOffset = contentComponent->getY();
  int x = component.getX();
  int y = component.getY();
  int w = component.getWidth();
  int h = component.getHeight();
  if( x < 0 && w < getWidth() )
    w = getWidth() - x;
  if( y < 0 && h < getHeight() )
    w = getHeight() - y;
  component.setSize(w, h);
  */

  // old:
  if( wasResized )
  {
    updateScrollBarBoundsAndVisibility();
    //contentComponent->setTopLeftPosition(0, 0);  
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void ComponentScrollContainer::updateScrollBarBoundsAndVisibility()
{
  // the logic to determine the available width/height and necessity for scrollbars are interrelated in a somewhat messy way...

  int availableHeight = getHeight();
  int requiredHeight  = contentComponent->getHeight();
  int availableWidth  = getWidth();
  int requiredWidth   = contentComponent->getWidth();

  bool upDownBarNeeded     = requiredHeight > availableHeight;
  bool leftRightBarNeeded  = requiredWidth  > availableWidth;
  bool leftRightBarNeeded2 = false;
  bool upDownBarNeeded2    = false;

  if( upDownBarNeeded )
  {
    availableWidth      -= scrollBarThickness;
    leftRightBarNeeded2  = requiredWidth  > availableWidth;
  }
  if( leftRightBarNeeded )
  {
    availableHeight  -= scrollBarThickness;
    upDownBarNeeded2  = requiredHeight > availableHeight;
  }

  if( (upDownBarNeeded == false) &&  (upDownBarNeeded2 == true) )
  {
    upDownBarNeeded  = true;
    availableWidth  -= scrollBarThickness;
  }
  if( (leftRightBarNeeded == false) && (leftRightBarNeeded2 == true) )
  {
    leftRightBarNeeded  = true;
    availableHeight    -= scrollBarThickness;
  }

  // OK, done with the mess - now set up the scrollbars:
  if( upDownBarNeeded )
  {
    upDownScrollBar->setVisible(true);
    upDownScrollBar->setBounds(getWidth()-scrollBarThickness, 0, scrollBarThickness, availableHeight);
    upDownScrollBar->setRangeLimits(0.0, requiredHeight);
    upDownScrollBar->setCurrentRange(-yOffset, availableHeight);
  }
  else
  {
    yOffset = 0;
    upDownScrollBar->setVisible(false);
  }

  if( leftRightBarNeeded )
  {
    leftRightScrollBar->setVisible(true);
    leftRightScrollBar->setBounds(0, getHeight()-scrollBarThickness, availableWidth, scrollBarThickness);
    leftRightScrollBar->setRangeLimits(0.0, requiredWidth);  // old
    //leftRightScrollBar->setRangeLimits(contentComponent->getX(), requiredWidth-contentComponent->getX());
    //xOffset = contentComponent->getX(); // new
    leftRightScrollBar->setCurrentRange(-xOffset, availableWidth);
  }
  else
  {
    xOffset = 0;
    leftRightScrollBar->setVisible(false);
  }

  if( upDownBarNeeded && leftRightBarNeeded )
  {
    upDownScrollBar->setBounds(getWidth()-scrollBarThickness, 0, scrollBarThickness, availableHeight+RWidget::outlineThickness);
    leftRightScrollBar->setBounds(0, getHeight()-scrollBarThickness, availableWidth+RWidget::outlineThickness, scrollBarThickness);
    bottomRightCoverage->setBounds(getWidth()-scrollBarThickness, getHeight()-scrollBarThickness, scrollBarThickness, scrollBarThickness);
    bottomRightCoverage->setVisible(true);
  }
  else
    bottomRightCoverage->setVisible(false);
}
