
//-------------------------------------------------------------------------------------------------
// construction/destruction and static member initialization:

const double Panel::minGridInterval = 0.000000001;

Panel::Panel(const String &componentName) : Component(componentName)
{
  scaleX                        =  1.0;
  scaleY                        =  1.0;
  horizontalCoarseGridInterval  =  1.0;
  horizontalFineGridInterval    =  0.1;
  verticalCoarseGridInterval    =  1.0;
  verticalFineGridInterval      =  0.1;
  minimumWidth                  =  0.001;
  minimumHeight                 =  0.001;

  maximumRange.setRangeX(-2.2, 2.2);
  maximumRange.setRangeY(-2.2, 2.2);
  currentRange.setRangeX(-2.2, 2.2);
  currentRange.setRangeY(-2.2, 2.2);
}

Panel::~Panel()
{
  //deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void Panel::setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  setMaximumRange( PanelRange(newMinX, newMaxX, newMinY, newMaxY) );
}

void Panel::setMaximumRange(PanelRange newMaximumRange)
{
  if( maximumRange != newMaximumRange )
  {
    maximumRange = newMaximumRange;
    constrainCurrentRange(true);
  }
}

void Panel::setMaximumRangeX(double newMinX, double newMaxX)
{
  if( newMinX != maximumRange.getMinX() || newMaxX != maximumRange.getMaxX() )
  {
    maximumRange.setRangeX(newMinX, newMaxX);
    constrainCurrentRange(true);
  }
}

void Panel::setMaximumRangeY(double newMinY, double newMaxY)
{
  if( newMinY != maximumRange.getMinY() || newMaxY != maximumRange.getMaxY() )
  {
    maximumRange.setRangeY(newMinY, newMaxY);
    constrainCurrentRange(true);
  }
}

void Panel::setMaximumRangeMinX(double newMinX)
{
  if( newMinX != maximumRange.getMinX() )
  {
    maximumRange.setMinX(newMinX);
    constrainCurrentRange(true);
  }
}

void Panel::setMaximumRangeMaxX(double newMaxX)
{
  if( newMaxX != maximumRange.getMaxX() )
  {
    maximumRange.setMaxX(newMaxX);
    constrainCurrentRange(true);
  }
}

void Panel::setMaximumRangeMinY(double newMinY)
{
  if( newMinY != maximumRange.getMinY() )
  {
    maximumRange.setMinY(newMinY);
    constrainCurrentRange(true);
  }
}

void Panel::setMaximumRangeMaxY(double newMaxY)
{
  if( newMaxY != maximumRange.getMaxY() )
  {
    maximumRange.setMaxY(newMaxY);
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  setCurrentRange( PanelRange(newMinX, newMaxX, newMinY, newMaxY) );
}

void Panel::setCurrentRange(PanelRange newRange)
{
  if( currentRange != newRange )
  {
    currentRange = newRange;
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRangeX(double newMinX, double newMaxX)
{
  if( newMinX != currentRange.getMinX() || newMaxX != currentRange.getMaxX() )
  {
    currentRange.setRangeX(newMinX, newMaxX);
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRangeY(double newMinY, double newMaxY)
{
  if( newMinY != currentRange.getMinY() || newMaxY != currentRange.getMaxY() )
  {
    currentRange.setRangeY(newMinY, newMaxY);
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRangeMinX(double newMinX)
{
  if( newMinX != currentRange.getMinX() )
  {
    currentRange.setMinX(newMinX);
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRangeMaxX(double newMaxX)
{
  if( newMaxX != currentRange.getMaxX() )
  {
    currentRange.setMaxX(newMaxX);
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRangeMinY(double newMinY)
{
  if( newMinY != currentRange.getMinY() )
  {
    currentRange.setMinY(newMinY);
    constrainCurrentRange(true);
  }
}

void Panel::setCurrentRangeMaxY(double newMaxY)
{
  if( newMaxY != currentRange.getMaxY() )
  {
    currentRange.setMaxY(newMaxY);
    constrainCurrentRange(true);
  }
}

void Panel::setMinimumWidth(double newMinWidth)
{
  if( newMinWidth != minimumWidth && newMinWidth > 0.0 )
  {
    minimumWidth = newMinWidth;
    constrainCurrentRange(true);
  }
}

void Panel::setMinimumHeight(double newMinHeight)
{
  if( newMinHeight != minimumHeight && newMinHeight > 0.0 )
  {
    minimumHeight = newMinHeight;
    constrainCurrentRange(true);
  }
}

void Panel::setHorizontalCoarseGridInterval(double newGridInterval)
{
  jassert(newGridInterval >= minGridInterval);   // too small interval
  if( newGridInterval >= minGridInterval && horizontalCoarseGridInterval != newGridInterval)
  {
    horizontalCoarseGridInterval  = newGridInterval;
    setDirty();
  }
}

void Panel::setHorizontalFineGridInterval(double newGridInterval)
{
  jassert(newGridInterval >= minGridInterval);   // too small interval
  if( newGridInterval >= minGridInterval && horizontalFineGridInterval != newGridInterval)
  {
    horizontalFineGridInterval  = newGridInterval;
    setDirty();
  }
}

void Panel::setVerticalCoarseGridInterval(double newGridInterval)
{
  jassert(newGridInterval >= minGridInterval);   // too small interval
  if( newGridInterval >= minGridInterval && verticalCoarseGridInterval != newGridInterval)
  {
    verticalCoarseGridInterval  = newGridInterval;
    setDirty();
  }
}

void Panel::setVerticalFineGridInterval(double newGridInterval)
{
  jassert(newGridInterval >= minGridInterval);   // too small interval
  if( newGridInterval >= minGridInterval && verticalFineGridInterval != newGridInterval)
  {
    verticalFineGridInterval  = newGridInterval;
    setDirty();
  }
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void Panel::resized()
{
  updateScaleFactors();
  //ThreadedDrawingComponent::resized();
}
  
/*
void Panel::drawComponent(Image* imageToDrawOnto)
{
  Graphics g(*imageToDrawOnto);
  int w = imageToDrawOnto->getWidth();
  int h = imageToDrawOnto->getHeight();
  g.setColour(Colours::black);
  //g.drawFittedText(String(T("Panel")), 0, 0, w, h, Justification::centred, 1);
    // triggers a JUCE-breakpoint when called early on app-startup
}
*/

//-------------------------------------------------------------------------------------------------
// coordinate transformations:

void Panel::transformToComponentsCoordinates(double &x, double &y) const
{
  x -= currentRange.getMinX();	 // shift origin left/right
  x *= scaleX;	                 // scale to fit width
  y -= currentRange.getMinY();   // shift origin up/down
  y *= scaleY;	                 // scale to fit height
  y  = getHeight()-y;	           // invert (pixels begin at top-left)
}

void Panel::transformToComponentsCoordinates(float &x, float &y) const
{
  double xd = (double) x;
  double yd = (double) y;
  transformToComponentsCoordinates(xd, yd);
  x = (float) xd;
  y = (float) yd;
}

void Panel::transformFromComponentsCoordinates(double &x, double &y) const
{
  x /= scaleX;                    // scale to fit width
  x += currentRange.getMinX();    // shift origin left/right
  y  = getHeight()-y;             // shift origin up/down
  y /= scaleY;                    // scale to fit height
  y += currentRange.getMinY();    // shift origin up/down
}

void Panel::transformFromComponentsCoordinates(float &x, float &y) const
{
  double xd = (double) x;
  double yd = (double) y;
  transformFromComponentsCoordinates(xd, yd);
  x = (float) xd;
  y = (float) yd;
}

//-------------------------------------------------------------------------------------------------
// others:

void Panel::setDirty(bool shouldSetToDirty)
{
  updateScaleFactors();
}

void Panel::constrainCurrentRange(bool callSetDirty)
{ 
  // ensure a minimum range:
  if( currentRange.getMaxX()-currentRange.getMinX() < minimumWidth )
    currentRange.setMaxX(currentRange.getMinX() + minimumWidth);
  if( currentRange.getMaxY()-currentRange.getMinY() < minimumHeight )
    currentRange.setMaxY(currentRange.getMinY() + minimumHeight);

  // ensure to not exceed the maxmimum range:
  currentRange.clipRange(maximumRange);

  updateScaleFactors();
  if( callSetDirty == true )
    setDirty();
}

void Panel::updateScaleFactors()
{  
  scaleX = getWidth()  / (currentRange.getMaxX()-currentRange.getMinX()); 
  scaleY = getHeight() / (currentRange.getMaxY()-currentRange.getMinY());
}











