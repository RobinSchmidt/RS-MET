
//-------------------------------------------------------------------------------------------------
// construction/destruction:

PanelRange::PanelRange(double initMinX, double initMaxX, double initMinY, double initMaxY)
{
  minX = initMinX;
  maxX = initMaxX;
  minY = initMinY;
  maxY = initMaxY;
}

PanelRange::~PanelRange()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void PanelRange::setMinX(double newMinX) 
{ 
  jassert( newMinX < maxX );
  if( newMinX < maxX )
    minX = newMinX;
}

void PanelRange::setMaxX(double newMaxX) 
{ 
  jassert( newMaxX > minX );
  if( newMaxX > minX )
    maxX = newMaxX;
}

void PanelRange::setRangeX(double newMinX, double newMaxX)
{
  jassert( newMaxX > newMinX );
  if( newMaxX > newMinX )
  {
    minX = newMinX;
    maxX = newMaxX;
  }
}

void PanelRange::setMinY(double newMinY) 
{ 
  jassert( newMinY < maxY );
  if( newMinY < maxY )
    minY = newMinY;
}

void PanelRange::setMaxY(double newMaxY) 
{ 
  jassert( newMaxY > minY );
  if( newMaxY > minY )
    maxY = newMaxY;
}

void PanelRange::setRangeY(double newMinY, double newMaxY)
{
  jassert( newMaxY > newMinY );
  if( newMaxY > newMinY )
  {
    minY = newMinY;
    maxY = newMaxY;
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void PanelRange::clipRange(PanelRange rangeToClipTo)
{
  if( minX < rangeToClipTo.getMinX() )
    minX = rangeToClipTo.getMinX();
  if( maxX > rangeToClipTo.getMaxX() )
    maxX = rangeToClipTo.getMaxX();
  if( minY < rangeToClipTo.getMinY() )
    minY = rangeToClipTo.getMinY();
  if( maxY > rangeToClipTo.getMaxY() )
    maxY = rangeToClipTo.getMaxY();
}


