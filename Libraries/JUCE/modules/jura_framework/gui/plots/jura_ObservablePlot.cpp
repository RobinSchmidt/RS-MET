

rsObservablePlot::rsObservablePlot(const String &name) 
: rsPlot(name)
{

}

rsObservablePlot::~rsObservablePlot()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// range management:

void rsObservablePlot::setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  rsPlotRange r = getMaximumRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() || newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    rsPlot::setMaximumRange(newMinX, newMaxX, newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRange(rsPlotRange newMaximumRange)
{
  if( newMaximumRange != getMaximumRange() )
  {
    rsPlot::setMaximumRange(newMaximumRange);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRangeX(double newMinX, double newMaxX)
{
  rsPlotRange r = getMaximumRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() )
  {
    rsPlot::setMaximumRangeX(newMinX, newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRangeY(double newMinY, double newMaxY)
{
  rsPlotRange r = getMaximumRange();
  if( newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    rsPlot::setMaximumRangeY(newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRangeMinX(double newMinX)
{
  if( newMinX != getMaximumRange().getMinX() )
  {
    rsPlot::setMaximumRangeMinX(newMinX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRangeMaxX(double newMaxX)
{
  if( newMaxX != getMaximumRange().getMaxX() )
  {
    rsPlot::setMaximumRangeMaxX(newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRangeMinY(double newMinY)
{
  if( newMinY != getMaximumRange().getMinY() )
  {
    rsPlot::setMaximumRangeMinY(newMinY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setMaximumRangeMaxY(double newMaxY)
{
  if( newMaxY != getMaximumRange().getMaxY() )
  {
    rsPlot::setMaximumRangeMaxY(newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  rsPlotRange r = getCurrentRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() || newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    rsPlot::setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRange(rsPlotRange newCurrentRange)
{
  if( newCurrentRange != getCurrentRange() )
  {
    rsPlot::setCurrentRange(newCurrentRange);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRangeX(double newMinX, double newMaxX)
{
  rsPlotRange r = getCurrentRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() )
  {
    rsPlot::setCurrentRangeX(newMinX, newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRangeY(double newMinY, double newMaxY)
{
  rsPlotRange r = getCurrentRange();
  if( newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    rsPlot::setCurrentRangeY(newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRangeMinX(double newMinX)
{
  if( newMinX != getCurrentRange().getMinX() )
  {
    rsPlot::setCurrentRangeMinX(newMinX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRangeMaxX(double newMaxX)
{
  if( newMaxX != getCurrentRange().getMaxX() )
  {
    rsPlot::setCurrentRangeMaxX(newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRangeMinY(double newMinY)
{
  if( newMinY != getCurrentRange().getMinY() )
  {
    rsPlot::setCurrentRangeMinY(newMinY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void rsObservablePlot::setCurrentRangeMaxY(double newMaxY)
{
  if( newMaxY != getCurrentRange().getMaxY() )
  {
    rsPlot::setCurrentRangeMaxY(newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}





/*
void rsObservablePlot::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  rsPlot::setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRange(CoordinateSystemRangeOld newRange)
{
  rsPlot::setCurrentRange(newRange);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRangeX(double newMinX, double newMaxX)
{
  rsPlot::setCurrentRangeX(newMinX, newMaxX);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRangeY(double newMinY, double newMaxY)
{
  rsPlot::setCurrentRangeY(newMinY, newMaxY);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRangeMinX(double newMinX)
{
  rsPlot::setCurrentRangeMinX(newMinX);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRangeMaxX(double newMaxX)
{
  rsPlot::setCurrentRangeMaxX(newMaxX);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRangeMinY(double newMinY)
{
  rsPlot::setCurrentRangeMinY(newMinY);
  sendCoordinateSystemChangedMessage(this);
}

void rsObservablePlot::setCurrentRangeMaxY(double newMaxY)
{
  rsPlot::setCurrentRangeMaxY(newMaxY);
  sendCoordinateSystemChangedMessage(this);
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------
// others;

void rsObservablePlot::addCoordinateSystemOldObserver(CoordinateSystemOldObserver* observerToAdd)
{
  observers.getLock().enter();
  observers.addIfNotAlreadyThere(observerToAdd);
  observers.getLock().exit();
}

void rsObservablePlot::removeCoordinateSystemOldObserver(CoordinateSystemOldObserver* observerToRemove)
{
  observers.getLock().enter();
  observers.removeFirstMatchingValue(observerToRemove);
  observers.getLock().exit();
}

void rsObservablePlot::removeAllCoordinateSystemOldObservers()
{
  observers.getLock().enter();
  observers.clear();
  observers.getLock().exit();
}

void rsObservablePlot::sendCoordinateSystemChangedMessage(rsObservablePlot *coordinateSystemThatHasChanged)
{
  observers.getLock().enter();
  for(int i=0; i<observers.size(); i++)
    observers[i]->coordinateSystemChanged(coordinateSystemThatHasChanged);
  observers.getLock().exit();
}








