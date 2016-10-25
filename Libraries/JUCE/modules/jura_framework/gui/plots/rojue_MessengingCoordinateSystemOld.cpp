#include "rojue_MessengingCoordinateSystemOld.h"
using namespace rojue;

MessengingCoordinateSystemOld::MessengingCoordinateSystemOld(const String &name) 
: CoordinateSystemOld(name)
{

}

MessengingCoordinateSystemOld::~MessengingCoordinateSystemOld()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// range management:

void MessengingCoordinateSystemOld::setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  CoordinateSystemRangeOld r = getMaximumRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() || newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    CoordinateSystemOld::setMaximumRange(newMinX, newMaxX, newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRange(CoordinateSystemRangeOld newMaximumRange)
{
  if( newMaximumRange != getMaximumRange() )
  {
    CoordinateSystemOld::setMaximumRange(newMaximumRange);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRangeX(double newMinX, double newMaxX)
{
  CoordinateSystemRangeOld r = getMaximumRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() )
  {
    CoordinateSystemOld::setMaximumRangeX(newMinX, newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRangeY(double newMinY, double newMaxY)
{
  CoordinateSystemRangeOld r = getMaximumRange();
  if( newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    CoordinateSystemOld::setMaximumRangeY(newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRangeMinX(double newMinX)
{
  if( newMinX != getMaximumRange().getMinX() )
  {
    CoordinateSystemOld::setMaximumRangeMinX(newMinX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRangeMaxX(double newMaxX)
{
  if( newMaxX != getMaximumRange().getMaxX() )
  {
    CoordinateSystemOld::setMaximumRangeMaxX(newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRangeMinY(double newMinY)
{
  if( newMinY != getMaximumRange().getMinY() )
  {
    CoordinateSystemOld::setMaximumRangeMinY(newMinY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setMaximumRangeMaxY(double newMaxY)
{
  if( newMaxY != getMaximumRange().getMaxY() )
  {
    CoordinateSystemOld::setMaximumRangeMaxY(newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  CoordinateSystemRangeOld r = getCurrentRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() || newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    CoordinateSystemOld::setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRange(CoordinateSystemRangeOld newCurrentRange)
{
  if( newCurrentRange != getCurrentRange() )
  {
    CoordinateSystemOld::setCurrentRange(newCurrentRange);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRangeX(double newMinX, double newMaxX)
{
  CoordinateSystemRangeOld r = getCurrentRange();
  if( newMinX != r.getMinX() || newMaxX != r.getMaxX() )
  {
    CoordinateSystemOld::setCurrentRangeX(newMinX, newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRangeY(double newMinY, double newMaxY)
{
  CoordinateSystemRangeOld r = getCurrentRange();
  if( newMinY != r.getMinY() || newMaxY != r.getMaxY() )
  {
    CoordinateSystemOld::setCurrentRangeY(newMinY, newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRangeMinX(double newMinX)
{
  if( newMinX != getCurrentRange().getMinX() )
  {
    CoordinateSystemOld::setCurrentRangeMinX(newMinX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRangeMaxX(double newMaxX)
{
  if( newMaxX != getCurrentRange().getMaxX() )
  {
    CoordinateSystemOld::setCurrentRangeMaxX(newMaxX);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRangeMinY(double newMinY)
{
  if( newMinY != getCurrentRange().getMinY() )
  {
    CoordinateSystemOld::setCurrentRangeMinY(newMinY);
    sendCoordinateSystemChangedMessage(this);
  }
}

void MessengingCoordinateSystemOld::setCurrentRangeMaxY(double newMaxY)
{
  if( newMaxY != getCurrentRange().getMaxY() )
  {
    CoordinateSystemOld::setCurrentRangeMaxY(newMaxY);
    sendCoordinateSystemChangedMessage(this);
  }
}





/*
void MessengingCoordinateSystemOld::setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  CoordinateSystemOld::setCurrentRange(newMinX, newMaxX, newMinY, newMaxY);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRange(CoordinateSystemRangeOld newRange)
{
  CoordinateSystemOld::setCurrentRange(newRange);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRangeX(double newMinX, double newMaxX)
{
  CoordinateSystemOld::setCurrentRangeX(newMinX, newMaxX);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRangeY(double newMinY, double newMaxY)
{
  CoordinateSystemOld::setCurrentRangeY(newMinY, newMaxY);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRangeMinX(double newMinX)
{
  CoordinateSystemOld::setCurrentRangeMinX(newMinX);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRangeMaxX(double newMaxX)
{
  CoordinateSystemOld::setCurrentRangeMaxX(newMaxX);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRangeMinY(double newMinY)
{
  CoordinateSystemOld::setCurrentRangeMinY(newMinY);
  sendCoordinateSystemChangedMessage(this);
}

void MessengingCoordinateSystemOld::setCurrentRangeMaxY(double newMaxY)
{
  CoordinateSystemOld::setCurrentRangeMaxY(newMaxY);
  sendCoordinateSystemChangedMessage(this);
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------
// others;

void MessengingCoordinateSystemOld::addCoordinateSystemOldObserver(CoordinateSystemOldObserver* observerToAdd)
{
  observers.getLock().enter();
  observers.addIfNotAlreadyThere(observerToAdd);
  observers.getLock().exit();
}

void MessengingCoordinateSystemOld::removeCoordinateSystemOldObserver(CoordinateSystemOldObserver* observerToRemove)
{
  observers.getLock().enter();
  observers.removeValue(observerToRemove);
  observers.getLock().exit();
}

void MessengingCoordinateSystemOld::removeAllCoordinateSystemOldObservers()
{
  observers.getLock().enter();
  observers.clear();
  observers.getLock().exit();
}

void MessengingCoordinateSystemOld::sendCoordinateSystemChangedMessage(MessengingCoordinateSystemOld *coordinateSystemThatHasChanged)
{
  observers.getLock().enter();
  for(int i=0; i<observers.size(); i++)
    observers[i]->coordinateSystemChanged(coordinateSystemThatHasChanged);
  observers.getLock().exit();
}








