//#include "romos_NoteEvent.h"
//using namespace romos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

NoteEvent::NoteEvent(unsigned int _deltaFrames, unsigned int _key, unsigned int _velocity)
{
  this->deltaFrames = _deltaFrames;
  this->key         = _key;
  this->velocity    = _velocity;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// related non-member functions:

bool noteEventLessByDeltaFrames(NoteEvent left, NoteEvent right)
{
  if( left.getDeltaFrames() < right.getDeltaFrames() )
    return true;
  else if( left.getDeltaFrames() > right.getDeltaFrames() )
    return false;
  else
  {
    if( left.getVelocity() < right.getVelocity() )
      return true;
    else if( left.getVelocity() > right.getVelocity() )
      return false;
    else
    {
      if( left.getKey() < right.getKey() )
        return true;
      else if( left.getKey() > right.getKey() )
        return false;
    }
  }
  return false;
}
