#include "rosic_BreakpointModulator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BreakpointModulator::BreakpointModulator()
{
  data         = new BreakpointModulatorData;
  isMaster     = true;
  currentShape =  ModBreakpoint::LINEAR;
  initialize();
  setToDefaultValues();
}

BreakpointModulator::~BreakpointModulator()
{
  if( isMaster && data != NULL )
  {
    data->mutex.lock();
    data->breakpoints.clear();
    data->mutex.unlock();  
    delete data;
  }
}

void BreakpointModulator::copyDataFrom(const BreakpointModulator &source)
{
  // copy everything except the mutex-member:
  state1                      = source.state1;
  state2                      = source.state2;
  state1_change               = source.state1_change;
  state2_change               = source.state2_change;
  leftLevel                   = source.leftLevel;
  rightLevel                  = source.rightLevel;
  levelDelta                  = source.levelDelta;
  state1_min                  = source.state1_min;
  state2_min                  = source.state2_min;
  state1_max                  = source.state1_max;
  state2_max                  = source.state2_max;
  scaler1                     = source.scaler1;
  scaler2                     = source.scaler2;
  currentShape                = source.currentShape;
  samplesToNextBreakpoint     = source.samplesToNextBreakpoint;
  accumulatedTimingError      = source.accumulatedTimingError;
  leftIndex                   = source.leftIndex;
  rightIndex                  = source.rightIndex;
  data->loopStartIndex        = source.data->loopStartIndex;
  data->loopEndIndex          = source.data->loopEndIndex;
  data->loopIsOn              = source.data->loopIsOn;
  noteIsOn                    = source.noteIsOn;
  outLevelIsConstant          = source.outLevelIsConstant;
  previousOut                 = source.previousOut;
  data->editMode              = source.data->editMode;
  data->syncMode              = source.data->syncMode;
  data->bpm                   = source.data->bpm;
  timeScaleFactor             = source.timeScaleFactor;
  data->breakpoints           = source.data->breakpoints;
  data->sampleRate            = source.data->sampleRate;
  data->minimumAllowedLevel   = source.data->minimumAllowedLevel;
  data->maximumAllowedLevel   = source.data->maximumAllowedLevel;
  data->endLevel              = source.data->endLevel;
  data->endLevelFixedAtZero   = source.data->endLevelFixedAtZero;
  data->minBreakpointDistance = source.data->minBreakpointDistance;
  endIsReached                = source.endIsReached;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void BreakpointModulator::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    data->sampleRate = newSampleRate;
}

void BreakpointModulator::setScaleFactor(double newScaleFactor)
{
  data->scaleFactor = newScaleFactor;
  //markPresetAsDirty();
}

void BreakpointModulator::setOffset(double newOffset)
{
  data->offset = newOffset;
  //markPresetAsDirty();
}

void BreakpointModulator::setMinimumAllowedLevel(double newMinimum)
{
  if( newMinimum < data->maximumAllowedLevel-0.001 )
    data->minimumAllowedLevel = newMinimum;

  for(int i=0; i <= lastBreakpointIndex(); i++)
  {
    if( data->breakpoints[i].level < data->minimumAllowedLevel )
      data->breakpoints[i].level = data->minimumAllowedLevel;
  }
}

void BreakpointModulator::setMaximumAllowedLevel(double newMaximum)
{
  if( newMaximum > data->minimumAllowedLevel+0.001 )
    data->maximumAllowedLevel = newMaximum;

  for(int i=0; i <= lastBreakpointIndex(); i++)
  {
    if( data->breakpoints[i].level > data->maximumAllowedLevel )
      data->breakpoints[i].level = data->maximumAllowedLevel;
  }
}

void BreakpointModulator::fixEndLevelAtZero(bool shouldBeFixed)
{
  data->endLevelFixedAtZero = shouldBeFixed;
}

int BreakpointModulator::insertBreakpoint(double newTimeStamp, 
                                          double newLevel, 
                                          int    newShape, 
                                          double newShapeAmount)
{
  if( newTimeStamp >= getStartTime() )
  {
    //markPresetAsDirty();

    // Loop through the existing data->breakpoints and inspect their time-stamps. As 
    // soon as this time-stamp gets (strictly) larger than the time-stamp of the new 
    // breakpoint to be inserted, we know that this is the index, right before 
    // which the new breakpoint has to be inserted:
    int    bpIndex = 1;  // omit the first breakpoint(index 0)
    double bpTime  = 0.0;
    while( bpIndex < (int) data->breakpoints.size() )
    {
      bpTime = data->breakpoints[bpIndex].timeStamp;
      if( bpTime > newTimeStamp )
      {
        // make sure that the new  breakpoint to be inserted is not too close to 
        // another existing breakpoint:
        if( newTimeStamp-data->breakpoints[bpIndex-1].timeStamp < data->minBreakpointDistance 
          || bpTime-newTimeStamp < data->minBreakpointDistance )
        {
          return -1;
        }

        // create a new breakpoint for insertion right before the breakpoint with 
        // the current index:
        ModBreakpoint newBreakpoint;
        newBreakpoint.timeStamp   = newTimeStamp;
        newBreakpoint.level       = clipLevelToRange(newLevel);

        // copy the shape and shapeAmount from the succeeding breakpoint into,
        // the new one if no shape and shapeAmount were specified (i.e. when they 
        // are equal to their default-arguments which are zero), otherwise use the 
        // passed arguments:
        if( newShape == 0 )
          newBreakpoint.shape = data->breakpoints[bpIndex].shape;
        else
          newBreakpoint.shape = newShape;
        if( newShapeAmount == 0.0 )
          newBreakpoint.shapeAmount = data->breakpoints[bpIndex].shapeAmount;
        else
          newBreakpoint.shapeAmount = newShapeAmount;

        // ! enter critical section !
        data->mutex.lock();

        // insert the new breakpoint right before the breakpoint with the
        // current index:
        data->breakpoints.insert(data->breakpoints.begin()+bpIndex, 1, newBreakpoint);

        // advance the loopStart- and loopEnd-index by one, if the new breakpoint
        // was inserted before the current loop:
        if( bpIndex <= data->loopStartIndex )
        {
          data->loopStartIndex++;
          data->loopEndIndex++;
        }
        // if the breakpoint was inserted inside the current loop, increment the 
        // loopEnd-index only:
        else if( bpIndex <= data->loopEndIndex )
          data->loopEndIndex++;

        // a breakpoint may also be inserted during calls to getSample(), so we 
        // need to keep the indices which are dereferenced there up to date:
        if( bpIndex <= rightIndex )
        {
          leftIndex++;
          rightIndex++;
        }

        // ! leave critical section !
        data->mutex.unlock();


        // nothing more to do here. jump out of the function:
        return bpIndex;
      } // end of  if( bpTime >= newTimeStamp )
      else
        bpIndex++; // lets look at the next breakpoint...

    } // end of while( bpIndex < (int) data->breakpoints.size() )

    return -1; // this command should actually never be executed due to the 
    // logic of the stuff above, but we need it to suppress a 
    // compiler-warning

  } // end of  if( newTimeStamp >= 0.0 )
  else
    return -1; // data->breakpoints with newTimeStamp<0 are illegal
}

bool BreakpointModulator::removeBreakpoint(int index)
{
  // do not allow removal of the first and last breakpoint:
  if( index > 0 && index < lastBreakpointIndex() )
  {
    //markPresetAsDirty();

    // ! enter critical section !
    data->mutex.lock();

    data->breakpoints.erase(data->breakpoints.begin()+index);

    // decrement the loopStart- and loopEnd-index by one, if the new breakpoint
    // was removed before the current loop:
    if( index <= data->loopStartIndex )
    {
      data->loopStartIndex--;
      data->loopEndIndex--;
    }
    // if the breakpoint was removed inside the current loop, decrement the 
    // loopEnd-index only:
    else if( index <= data->loopEndIndex )
      data->loopEndIndex--;

    if( data->loopEndIndex <= data->loopStartIndex )
      data->loopEndIndex++;

    // a breakpoint may also be removed during calls to getSample(), so we 
    // need to keep the indices which are dereferenced there up to date:
    if( index <= rightIndex )
    {
      leftIndex--;
      rightIndex--;
    }

    // ! leave critical section !
    data->mutex.unlock();

    return true;
  }
  else
    return false;
}

bool BreakpointModulator::modifyBreakpoint(int    index, 
                                           double newTimeStamp, 
                                           double newLevel, 
                                           int    newShape, 
                                           double newShapeAmount)
{
  if( index >= 0 && index < (int) data->breakpoints.size()  )   
  {
    //markPresetAsDirty();

    // check if the first and last data->breakpoints are being modified and impose
    // some restrictions on them....
    /*
    if( index == 0 )
    {
      //data->breakpoints[index].timeStamp = 0.0;
      data->breakpoints[index].level     = clipLevelToRange(newLevel);
      if( newShape != 0 )
        data->breakpoints[index].shape = newShape;
      if( newShapeAmount >= 0.01 )
        data->breakpoints[index].shapeAmount = newShapeAmount;

      return true;
    }
    */
    if( index == lastBreakpointIndex() )
    {
      // make sure that the breakpoint does not come too close to its neighbours:
      if( newTimeStamp-data->breakpoints[index-1].timeStamp < data->minBreakpointDistance )
        newTimeStamp = data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;

      data->breakpoints[index].timeStamp = newTimeStamp;

      if( data->endLevelFixedAtZero )
      {
        data->endLevel                 = 0.0;
        data->breakpoints[index].level = 0.0;
      }
      else
      {
        data->endLevel                 = clipLevelToRange(newLevel);
        data->breakpoints[index].level = data->endLevel;
      }

      if( newShape != 0 )
        data->breakpoints[index].shape = newShape;
      if( newShapeAmount >= 0.01 )
        data->breakpoints[index].shapeAmount = newShapeAmount;

      return true;
    }
    else // we deal with an intermediate breakpoint
    {
      // edit the breakpoint without affecting the others:
      if( data->editMode == EDIT_WITHOUT_SHIFT )
      {
        // make sure that the breakpoint does not come too close to its neighbours:
        if( index > 0 )
        {
          if( newTimeStamp-data->breakpoints[index-1].timeStamp < data->minBreakpointDistance )
          {
            newTimeStamp = data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;
          }
        }
        if( data->breakpoints[index+1].timeStamp-newTimeStamp < data->minBreakpointDistance )
          newTimeStamp = data->breakpoints[index+1].timeStamp - data->minBreakpointDistance;

        data->breakpoints[index].timeStamp = newTimeStamp;
        data->breakpoints[index].level     = clipLevelToRange(newLevel);
        if( newShape != 0 )
          data->breakpoints[index].shape = newShape;
        if( newShapeAmount >= 0.01 )
          data->breakpoints[index].shapeAmount = newShapeAmount;
      }
      // edit the breakpoint and shift its successors in time:
      else 
      {
        // make sure that the breakpoint does not come too close to its 
        // predecessor:
        if( index > 0 )
        {
          if( newTimeStamp-data->breakpoints[index-1].timeStamp < data->minBreakpointDistance )
          {
            newTimeStamp = data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;
          }
        }

        // move all the succeeding data->breakpoints in time by the difference between
        // the new and the old value:
        double timeShift = newTimeStamp - data->breakpoints[index].timeStamp;
        for(int bp=index+1; bp <= lastBreakpointIndex(); bp++)
          data->breakpoints[bp].timeStamp += timeShift;

        // modify the breakpoint-data at "index":
        data->breakpoints[index].timeStamp = newTimeStamp;
        data->breakpoints[index].level     = clipLevelToRange(newLevel);
        if( newShape != 0 )
          data->breakpoints[index].shape = newShape;
        if( newShapeAmount >= 0.01 )
          data->breakpoints[index].shapeAmount = newShapeAmount;
      }
      return true;
    } // end of if(index==0), else if(index==lastBreakpointIndex()), else

  } // end of if( index >= 0 && index < (int) data->breakpoints.size()  )  
  else
    return false;
}

void BreakpointModulator::setLoopMode(bool shouldBeLooped)
{
  data->loopIsOn = shouldBeLooped;
}

bool BreakpointModulator::setLoopStartIndex(int newLoopStartIndex)
{
  // make sure, that the new loop start makes sense and then update our 
  // member accordingly:
  if( newLoopStartIndex >=  0 && 
    newLoopStartIndex <  lastBreakpointIndex() &&
    newLoopStartIndex <  data->loopEndIndex )
  {
    data->loopStartIndex = newLoopStartIndex;
    //markPresetAsDirty();
    return true;
  }
  else
    return false;
}

bool BreakpointModulator::setLoopEndIndex(int newLoopEndIndex)
{
  // make sure, that the new loop end makes sense and then update our 
  // member accordingly:
  if( newLoopEndIndex >  0 && 
    newLoopEndIndex <= lastBreakpointIndex() &&
    newLoopEndIndex >  data->loopStartIndex )
  {
    data->loopEndIndex = newLoopEndIndex;
    //markPresetAsDirty();
    return true;
  }
  else
    return false;
}

void BreakpointModulator::setEditMode(int newEditMode)
{
  data->editMode = newEditMode;
}

void BreakpointModulator::setNumCyclesInLoop(int newNumberOfCyclesInLoop)
{
  if( newNumberOfCyclesInLoop >= 1 )
    data->numCyclesInLoop = newNumberOfCyclesInLoop;
  //markPresetAsDirty();
}

void BreakpointModulator::setSyncMode(bool shouldBeSynced)
{
  data->syncMode = shouldBeSynced;
  //markPresetAsDirty();
}

void BreakpointModulator::setBeatsPerMinute(double newBpm)
{
  if( newBpm > 0.0 )
    data->bpm = newBpm;
}

void BreakpointModulator::setTimeScale(double newTimeScale)
{ 
  if( newTimeScale >= 0.0001 )
  {
    data->timeScale  = newTimeScale;
    updateTimeScaleFactor();
    //markPresetAsDirty();
  }
}

void BreakpointModulator::setTimeScaleByKey(double newTimeScaleByKey)
{ 
  data->timeScaleByKey  = newTimeScaleByKey;
  updateTimeScaleFactor();
  //markPresetAsDirty();
}

void BreakpointModulator::setTimeScaleByVel(double newTimeScaleByVel)
{ 
  data->timeScaleByVel  = newTimeScaleByVel;
  updateTimeScaleFactor();
  //markPresetAsDirty();
}

void BreakpointModulator::setDepth(double newDepth)
{ 
  data->depth = newDepth;
  //markPresetAsDirty();
}

void BreakpointModulator::setDepthByKey(double newDepthByKey)
{ 
  data->depthByKey = newDepthByKey;
  //markPresetAsDirty();
}

void BreakpointModulator::setDepthByVel(double newDepthByVel)
{ 
  data->depthByVel = newDepthByVel;
  //markPresetAsDirty();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int BreakpointModulator::getNumBreakpoints() const
{
  return (int) data->breakpoints.size();
}

int BreakpointModulator::lastBreakpointIndex() const
{
  return (int) data->breakpoints.size()-1;
}

double BreakpointModulator::getScaleFactor() const
{
  return data->scaleFactor;
}

double BreakpointModulator::getOffset() const
{
  return data->offset;
}

double BreakpointModulator::getStartTime() const
{
  return data->breakpoints[0].timeStamp;
}

double BreakpointModulator::getEndTime() const
{
  if( data->breakpoints.size() > 1)
    return data->breakpoints[data->breakpoints.size()-1].timeStamp;
  else
    return 0.0;
}

double BreakpointModulator::getMinLevel() const
{
  double min = data->breakpoints[0].level;
  for(int p=0; p < (int) data->breakpoints.size(); p++)
  {
    if( data->breakpoints[p].level < min )
      min = data->breakpoints[p].level;
  }
  return min;
}

double BreakpointModulator::getMaxLevel() const
{
  double max = data->breakpoints[0].level;
  for(int p=0; p < (int) data->breakpoints.size(); p++)
  {
    if( data->breakpoints[p].level > max )
      max = data->breakpoints[p].level;
  }
  return max;
}

double BreakpointModulator::getBreakpointTime(int index) const
{
  if( index < 0 || index > (int) data->breakpoints.size()-1 )
    return -1.0;
  else
    return data->breakpoints[index].timeStamp;
}

double BreakpointModulator::getBreakpointMinTime(int index) const
{
  if( index < 0 || index > (int) data->breakpoints.size()-1 )
    return 0.0;
  else
  {
    if( index == 0 )
      return data->breakpoints[index].timeStamp - 1.0;
      //return 0.0;
    else
      return data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;
  }
}

double BreakpointModulator::getBreakpointMaxTime(int index) const
{
  if( index < 0 || index > (int) data->breakpoints.size()-1 )
    return 0.0;
  else
  {
    //if( index == 0 )
    //  return 0.0;
    if( index == lastBreakpointIndex() || data->editMode == EDIT_WITH_SHIFT )
      return 30.0; // this is somewhat arbitrary
    else
      return data->breakpoints[index+1].timeStamp - data->minBreakpointDistance;
  }
}

bool BreakpointModulator::setBreakpointTime(int index, double newTimeStamp)
{
  // check if the first and last data->breakpoints are being modified and impose
  // some restrictions on them....
  if( index == 0 )
  {
    /*
    if( newTimeStamp == 0.0 )
      return true;
    else
      return false;
    // the first breakpoint is fixed in time at time zero
    */
    if( data->breakpoints[index+1].timeStamp-newTimeStamp < data->minBreakpointDistance )
      newTimeStamp = data->breakpoints[index+1].timeStamp - data->minBreakpointDistance;
    data->breakpoints[index].timeStamp = newTimeStamp;
    return true;
  }
  else if( index == lastBreakpointIndex() )
  {
    if( index<1 ) return false; // should never happen -> would result in an access violation!

    // make sure that the breakpoint does not come too close to its predecessor:
    if( newTimeStamp-data->breakpoints[index-1].timeStamp < data->minBreakpointDistance )
      newTimeStamp = data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;
    data->breakpoints[index].timeStamp = newTimeStamp;
    return true;
  }
  else if( index > 0 && index < lastBreakpointIndex() )
    // we deal with an intermediate breakpoint
  {
    // edit the breakpoint without affecting the others:
    if( data->editMode == EDIT_WITHOUT_SHIFT )
    {
      // make sure that the breakpoint does not come too close to its neighbours:
      if( newTimeStamp-data->breakpoints[index-1].timeStamp < data->minBreakpointDistance )
        newTimeStamp = data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;
      if( data->breakpoints[index+1].timeStamp-newTimeStamp < data->minBreakpointDistance )
        newTimeStamp = data->breakpoints[index+1].timeStamp - data->minBreakpointDistance;
      data->breakpoints[index].timeStamp = newTimeStamp;
      return true;
    }
    // edit the breakpoint and shift its successors in time:
    else 
    {
      // make sure that the breakpoint does not come too close to its 
      // predecessor:
      if( newTimeStamp-data->breakpoints[index-1].timeStamp < data->minBreakpointDistance )
        newTimeStamp = data->breakpoints[index-1].timeStamp + data->minBreakpointDistance;

      // move all the succeeding data->breakpoints in time by the difference between
      // the new and the old value:
      double timeShift = newTimeStamp - data->breakpoints[index].timeStamp;
      for(int bp=index+1; bp <= lastBreakpointIndex(); bp++)
        data->breakpoints[bp].timeStamp += timeShift;

      // modify the breakpoint-data at "index":
      data->breakpoints[index].timeStamp = newTimeStamp;
    }
    return true;
  } // end of if(index==0), else if(index==lastBreakpointIndex()), else
  else
    return false;
}

double BreakpointModulator::getBreakpointLevel(int index) const
{
  if( index < 0 || index > (int) data->breakpoints.size()-1 )
    return 0.0;
  else
    return data->breakpoints[index].level;
}

bool BreakpointModulator::setBreakpointLevel(int index, double newLevel)
{
  if( index == lastBreakpointIndex() )
  {
    if( data->endLevelFixedAtZero )
      data->breakpoints[index].level = 0.0;
    else
      data->breakpoints[index].level = clipLevelToRange(newLevel);
    return true;
  }
  else if( index >= 0 && index <= lastBreakpointIndex() )
  {
    data->breakpoints[index].level = clipLevelToRange(newLevel);
    return true;
  }
  else
    return false;
}

int BreakpointModulator::getBreakpointShape(int index) const
{
  if( index < 0 || index > (int) data->breakpoints.size()-1 )
    return -1;
  else
    return data->breakpoints[index].shape;
}

bool BreakpointModulator::setBreakpointShape(int index, int newShape)
{
  if( index >= 0 && index <= lastBreakpointIndex() )
  {
    if( newShape >= ModBreakpoint::STAIRSTEP &&
      newShape <= ModBreakpoint::SINE_2 )
    {
      data->breakpoints[index].shape = newShape;
      return true;
    }
    else
      return false;
  }
  else
    return false;
}

double BreakpointModulator::getBreakpointShapeAmount(int index) const
{
  if( index < 0 || index > (int) data->breakpoints.size()-1 )
    return 0.0;
  else
  {
    return data->breakpoints[index].shapeAmount;
  }
}

bool BreakpointModulator::setBreakpointShapeAmount(int index, double newShapeAmount)
{
  if( index >= 0 && index <= lastBreakpointIndex() )
  {
    if( newShapeAmount >= 0.01 )
    {
      data->breakpoints[index].shapeAmount = newShapeAmount;
      return true;
    }
    else return false;
  }
  else
    return false;
}

int BreakpointModulator::getLoopMode() const
{
  return data->loopIsOn;
}

int BreakpointModulator::getLoopStartIndex() const
{
  return data->loopStartIndex;
}

int BreakpointModulator::getLoopEndIndex() const
{
  return data->loopEndIndex;
}

int BreakpointModulator::getNumCyclesInLoop() const
{
  return data->numCyclesInLoop;
}

bool BreakpointModulator::isInSyncMode() const
{
  return data->syncMode;
}

double BreakpointModulator::getBeatsPerMinute() const
{
  return data->bpm;
}

double BreakpointModulator::getTimeScale() const
{
  return data->timeScale;
}

double BreakpointModulator::getTimeScaleByKey() const
{
  return data->timeScaleByKey;
}

double BreakpointModulator::getTimeScaleByVel() const
{
  return data->timeScaleByVel;
}

double BreakpointModulator::getDepth() const
{
  return data->depth;
}

double BreakpointModulator::getDepthByKey() const
{
  return data->depthByKey;
}

double BreakpointModulator::getDepthByVel() const
{
  return data->depthByVel;
}

//-------------------------------------------------------------------------------------------------
// event handling:

void BreakpointModulator::noteOn(bool startFromCurrentLevel, int newKey, int newVel)
{
  noteIsOn   = true;
  leftIndex  = 0;
  rightIndex = 1;
  currentKey = newKey;
  currentVel = newVel;
  updateTimeScaleFactor();
  accumulatedTimingError = 0.0;

  /*
  timeScaleFactor  = data->timeScale;
  timeScaleFactor *= pow(2.0, (0.01*data->timeScaleByKey/12.0) * (currentKey-64) );
  timeScaleFactor *= pow(2.0, (0.01*data->timeScaleByVel/63.0) * (currentVel-64) ); 
  */

  data->mutex.lock();

  // set up the 'countdown' variable:
  updateSamplesToNextBreakpoint();

  // get the start-level and calculate the difference to the target-level:
  if( startFromCurrentLevel == true && !endIsReached && data->breakpoints[0].level == 0.0)  
  //if( startFromCurrentLevel == true && !endIsReached )
  {
    leftLevel = previousOut;
  }
  else
  {
    leftLevel = data->breakpoints[leftIndex].level;
    leftLevel = scaleLevelByKeyAndVelocity(leftLevel);
  }
  rightLevel = data->breakpoints[rightIndex].level;
  rightLevel = scaleLevelByKeyAndVelocity(rightLevel);
  levelDelta = rightLevel - leftLevel;

  data->mutex.unlock();

  endIsReached       = false;
  outLevelIsConstant = false;

  // set up the state variables:
  setupStateVariables();
}

void BreakpointModulator::noteOnAndAdvanceTime(int sampleIndexToStartFrom)
{
  noteOn(false);
  for(int i=0; i<sampleIndexToStartFrom; i++)
    getSample(); 
  // crude implementation - runs through all breakpoints which possibly may have been skipped
}

void BreakpointModulator::noteOff(bool startFromCurrentLevel)
{
  noteIsOn   = false;
  leftIndex  = data->loopEndIndex;
  leftIndex  = getLastSimultaneousIndex(leftIndex);
  rightIndex = leftIndex + 1;

  if( leftIndex > lastBreakpointIndex()-1 )
  {
    leftIndex  = lastBreakpointIndex()-1;
    rightIndex = lastBreakpointIndex();
  }

  data->mutex.lock();

  // set up the 'countdown' variable:
  updateSamplesToNextBreakpoint();

  // get the start-level and calculate the difference to the target-level:
  if( startFromCurrentLevel == true )
    leftLevel = previousOut;
  else
  {
    leftLevel = data->breakpoints[leftIndex].level;
    leftLevel = scaleLevelByKeyAndVelocity(leftLevel);
  }
  rightLevel = data->breakpoints[rightIndex].level;
  rightLevel = scaleLevelByKeyAndVelocity(rightLevel);
  levelDelta = rightLevel - leftLevel;

  data->mutex.unlock();


  outLevelIsConstant = false;

  // set up the state variables:
  setupStateVariables();
}

void BreakpointModulator::handleBreakpointArrival()
{
  // aquire mutual exclusive variable access to ensure that the vector is not modified during the 
  // following actions:
  data->mutex.lock();

  // increment the breakpoint indices to the next pair of breakpoints:
  //leftIndex  = getNextNonSimultaneousIndex(leftIndex);
  leftIndex++;
  leftIndex  = getLastSimultaneousIndex(leftIndex);
  rightIndex = leftIndex + 1;
  //leftIndex++;
  //rightIndex++;

  // if we have incremented our right index beyond the end of the breakpoint-array, we have reached
  // the end:
  if( leftIndex >= lastBreakpointIndex() )
  {
    // remain looping when the loop end happens to be the very last breakpoint:
    if( lastBreakpointIndex() == getLoopEndIndex() && noteIsOn && data->loopIsOn )
    {
      wrapAroundIndicesForLoop();
    }
    else
    {
      rightIndex         = leftIndex;
      leftLevel          = scaleLevelByKeyAndVelocity(data->breakpoints[leftIndex].level);
      rightLevel         = leftLevel;
      levelDelta         = 0.0;
      previousOut        = leftLevel;
      endIsReached       = true;
      outLevelIsConstant = true;
      data->mutex.unlock();
      return;
    }
  }
  else if( leftIndex >= data->loopEndIndex && noteIsOn && data->loopIsOn )
  {
    wrapAroundIndicesForLoop();
  }
  
  // these wraparounds may result in leftIndex==rightIndex -> this condition indicates a loop of
  // length zero - in this case we just pop out the constant value at the left/right index:
  if( leftIndex == rightIndex )
  {
    leftLevel          = scaleLevelByKeyAndVelocity(data->breakpoints[leftIndex].level);
    rightLevel         = leftLevel;
    levelDelta         = 0.0;
    previousOut        = leftLevel;
    outLevelIsConstant = true;
    data->mutex.unlock();
    return;
  }

  // set up the 'countdown' variable:
  updateSamplesToNextBreakpoint();

  // get the current level on the left side and its difference to the level 
  // on the right side:
  leftLevel  = data->breakpoints[leftIndex].level;
  leftLevel  = scaleLevelByKeyAndVelocity(leftLevel);
  rightLevel = data->breakpoints[rightIndex].level;
  rightLevel = scaleLevelByKeyAndVelocity(rightLevel);
  levelDelta = rightLevel - leftLevel;

  // set up the internal state variables for the recursive formulas:
  setupStateVariables();

  // the output will not be constant when we made to this point in the function (without leaving 
  // it prematurely):
  outLevelIsConstant = false;

  // release the mutex-lock:
  data->mutex.unlock();
}

void BreakpointModulator::wrapAroundIndicesForLoop()
{
  data->mutex.lock();

  if( leftIndex >= data->loopEndIndex )
    leftIndex  = data->loopStartIndex;

  if( data->loopEndIndex > data->loopStartIndex )
    rightIndex = leftIndex+1;
  else 
    rightIndex = leftIndex; // loop of zero length

  data->mutex.unlock();
}

int BreakpointModulator::getNextNonSimultaneousIndex(int startIndex)
{
  //data->mutex.lock();
  // mutex is not necesarry - this function is called only internally from functions which 
  // already aquire the lock

  // when there are no breakpoints at the same time instant, this will be the default value:
  int foundIndex = startIndex + 1;

  // make sure to not access out-of-range indices: 
  if( foundIndex > lastBreakpointIndex() )
    return 0;

  while( data->breakpoints[startIndex].timeStamp >= data->breakpoints[foundIndex].timeStamp )
  {
    // timeStamp at foundIndex is not strictly larger than at startIndex - skip to next breakpoint:
    foundIndex++;

    // make sure to not access out-of-range indices: 
    if( foundIndex > lastBreakpointIndex() )
      return 0;
  }

  // O.K. we now have skipped to the first index for which the timeStamp is actually larger than 
  // at the startIndex (and returned 0, if that would be out-of-range) - now we can return our 
  // finding:
  return foundIndex;
}

int BreakpointModulator::getLastSimultaneousIndex(int index)
{
  // make sure to not access out-of-range indices: 
  if( index+1 > lastBreakpointIndex() )
    return lastBreakpointIndex();

  while( data->breakpoints[index].timeStamp >= data->breakpoints[index+1].timeStamp )
  {
    // timeStamp at foundIndex is not strictly larger than at startIndex - skip to next breakpoint:
    index++;

    // make sure to not access out-of-range indices: 
    if( index+1 > lastBreakpointIndex() )
      return lastBreakpointIndex();
  }

  return index;
}

void BreakpointModulator::updateSamplesToNextBreakpoint()
{
  // get the length of the envelope-segment to be generated (in seconds):
  double timeDelta = data->breakpoints[rightIndex].timeStamp 
    - data->breakpoints[leftIndex].timeStamp;
  if( data->syncMode == true )
    timeDelta = beatsToSeconds(timeDelta, data->bpm);

  // convert to samples and use this length as initial value for our countdown-variable:
  double samplesExact     = timeScaleFactor * timeDelta * data->sampleRate;
  samplesToNextBreakpoint = roundToInt(samplesExact);
  double error            = (double) samplesToNextBreakpoint - samplesExact;
  accumulatedTimingError += error;
  if( accumulatedTimingError > 0.5 )
  {
    samplesToNextBreakpoint -= 1;
    accumulatedTimingError  -= 1.0;
  }
  else if( accumulatedTimingError < -0.5 )
  {
    samplesToNextBreakpoint += 1;
    accumulatedTimingError  += 1.0;
  }
}

void BreakpointModulator::setupStateVariables()
{
  // do the specific initializations for the different envelope shapes (for 
  // details about what's going on, refer to comments in the MatLab
  // implementation):
  int segmentLength = samplesToNextBreakpoint;

  if( segmentLength < 1 )
  {
    // TODO: catch this special condition (can occur when 
    // loopEndIndex == loopStartIndex which is allowed).....
    return;
  }

  currentShape = data->breakpoints[rightIndex].shape;
  switch( currentShape )
  {
  case ModBreakpoint::STAIRSTEP:
    {
      state1        = leftLevel;
      //state1_change = 0.0;
    }
    break;
  case ModBreakpoint::LINEAR:
    {
      state1        = leftLevel;
      state1_change = levelDelta / (double) segmentLength;
    }
    break;
  case ModBreakpoint::SMOOTH:
    {
      double omega  = PI / (double) segmentLength;
      state1_change = 2.0*cos(omega);
      state1        = sin( -(0.5*PI) - omega );
      state2        = sin( -(0.5*PI) - 2.0*omega );
    }
    break;
  case ModBreakpoint::ANALOG:
    {
      state1_min    = pow(0.01, data->breakpoints[rightIndex].shapeAmount);
      state1_max    = pow(state1_min, 1.0 / (double) (segmentLength+1));
      scaler1       = levelDelta / (state1_max-state1_min);
      state1        = state1_max;
      state1_change = state1_max;
    }
    break;
  case ModBreakpoint::GROWING:
    {
      state1_min    = pow(0.01, data->breakpoints[rightIndex].shapeAmount);
      state1_max    = pow(state1_min, 1.0 / (double) (segmentLength+1));
      scaler1       = levelDelta / (state1_max-state1_min);
      state1        = state1_min;
      state1_change = 1.0/state1_max;
    }
    break;
  case ModBreakpoint::SIGMOID:
    {
      state1_min    = pow(0.01, data->breakpoints[rightIndex].shapeAmount);
      state1_max    = pow(state1_min, 1.0 / (double) (segmentLength+1));
      scaler1       = levelDelta / (state1_max-state1_min);
      state1        = state1_max;
      state1_change = state1_max;
      //state2_min    = state1_min; // these two variables are actually not used
      //state2_max    = state1_max; // therefore their assignment is omitted
      state2        = state1_min;
      state2_change = 1.0 / state1_max;
    }
    break;
  case ModBreakpoint::SPIKEY:
    {
      state1_min    = pow(0.01, data->breakpoints[rightIndex].shapeAmount);
      state1_max    = pow(state1_min, 1.0 / (double) (segmentLength+1));
      scaler1       = 1.0 / (state1_max-state1_min);
      state1        = state1_max;
      state1_change = state1_max;

      //state2_min    = state1_min; // these two variables are actually not used
      //state2_max    = state1_max; // therefore their assignment is omitted
      state2        = state1_min;
      state2_change = 1.0 / state1_max;
    }
    break;
  case ModBreakpoint::SINE_1:
    {
      double omega  = 0.5*PI / (double) segmentLength;
      state1_change = 2.0*cos(omega);
      state1        = sin( -(0.0*PI) - omega );
      state2        = sin( -(0.0*PI) - 2.0*omega );
    }
    break;
  case ModBreakpoint::SINE_2:
    {
      double omega  = 0.5*PI / (double) segmentLength;
      state1_change = 2.0*cos(omega);
      state1        = sin( -(0.5*PI) - omega );
      state2        = sin( -(0.5*PI) - 2.0*omega );
    }
    break;
  default:  // linear by default
    {
      state1        = leftLevel;
      state1_change = levelDelta / (double) segmentLength;
    }

  } // end of switch( shape )

  // a very unelegant method to skip the samples which are to be skipped:
  /*
  double dummy;
  for(int i = 0; i < sampleIndexToStart; i++)
    dummy = getSample();
  */
}

//-------------------------------------------------------------------------------------------------
// initialization:

void BreakpointModulator::setToDefaultValues()
{
  data->editMode            =  EDIT_WITHOUT_SHIFT;
  //loopMode                  =  NO_LOOP;
  samplesToNextBreakpoint   =  0;
  accumulatedTimingError    =  0.0;
  leftIndex                 = -1;
  rightIndex                =  0;
  data->loopStartIndex      =  2;
  data->loopEndIndex        =  3;
  data->loopIsOn            =  true;
  endIsReached              =  true;
  outLevelIsConstant        =  true;
  leftLevel                 =  0.0;
  rightLevel                =  0.0;
  levelDelta                =  0.0;
  timeScaleFactor           =  1.0;

  state1                    = 0.0;
  state2                    = 0.0;
  state1_change             = 0.0;
  state2_change             = 0.0;
  state1_min                = 0.0;
  state2_min                = 0.0;
  state1_max                = 0.0;
  state2_max                = 0.0;
  scaler1                   = 1.0;
  scaler2                   = 1.0;

  previousOut               = 0.0;

  // ! enter critical section !
  data->mutex.lock();

  // initialize the breakpoint-vector with two entries, these two will always 
  // be there (their data can be modified, though), additional entries can be 
  // inserted and removed at will in between:
  data->breakpoints.clear();
  ModBreakpoint newBreakpoint;

  newBreakpoint.timeStamp   = 0.0;
  newBreakpoint.level       = 0.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  newBreakpoint.timeStamp   = 0.5;
  newBreakpoint.level       = 1.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  newBreakpoint.timeStamp   = 1.0;
  newBreakpoint.level       = 0.5;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  newBreakpoint.timeStamp   = 2.0;
  newBreakpoint.level       = 0.5;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  newBreakpoint.timeStamp   = 3.0;
  newBreakpoint.level       = 0.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  // ! leave critical section !
  data->mutex.unlock();
}

void BreakpointModulator::initialize()
{
  data->editMode            =  EDIT_WITHOUT_SHIFT;
  samplesToNextBreakpoint   =  0;
  accumulatedTimingError    =  0.0;
  leftIndex                 = -1;
  rightIndex                =  0;
  data->loopStartIndex      =  0;
  data->loopEndIndex        =  1;
  data->loopIsOn            =  false;
  endIsReached              =  true;
  outLevelIsConstant        =  true;
  leftLevel                 =  0.0;
  rightLevel                =  0.0;
  levelDelta                =  0.0;
  timeScaleFactor           =  1.0;
  currentKey                =  0;
  currentVel                =  0;

  state1                    = 0.0;
  state2                    = 0.0;
  state1_change             = 0.0;
  state2_change             = 0.0;
  state1_min                = 0.0;
  state2_min                = 0.0;
  state1_max                = 0.0;
  state2_max                = 0.0;
  scaler1                   = 1.0;
  scaler2                   = 1.0;

  previousOut               = 0.0;

  // ! enter critical section !
  data->mutex.lock();

  // initialize the breakpoint-vector with two entries, these two will always 
  // be there (their data can be modified, though), additional entries can be 
  // inserted and removed at will in between:
  data->breakpoints.clear();
  ModBreakpoint newBreakpoint;

  newBreakpoint.timeStamp   = 0.0;
  newBreakpoint.level       = 1.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  newBreakpoint.timeStamp   = 1.0;
  newBreakpoint.level       = 1.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);


  // ! leave critical section !
  data->mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// audio processing:

void BreakpointModulator::fillBufferWithEnvelope(double *buffer, int length, bool useOnlyLoop)
{
  // retrieve a couple of parameters that must be modified for later restore:
  double tmpSampleRate = data->sampleRate;
  bool   tmpSync       = data->syncMode;
  bool   tmpLoop       = data->loopIsOn;

  // fill the buffer:
  if( useOnlyLoop == false )
  {
    double lengthInSeconds = getEndTime() - getStartTime();
    data->sampleRate = (double) length / lengthInSeconds;
    data->syncMode   = false;
    data->loopIsOn   = false;
    noteOn();
    for(int n=0; n<length; n++)
      buffer[n] = getSample();
  }
  else
  {
    double end             = getBreakpointTime(getLoopEndIndex());
    double start           = getBreakpointTime(getLoopStartIndex());
    double lengthInSeconds = end - start;
    data->sampleRate = (double) length / lengthInSeconds;
    data->syncMode   = false;
    data->loopIsOn   = false;
    double startSample = start*data->sampleRate; // for test only
    noteOnAndAdvanceTime(roundToInt(start*data->sampleRate));
    for(int n=0; n<length; n++)
      buffer[n] = getSample();
  }

  // restore the original parameters:
  data->sampleRate = tmpSampleRate;
  data->syncMode   = tmpSync;
  data->loopIsOn   = tmpLoop;
}

//-------------------------------------------------------------------------------------------------
// internal helper functions:

void BreakpointModulator::updateTimeScaleFactor()
{
  timeScaleFactor  = data->timeScale;  
  timeScaleFactor *= pow(2.0, (0.01*data->timeScaleByKey/12.0) * (currentKey-64) );
  timeScaleFactor *= pow(2.0, (0.01*data->timeScaleByVel/63.0) * (currentVel-64) ); 
      
  for(int s=0; s < (int) slaves.size(); s++)  
    slaves[s]->updateTimeScaleFactor();
}

double BreakpointModulator::clipLevelToRange(double inLevel)
{
  double outLevel = inLevel;

  if( outLevel < data->minimumAllowedLevel )
    outLevel = data->minimumAllowedLevel;
  if( outLevel > data->maximumAllowedLevel )
    outLevel = data->maximumAllowedLevel;

  return outLevel;
}

double BreakpointModulator::scaleLevelByKeyAndVelocity(double unscaledLevel)
{
  double scaledLevel;
  scaledLevel = powBipolar(unscaledLevel, data->depth);  
  scaledLevel = powBipolar(scaledLevel,   pow(2.0, (0.01*data->depthByKey/12.0)*(currentKey-64)) );  
  scaledLevel = powBipolar(scaledLevel,   pow(2.0, (0.01*data->depthByVel/63.0)*(currentVel-64)) );
  return scaledLevel;
}

//-------------------------------------------------------------------------------------------------
// master/slave config:

void BreakpointModulator::addSlave(BreakpointModulator* newSlave)
{
  // add the new slave to the vector of slaves:
  slaves.push_back(newSlave);

  // delete the original parameter-set of the new slave and redirect it to ours (with some safety 
  // checks):
  if( newSlave->data != NULL && newSlave->data != this->data )
  {
    newSlave->data->mutex.lock();
    newSlave->data->breakpoints.clear();
    newSlave->data->mutex.unlock();  
    delete newSlave->data;
    newSlave->data = this->data;
  }
  else
  {
    DEBUG_BREAK; 
    // the object to be added as slave did not contain a valid parameter-pointer - maybe it has 
    // been already added as slave to another master?
  }

  // set the isMaster-flag of the new slave to false: 
  newSlave->isMaster = false;

  // this flag will prevent the destructor of the slave from trying to delete the parameter-set 
  // which is now shared - only masters delete their parameter-set on destruction
}