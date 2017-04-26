#include "BreakpointModulator.h"

//----------------------------------------------------------------------------
// construction/destruction:
BreakpointModulator::BreakpointModulator()
{
	sampleRate = 44100.0;
 initialize();
}

BreakpointModulator::~BreakpointModulator()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void BreakpointModulator::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;
}

int BreakpointModulator::insertBreakpoint(double newTimeStamp, double newLevel, 
                                           int newShape, double newShapeAmount)
{
 if( newTimeStamp >= 0.0 )
 {
  // Loop through the existing breakpoints and accumulate all the time deltas
  // to obtain an absolute time-value. As soon as this absolute time get 
  // larger than the time-stamp of the new breakpoint to be inserted, we know
  // where we need to insert the new breakpoint:
  int    bpIndex = 1;  // omit the timeDelta of the first breapoint(index 0)
  double bpTime  = 0.0;
  while( bpIndex < (int) breakpoints.size() )
  {
   bpTime = breakpoints[bpIndex].timeStamp;
   if( bpTime >= newTimeStamp )
   {
    // make sure that the new  breakpoint to be inserted is not too close to 
    // another existing breakpoint:
    if( newTimeStamp-breakpoints[bpIndex-1].timeStamp < minBreakpointDistance 
        || bpTime-newTimeStamp < minBreakpointDistance )
    {
     return -1;
    }

    // create a new breakpoint for insertion right before the breakpoint with 
    // the current index:
    ModBreakpoint newBreakpoint;
    newBreakpoint.timeStamp   = newTimeStamp;
    newBreakpoint.level       = newLevel;

    // copy the shape and shapeAmount from the succeeding breakpoint into,
    // the new one if no shape and shapeAmount were specified (i.e. when they 
    // are equal to their default-arguments which are zero), otherwise use the 
    // passed arguments:
    if( newShape == 0 )
     newBreakpoint.shape = breakpoints[bpIndex].shape;
    else
     newBreakpoint.shape = newShape;
    if( newShapeAmount == 0.0 )
     newBreakpoint.shapeAmount = breakpoints[bpIndex].shapeAmount;
    else
     newBreakpoint.shapeAmount = newShapeAmount;

    // insert the new breakpoint right before the breakpoint with the
    // current index:
    breakpoints.insert(breakpoints.begin()+bpIndex, 1, newBreakpoint);

    // advance the loopStart- and loopEnd-index by one, if the new breakpoint
    // was inserted before the current loop:
    if( bpIndex <= loopStartIndex )
    {
     loopStartIndex++;
     loopEndIndex++;
    }
    // if the breakpoint was inserted inside the curren loop, increment the 
    // loopEnd-index only:
    else if( bpIndex <= loopEndIndex )
     loopEndIndex++;

    // nothing more to do here. jump out of the function:
    return bpIndex;
   } // end of  if( bpTime >= newTimeStamp )
   else
    bpIndex++; // lets look at the next breakpoint...

  } // end of while( bpIndex < (int) breakpoints.size() )

  return -1; // this command should actually never be executed due to the 
             // logic of the stuff above, but we need it to suppress a 
             // compiler-warning

 } // end of  if( newTimeStamp >= 0.0 )
 else
  return -1; // breakpoints with newTimeStamp<0 are illegal
}

bool BreakpointModulator::modifyBreakpoint(int index, double newTimeStamp, 
                                           double newLevel, int newShape, 
                                           double newShapeAmount)
{
 if( index >= 0 && index < (int) breakpoints.size()  )   
 {
  // check if the first and last breakpoints are being modified and impose
  // some restrictions on them....
  if( index == 0 )
  {
    breakpoints[index].timeStamp = 0.0;
    breakpoints[index].level     = newLevel;
    if( newShape != 0 )
     breakpoints[index].shape = newShape;
    if( newShapeAmount >= 0.01 )
     breakpoints[index].shapeAmount = newShapeAmount;

   return true;
  }
  else if( index == lastBreakpointIndex() )
  {
    // make sure that the breakpoint does not come too close to its neighbours:
    if( newTimeStamp-breakpoints[index-1].timeStamp < minBreakpointDistance )
     newTimeStamp = breakpoints[index-1].timeStamp + minBreakpointDistance;

    breakpoints[index].timeStamp = newTimeStamp;
    breakpoints[index].level     = 0.0;
    if( newShape != 0 )
     breakpoints[index].shape = newShape;
    if( newShapeAmount >= 0.01 )
     breakpoints[index].shapeAmount = newShapeAmount;

   return true;
  }
  else // we deal with an intermediate breakpoint
  {
   // edit the breakpoint without affecting the others:
   if( editMode == EDIT_WITHOUT_SHIFT )
   {
    // make sure that the breakpoint does not come too close to its neighbours:
    if( newTimeStamp-breakpoints[index-1].timeStamp < minBreakpointDistance )
     newTimeStamp = breakpoints[index-1].timeStamp + minBreakpointDistance;
    if( breakpoints[index+1].timeStamp-newTimeStamp < minBreakpointDistance )
     newTimeStamp = breakpoints[index+1].timeStamp - minBreakpointDistance;

    breakpoints[index].timeStamp = newTimeStamp;
    breakpoints[index].level     = newLevel;
    if( newShape != 0 )
     breakpoints[index].shape = newShape;
    if( newShapeAmount >= 0.01 )
     breakpoints[index].shapeAmount = newShapeAmount;
   }
   // edit the breakpoint and shift its successors in time:
   else 
   {
    // make sure that the breakpoint does not come too close to its 
    // predecessor:
    if( newTimeStamp-breakpoints[index-1].timeStamp < minBreakpointDistance )
     newTimeStamp = breakpoints[index-1].timeStamp + minBreakpointDistance;

    // move all the succeeding breakpoints in time by the difference between
    // the new and the old value:
    double timeShift = newTimeStamp - breakpoints[index].timeStamp;
    for(int bp=index+1; bp <= lastBreakpointIndex(); bp++)
     breakpoints[bp].timeStamp += timeShift;

    // modify the breakpoint-data at "index":
    breakpoints[index].timeStamp = newTimeStamp;
    breakpoints[index].level     = newLevel;
    if( newShape != 0 )
     breakpoints[index].shape = newShape;
    if( newShapeAmount >= 0.01 )
     breakpoints[index].shapeAmount = newShapeAmount;
   }
   return true;
  } // end of if(index==0), else if(index==lastBreakpointIndex()), else

 } // end of if( index >= 0 && index < (int) breakpoints.size()  )  
 else
  return false;
}

bool BreakpointModulator::removeBreakpoint(int index)
{
 // do not allow removal of the first and last breakpoint:
 if( index > 0 && index < lastBreakpointIndex() )
 {
  breakpoints.erase(breakpoints.begin()+index);

  // decrement the loopStart- and loopEnd-index by one, if the new breakpoint
  // was inserted before the current loop:
  if( index <= loopStartIndex )
  {
   loopStartIndex--;
   loopEndIndex--;
  }
  // if the breakpoint was inserted inside the curren loop, increment the 
  // loopEnd-index only:
  else if( index <= loopEndIndex )
   loopEndIndex--;

  return true;
 }
 else
  return false;
}

void BreakpointModulator::setLoopMode(bool shouldBeLooped)
{
 loopIsOn = shouldBeLooped;
}

bool BreakpointModulator::setLoopStartIndex(int newLoopStartIndex)
{
 // make sure, that the new loop start makes senese and then update our 
 // member accordingly:
 if( newLoopStartIndex >  0 && 
     newLoopStartIndex <  (int) breakpoints.size()-1 &&
     newLoopStartIndex <= loopEndIndex )
 {
  loopStartIndex = newLoopStartIndex;
  return true;
 }
 else
  return false;
}

bool BreakpointModulator::setLoopEndIndex(int newLoopEndIndex)
{
 // make sure, that the new loop end makes senese and then update our 
 // member accordingly:
 if( newLoopEndIndex >  0 && 
     newLoopEndIndex <  (int) breakpoints.size()-1 &&
     newLoopEndIndex >= loopStartIndex )
 {
  loopEndIndex = newLoopEndIndex;
  return true;
 }
 else
  return false;
}

//----------------------------------------------------------------------------
// inquiry:

int BreakpointModulator::lastBreakpointIndex()
{
 return breakpoints.size()-1;
}

double BreakpointModulator::getStartTime()
{
 return breakpoints[0].timeStamp;
}

double BreakpointModulator::getEndTime()
{
 if( breakpoints.size() > 1)
  return breakpoints[breakpoints.size()-1].timeStamp;
 else
  return 0.0;
}

double BreakpointModulator::getMinLevel()
{
 double min = breakpoints[0].level;
 for(int p=0; p < (int) breakpoints.size(); p++)
 {
  if( breakpoints[p].level < min )
   min = breakpoints[p].level;
 }
 return min;
}

double BreakpointModulator::getMaxLevel()
{
 double max = breakpoints[0].level;
 for(int p=0; p < (int) breakpoints.size(); p++)
 {
  if( breakpoints[p].level > max )
   max = breakpoints[p].level;
 }
 return max;
}

double BreakpointModulator::getBreakpointTime(int index)
{
 if( index < 0 || index > (int) breakpoints.size()-1 )
  return -1.0;
 else
  return breakpoints[index].timeStamp;
}

double BreakpointModulator::getBreakpointLevel(int index)
{
 if( index < 0 || index > (int) breakpoints.size()-1 )
  return 0.0;
 else
  return breakpoints[index].level;
}

int BreakpointModulator::getLoopMode()
{
 return loopIsOn;
}

int BreakpointModulator::getLoopStartIndex()
{
 return loopStartIndex;
}

int BreakpointModulator::getLoopEndIndex()
{
 return loopEndIndex;
}


//----------------------------------------------------------------------------
// others:

void BreakpointModulator::noteOn(bool startFromCurrentLevel)
{
 noteIsOn     = true;
 endIsReached = false;
 leftIndex    = 0;
 rightIndex   = 1;

 // get the length of the envelope-segment to be generated (in samples):
 double timeDelta  =   breakpoints[rightIndex].timeStamp 
                     - breakpoints[leftIndex].timeStamp;
 int segmentLength = MoreMath::roundToInt(timeScale * timeDelta * sampleRate);

 // use this length as initial value for our countdown-variable:
 samplesToNextBreakpoint = segmentLength;

 // get the start-level and calculate the difference to the target-level:
 if( startFromCurrentLevel == true )
  leftLevel = previousOut;
 else
  leftLevel = breakpoints[leftIndex].level;
 levelDelta = breakpoints[rightIndex].level - leftLevel;

 // set up the state variables:
 setupStateVariables();
}

void BreakpointModulator::noteOff(bool startFromCurrentLevel)
{
 noteIsOn   = false;
 leftIndex  = loopEndIndex;
 rightIndex = leftIndex + 1;

 // get the length of the envelope-segment to be generated (in samples):
 double timeDelta  =   breakpoints[rightIndex].timeStamp 
                     - breakpoints[leftIndex].timeStamp;
 int segmentLength = MoreMath::roundToInt(timeScale * timeDelta * sampleRate);

 // use this length as initial value for our countdown-variable:
 samplesToNextBreakpoint = segmentLength;

 // get the start-level and calculate the difference to the target-level:
 if( startFromCurrentLevel == true )
  leftLevel = previousOut;
 else
  leftLevel = breakpoints[leftIndex].level;
 levelDelta = breakpoints[rightIndex].level - leftLevel;

 // set up the state variables:
 setupStateVariables();
}

void BreakpointModulator::setupStateVariables()
{
 // do the specific initializations for the different envelope shapes (for 
 // details about what's going on, refer to comments in the MatLab
 // implementation):
 int segmentLength = samplesToNextBreakpoint;
 switch( breakpoints[rightIndex].shape )
 {
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
   state1_min    = pow(0.01, breakpoints[rightIndex].shapeAmount);
   state1_max    = pow(state1_min, 1.0 / (double) (segmentLength+1));
   scaler1       = levelDelta / (state1_max-state1_min);
   state1        = state1_max;
   state1_change = state1_max;
  }
  break;
 case ModBreakpoint::GROWING:
  {
   state1_min    = pow(0.01, breakpoints[rightIndex].shapeAmount);
   state1_max    = pow(state1_min, 1.0 / (double) (segmentLength+1));
   scaler1       = levelDelta / (state1_max-state1_min);
   state1        = state1_min;
   state1_change = 1.0/state1_max;
  }
  break;
 case ModBreakpoint::SIGMOID:
  {
   state1_min    = pow(0.01, breakpoints[rightIndex].shapeAmount);
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
   state1_min    = pow(0.01, breakpoints[rightIndex].shapeAmount);
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
 } // end of switch( shape )
}

void BreakpointModulator::initialize()
{
 editMode                  =  EDIT_WITHOUT_SHIFT;
 //loopMode                  =  NO_LOOP;
 minBreakpointDistance     =  0.01;
 samplesToNextBreakpoint   =  0;
 leftIndex                 = -1;
 rightIndex                =  0;
 loopStartIndex            =  2;
 loopEndIndex              =  3;
 loopIsOn                  =  false;
 endIsReached              =  false;
 leftLevel                 =  0.0;
 levelDelta                =  0.0;
 timeScale                 =  1.0;

 state1                    = 0.0;
 state2                    = 0.0;
 state1_change             = 0.0;
 state2_change             = 0.0;
 state1_min                = 0.0;
 state2_min                = 0.0;
 state1_max                = 0.0;
 state2_max                = 0.0;

 previousOut               = 0.0;

 // initialize the breakpoint-vector with two entries, these two will always 
 // be there (their data can be modified, though), additional entries can be 
 // inserted and removed at will in between:
 breakpoints.clear();
 ModBreakpoint newBreakpoint;

 newBreakpoint.timeStamp   = 0.0;
 newBreakpoint.level       = 0.0;
 newBreakpoint.shape       = 3;
 newBreakpoint.shapeAmount = 1.0;
 breakpoints.push_back(newBreakpoint);

 newBreakpoint.timeStamp   = 0.5;
 newBreakpoint.level       = 1.0;
 newBreakpoint.shape       = 3;
 newBreakpoint.shapeAmount = 1.0;
 breakpoints.push_back(newBreakpoint);

 newBreakpoint.timeStamp   = 1.0;
 newBreakpoint.level       = 0.5;
 newBreakpoint.shape       = 3;
 newBreakpoint.shapeAmount = 1.0;
 breakpoints.push_back(newBreakpoint);

 newBreakpoint.timeStamp   = 2.0;
 newBreakpoint.level       = 0.5;
 newBreakpoint.shape       = 3;
 newBreakpoint.shapeAmount = 1.0;
 breakpoints.push_back(newBreakpoint);

 newBreakpoint.timeStamp   = 3.0;
 newBreakpoint.level       = 0.0;
 newBreakpoint.shape       = 2;
 newBreakpoint.shapeAmount = 1.0;
 breakpoints.push_back(newBreakpoint);
}


