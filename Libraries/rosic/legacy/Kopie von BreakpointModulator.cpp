#include "BreakpointModulator.h"

//----------------------------------------------------------------------------
// construction/destruction:
BreakpointModulator::BreakpointModulator()
{
	sampleRate                =  44100.0;
 numBreakpoints            =  8;
 countdownToNextBreakpoint =  0;
 leftIndex                 = -1;
 rightIndex                =  0;
 loopStartIndex            =  0;
 loopEndIndex              =  0;
 leftLevel                 =  0.0;
 levelDelta                =  0.0;

 state1                    = 0.0;
 state2                    = 0.0;
 state1_change             = 0.0;
 state2_change             = 0.0;
 state1_min                = 0.0;
 state2_min                = 0.0;
 state1_max                = 0.0;
 state2_max                = 0.0;
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

void BreakpointModulator::insertBreakpoint(double newTimeStamp, double newLevel, 
                                           int newShape, double newShapeAmount)
{
 if( newTimeStamp >= 0.0 )
 {
  // create a new breakpoint:
  ModBreakpoint newBreakpoint;
  newBreakpoint.timeStamp   = newTimeStamp;
  newBreakpoint.level       = newLevel;
  newBreakpoint.shape       = newShape;
  newBreakpoint.shapeAmount = newShapeAmount;

  // append it to out breakpoints-vector:
  breakpoints.push_back(newBreakpoint);

  // sort the breakpoints-vector:

 }


}

void BreakpointModulator::removeBreakpoint(int index)
{

}

void BreakpointModulator::modifyBreakpoint(int index, double newTimeStamp, 
                                           double newLevel, int newShape, 
                                           double newShapeAmount)
{
 if( index < (int) breakpoints.size() && newTimeStamp >= 0.0 )
 {
  breakpoints[index].timeStamp    = newTimeStamp;
  breakpoints[index].level        = newLevel;
  breakpoints[index].shape        = newShape;
  if( newShapeAmount >= 0.01 )
   breakpoints[index].shapeAmount = newShapeAmount;
 }
}

//----------------------------------------------------------------------------
// others:

void BreakpointModulator::reset()
{
 countdownToNextBreakpoint =  0;
 leftIndex                 = -1;
 rightIndex                =  0;
 leftLevel                 =  0.0;
 levelDelta                =  0.0;

 state1                    =  0.0;
 state2                    =  0.0;
 state1_change             =  0.0;
 state2_change             =  0.0;
 state1_min                =  0.0;
 state2_min                =  0.0;
 state1_max                =  0.0;
 state2_max                =  0.0;
}

void BreakpointModulator::trigger(bool startFromCurrentValue)
{
 countdownToNextBreakpoint =  0;
 leftIndex                 = -1;
 rightIndex                =  0;
 leftLevel                 =  0.0;
 levelDelta                =  0.0;
}

void BreakpointModulator::noteOff()
{

}

bool BreakpointModulator::endIsReached()
{
 return false;
}