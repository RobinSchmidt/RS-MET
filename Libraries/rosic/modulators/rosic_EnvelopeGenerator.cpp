#include "rosic_EnvelopeGenerator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

EnvelopeGenerator::EnvelopeGenerator()
{
  data->breakpoints.clear();
  ModBreakpoint newBreakpoint;

  // start:
  newBreakpoint.timeStamp   = 0.0;
  newBreakpoint.level       = 0.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  // peak:
  newBreakpoint.timeStamp   = 0.25;
  newBreakpoint.level       = 1.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  // sustain-loop start:
  newBreakpoint.timeStamp   = 1.0;
  newBreakpoint.level       = 0.5;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  // sustain-loop end:
  newBreakpoint.timeStamp   = 2.0;
  newBreakpoint.level       = 0.5;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  // end:
  newBreakpoint.timeStamp   = 3.0;
  newBreakpoint.level       = 0.0;
  newBreakpoint.shape       = ModBreakpoint::ANALOG;
  newBreakpoint.shapeAmount = 1.0;
  data->breakpoints.push_back(newBreakpoint);

  data->loopStartIndex =  2;
  data->loopEndIndex   =  3;
  data->loopIsOn       =  true;
}


