#ifndef BreakpointModulator_h
#define BreakpointModulator_h

#include "Definitions.h"
#include "MoreMath.h"

#include <vector>
using namespace std; // for the Breakpoint-sequence

/**

This is a class which generates a modulation trajectory based on a set of
breakpoints......

*/

// a small helper class, combining all the required data for one breakpoint 
// into one class:
class ModBreakpoint
{
public:

 /** This is an enumeration of the available breakpoint shapes. */
 enum shapes
 {
  LINEAR = 1,
  SMOOTH,
  ANALOG,
  EXPLOSIVE,
  SIGMOID,
  SPIKEY
 };

 // member variables:
 //int    index;    
 double timeStamp;
 double level;
 int    shape;
 double shapeAmount;

 /** Standard-Contructor. */
 ModBreakpoint()
 {
  //index       = -1;   // the index is invalid after construction
  timeStamp   = 0.0;
  level       = 0.0;
  shape       = 1;
  shapeAmount = 1.0;
 }

 ~ModBreakpoint() {}; ///< Destructor.
};


//============================================================================
// the actual BreakpointModulator-class:

class BreakpointModulator
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 BreakpointModulator();  ///< Constructor.
 ~BreakpointModulator();  ///< Destructor.


 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Sets the sample-rate for this module.

 void insertBreakpoint(double newTimeStamp, double newLevel, int newShape,
                       double newShapeAmount);
 ///< Inserts a new breakpoint into the list.

 void removeBreakpoint(int index);
 ///< Removes a breakpoint from the list.

 void modifyBreakpoint(int index, double newTimeStamp, double newLevel, 
                       int newShape, double newShapeAmount);
 ///< Modifies the data of an existing breakpoint.

 void setLoopStartIndex(int newLoopStartIndex);
 /**< Sets the start index for a loop - can be used for a generalized 
      'sustain'-phase. */

 void setLoopEndIndex(int newLoopEndIndex);
 /**< Sets the end index for a loop - can be used for a generalized 
      'sustain'-phase. */

 //---------------------------------------------------------------------------
 // audio processing:
 INLINE double getSample();    
 ///< Calculates one output sample at a time.

 //---------------------------------------------------------------------------
 // event-handling:
 void trigger(bool startFromCurrentValue = false);  
 /**< Causes the envelope to start with its attack-phase. When the parameter
      "startFromCurrentValue" is true, the internal state will not be reset 
      to startAmp, such that the curve begins at the level, where the envelope
      currently is. */  

 void noteOff();  
 ///< Causes the envelope to start with its release-phase.

 void reset();    
 ///< Resets the time variable.

 //---------------------------------------------------------------------------
 // inquiry:
 bool endIsReached();  
 ///< True, if output is below 40 dB.

 bool noteIsOn;   // indicates if note is on, if not, the envelope starts 
                    // with its release phase

 //===========================================================================

protected:

 // buffering-variables for the variuos recursions: 
 double state1, state2;

 // ... their associated increment/decrement (which can have multiplictive or
 // additive nature depending on the chosen shape, so we call it genreally
 // 'change'):
 double state1_change, state2_change;

 // ... their respective minimum-values:
 double state1_min, state2_min;

 // ... their respective maximum-values:
 double state1_max, state2_max;

 // some scale-factors:
 double scaler1, scaler2;
 

 // an array with the breakpoints:
 //static const int maxNumBreakpoints = 8;
 //ModBreakpoint breakpoints[maxNumBreakpoints];

 vector<ModBreakpoint> breakpoints;

 double sampleRate;

 int numBreakpoints; // number of breakpoints

 int countdownToNextBreakpoint;
  // a counter which counts down to the next breakpoint. when the counter 
  // reaches zero, the envelope has reached a new breakpoint

 int leftIndex;  
  // index of the breakpoint from which we come, this is  to the left on
  // the time-axis unless we sit exactly on it

 int rightIndex; 
  // index of the breakpoint to which we go

 int loopStartIndex, loopEndIndex;
  // indices of the start- and end-point of a loop in the modulator. The loop
  // realizes a generalized sustain phase. For an usual sustain, it is 
  // possible to have the loops start- and end-breakpoit equal.


 double leftLevel;  
  // level of the breakpoint from which we come

 double levelDelta;
 //double levelDeltaHalf;
  // the difference between the levels of the two breakpoints to the left
  // and to the right (and half its value)

};

//----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE double BreakpointModulator::getSample()
{
 double tmp1, tmp2; // temporary variables for intermediate results
 double out = 0.0;  // output variable

 // do the necessary initializations for the various variables when we are 
 // at a breakpoint:
 if( countdownToNextBreakpoint == 0 )
 {
  // increment the breakpoint indices to the next pair of breakpoints:
  leftIndex++;
  rightIndex++; 

  // wrap-around (introduce arbitrary loop-start and -end points here later):
  if( leftIndex >= (int) breakpoints.size() )
   leftIndex = 0;
  if( rightIndex >= (int) breakpoints.size() )
   rightIndex = 0;

  // get the time-difference between the two breakpoints (in seconds):
  double timeDelta =   breakpoints[rightIndex].timeStamp 
                     - breakpoints[leftIndex].timeStamp;

  // get the length of the envelope-segment to be generated in samples:
  int segmentLength = MoreMath::roundToInt(timeDelta*sampleRate);

  // use this length as initial value for our countdown-variable:
  countdownToNextBreakpoint = segmentLength;

  // get the current level on the left side and its difference to the level 
  // on the right side:
  leftLevel  = breakpoints[leftIndex].level;
  levelDelta = breakpoints[rightIndex].level - leftLevel;

  // do the specific initializations for the different envelope shapes (for 
  // details about what's going on, refer to comments in the MatLab
  // implementation):
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
  case ModBreakpoint::EXPLOSIVE:
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
  } // end of switch( breakpoints[rightIndex].shape )
 } // end of if( countdownToNextBreakpoint == 0 )



 // do the specific per-sample calculations for the different envelope shapes
 // (for  details about what's going on, refer to the comments in the MatLab
 // implementation):
 switch( breakpoints[rightIndex].shape )
 {
 case ModBreakpoint::LINEAR:
  {
   out     = state1;
   state1 += state1_change;
  }
  break;
 case ModBreakpoint::SMOOTH:
  {
   tmp1   = state1_change*state1 - state2;
   state2 = state1;
   state1 = tmp1;
   out    = leftLevel + levelDelta*(0.5*tmp1+0.5);
  }
  break;
 case ModBreakpoint::ANALOG:
  {
   out     = leftLevel + levelDelta - scaler1*(state1-state1_min);
   state1 *= state1_change;
  }
  break;
 case ModBreakpoint::EXPLOSIVE:
  {
   out     = leftLevel + scaler1*(state1-state1_min);
   state1 *= state1_change;
  }
  break;
 case ModBreakpoint::SIGMOID:
  {
   tmp1 = scaler1*(state1-state1_min);
   tmp2 = scaler1*(state2-state1_min); // state2_min == state1_min
   if( tmp1+tmp2 == 0.0 )
    tmp1 = 0.0;
   else
    tmp1 = (tmp2-tmp1) / (tmp2+tmp1);
   tmp1    = 0.5 * tmp1 + 0.5; // may be optimized away later
   out     = leftLevel + levelDelta*tmp1; 
   state1 *= state1_change;
   state2 *= state2_change;
  }
  break;
 case ModBreakpoint::SPIKEY:
  {
   tmp1    = scaler1*(state1-state1_min);
   tmp2    = scaler1*(state2-state1_min); // state2_min == state1_min
   tmp1    = 0.5*(tmp2-tmp1) + 0.5;
   out     = leftLevel + levelDelta*tmp1; 
   state1 *= state1_change;
   state2 *= state2_change;
  }
  break;
 } // end of switch( breakpoints[rightIndex].shape )

 // decrement countdown to next breakpoint:
 countdownToNextBreakpoint--; 

 // and return the generated output:
 return out;
}

#endif // BreakpointModulator_h
