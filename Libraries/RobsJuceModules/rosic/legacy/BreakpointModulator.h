#ifndef BreakpointModulator_h
#define BreakpointModulator_h

#include "Definitions.h"
#include "MoreMath.h"

#include <float.h>
#include <iterator>
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
  GROWING,
  SIGMOID,
  SPIKEY
 };

 // member variables:
 double timeStamp;
 double level;
 double shapeAmount;
 int    shape;

 /** Standard-Contructor. */
 ModBreakpoint()
 {
  timeStamp   = 0.0;  
  level       = 0.0;
  shapeAmount = 1.0;
  shape       = 1;
 }

 ~ModBreakpoint() {}; ///< Destructor.
};


//============================================================================
// the actual BreakpointModulator-class:

class BreakpointModulator
{
public:

 /** This is an enumeration of the availabel edit-modes EDIT_WITH_SHIFT means
     that all succeeding breakpoints will be time-shifted when one breakpoint
     is being inserted, moved or removed. */
 enum editModes
 {
  EDIT_WITHOUT_SHIFT = 1,
  EDIT_WITH_SHIFT
 };

 enum loopModes
 {
  NO_LOOP = 0,
  FORWARD_LOOP
 };

 //---------------------------------------------------------------------------
 // construction/destruction:

 BreakpointModulator();  ///< Constructor.
 ~BreakpointModulator();  ///< Destructor.


 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Sets the sample-rate for this module.

 int  insertBreakpoint(double newTimeStamp, double newLevel, int newShape = 0,
                       double newShapeAmount = 0.0);
 /**< Inserts a new breakpoint into the vector. The new breakpoint must 
      satisfy some preconditions in order to be successfully inserted (it is 
      not allowed to be too close to an already existing breakpoint, for 
      example). The integer return-value informs, at which index if the new 
      breakpoint was inserted. It will return -1, when the breakpoint could not be 
      inserted. */

 bool modifyBreakpoint(int index, double newTimeStamp, double newLevel, 
                       int newShape = 0, double newShapeAmount = 0.0);
 /**< Modifies the data of an existing breakpoint. If you try to modify a 
      non-existent breakpoint or try to modify the data in such a way which is 
      not allowed, it will return false. */

 bool removeBreakpoint(int index);
 /**< Removes a breakpoint from the vector. The return-value informs, if 
      there was actually a breakpoint removed (if you try to remove a 
      non-existing breakpoint, or a breakpoint which cannot be removed, it 
      will return false */

 void setLoopMode(bool shouldBeLooped);
 ///< Turns looping on or off.

 bool setLoopStartIndex(int newLoopStartIndex);
 /**< Sets the start index for a loop - can be used for a generalized 
      'sustain'-phase. The loop-start index must be at least 1 and at most 
      the numberOfBreakpoints-2 ...... */

 bool setLoopEndIndex(int newLoopEndIndex);
 /**< Sets the end index for a loop - can be used for a generalized 
      'sustain'-phase. .........*/

 //---------------------------------------------------------------------------
 // inquiry:

 int lastBreakpointIndex();
 ///< Returns the index of the last breakpoint (the end-point).

 double getStartTime();
 /**< Returns the starting time (in seconds), is usually zero. */

 double getEndTime();
 /**< Returns the end time (in seconds), i.e. the time stamp of the last 
      breakpoint. */

 double getMinLevel();
 /**< Returns the minimum level of the curve. */

 double getMaxLevel();
 /**< Returns the maximum level of the curve. */

 double getBreakpointTime(int index);
 /**< Returns the absolute time (in seconds) of one breakpoint. If the index 
      is out of range, it will return -1.0. */

 double getBreakpointLevel(int index);
 /**< Returns the level of one breakpoint. If the index is out of range, it 
      will return 0.0. */

 int getLoopMode();
 /**< Returns the current loop mode. */

 int getLoopStartIndex();
 /**< Returns the index of the breakpoint at which the loop starts. */

 int getLoopEndIndex();
 /**< Returns the index of the breakpoint at which the loop ends. */

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE double getSample();    
 ///< Calculates one output sample at a time.

 //---------------------------------------------------------------------------
 // event-handling:

 void noteOn(bool startFromCurrentLevel = false);  
 /**< Causes the envelope to start with its attack-phase. When the parameter
      "startFromCurrentLevel" is true, the internal state will not be reset 
      to startAmp, such that the curve begins at the level, where the envelope
      currently is. */  

 void noteOff(bool startFromCurrentLevel = true);  
 ///< Causes the envelope to start with its release-phase.

 void initialize();    
 ///< Sets everything except the sampleRate to default-values.


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

 double timeScale;
  // scale the overall time-length of the modulator.
 
 vector<ModBreakpoint> breakpoints;

 double sampleRate;

 int samplesToNextBreakpoint;
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

 bool loopIsOn;
  // flag to indicate, if the modulator should run a loop

 double leftLevel;  
  // level of the breakpoint from which we come

 double levelDelta;
 //double levelDeltaHalf;
  // the difference between the levels of the two breakpoints to the left
  // and to the right (and half its value)

 bool noteIsOn;   
  // indicates if note is on, if not, the envelope will not trap into the 
  // loop

 double previousOut;

 int    editMode;

 double minBreakpointDistance;
  // a minimum distance between two successive breakpoints

 void setupStateVariables();
  // this function calculates and initializes the required state-variables 
  // for the recursions.

public:

 bool endIsReached;
  // flag to indicate, if the modulator has reched its end-point

};

//----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE double BreakpointModulator::getSample()
{
 double tmp1, tmp2; // temporary variables for intermediate results
 double out = 0.0;  // output variable

 // return the value of the last breakpoint, when our endisReached-flag is 
 // true:
 if( endIsReached )
  return breakpoints[breakpoints.size()-1].level;

 // do the necessary initializations for the various variables when we are 
 // at a breakpoint:
 if( samplesToNextBreakpoint == 0 )
 {
  // increment the breakpoint indices to the next pair of breakpoints:
  leftIndex++;
  rightIndex++; 

  // wrap-around for looping:
  if( noteIsOn && loopIsOn )
  {
   if( leftIndex >= loopEndIndex )
    leftIndex = loopStartIndex;
   if( rightIndex > loopEndIndex )
    rightIndex = loopStartIndex+1;
  }

  // if we have incremented our right index beyond the end of the 
  // breakpoint-array, we have reached the end:
  if( rightIndex >= (int) breakpoints.size() )
  {
   endIsReached = true;
   return breakpoints[breakpoints.size()-1].level;
  }

  // get the length of the envelope-segment to be generated (in samples):
  double timeDelta  =   breakpoints[rightIndex].timeStamp 
                      - breakpoints[leftIndex].timeStamp;
  int segmentLength = MoreMath::roundToInt(timeScale*timeDelta*sampleRate);

  // use this length as initial value for our countdown-variable:
  samplesToNextBreakpoint = segmentLength;

  // get the current level on the left side and its difference to the level 
  // on the right side:
  leftLevel  = breakpoints[leftIndex].level;
  levelDelta = breakpoints[rightIndex].level - leftLevel;

  // set up the internal state variables for the recursive formulas:
  setupStateVariables();

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
 case ModBreakpoint::GROWING:
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
 samplesToNextBreakpoint--; 

 // remember the output for the next call (might be needed, if we enter the 
 // release-phase (the segments(s) behind the loop) due to a note-off):
 previousOut = out;

 // and return the generated output:
 return out;
}

#endif // BreakpointModulator_h
