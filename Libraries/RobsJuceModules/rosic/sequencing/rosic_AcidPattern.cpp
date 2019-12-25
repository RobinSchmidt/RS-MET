//#include "rosic_AcidPattern.h"
//using namespace rosic;

AcidPattern::AcidPattern()
{
  numSteps   = 16;
  stepLength = 0.5;
}

//-------------------------------------------------------------------------------------------------   
// setup:

void AcidPattern::clear()
{
  for(int i=0; i<maxNumSteps; i++)
  {
    notes[i].key    = 0;
    notes[i].octave = 0;
    notes[i].accent = false;
    notes[i].slide  = false;
    notes[i].gate   = false;
  }
}

void AcidPattern::randomize()
{
  for(int i=0; i<maxNumSteps; i++)
  {
    notes[i].key    = roundToInt(RAPT::rsRandomUniform( 0, 11));
    notes[i].octave = roundToInt(RAPT::rsRandomUniform(-2,  2));
    notes[i].accent = roundToInt(RAPT::rsRandomUniform( 0,  1)) == 1;
    notes[i].slide  = roundToInt(RAPT::rsRandomUniform( 0,  1)) == 1;
    notes[i].gate   = roundToInt(RAPT::rsRandomUniform( 0,  1)) == 1;
  }
}

void AcidPattern::circularShift(int numStepsToShift)
{
  RAPT::rsArrayTools::circularShift(notes, maxNumSteps, numStepsToShift);
}

//-------------------------------------------------------------------------------------------------   
// inquiry:

bool AcidPattern::isEmpty() const
{
  for(int i=0; i<maxNumSteps; i++)
  {
    if( notes[i].gate == true )
      return false;
  }
  return true;
}