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

void AcidPattern::circularShiftAll(int numStepsToShift)
{
  RAPT::rsArrayTools::circularShift(notes, maxNumSteps, numStepsToShift);
}


int rsMod(int val, int modulus)
{
  while(val > modulus) val -= modulus;
  while(val < 0      ) val += modulus;
  return val;
}
// move somewhere else

void AcidPattern::circularShiftAccents(int shift)
{
  AcidPattern tmp = *this;
  for(int i = 0; i < numSteps; i++)
    notes[i].accent = tmp.notes[rsMod(i-shift, numSteps)].accent;
}
// seems not to work - maybe the mod operation doe not work for negative numbers

void AcidPattern::circularShiftSlides(int shift)
{
  AcidPattern tmp = *this;
  for(int i = 0; i < numSteps; i++)
    notes[i].slide = tmp.notes[rsMod(i-shift, numSteps)].slide;
}

void AcidPattern::circularShiftOctaves(int shift)
{
  AcidPattern tmp = *this;
  for(int i = 0; i < numSteps; i++)
    notes[i].octave = tmp.notes[rsMod(i-shift, numSteps)].octave;
}

void AcidPattern::circularShiftNotes(int shift)
{
  AcidPattern tmp = *this;
  for(int i = 0; i < numSteps; i++) {
    notes[i].key  = tmp.notes[rsMod(i-shift, numSteps)].key;
    notes[i].gate = tmp.notes[rsMod(i-shift, numSteps)].gate; }
}

// maybe have functions to swap slides and accents etc - maybe this functionality would be easier
// to implement with parallel arrays


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