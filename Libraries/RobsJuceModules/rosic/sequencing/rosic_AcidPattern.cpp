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
  // todo: use PRNG member, provide the options to randomize only accents, slides, etc.
}

void AcidPattern::setNumSteps(int newNumber)
{
  RAPT::rsAssert(newNumber <= maxNumSteps);
  numSteps = maxNumSteps;
}

void AcidPattern::circularShiftAll(int numStepsToShift)
{
  RAPT::rsArrayTools::circularShift(notes, maxNumSteps, numStepsToShift);
}


int rsMod(int val, int modulus)
{
  while(val >= modulus) val -= modulus;
  while(val < 0       ) val += modulus;
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

void AcidPattern::reverseAll()
{
  for(int i = 0; i < numSteps/2; i++)
    RAPT::rsSwap(notes[i], notes[numSteps-1-i]);
}

void AcidPattern::reverseAccents()
{
  for(int i = 0; i < numSteps/2; i++)
    RAPT::rsSwap(notes[i].accent, notes[numSteps-1-i].accent);
}

void AcidPattern::reverseSlides()
{
  for(int i = 0; i < numSteps/2; i++)
    RAPT::rsSwap(notes[i].slide, notes[numSteps-1-i].slide);
}

void AcidPattern::reverseOctaves()
{
  for(int i = 0; i < numSteps/2; i++)
    RAPT::rsSwap(notes[i].octave, notes[numSteps-1-i].octave);
}

void AcidPattern::reverseNotes()
{
  for(int i = 0; i < numSteps/2; i++)
    RAPT::rsSwap(notes[i].key, notes[numSteps-1-i].key);
}

void AcidPattern::invertAccents()
{
  for(int i = 0; i < numSteps; i++)
    notes[i].accent = !notes[i].accent;
}

void AcidPattern::invertSlides()
{
  for(int i = 0; i < numSteps; i++)
    notes[i].slide = !notes[i].slide;
}

void AcidPattern::invertOctaves()
{
  for(int i = 0; i < numSteps; i++)
    notes[i].octave = -notes[i].octave;
}

void AcidPattern::swapAccentsWithSlides()
{
  for(int i = 0; i < numSteps; i++)
    RAPT::rsSwap(notes[i].accent, notes[i].slide);
}

void AcidPattern::xorAccentsWithSlides()
{
  for(int i = 0; i < numSteps; i++)
    notes[i].accent = rsXor(notes[i].accent, notes[i].slide);
}

void AcidPattern::xorSlidesWithAccents()
{
  for(int i = 0; i < numSteps; i++)
    notes[i].accent = rsXor(notes[i].accent, notes[i].slide);
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