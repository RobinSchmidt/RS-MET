#include "rosic_PitchShifterGrainAdaptive.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

PitchShifterGrainAdaptive::PitchShifterGrainAdaptive() // : formantPreserver(64)
{
  grainLengthInMilliseconds = 20.0;
  grainLengthInPitchCycles  = 8.0;
  grainLengthInBeats        = 1.0/2.0;
  tempoInBpm                = 120.0;
  grainLengthUnit           = MILLISECONDS;
  formantPreserve           = false;

  smoother.setMode(FourPoleFilterParameters::LOWPASS_6);
  smoother.useTwoStages(true);
  smoother.setFrequency(60.0);
  smoother.setSampleRate(sampleRate);

  pitchDetector.setSampleRate(sampleRate);
}

PitchShifterGrainAdaptive::~PitchShifterGrainAdaptive()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void PitchShifterGrainAdaptive::setSampleRate(double newSampleRate)
{
  PitchShifter::setSampleRate(newSampleRate);
  smoother.setSampleRate(newSampleRate);
  pitchDetector.setSampleRate(newSampleRate);
}

void PitchShifterGrainAdaptive::setBeatsPerMinute(double newBpm)
{
  if( newBpm > 0.0 )
  {
    tempoInBpm = newBpm;
    setGrainLengthUnit(grainLengthUnit); // takes care for updating the grainlength
    /*
    if( grainLengthUnit == BEATS )
      PitchShifter::setGrainLength(1000.0*beatsToSeconds(grainLengthInBeats, tempoInBpm));
    */
  }
}

void PitchShifterGrainAdaptive::setGrainLengthInMilliseconds(double newGrainLength)
{
  if( newGrainLength >= 1.0 )
  {
    if( 0.001*newGrainLength*sampleRate < delayLineLength-10 ) // some arbitrary safety-margin
    {
      grainLengthInMilliseconds = newGrainLength;
      if( grainLengthUnit == MILLISECONDS )
      {
        PitchShifter::setGrainLength(grainLengthInMilliseconds);
      }
    }
  }
}

void PitchShifterGrainAdaptive::setGrainLengthInPitchCycles(double newNumCycles)
{
  if( newNumCycles >= 0.25 )
    grainLengthInPitchCycles = newNumCycles;
  else
    DEBUG_BREAK; // values <= 0.25 are invalid
}

void PitchShifterGrainAdaptive::setGrainLengthInBeats(double newNumBeats)
{
  if( newNumBeats >= 0.0 )
  {
    grainLengthInBeats = newNumBeats;
    setGrainLengthUnit(grainLengthUnit); // takes care for updating the grainlength
    //PitchShifter::setGrainLength(1000.0*beatsToSeconds(grainLengthInBeats, tempoInBpm));
  }
  else
    DEBUG_BREAK; // values <= 0 are invalid
}

void PitchShifterGrainAdaptive::setGrainLengthUnit(int newUnit)
{
  if( newUnit < MILLISECONDS || newUnit > BEATS )
    DEBUG_BREAK;

  grainLengthUnit = newUnit;
  if( grainLengthUnit == MILLISECONDS )
    PitchShifter::setGrainLength(grainLengthInMilliseconds);
  else if( grainLengthUnit == BEATS )
    PitchShifter::setGrainLength(1000.0*beatsToSeconds(grainLengthInBeats, tempoInBpm));

  // else... unit is PITCH_CYCLES -> setGrainLength will be called in getSample()
}

void PitchShifterGrainAdaptive::setFormantPreserve(bool shouldPreserveFormants)
{
  formantPreserve = shouldPreserveFormants;
  formantPreserver.reset();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

/*
double PitchShifterGrainAdaptive::getGrainLengthInMilliseconds()
{
  return grainLengthInMilliseconds;
}

double PitchShifterGrainAdaptive::getGrainLengthInPitchCycles()
{
  return grainLengthInPitchCycles;
}

double PitchShifterGrainAdaptive::getGrainLengthInBeats()
{
  return grainLengthInBeats;
}

int PitchShifterGrainAdaptive::getGrainLengthUnit()
{
  return grainLengthUnit;
}

bool PitchShifterGrainAdaptive::getFormantPreserve()
{
  return formantPreserve;
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------
// internal functions:





