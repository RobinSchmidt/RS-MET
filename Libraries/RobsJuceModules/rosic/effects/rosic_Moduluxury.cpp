//#include "rosic_Moduluxury.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Moduluxury::Moduluxury()
{
  modulationSource    = BREAKPOINT_MODULATOR;
  modulationTarget    = AMPLITUDE;
  panRule             = CONSTANT_POWER_PAN;
  filterFrequency     = 1000.0;
  delayModulationMode = ADDITIVE;
  delay               = 0.01;
  midiControlNumber   = 7;
  scaleFactor         = 1.0;
  offset              = 0.0;
}

Moduluxury::~Moduluxury()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void Moduluxury::setSampleRate(double newSampleRate)
{
  Modulator::setSampleRate(newSampleRate);
  filter.setSampleRate(newSampleRate);
}

void Moduluxury::setModulationTarget(int newModulationTarget)
{
  if( newModulationTarget >= AMPLITUDE       &&
    newModulationTarget <= MIDI_CONTROLLER      )
  {
    modulationTarget = newModulationTarget;
  }
}

void Moduluxury::setPanRule(int newPanRule)
{
  if( newPanRule >= LINEAR_PAN &&
      newPanRule <= CONSTANT_POWER_PAN                )
  {
    panRule = newPanRule;
  }
}

void Moduluxury::setFilterFrequency(double newFilterFrequency)
{
  if( newFilterFrequency <= 20.0 )
    filterFrequency = 20.0;
  else if( newFilterFrequency >= 20000.0 )
    filterFrequency = 20000.0;
  else
    filterFrequency = newFilterFrequency;
}

void Moduluxury::setDelay(double newDelay)
{
  if( newDelay >= 0.0 )
    delay = newDelay;
}

void Moduluxury::setDelayModulationMode(int newDelayModulationMode)
{
  if( newDelayModulationMode >= ADDITIVE &&
      newDelayModulationMode <= MULTIPLICATIVE  )
  {
    delayModulationMode = newDelayModulationMode;
  }
}

void Moduluxury::setMidiControllerNumber(int newControllerNumber)
{
  if( newControllerNumber >= 0 && newControllerNumber <= 127 )
    midiControlNumber = newControllerNumber;
}

/*
void Moduluxury::setScaleFactor(double newScaleFactor)
{
  if( newScaleFactor >= 0.0 )
    scaleFactor = newScaleFactor;
}

void Moduluxury::setOffset(double newOffset)
{
  offset = newOffset;
}
*/

//-------------------------------------------------------------------------------------------------
// inquiry (get-, is-, etc. functions):

int Moduluxury::getModulationTarget()
{
  return modulationTarget;
}

int Moduluxury::getPanRule()
{
  return panRule;
}

double Moduluxury::getFilterFrequency()
{
  return filterFrequency;
}

int Moduluxury::getDelayModulationMode()
{
  return delayModulationMode;
}

double Moduluxury::getDelay()
{
  return delay;
}

int Moduluxury::getMidiControllerNumber()
{
  return midiControlNumber;
}

/*
double Moduluxury::getScaleFactor()
{
  return scaleFactor;
}

double Moduluxury::getOffset()
{
  return offset;
}
*/