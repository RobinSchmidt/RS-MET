#include "rosic_ModulationEffect.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulationEffect::ModulationEffect()
{
  sampleRate        = 44100.0;
  bpm               = 120.0;
  depth             = 0.5;
  stereoPhaseOffset = 0.0;
  resetOscillatorPhases();
  setCycleLength(0.5);     // 0.5 beats
}

ModulationEffect::~ModulationEffect()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void ModulationEffect::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    lfo.setSampleRate(sampleRate);
  }
}

void ModulationEffect::setCycleLength(double newCycleLength)
{ 
  lfo.setCycleLength(newCycleLength);
}

void ModulationEffect::setTempoSync(bool shouldTempoSync)
{
  lfo.setTempoSync(shouldTempoSync);
}

void ModulationEffect::setTempoInBPM(double newTempoInBPM)
{
  lfo.setBeatsPerMinute(newTempoInBPM);
}

//-------------------------------------------------------------------------------------------------
// others:

void ModulationEffect::resetOscillatorPhases(double phaseShift)
{
  lfo.triggerWithPhase(phaseShift);
  /*
  double phiL = lfoL.getStartPhase()+(phaseShift-0.5*stereoPhaseOffset); 
  double phiR = lfoR.getStartPhase()+(phaseShift+0.5*stereoPhaseOffset); 
  lfoL.triggerWithPhase(phiL);
  lfoR.triggerWithPhase(phiR);
  */
}



