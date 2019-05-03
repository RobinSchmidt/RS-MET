//#include "rosic_LorentzSystem.h"
//using namespace rosic;

// Construction/Destruction:

LorentzSystem::LorentzSystem()
{
  x = 0.5;
  y = 0.0;
  z = 0.0;

  sigma = 10.0;
  rho   = 28.0;
  beta  = 8.0/3.0;

  sampleRate      = 44100.0;
  pseudoFrequency = 1000.0;
  updateStepSize();
}

// Setup:

void LorentzSystem::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  updateStepSize();
}

void LorentzSystem::setPseudoFrequency(double newPseudoFrequency)
{
  pseudoFrequency = newPseudoFrequency;
  updateStepSize();
}
