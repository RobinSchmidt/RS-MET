#include "rosic_Phaser.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

Phaser::Phaser()
{
  pitch          = freqToPitch(1000.0);
  dryWet         = 1.0;
  feedbackFactor = 0.0;
  bypass         = false;
  reset();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void Phaser::setSampleRate(double newSampleRate)
{
  ModulationEffect::setSampleRate(newSampleRate);
  allpassL.setSampleRate(newSampleRate);
  allpassR.setSampleRate(newSampleRate);
}
   
void Phaser::setFilterMode(int newMode)
{
  allpassL.setMode(newMode);
  allpassR.setMode(newMode);
}
   
void Phaser::setFrequency(double newFrequency)
{
  pitch = freqToPitch(newFrequency);
}

void Phaser::setQ(double newQ)
{
  allpassL.setQ(newQ);
  allpassR.setQ(newQ);
}

void Phaser::setFeedbackFactor(double newFeedbackFactor)
{
  feedbackFactor = newFeedbackFactor;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void Phaser::reset()
{
  allpassL.reset();
  allpassR.reset();
  yL = 0.0;
  yR = 0.0;
}

