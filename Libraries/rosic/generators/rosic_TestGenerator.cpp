#include "rosic_TestGenerator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

TestGenerator::TestGenerator()
{
  // initialize parameters:
  //sampleRate = 44100.0;
  //frequency  = 1000.0;
  //level      = 0.0;

  sineOscillator.trigger();
}

TestGenerator::~TestGenerator()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void TestGenerator::setSampleRate(double newSampleRate)
{
  sineOscillator.setSampleRate(newSampleRate);
  sineOscillator.trigger();
}

void TestGenerator::setFrequency(double newFrequency)
{
  sineOscillator.setFrequency(newFrequency);
  sineOscillator.trigger();
}

void TestGenerator::setLevel(double newLevel)
{
  sineOscillator.setAmplitude(dB2amp(newLevel));
  sineOscillator.trigger();
}

