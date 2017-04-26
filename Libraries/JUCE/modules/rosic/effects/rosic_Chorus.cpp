#include "rosic_Chorus.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Chorus::Chorus(int bufferLengthToAllocate) : Vibrato(bufferLengthToAllocate)
{
  dryWetRatio  = 0.5;
  setAverageDelayTime(20.0);
  setStereoPhaseOffsetInDegrees(180.0);
  setCycleLength(4.0);
  setDepth(0.25);
  setGlobalFeedback(0.0);
  for(int i=0; i<numVoices; i++)
  {
    voiceIsActive[i] = true;
    delayScales[i]   = 1.f;
    depthScales[i]   = 1.f;
    ampScales[i]     = 1.f;
  }
  calculateGainCompensation();

  //feedbackFilter.setMode(OnePoleFilter::HIGHPASS);
  //feedbackFilter.setCutoff(50.0);

  reset();
}

Chorus::~Chorus()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Chorus::setSampleRate(double newSampleRate)
{
  Vibrato::setSampleRate(newSampleRate);
  //feedbackFilter.setSampleRate(newSampleRate);
}

//-------------------------------------------------------------------------------------------------
// others:

void Chorus::calculateGainCompensation()
{
  double tmp = 0.0;
  for(int i=0; i<numVoices; i++)
  {
    if( voiceIsActive[i] )
      tmp += ampScales[i]*ampScales[i];
  }
  tmp = sqrt(tmp);
  tmp = rmax(tmp, 0.5);
  voiceGainCompensator = 1.0/tmp;
}

void Chorus::reset()
{
  Vibrato::reset();
  //feedbackFilter.resetBuffers();
  yOldL = 0.0;
  yOldR = 0.0;
}