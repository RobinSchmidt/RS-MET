//#include "rosic_Chorus.h"
//using namespace rosic;

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
  tmp = RAPT::rsMax(tmp, 0.5);
  voiceGainCompensator = 1.0/tmp;
}

void Chorus::reset()
{
  Vibrato::reset();
  //feedbackFilter.resetBuffers();
  yOldL = 0.0;
  yOldR = 0.0;
}

/*

Ideas:
-Combine regular chorus techniques based on delayline modulation (i.e. slow frequency modulation) 
 with ringmodulation, amplitude modulation and single-sideband modulation. Maybe like: RM a signal
 by 5 Hz and then FM the output of that by 7 Hz or the other way around. And/or maybe try parallel
 connections.
-Add a fade-in parameter. The goal is to leave the transients dry. Or maybe do a 
 transient/release separation before the chorus unit ad feed a use adjustable mix of it into the 
 chorus. At default settings, both components have a gain of 1. That could be useful for other
 effects, too
-Maybe before applying any sort of effect we should first:
 -convert L/R to M/S
 -split low from high frequencies
 -split both into attack and sustain/release (with longer settings for the bass part)

*/