#include "rosic_PingPongEcho.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

PingPongEcho::PingPongEcho(int maximumDelayInSamples)
{
  if( maximumDelayInSamples < 512 )
    DEBUG_BREAK;  // you need allocate some reasonable amount of memory

  length         = maximumDelayInSamples;
  buffer1        = new(std::nothrow) double[length];
  buffer2        = new(std::nothrow) double[length]; 
  tapIn          = 0;
  tapOut         = 0;
  delayTime      = 0.25;
  sampleRate     = 44100.0;
  bpm            = 120.0;
  feedback       = 0.0;
  pan            = 0.0;
  gLL            = sqrt(0.5);
  gRR            = sqrt(0.5);
  gLR            = 0.0;
  gRL            = 0.0;
  dryWetRatio    = 0.5;
  dry            = sqrt(0.5);
  wet            = sqrt(0.5);
  g              = 1.0;
  tempoSync      = true;
  pingPongMode   = true;
  trueStereoMode = false;
  mute           = false;
  bypass         = false;
  stereoSwap     = false;
  reset();
  setDelayTime(0.5);     // 0.5 beats
}

PingPongEcho::~PingPongEcho()
{
  delete[] buffer1;
  delete[] buffer2;
}

//-------------------------------------------------------------------------------------------------
// setup:

void PingPongEcho::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.01 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
    filter1.setSampleRate(sampleRate);
    filter2.setSampleRate(sampleRate);
  }
}

void PingPongEcho::setDelayTime(double newDelayTime)
{
  if( newDelayTime != delayTime )
  {
    delayTime = clip(newDelayTime, 0.0, 4.25);  
    setupDelayInSamples();
  }
}

void PingPongEcho::setSyncMode(bool shouldTempoSync)
{
  tempoSync = shouldTempoSync;
  setupDelayInSamples();
}

void PingPongEcho::setTempoInBPM(double newTempoInBPM)
{
  if( newTempoInBPM < 0.0 )
  {
    //DEBUG_BREAK;
    return;
  }
  if( newTempoInBPM > 0.0 && newTempoInBPM != bpm )
  {
    bpm = newTempoInBPM;
    setupDelayInSamples();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void PingPongEcho::reset()
{
  for(int i=0; i<length; i++)
  {
    buffer1[i] = 0.0;
    buffer2[i] = 0.0;
  }
  filter1.reset();
  filter2.reset();
}

void PingPongEcho::setupDelayInSamples()
{
  double delayInSeconds;
  if( tempoSync )
    delayInSeconds = beatsToSeconds(delayTime, bpm);
  else
    delayInSeconds = delayTime;

  double delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = clip(delayInSamples, 0.0, (double)length);

  /*
  // update member delayTime to the clipped value:
  delayInSeconds = delayInSamples / sampleRate;   
  if( tempoSync )
    delayTime = secondsToBeats(delayInSeconds, bpm);
  else
    delayTime = delayInSeconds;
  */

  int dInt = (int) round(delayInSamples);
  tapOut   = tapIn - dInt;
  tapOut   = wrapAround(tapOut, length);
      
  reset();
}

void PingPongEcho::calculateGainFactors()
{
  if( trueStereoMode == true )
    StereoPan::panLawLinearCrossMixSquareNormalized(pan, gLL, gRL, gLR, gRR);
  else
    StereoPan::panLawSinCos(pan, gLL, gRL, gLR, gRR);
  gLL *= g;
  gLR *= g;
  gRL *= g;
  gRR *= g;
}