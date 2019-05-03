//#include "rosic_EchoLabDelayLine.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

EchoLabDelayLine::EchoLabDelayLine(int maximumDelayInSamples)
{
  length       = maximumDelayInSamples + 1;
  delayBuffer1 = new(std::nothrow) float[length+1];
  delayBuffer2 = new(std::nothrow) float[length+1];  // are all these '+1's really needed?
  if( delayBuffer1 == NULL || delayBuffer2 == 0 )
    pointersInvalid = true;
  else
    pointersInvalid = false;

  tapIn        = 0;
  tapOut       = 0;
  delayTime    = 0.25;
  sampleRate   = 44100.0;
  bpm          = 120.0;
  feedback     = 0.0;
  pan          = 0.0;
  gL           = sqrt(0.5);
  gR           = sqrt(0.5);
  g            = 1.0;
  tempoSync    = false;
  pingPongMode = false;
  mute         = false;

  clearBuffers();
  setDelayTime(0.5);     // 20 ms
}

EchoLabDelayLine::~EchoLabDelayLine()
{
  if( delayBuffer1 != NULL )
  {
    delete[] delayBuffer1;
    delayBuffer1 = NULL;
  }
  if( delayBuffer2 != NULL )
  {
    delete[] delayBuffer2;
    delayBuffer2 = NULL;
  }
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void EchoLabDelayLine::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
    inputEqualizer.setSampleRate(sampleRate);
    feedbackEqualizer.setSampleRate(sampleRate);
  }
}

void EchoLabDelayLine::setDelayTime(double newDelayTime)
{
  if( newDelayTime != delayTime )
  {
    delayTime = RAPT::rsClip(newDelayTime, 0.0, 4.25);  
    setupDelayInSamples();
  }
}

void EchoLabDelayLine::setSyncMode(bool shouldTempoSync)
{
  if( shouldTempoSync != tempoSync )
  {
    tempoSync = shouldTempoSync;
    setupDelayInSamples();
  }
}

void EchoLabDelayLine::setTempoInBPM(double newTempoInBPM)
{
  if( newTempoInBPM >= 0.0 && newTempoInBPM != bpm )
  {
    bpm = newTempoInBPM;
    setupDelayInSamples();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void EchoLabDelayLine::clearBuffers()
{
  if( pointersInvalid )
    return;
  for(int i=0; i<length+1; i++)
  {
    delayBuffer1[i] = 0.0;
    delayBuffer2[i] = 0.0;
  }
  inputEqualizer.reset();
  feedbackEqualizer.reset();
}

void EchoLabDelayLine::setupDelayInSamples()
{
  double delayInSeconds;
  if( tempoSync )
    delayInSeconds = RAPT::rsBeatsToSeconds(delayTime, bpm);
  else
    delayInSeconds = delayTime;

  delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = RAPT::rsClip(delayInSamples, 0.0, (double)(length-1-1));

  // update member delayTime to the clipped value:
  delayInSeconds = delayInSamples / sampleRate;   
  if( tempoSync )
    delayTime = RAPT::rsSecondsToBeats(delayInSeconds, bpm);
  else
    delayTime = delayInSeconds;

  int dInt = (int) round(delayInSamples);
  tapOut   = tapIn - dInt;
  tapOut   = wrapAround(tapOut);
}