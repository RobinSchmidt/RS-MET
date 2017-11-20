#include "TimeStretcher.h"

TimeStretcher::TimeStretcher(int outputBufferSize)
{
  stretchFactor     = 1.f;
  pitchFactor       = 1.f;
  readIndexActual   = 0;
  readIndexApparent = 0;
  outFrameIndex     = 0;
  outBufferLength   = outputBufferSize;
  outBuffer         = new float[outBufferLength];
  tmpBuffer         = new float[outBufferLength];
  clearBuffer(outBuffer, outBufferLength);
  clearBuffer(tmpBuffer, outBufferLength);
}

TimeStretcher::~TimeStretcher()
{
  delete[] outBuffer;
  delete[] tmpBuffer;
}

void TimeStretcher::setInputStream(float *newInputStream, int newLength)
{
  inStream.setInputDataAddress(newInputStream, newLength);
  reset();
}

void TimeStretcher::setStretchFactor(float newStretchFactor) 
{ 
  setStretchAndPitchFactor(newStretchFactor, pitchFactor); 
}

void TimeStretcher::setPitchFactor(float newPitchFactor) 
{ 
  setStretchAndPitchFactor(stretchFactor, newPitchFactor); 
}

void TimeStretcher::getBuffer(float *buffer, int numFrames)
{
  int numUnconsumedOutputFrames = outBufferLength - outFrameIndex;
  if( numFrames <= numUnconsumedOutputFrames )
  {
    copyBuffer(&outBuffer[outFrameIndex], buffer, numFrames); 
    outFrameIndex += numFrames;
  }
  else
  {
    int numRemainingFrames = numFrames;  // # frames that we still have to write into "buffer"
    int bufferIndex        = 0;          // current frame index in "buffer"
    int copyLength         = 0;
    while( true )
    {
      numUnconsumedOutputFrames = outBufferLength - outFrameIndex;
      copyLength = min(numUnconsumedOutputFrames, numRemainingFrames);
      copyBuffer(&outBuffer[outFrameIndex], &buffer[bufferIndex], copyLength); 
      bufferIndex        += copyLength;
      outFrameIndex      += copyLength;
      numRemainingFrames -= copyLength;
      if( numRemainingFrames == 0 )
        break;
      fillOutputBuffer();                // we need fresh data in outBuffer
    }
  }
  readIndexApparent += numFrames / stretchFactor;
}

void TimeStretcher::reset()
{
  clearBuffer(outBuffer, outBufferLength);
  clearBuffer(tmpBuffer, outBufferLength);
  readIndexActual   = 0;
  readIndexApparent = 0;
  outFrameIndex     = 0;
  resetSubclassVariables();
}

void TimeStretcher::prepareForCrossfade()
{
  int oldReadIndexActual = readIndexActual;
  long double oldReadIndexApparent = readIndexApparent;
  getBuffer(tmpBuffer, outBufferLength); 
  readIndexActual   = oldReadIndexActual ;
  readIndexApparent = oldReadIndexApparent;
}

void TimeStretcher::applyCrossfade()
{
  for(int n = 0; n < outBufferLength; n++)
  {
    float x = (float) n / (float) (outBufferLength-1);
    float fadeOut = (float) (cos(PI*x)+1)/2;
    float fadeIn  = 1 - fadeOut;
    outBuffer[n] = fadeOut * tmpBuffer[n] + fadeIn * outBuffer[n];
  }
}

//-------------------------------------------------------------------------------------------------

TimeStretcherElastique::TimeStretcherElastique(int outputBufferSize)
  : TimeStretcher(outputBufferSize)
{
  if(outputBufferSize > 1024)
    printf("%s", "Elastique doesn't support output buffer sizes > 1024.\n");

  pcElastiqueHandle = 0;
  int error = CElastiqueProIf::CreateInstance(pcElastiqueHandle, outBufferLength,
    1,                  // 1 channel
    44100.f,            // sampleRate
    CElastiqueProIf::_elastiquePro_proc_mode::kProDefaultMode);
  if(error)
    printf("%s", "Error creating Elastique instance.\n");

  maxInBufferLength = pcElastiqueHandle->GetMaxFramesNeeded();
  inBuffer  = new float[maxInBufferLength];
  ppIn  = &inBuffer;
  ppOut = &outBuffer;

  resetSubclassVariables();
}

TimeStretcherElastique::~TimeStretcherElastique()
{
  CElastiqueProIf::DestroyInstance(pcElastiqueHandle);
  delete[] inBuffer;
}

void TimeStretcherElastique::setStretchAndPitchFactor(float newStretchFactor, float newPitchFactor)
{
  if( stretchFactor == newStretchFactor && pitchFactor == newPitchFactor )
    return;
  prepareForCrossfade();
  stretchFactor = newStretchFactor;
  pitchFactor   = newPitchFactor;
  pcElastiqueHandle->SetStretchPitchQFactor(stretchFactor, pitchFactor, false);    
  pcElastiqueHandle->Reset();
  warmUpElastique();
  applyCrossfade();
}

void TimeStretcherElastique::jumpToInputFrame(int frameIndex)
{
  prepareForCrossfade();
  readIndexApparent = frameIndex;
  pcElastiqueHandle->Reset();
  warmUpElastique();
  applyCrossfade();
}

void TimeStretcherElastique::resetSubclassVariables()
{
  clearBuffer(inBuffer,  maxInBufferLength);
  pcElastiqueHandle->Reset();
  warmUpElastique();
}

void TimeStretcherElastique::fillOutputBuffer()
{
  int numInFrames = pcElastiqueHandle->GetFramesNeeded();
  inStream.getBuffer(inBuffer, readIndexActual, numInFrames);
  int error = pcElastiqueHandle->ProcessData(ppIn, numInFrames, ppOut);
  readIndexActual += numInFrames;
  outFrameIndex = 0;
}

void TimeStretcherElastique::warmUpElastique()
{
  readIndexApparent         = floor(readIndexApparent);
  int numExtraWarmUpBuffers = 3; // 3 seems to be required to avoid artifacts with a sinusoid 
  int numFramesNeeded       = pcElastiqueHandle->GetFramesNeeded();
  int readIndexApparentInt  = (int) readIndexApparent;
  readIndexActual           = readIndexApparentInt - numExtraWarmUpBuffers * outBufferLength;
  inStream.getBuffer(inBuffer, readIndexActual, numFramesNeeded);
  int error = pcElastiqueHandle->ProcessData(ppIn, numFramesNeeded, ppOut);
  readIndexActual += numFramesNeeded;
  outFrameIndex = 0;
  for(int i = 0; i < numExtraWarmUpBuffers; i++)
  {
    numFramesNeeded = pcElastiqueHandle->GetFramesNeeded();
    inStream.getBuffer(inBuffer, readIndexActual, numFramesNeeded);
    int error = pcElastiqueHandle->ProcessData(ppIn, numFramesNeeded, ppOut);
    readIndexActual += numFramesNeeded;
  }
}

//-------------------------------------------------------------------------------------------------

TimeStretcherSoundTouch::TimeStretcherSoundTouch(int outputBufferSize)
  : TimeStretcher(outputBufferSize)
{
  soundTouch.setChannels(1);
  soundTouch.setSampleRate(44100);
  soundTouch.setSetting(SETTING_USE_AA_FILTER,  1); 
  soundTouch.setSetting(SETTING_SEEKWINDOW_MS, 10);  // default:  22
  soundTouch.setSetting(SETTING_OVERLAP_MS,    10);  // default:   8

  inBufferLength = 512; 
  inBuffer  = new float[inBufferLength];

  setStretchAndPitchFactor(1.f, 1.f);
  resetSubclassVariables();
}

TimeStretcherSoundTouch::~TimeStretcherSoundTouch()
{
  delete[] inBuffer;
}

void TimeStretcherSoundTouch::setStretchAndPitchFactor(float newStretchFactor, float newPitchFactor)
{
  if( stretchFactor == newStretchFactor && pitchFactor == newPitchFactor )
    return;
  prepareForCrossfade();
  stretchFactor = newStretchFactor;
  pitchFactor   = newPitchFactor;

  soundTouch.setPitch(pitchFactor);
  soundTouch.setTempo(1.f/stretchFactor);

  // set SoundTouch's "sequence-length" according to the pitch/stretch settings:  
  float factor    = stretchFactor * pitchFactor;
  int   seqLength = (int) ceil((factor*25.f));
  soundTouch.setSetting(SETTING_SEQUENCE_MS, seqLength);

  soundTouch.clear();
  warmUpSoundTouch();
  applyCrossfade();
}

void TimeStretcherSoundTouch::jumpToInputFrame(int frameIndex)
{
  prepareForCrossfade();
  readIndexApparent = frameIndex;
  soundTouch.clear();
  warmUpSoundTouch();
  applyCrossfade();
}

void TimeStretcherSoundTouch::resetSubclassVariables()
{
  clearBuffer(inBuffer, inBufferLength);
  soundTouch.clear();
  warmUpSoundTouch();
}

void TimeStretcherSoundTouch::fillOutputBuffer()
{
  int framesRequired = outBufferLength;
  int framesReceived = 0;
  while( framesRequired > 0 )
  {
    framesReceived  = soundTouch.receiveSamples(&outBuffer[outBufferLength-framesRequired], 
      framesRequired);
    framesRequired -= framesReceived;
    if( framesRequired > 0 )
    {
      inStream.getBuffer(inBuffer, readIndexActual, inBufferLength);
      soundTouch.putSamples(inBuffer, inBufferLength);
      readIndexActual += inBufferLength;
    }
  }
  outFrameIndex = 0;
}

void TimeStretcherSoundTouch::warmUpSoundTouch()
{
  readIndexApparent        = floor(readIndexApparent);
  int readIndexApparentInt = (int) readIndexApparent;
  readIndexActual          = readIndexApparentInt;
  fillOutputBuffer();
}
