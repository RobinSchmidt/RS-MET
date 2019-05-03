//#include "rosic_SamplePlaybackParameters.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SamplePlaybackParameters::SamplePlaybackParameters()
{
  recordingSampleRate  = 44100.0;
  //fundamentalFrequency = 440.0;
  sampleName           = NULL;
  //sampleNameLength     = 0;
  numSamples           = 0;
  numChannels          = 0;

  // the rest of the member variables can be initialized by a member-function:
  initSettings();
}

SamplePlaybackParameters::~SamplePlaybackParameters()
{
  if( sampleName != NULL )
    delete[] sampleName;
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void SamplePlaybackParameters::setSampleName(char *newSampleName)
{
  // free old and allocate new memory for the name:
  if( sampleName != NULL )
  {
    delete[] sampleName;
    sampleName = NULL;
  }
  if( newSampleName != NULL )
  {
    int newLength    = (int) strlen(newSampleName);
    sampleName       = new char[newLength+1];
    for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
      sampleName[c] = newSampleName[c];
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void SamplePlaybackParameters::ensureDataValidity()
{
  if( numSamples < 0 )
    numSamples = 0;
  if( numChannels < 0 )
    numChannels = 0;
  /*
  if( fundamentalFrequency < 0.0 )
    fundamentalFrequency = 0.0;
  */
  if( recordingSampleRate <= 0.1 )
    recordingSampleRate = 0.1;

  ensureLoopValidity();
  ensureKeyRangeValidity();
  ensureVelRangeValidity();
}

void SamplePlaybackParameters::ensureLoopValidity()
{
  if( playbackStart < 0 )
    playbackStart = 0;
  else if( playbackStart > numSamples-1 )
    playbackStart = numSamples-1;

  if( playbackEnd < 0 )
    playbackEnd = 0;
  else if( playbackEnd > numSamples-1 )
    playbackEnd = numSamples-1;

  if( loopStart >= loopEnd )
    loopStart = loopEnd-1;
  if( loopStart <= 0 )
    loopStart = 0;
  else if( loopStart > playbackEnd )
    loopStart = playbackEnd;

  if( loopEnd <= loopStart )
    loopEnd = loopStart+1;
  if( loopEnd <= 0 )
    loopEnd = 0;
  else if( loopEnd > playbackEnd )
    loopEnd = playbackEnd;

  if( numSamples == 0 )
  {
    playbackStart = loopStart = loopEnd = 0;
  }

  double loopLength = loopEnd - loopStart;
  if( loopCrossfadeLength < 0 )
    loopCrossfadeLength = 0;
  else if( loopCrossfadeLength > loopLength/2 )
    loopCrossfadeLength = loopLength/2;

  if( loopCrossfadeShape < 0.0 )
    loopCrossfadeShape = 0.0;
  if( loopCrossfadeShape > 1.0 )
    loopCrossfadeShape = 1.0;

  if( pitchCyclesInLoop < 1.0 )
    pitchCyclesInLoop = 1.0;
}

void SamplePlaybackParameters::ensureKeyRangeValidity()
{
  if( loKey < 0 )
    loKey = 0;
  if( loKey > 127 )
    loKey = 127;
  if( hiKey < 0 )
    hiKey = 0;
  if( hiKey > 127 )
    hiKey = 127;
  if( loKey > hiKey )
    loKey = hiKey;
  if( fadeInToKey < loKey )
    fadeInToKey = loKey;
  if( fadeOutFromKey > hiKey )
    fadeOutFromKey = hiKey;
}

void SamplePlaybackParameters::ensureVelRangeValidity()
{
  if( loVel < 0 )
    loVel = 0;
  if( loVel > 127 )
    loVel = 127;
  if( hiVel < 0 )
    hiVel = 0;
  if( hiVel > 127 )
    hiVel = 127;
  if( loVel > hiVel )
    loVel = hiVel;
  if( fadeInToVel < loVel )
    fadeInToVel = loVel;
  if( fadeOutFromVel > hiVel )
    fadeOutFromVel = hiVel;
}

void SamplePlaybackParameters::initSettings()
{
  initLoopSettings();
  initKeyAndVelRanges();
  initLevelSettings();
  initPitchSettings();
  initPanSettings();
  initFilterSettings();
}

void SamplePlaybackParameters::initLoopSettings()
{
  loopMode            = 0;
  pitchCyclesInLoop   = 1.0;
  playbackStart       = 0.0;
  playbackStartByVel  = 0.0;
  playbackStartByKey  = 0.0;
  loopStart           = 0.0;
  loopCrossfadeLength = 0.0;
  loopCrossfadeShape  = 0.0;
  loopLengthLock      = false;

  if( numSamples > 0 )
  {
    playbackEnd = loopEnd = (double) (numSamples-1);
  }
  else
  {
    playbackEnd = loopEnd = 0;
  }
}

void SamplePlaybackParameters::initKeyAndVelRanges()
{
  rootKey        = 64;
  //rootDetune     = 0.0;
  loKey          = 0;
  hiKey          = 127;
  loVel          = 0;
  hiVel          = 127;
  fadeInToKey    = 0;
  fadeOutFromKey = 127;
  fadeInToVel    = 0;
  fadeOutFromVel = 127;
}

void SamplePlaybackParameters::initLevelSettings()
{
  mute          = false;
  solo          = false;
  level         = 0.0;
  levelByKey    = 0.0;
  levelByVel    = 0.0;
  levelStart    = 0.0;
  levelTime     = 0.0;
  levelRelease  = 0.0;
}

void SamplePlaybackParameters::initPitchSettings()
{
  tune      = 0.0;
  detuneHz  = 0.0;
  tuneByKey = 100.0;
  tuneByVel = 0.0;
}

void SamplePlaybackParameters::initPanSettings()
{
  pan      = 0.0;
  panByKey = 0.0;
  panByVel = 0.0;
  panCurve = 0.0;
  midSide  = 0.0;
}

void SamplePlaybackParameters::initFilterSettings()
{
  filterOnOff    = false;
  lowpassCutoff  = 20000.0;
  lowpassByKey   = 0.0;
  lowpassByVel   = 0.0;
  highpassCutoff = 20.0;
  highpassByKey  = 0.0;
  highpassByVel  = 0.0;
}




