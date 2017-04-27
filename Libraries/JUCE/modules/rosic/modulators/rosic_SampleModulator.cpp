#include "rosic_SampleModulator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SampleModulator::SampleModulator()
{
  // parameter initializition:
  sampleRate                 = 44100.0; // default sample rate
  position                   = 0.0;
  increment                  = 1.0;
  bpm                        = 120.0;
  numCyclesPerTimeUnit       = 1.0;
  numCyclesInLoop            = 1;
  numSamples                 = 0;
  start                      = 0;
  //startR                     = 0;
  loopStart                  = 0;
  loopEnd                    = 0;
  loopLength                 = 0;
  sampleData                 = NULL;
  loopIsOn                   = false;
  syncMode                   = true;
  endIsReached               = false;

  sampleName                 = NULL;
  sampleNameLength           = 0;

  audioProcessingIsSuspended = false;

  sampleDataEnd        = NULL;
  //int dummy            = 1;
}

SampleModulator::~SampleModulator()
{

}

void SampleModulator::copyDataFrom(const SampleModulator &source)
{
  audioProcessingIsSuspended = source.audioProcessingIsSuspended;
  bpm                        = source.bpm;
  increment                  = source.increment;
  loopEnd                    = source.loopEnd;
  loopIsOn                   = source.loopIsOn;
  loopLength                 = source.loopLength;
  loopStart                  = source.loopStart;
  numCyclesInLoop            = source.numCyclesInLoop;
  numCyclesPerTimeUnit       = source.numCyclesPerTimeUnit;
  numSamples                 = source.numSamples;
  position                   = source.position;
  sampleRate                 = source.sampleRate;
  start                      = source.start;
  //startR                     = source.startR;
  syncMode                   = source.syncMode;

  // create deep copies of the member arrays (WARNING: this will create memory leaks when this
  // object already has it's pointers assigned - i.e. the copyDataFrom was not called immediately
  // after construction):
  if( source.sampleName != NULL )
  {
    sampleNameLength     = source.sampleNameLength;
    sampleName           = new char[sampleNameLength+1];
    for(int c=0; c<=sampleNameLength; c++) // the <= is valid here, because we have one more cell allocated
      sampleName[c] = source.sampleName[c];
  }

  if( source.sampleData != NULL )
  {
    mutex.lock();
    sampleData = new double[numSamples+1];
    for(int n=0; n<numSamples; n++)
      sampleData[n] = source.sampleData[n];
    sampleData[numSamples] = sampleData[0];
    mutex.unlock();
  }
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void SampleModulator::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
    sampleRate = newSampleRate;


  //.....

  updateIncrement();
}

void SampleModulator::setSampleData(double *newSampleData, int newNumSamples)
{
  if( newSampleData == NULL || newNumSamples == 0 )
    return;

  // we need to aquire a mutex-lock because we re-allocate memory which is dereferenced in the
  // getSample()-function:
  suspendAudioProcessing();
  mutex.lock();

  // free old and allocate new memory if necesarry:
  if( newNumSamples != numSamples )
  {
    if( sampleData != NULL )
    {
      delete[] sampleData;
      sampleData = NULL;
    }
    numSamples = newNumSamples;
    sampleData = new double[numSamples+1];
  }

  // copy the sample-data into our own local buffer now, thereby repeat the first sample at the end
  // of the buffer for the linear interpolator:
  for(int n=0; n<numSamples; n++)
    sampleData[n] = newSampleData[n];
  sampleData[numSamples] = sampleData[0];

  // initialize the pointers:
  start      = 0;
  //startR     = 0;
  loopStart  = 0;
  loopEnd    = numSamples-1;
  loopLength = loopEnd - loopStart + 1;
  position   = 0.0;

  updateIncrement();

  // for debug:
  sampleDataEnd = &(sampleData[numSamples-1]);

  // release the mutex-lock:
  mutex.unlock();
  resumeAudioProcessing();
}

void SampleModulator::setSampleData(float *newSampleData, int newNumSamples)
{
  if( newSampleData == NULL || newNumSamples == 0 )
    return;

  // we need to aquire a mutex-lock because we re-allocate memory which is dereferenced in the
  // getSample()-function:
  suspendAudioProcessing();
  mutex.lock();

  // free old and allocate new memory if necesarry:
  if( newNumSamples != numSamples )
  {
    if( sampleData != NULL )
    {
      delete[] sampleData;
      sampleData = NULL;
    }
    numSamples = newNumSamples;
    sampleData = new double[numSamples+1];
  }

  // copy the sample-data into our own local buffer now, thereby repeat the first sample at the end
  // of the buffer for the linear interpolator:
  for(int n=0; n<numSamples; n++)
    sampleData[n] = (double) newSampleData[n];
  sampleData[numSamples] = sampleData[0];

  // initialize the pointers:
  start     = 0;
  //startR     = 0;
  loopStart  = 0;
  loopEnd    = numSamples-1;
  loopLength = loopEnd - loopStart + 1;
  position   = 0.0;

  updateIncrement();

  // for debug:
  sampleDataEnd = &(sampleData[numSamples-1]);

  // release the mutex-lock:
  mutex.unlock();
  resumeAudioProcessing();
}

void SampleModulator::setSampleName(char *newSampleName, int newLength)
{
  // free old and allocate new memory for the name:
  if( sampleName != NULL )
  {
    delete[] sampleName;
    sampleName = NULL;
  }
  sampleName       = new char[newLength+1];
  sampleNameLength = newLength;
  for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
    sampleName[c] = newSampleName[c];
}

void SampleModulator::setStartSample(int newStartSample)
{
  start = newStartSample;

  // make sure that the start-sample is in the valid range:
  if( start <= 0 )
    start = 0;
  else if( start >= numSamples )
    start = numSamples-1;
}

/*
void SampleModulator::setStartSampleLeft(int newStartSampleLeft)
{
  startL = newStartSampleLeft;

  // make sure that the start-sample is in the valid range:
  if( startL <= 0 )
    startL = 0;
  else if( startL >= numSamples )
    startL = numSamples-1;
}

void SampleModulator::setStartSampleRight(int newStartSampleRight)
{
  startR = newStartSampleRight;

  // make sure that the start-sample is in the valid range:
  if( startR <= 0 )
    startR = 0;
  else if( startR >= numSamples )
    startR = numSamples-1;
}
*/

void SampleModulator::setLoopMode(int newLoopMode)
{
  loopIsOn = (newLoopMode != NO_LOOP);
}

void SampleModulator::setLoopStartSample(int newLoopStartSample)
{
  loopStart = newLoopStartSample;

  // make sure that the loop start is before the loop end:
  if( loopStart >= loopEnd )
    loopStart = loopEnd-1;

  // make sure that the loop start is in the valid range:
  if( loopStart <= 0 )
    loopStart = 0;
  else if( loopStart >= numSamples )
    loopStart = numSamples-1;

  loopLength = loopEnd - loopStart + 1;

  updateIncrement();
}

void SampleModulator::setLoopEndSample(int newLoopEndSample)
{
  loopEnd = newLoopEndSample;

  // make sure that the loop end is behind the loop start:
  if( loopEnd <= loopStart )
    loopEnd = loopStart+1;

  // make sure that the loop start is in the valid range:
  if( loopEnd <= 0 )
    loopEnd = 0;
  else if( loopEnd >= numSamples )
    loopEnd = numSamples-1;

  loopLength = loopEnd - loopStart + 1;

  updateIncrement();
}

void SampleModulator::setNumCyclesInLoop(int newNumberOfCyclesInLoop)
{
  if( newNumberOfCyclesInLoop >= 1 )
    numCyclesInLoop = newNumberOfCyclesInLoop;
  updateIncrement();
}

void SampleModulator::setSyncMode(bool shouldBeSynced)
{
  syncMode = shouldBeSynced;
  updateIncrement();
}

void SampleModulator::setBeatsPerMinute(double newBpm)
{
  if( newBpm > 0.0 )
    bpm = newBpm;

  ///< \todo uncomment this, it's commented only for debug purposes:
  //updateIncrement();
}

void SampleModulator::setNumCyclesPerTimeUnit(double newNumCycles)
{
  if( newNumCycles > 0.0 )
    numCyclesPerTimeUnit = newNumCycles;

  updateIncrement();
}


//-------------------------------------------------------------------------------------------------
// inquiry (get-, is-, etc. functions):

double SampleModulator::getSampleRate()
{
  return sampleRate;
}

char* SampleModulator::getSampleName()
{
  return sampleName;
}

int SampleModulator::getNumSamples()
{
  return numSamples;
}

double* SampleModulator::getSampleData()
{
  return sampleData;
}

int SampleModulator::getStartSample()
{
  return start;
}

/*
int SampleModulator::getStartSampleRight()
{
  return startR;
}
*/

int SampleModulator::getLoopMode()
{
  return (int) loopIsOn;
}

int SampleModulator::getLoopStartSample()
{
  return loopStart;
}

int SampleModulator::getLoopEndSample()
{
  return loopEnd;
}

int SampleModulator::getNumCyclesInLoop()
{
  return numCyclesInLoop;
}

bool SampleModulator::isInSyncMode()
{
  return syncMode;
}

double SampleModulator::getNumCyclesPerTimeUnit()
{
  return numCyclesPerTimeUnit;
}

//-------------------------------------------------------------------------------------------------
// event-handling:

void SampleModulator::trigger()
{
  position = (double) start;
}

//-------------------------------------------------------------------------------------------------
// others:

void SampleModulator::updateIncrement()
{
  // convert the number of cycles per whole note into a number of cycles per second if necesarry:
  double cyclesPerSecond;
  if( syncMode == true )
    cyclesPerSecond = numCyclesPerTimeUnit / wholeNotesToSeconds(1.0, bpm);
  else
    cyclesPerSecond = numCyclesPerTimeUnit;

  increment = (double) loopLength*cyclesPerSecond / (sampleRate*numCyclesInLoop);
}

void SampleModulator::suspendAudioProcessing()
{
  audioProcessingIsSuspended = true;
}

void SampleModulator::resumeAudioProcessing()
{
  audioProcessingIsSuspended = false;
}

void SampleModulator::lockSampleData()
{
  mutex.lock();
}

void SampleModulator::unlockSampleData()
{
  mutex.unlock();
}



