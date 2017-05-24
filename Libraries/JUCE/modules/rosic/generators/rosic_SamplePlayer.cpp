//#include "rosic_SamplePlayer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SamplePlayer::SamplePlayer()
{
  // parameter initializition:
  sampleRate       = 44100.0; // default sample rate
  position         = 0.0;
  increment        = 1.0;
  frequencyNominal = 440.0;         
  frequencyDetuned = 440.0;      
  gainFactor       = 1.0;
  key              = 64;
  vel              = 64;
  loopSnap         = true;
  theBuffer        = NULL;
  parameters       = new SamplePlaybackParameters;
  isMaster         = true;
}

SamplePlayer::~SamplePlayer()
{
  if( isMaster ) 
  { 
    if( parameters != NULL )
      delete parameters;
    /*
    if( theBuffer != NULL )
      delete theBuffer;
    */
  }
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void SamplePlayer::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    calculateIncrement();
    filterL.setSampleRate(sampleRate);
    filterR.setSampleRate(sampleRate);
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->setSampleRate(sampleRate);
  }
}

void SamplePlayer::setSampleBufferToUse(SampleBuffer *newSampleBufferToUse)
{
  theBuffer = newSampleBufferToUse;
}

bool SamplePlayer::setSampleData(float **newSampleData, int newNumSamples, int newNumChannels)
{
  bool success;
  if( theBuffer != NULL )
  {
    success = theBuffer->setSampleData(newSampleData, newNumSamples, newNumChannels);

    // update the playback-info:
    parameters->setNumSamples(theBuffer->getNumSamples());
    parameters->setNumChannels(theBuffer->getNumChannels());
    parameters->initSettings();
  }
  else
    success = false;

  return success;
}

void SamplePlayer::setRootKey(double newRootKey)
{
  parameters->setRootKey(newRootKey);
  calculateIncrement();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->calculateIncrement();
}

/*
void SamplePlayer::setRootDetune(double newRootDetune)
{
  parameters->setRootDetune(newRootDetune);
  calculateIncrement();
  for(int s=0; s<slaves.size(); s++)
    slaves[s]->calculateIncrement();
}
*/

void SamplePlayer::setLoopStart(double newLoopStart)
{
  if( loopSnap == true )
  {
    parameters->setLoopStart(theBuffer->findNearestUpwardZeroCrossing(newLoopStart));
    //parameters->setLoopStart( findNearestUpwardZeroCrossing(
    //  theBuffer->getChannelSampleData(0), theBuffer->getNumSamples(), newLoopStart) );
  }
  else
    parameters->setLoopStart(newLoopStart);
}

void SamplePlayer::setLoopEnd(double newLoopEnd)
{
  if( loopSnap == true )
  {
    parameters->setLoopEnd(theBuffer->findNearestUpwardZeroCrossing(newLoopEnd));
    //parameters->setLoopEnd( findNearestUpwardZeroCrossing(
    //  theBuffer->getChannelSampleData(0), theBuffer->getNumSamples(), newLoopEnd) );
  }
  else
    parameters->setLoopEnd(newLoopEnd);
}

void SamplePlayer::setLoopLength(double newLength)
{
  setLoopEnd(parameters->getLoopStart() + newLength);
}

void SamplePlayer::setSnapLoopToZeros(bool shouldSnap)
{
  loopSnap = shouldSnap;
}

void SamplePlayer::setSampleName(char *newSampleName)
{
  parameters->setSampleName(newSampleName);
}

//-------------------------------------------------------------------------------------------------
// inquiry (get-, is-, etc. functions):

double SamplePlayer::getSampleRate()
{
  return sampleRate;
}

char* SamplePlayer::getSampleName()
{
  return parameters->getSampleName();
}

/*
float** SamplePlayer::getSampleData()
{
  if( theBuffer != NULL )
    return theBuffer->getSampleData();
  else
    return NULL;
}
*/

int SamplePlayer::getNumChannels()
{
  if( theBuffer != NULL )
    return theBuffer->getNumChannels();
  else
    return 0;
}

int SamplePlayer::getNumSamples()
{
  if( theBuffer != NULL )
    return theBuffer->getNumSamples();
  else
    return 0;
}

bool SamplePlayer::getLoopSnapToZeroMode()
{
  return loopSnap;
}

//-------------------------------------------------------------------------------------------------
// master/slave config:

void SamplePlayer::addSlave(SamplePlayer* newSlave)
{
  // add the new slave to the vector of slaves:
  slaves.push_back(newSlave);

  // delete the original parameter-set of the new slave and redirect it to ours (with some safety 
  // checks):
  if( newSlave->parameters != NULL && newSlave->parameters != this->parameters )
  {
    delete newSlave->parameters;
    newSlave->parameters = this->parameters;
  }
  else
  {
    DEBUG_BREAK; 
    // the object to be added as slave did not contain a valid parameter-pointer - maybe it has 
    // been already added as slave to another master?
  }

  // delete the original buffer of the new slave and redirect it to ours (with some safety 
  // checks):
  if( newSlave->theBuffer != NULL && newSlave->theBuffer != this->theBuffer )
    delete newSlave->theBuffer;
  newSlave->theBuffer = this->theBuffer;

  // set the isMaster-flag of the new slave to false: 
  newSlave->isMaster = false;

  // this flag will prevent the destructor of the slave from trying to delete the parameter-set 
  // which is now shared - only masters delete their parameter-set on destruction
}

//-------------------------------------------------------------------------------------------------
// others:

void SamplePlayer::calculateKeyAndVelocityDependentParameters()
{
  gainFactor  = dB2amp(parameters->getLevel());
  gainFactor *= dB2amp(parameters->getLevelByKey() * (double)(key-64)/63.0);
  gainFactor *= dB2amp(parameters->getLevelByVel() * (double)(vel-64)/63.0);

  /*
  // this can be optimized algebraically to call dB2amp only once:
  double level = parameters->getLevel() 
    + (  parameters->getLevelByKey() * (double)(key-64) 
       + parameters->getLevelByVel() * (double)(vel-64) ) / 63.0;
  gainFactor = dB2amp(level);
  */

  // todo: the pitch/frequency settings...
  double pitch = parameters->getRootKey() 
    + parameters->getTune()
    + parameters->getTuneByKey() * 0.01*(key-parameters->getRootKey())
    + parameters->getTuneByVel() * 0.01*(vel-64.0);
  double freq = pitchToFreq(pitch);
  setPlaybackFrequencyNominal(freq);

  //double pitchOffset = 
  //double freq =
  //freqWithKeyAndVel  = parameters->freqNominal;
  //freqWithKeyAndVel *= pow(2.0, (0.01*parameters->freqByKey/12.0) * (currentKey-64.0) );
  //freqWithKeyAndVel *= pow(2.0, (0.01*parameters->freqByVel/63.0) * (currentVel-64.0) );


  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->calculateKeyAndVelocityDependentParameters();
}

void SamplePlayer::reset()
{
  position = parameters->getPlaybackStart();
}

void SamplePlayer::lockSampleData()
{
  if( theBuffer != NULL )
    theBuffer->lockSampleData();
}

void SamplePlayer::unlockSampleData()
{
  if( theBuffer != NULL )
    theBuffer->unlockSampleData();
}





