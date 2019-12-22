//#include "rosic_OscillatorStereo.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

OscillatorStereo::OscillatorStereo()
{
  positionIntL       = 0;
  positionFracL      = 0.0;
  incrementIntL      = 1;
  incrementFracL     = 0.0;
  positionIntR       = 0;
  positionFracR      = 0.0;
  incrementIntR      = 1;
  incrementFracR     = 0.0;
  mipMapIndex        = 0;
  frequencyNominal   = 440.0;
  frequencyDetunedL  = 440.0;
  frequencyDetunedR  = 440.0;
  scaledAmplitude    = 1.0;
  modulatedAmplitude = 1.0;
  panFactorL         = 1.0;
  panFactorR         = 1.0;
  scaledStartPhase   = 0.0;
  key                = 64;
  vel                =  0;
  waveTable          = NULL;
  parameters         = new OscillatorStereoParameters;
  isMaster           = true;

  setSampleRate(44100.0);
}

OscillatorStereo::~OscillatorStereo()
{
  if( isMaster && parameters != NULL )
    delete parameters;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void OscillatorStereo::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.1 )
  {
    parameters->sampleRate    = newSampleRate;
    parameters->sampleRateRec = 1.0 / newSampleRate;
  }

  // trigger increment re-calculation:
  calculateIncrement();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->calculateIncrement();
}

void OscillatorStereo::setWaveTableToUse(MipMappedWaveTableStereo *newTableToUse)
{
  waveTable = newTableToUse;
  parameters->tableLength = waveTable->getTableLength();
  reset();
  calculateIncrement();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->calculateIncrement();
}

void OscillatorStereo::setTimeReverse(bool shouldBeReversed)
{
  if( waveTable != NULL )
  {
    waveTable->setTimeReverse(shouldBeReversed);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setPolarityInversion(bool shouldBeInversed)
{
  if( waveTable != NULL )
  {
    waveTable->setPolarityInversion(shouldBeInversed);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setFullWavePhaseWarp(double newWarpCoefficient)
{
  if( waveTable != NULL )
  {
    waveTable->setFullWavePhaseWarp(newWarpCoefficient);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setHalfWavePhaseWarp(double newWarpCoefficient)
{
  if( waveTable != NULL )
  {
    waveTable->setHalfWavePhaseWarp(newWarpCoefficient);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setKeyAndVel(int newKey, int newVel)
{
  key = newKey;
  vel = newVel;
  calculateKeyAndVelocityDependentParameters();
}

void OscillatorStereo::setMute(bool shouldBeMuted)
{
  parameters->mute = shouldBeMuted;
  //markPresetAsDirty();
}

void OscillatorStereo::setLevel(double newLevel)
{
  parameters->level     = newLevel;
  parameters->amplitude = RAPT::rsDbToAmp(newLevel);
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();

  // it's a bit wasteful to re-calculate all key and velocity dependent parameters when only one
  // of them has actually changed - however it is not supposed to be called at samplerate and makes
  // for cleaner code
}

void OscillatorStereo::setLevelByKey(double newLevelByKey)
{
  parameters->levelByKey = newLevelByKey;
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();
}

void OscillatorStereo::setLevelByVel(double newLevelByVel)
{
  parameters->levelByVel = newLevelByVel;
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();
}

void OscillatorStereo::setPan(double newPan)
{
  parameters->pan = newPan;
  double p = 0.25*PI*(newPan+1.0);
  RAPT::rsSinCos(p, &(parameters->panFactorR), &(parameters->panFactorL) );
  //markPresetAsDirty();
}

void OscillatorStereo::setPanByKey(double newPanByKey)
{
  parameters->panByKey = newPanByKey;
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();
}

void OscillatorStereo::setPanByVel(double newPanByVel)
{
  parameters->panByVel = newPanByVel;
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();
}

void OscillatorStereo::setMidSide(double newMidSide)
{
  parameters->midSide = newMidSide;
  double x = 0.5 * PI * newMidSide;
  RAPT::rsSinCos(x, &(parameters->sideScale), &(parameters->midScale) );
  //markPresetAsDirty();
}

void OscillatorStereo::setStartPhase(double newStartPhase)
{
  if( newStartPhase >= 0.0 && newStartPhase <= 360 )
  {
    parameters->startPhase    = newStartPhase;
    parameters->startPosition = parameters->tableLength * (newStartPhase/360.0);
  }
}

void OscillatorStereo::setStartPhaseByKey(double newStartPhaseByKey)
{
  parameters->startPhaseByKey = newStartPhaseByKey;
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();
}

void OscillatorStereo::setStartPhaseByVel(double newStartPhaseByVel)
{
  parameters->startPhaseByVel = newStartPhaseByVel;
  calculateKeyAndVelocityDependentParameters();
  //markPresetAsDirty();
}

void OscillatorStereo::setDetuneSemitones(double newDetuneSemitones)
{
  if( newDetuneSemitones != parameters->detuneSemitones )
  {
    parameters->detuneSemitones = newDetuneSemitones;
    parameters->detuneFactor    = RAPT::rsPitchOffsetToFreqFactor(newDetuneSemitones);
    parameters->detuneFactorL   = parameters->detuneFactor * parameters->stereoDetuneFactorL;
    parameters->detuneFactorR   = parameters->detuneFactor * parameters->stereoDetuneFactorR;

    // trigger increment re-calculation:
    calculateIncrement();
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->calculateIncrement();

    //markPresetAsDirty();
  }
}

void OscillatorStereo::setDetuneHz(double newDetuneHz)
{
  if( newDetuneHz != parameters->detuneHz )
  {
    parameters->detuneHz    = newDetuneHz;
    parameters->freqOffsetR = parameters->detuneHz + 0.5*parameters->stereoDetuneHz;
    parameters->freqOffsetL = parameters->detuneHz - 0.5*parameters->stereoDetuneHz;

    // trigger increment re-calculation:
    calculateIncrement();
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->calculateIncrement();

    //markPresetAsDirty();
  }
}

void OscillatorStereo::setStereoDetuneSemitones(double newStereoDetuneSemitones)
{
  if( newStereoDetuneSemitones != parameters->stereoDetuneSemitones )
  {
    parameters->stereoDetuneSemitones = newStereoDetuneSemitones;
    parameters->stereoDetuneFactorR   = RAPT::rsPitchOffsetToFreqFactor(0.5*newStereoDetuneSemitones);
    parameters->stereoDetuneFactorL   = 1.0 / parameters->stereoDetuneFactorR;
    parameters->detuneFactorL         = parameters->detuneFactor * parameters->stereoDetuneFactorL;
    parameters->detuneFactorR         = parameters->detuneFactor * parameters->stereoDetuneFactorR;

    // trigger increment re-calculation:
    calculateIncrement();
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->calculateIncrement();

    //markPresetAsDirty();
  }
}

void OscillatorStereo::setStereoDetuneHz(double newStereoDetuneHz)
{
  if( newStereoDetuneHz != parameters->stereoDetuneHz )
  {
    parameters->stereoDetuneHz = newStereoDetuneHz;
    parameters->freqOffsetR    = parameters->detuneHz + 0.5*parameters->stereoDetuneHz;
    parameters->freqOffsetL    = parameters->detuneHz - 0.5*parameters->stereoDetuneHz;

    // trigger increment re-calculation:
    calculateIncrement();
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->calculateIncrement();

    //markPresetAsDirty();
  }
}

void OscillatorStereo::setPitchEnvelopeDepth(double newDepth)
{
  parameters->pitchEnvDepth = newDepth;
}

void OscillatorStereo::setSpectralContrast(double newContrast)
{
  if( waveTable != NULL )
  {
    waveTable->setSpectralContrast(newContrast);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setSpectralSlope(double newSlope)
{
  if( waveTable != NULL )
  {
    waveTable->setSpectralSlope(newSlope);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setHighestHarmonicToKeep(int newHighestHarmonicToKeep)
{
  if( waveTable != NULL )
  {
    waveTable->setHighestHarmonicToKeep(newHighestHarmonicToKeep);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setLowestHarmonicToKeep(int newLowestHarmonicToKeep)
{
  if( waveTable != NULL )
  {
    waveTable->setLowestHarmonicToKeep(newLowestHarmonicToKeep);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setEvenOddRatio(double newRatio)
{
  if( waveTable != NULL )
  {
    waveTable->setEvenOddRatio(newRatio);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setEvenOddPhaseShift(double newPhaseShift)
{
  if( waveTable != NULL )
  {
    waveTable->setEvenOddPhaseShift(newPhaseShift);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setStereoPhaseShift(double newPhaseShift)
{
  if( waveTable != NULL )
  {
    waveTable->setStereoPhaseShift(newPhaseShift);
    //markPresetAsDirty();
  }
}

void OscillatorStereo::setEvenOddStereoPhaseShift(double newPhaseShift)
{
  if( waveTable != NULL )
  {
    waveTable->setEvenOddStereoPhaseShift(newPhaseShift);
    //markPresetAsDirty();
  }
}

/*
void OscillatorStereo::setTranspositionFactor(double newTranspositionFactor)
{
  if( newTranspositionFactor > 0.00001 )
    parameters->transpositionFactor = newTranspositionFactor;

  setFrequency(freq);
  for(int s=0; s<slaves.size(); s++)
    slaves[s]->setFrequency(slaves[s]->freq);
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

bool OscillatorStereo::isMuted()
{
  return parameters->mute;
}

double OscillatorStereo::getLevel()
{
  return parameters->level;
}

double OscillatorStereo::getLevelByKey()
{
  return parameters->levelByKey;
}

double OscillatorStereo::getLevelByVel()
{
  return parameters->levelByVel;
}

double OscillatorStereo::getPan()
{
  return parameters->pan;
}

double OscillatorStereo::getPanByKey()
{
  return parameters->panByKey;
}

double OscillatorStereo::getPanByVel()
{
  return parameters->panByVel;
}

double OscillatorStereo::getMidSide()
{
  return parameters->midSide;
}

double OscillatorStereo::getStartPhase()
{
  return parameters->startPhase;
}

double OscillatorStereo::getStartPhaseByKey()
{
  return parameters->startPhaseByKey;
}

double OscillatorStereo::getStartPhaseByVel()
{
  return parameters->startPhaseByVel;
}

double OscillatorStereo::getDetuneSemitones()
{
  return parameters->detuneSemitones;
}

double OscillatorStereo::getDetuneHz()
{
  return parameters->detuneHz;
}

double OscillatorStereo::getStereoDetuneSemitones()
{
  return parameters->stereoDetuneSemitones;
}

double OscillatorStereo::getStereoDetuneHz()
{
  return parameters->stereoDetuneHz;
}

void OscillatorStereo::getWaveformForDisplay(double **targetBuffer, int numSamplesToShow)
{
  if( waveTable == NULL )
    return;

  if( numSamplesToShow <= 0 || targetBuffer == NULL )
  {
    DEBUG_BREAK; // invalid size or address for the targetBuffer
    return;
  }

  // determine the ratio of the lengths of between the buffer to be read and the buffer to be
  // written - this is somewhat similar to a phase-increment for the read-buffer:
  double  readWriteRatio = (double) waveTable->getTableLength() / (double) numSamplesToShow;

  // we need some temporary arrays for the min/max extraction:
  int     tmpArrayLength = (int) (2.0*ceil(readWriteRatio));
  //int     tmpArrayLength = (int) floor(2 * readWriteRatio) - 1;
  double* tmpArrayL      = new double[tmpArrayLength];
  double* tmpArrayR      = new double[tmpArrayLength];
  double  tmpL, tmpR, tmpM, tmpS;

  // indices for read and write and length of the read-buffer chunk (to be decimated):
  int writeIndex = 0;
  int readIndex, nextReadIndex, readLength;
  int shift       = roundToInt(waveTable->getTableLength() * getStartPhase() / 360.0);

  //readLength = tmpArrayLength;

  while( (writeIndex+1) < numSamplesToShow )
  {
    readIndex     = (int) roundToInt( writeIndex    * readWriteRatio);
    nextReadIndex = (int) roundToInt((writeIndex+2) * readWriteRatio);
    readLength    = nextReadIndex - readIndex;

    // readout the table, thereby convert to L/R and apply M/S-scaling and pan:
    for(int n=0; n<tmpArrayLength; n++)
    {
      // this is a bit dirty because if there wouldn't be a whole table-set in the mip-map, we
      // could get access violations - the clean way would be to account for wrapaorunds here....

      tmpL         = waveTable->getValue(0, 0, RAPT::rsWrapAround(readIndex+n+shift, waveTable->getTableLength()), 0.0);
      tmpR         = waveTable->getValue(1, 0, RAPT::rsWrapAround(readIndex+n+shift, waveTable->getTableLength()), 0.0);
      tmpM         = parameters->midScale  * (tmpL+tmpR);
      tmpS         = parameters->sideScale * (tmpL-tmpR);
      tmpArrayL[n] = parameters->amplitude * parameters->panFactorL * (tmpM+tmpS);
      tmpArrayR[n] = parameters->amplitude * parameters->panFactorR * (tmpM-tmpS);
    }

    // extract min/max from the chunk and write into the target buffer:
    int minIndex = RAPT::rsArrayTools::minIndex(tmpArrayL, readLength);
    int maxIndex = RAPT::rsArrayTools::maxIndex(tmpArrayL, readLength);


    // debug:
    //double minValue = tmpArrayL[minIndex];
    //double maxValue = tmpArrayL[maxIndex];

    if( minIndex < maxIndex )
    {
      targetBuffer[0][writeIndex]   = tmpArrayL[minIndex];
      targetBuffer[0][writeIndex+1] = tmpArrayL[maxIndex];
    }
    else
    {
      targetBuffer[0][writeIndex]   = tmpArrayL[maxIndex];
      targetBuffer[0][writeIndex+1] = tmpArrayL[minIndex];
    }

    minIndex = RAPT::rsArrayTools::minIndex(tmpArrayR, tmpArrayLength);
    maxIndex = RAPT::rsArrayTools::maxIndex(tmpArrayR, tmpArrayLength);
    if( minIndex < maxIndex )
    {
      targetBuffer[1][writeIndex]   = tmpArrayR[minIndex];
      targetBuffer[1][writeIndex+1] = tmpArrayR[maxIndex];
    }
    else
    {
      targetBuffer[1][writeIndex]   = tmpArrayR[maxIndex];
      targetBuffer[1][writeIndex+1] = tmpArrayR[minIndex];
    }

    writeIndex += 2;
  }

  // in case of odd display-widths, we need to repeat the pre-final value (otherwise a zero will result as final value):
  if( RAPT::rsIsOdd(numSamplesToShow) )
  {
    targetBuffer[0][numSamplesToShow-1] = targetBuffer[0][numSamplesToShow-2];
    targetBuffer[1][numSamplesToShow-1] = targetBuffer[1][numSamplesToShow-2];
  }

  delete[] tmpArrayL;
  delete[] tmpArrayR;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// master/slave config:

void OscillatorStereo::addSlave(OscillatorStereo* newSlave)
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

  // set the isMaster-flag of the new slave to false:
  newSlave->isMaster = false;

  // this flag will prevent the destructor of the slave from trying to delete the parameter-set
  // which is now shared - only masters delete their parameter-set on destruction
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void OscillatorStereo::calculateIncrementForAllSlaves()
{
  calculateIncrement();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->calculateIncrement();
}

void OscillatorStereo::reset()
{
  double startPosition = (double)parameters->tableLength * (scaledStartPhase/360.0);
  positionIntL         = (int) floor(startPosition);
  positionFracL        = startPosition - positionIntL;
  positionIntR         = positionIntL;
  positionFracR        = positionFracL;
}

void OscillatorStereo::calculateKeyAndVelocityDependentParameters()
{
  scaledAmplitude  = RAPT::rsDbToAmp(parameters->level);
  scaledAmplitude *= RAPT::rsDbToAmp(parameters->levelByKey * (double)(key-64)/63.0);
  scaledAmplitude *= RAPT::rsDbToAmp(parameters->levelByVel * (double)(vel-64)/63.0);
    // this can be optimized algebraically to call dB2amp only once

  scaledStartPhase  = parameters->startPhase;
  scaledStartPhase += parameters->startPhaseByKey * (double)(key-64)/63.0;
  scaledStartPhase += parameters->startPhaseByVel * (double)(vel-64)/63.0;

  // todo: calculate pan factors...


  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->calculateKeyAndVelocityDependentParameters();
}
