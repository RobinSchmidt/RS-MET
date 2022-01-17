//#include "rosic_EchoLab.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

EchoLab::EchoLab()
{
  mutex.lock();
  sampleRate          = 44100.0;
  bpm                 = 120.0;
  soloedIndex         = -1;
  delayTimeSync       = false; 
  delayModulationSync = false; 
  ampModulationSync   = false;; 
  snapToTimeGrid      = true;
  wetLevel            = 0.0;
  setDryWet(0.5);
  mutex.unlock();
}

EchoLab::~EchoLab()
{
  mutex.lock();
  removeAllDelayLines();
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void EchoLab::copyDataFrom(EchoLab *other)
{
  mutex.lock();
  other->mutex.lock();

  this->removeAllDelayLines();
  for(int i=0; i<other->getNumDelayLines(); i++)
  {
    //this->addDelayLine(other->getDelayTime(i), other->getGainFactor(i));

    // \todo: bring the other per-delayline parameters in sync also (feedback, pan, eq-settings, whatever...)
    // maybe implement a copyDataFrom function in ModulatedDelayLine as well. ...and in Equalizer 
    // ...and use it then from here

  }

  other->mutex.unlock();
  mutex.unlock();
}

void EchoLab::setSampleRate(double newSampleRate)
{
  mutex.lock();
  if( newSampleRate <= 0.0 )
  {
    DEBUG_BREAK;
    mutex.unlock();
    return;
  }
  sampleRate = newSampleRate;
  for(unsigned int i = 0; i < delayLines.size(); i++)
    delayLines[i]->setSampleRate(sampleRate);
  mutex.unlock();
}

void EchoLab::setDryWet(double newDryWet)
{
  mutex.lock();
  dryWet = newDryWet;
  RAPT::rsEqualPowerGainFactors(dryWet, &dryFactor, &wetFactor, 0.0, 1.0);
  wetFactor *= RAPT::rsDbToAmp(wetLevel);
  mutex.unlock();
}

void EchoLab::setWetLevel(double newLevel)
{
  mutex.lock();
  wetLevel = newLevel;
  RAPT::rsEqualPowerGainFactors(dryWet, &dryFactor, &wetFactor, 0.0, 1.0);
  wetFactor *= RAPT::rsDbToAmp(wetLevel);
  mutex.unlock();
}

int EchoLab::addDelayLine(double newDelayTime, double newGainFactor)
{
  int result = -1;
  mutex.lock();

  EchoLabDelayLine* newDelayLine = new(std::nothrow) EchoLabDelayLine(524288);
  if( newDelayLine == NULL )
  {
    DEBUG_BREAK;  // memory allocation failed   
    result = (int) delayLines.size()-1;
    mutex.unlock();
    return result;
  }

  // setup the new delayline:
  newDelayLine->setDelayTime(newDelayTime);
  newDelayLine->setGlobalGainFactor(newGainFactor);
  newDelayLine->setSampleRate(sampleRate);
  newDelayLine->setTempoInBPM(bpm);
  newDelayLine->setSyncMode(delayTimeSync);
  //newDelayLine->setDelayModulationSyncMode(delayModulationSync);
  //newDelayLine->setAmplitudeModulationSyncMode(ampModulationSync);

  // add it to the vector and return the new number of delaylines:
  delayLines.push_back(newDelayLine);
  result = (int) delayLines.size()-1;
  mutex.unlock();
  return result;
}

bool EchoLab::removeDelayLine(int index)
{
  mutex.lock();

  if( index < 0 || index >= (int) delayLines.size() )
  {
    DEBUG_BREAK;
    mutex.unlock();
    return false;
  }
  delete delayLines[index];
  delayLines.erase(delayLines.begin()+index);

  mutex.unlock();
  return true;
}

bool EchoLab::removeAllDelayLines()
{
  mutex.lock();
  bool result = (delayLines.size() != 0);

  for(unsigned int i = 0; i < delayLines.size(); i++)
    delete delayLines[i];
  delayLines.clear();

  mutex.unlock();
  return result;
}

bool EchoLab::setDelayTime(int index, double newDelayTime)
{
  return setDelayLineParameterThreadSafe(DELAY_TIME, index, newDelayTime);
}

bool EchoLab::setGainFactor(int index, double newGainFactor)
{
  return setDelayLineParameterThreadSafe(GAIN_FACTOR, index, newGainFactor);
}

bool EchoLab::setDelayLineSolo(int index, bool shouldBeSolo)
{
  mutex.lock();
  if( index < 0 || index >= (int) delayLines.size() )
  {
    DEBUG_BREAK;
    soloedIndex = -1;
    mutex.unlock();
    return false;
  }

  if( shouldBeSolo == true )
    soloedIndex = index;
  else
    soloedIndex = -1;
  mutex.unlock();
  return true;
}

bool EchoLab::setDelayLineSolo(int index)
{
  mutex.lock();
  if( index < 0 || index >= (int) delayLines.size() )
  {
    soloedIndex = -1;
    mutex.unlock();
    return false;
  }
  soloedIndex = index;
  mutex.unlock();
  return true;
}

void EchoLab::setTempoInBPM(double newTempoInBPM)
{
  mutex.lock();
  if( newTempoInBPM > 0.0 )
  {
    bpm = newTempoInBPM;
    for(unsigned int i = 0; i < delayLines.size(); i++)
      delayLines[i]->setTempoInBPM(bpm);
  }
  mutex.unlock();
}

void EchoLab::setSyncForDelayTimes(bool shouldBeSynced)
{
  mutex.lock();  
  delayTimeSync = shouldBeSynced;
  for(unsigned int i = 0; i < delayLines.size(); i++)
    delayLines[i]->setSyncMode(delayTimeSync);
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int EchoLab::getNumDelayLines()
{
  int result = 0;
  mutex.lock();
  result = (int) delayLines.size();
  mutex.unlock();
  return result;
}

double EchoLab::getDelayTime(int index)
{
  return getDelayLineParameterThreadSafe(DELAY_TIME, index);
}

double EchoLab::getGainFactor(int index)
{
  return getDelayLineParameterThreadSafe(GAIN_FACTOR, index);
}

bool EchoLab::isDelayLineActive(int index)
{
  mutex.lock();
  if( index < 0 || index >= (int) delayLines.size() )
  {
    DEBUG_BREAK;
    mutex.unlock();
    return false;
  }
  bool m      = getDelayLineParameterThreadSafe(MUTE, index) != 0.0;     
  bool s      = soloedIndex != -1 && soloedIndex != index; 
                // some delayline is solo and it is not the one with 'index'
  bool result = !m && !s; // not muted and no other delayline is solo
  mutex.unlock();
  return result;
}

EchoLabDelayLine* EchoLab::getDelayLine(int index)
{
  if( index < 0 || index >= (int) delayLines.size() )
    return NULL;
  else
    return delayLines[index];
  // we don't wrap that into mutex-locks, because thread-safety has to be ensured by the caller 
  // anyway
}

EqualizerStereo* EchoLab::getFeedbackEqualizer(int index)
{
  if( index < 0 || index >= (int) delayLines.size() )
    return NULL;
  else
    return &(delayLines[index]->feedbackEqualizer);
}

EqualizerStereo* EchoLab::getInputEqualizer(int index)
{
  if( index < 0 || index >= (int) delayLines.size() )
    return NULL;
  else
    return &(delayLines[index]->inputEqualizer);
}

//-------------------------------------------------------------------------------------------------
// others:

void EchoLab::resetDelayLines()
{
  mutex.lock();
  for(unsigned int s = 0; s < delayLines.size(); s++)
    delayLines[s]->clearBuffers();
  mutex.unlock();
}

bool EchoLab::setDelayLineParameterThreadSafe(int parameterIndex, int delayLineIndex, double newValue)
{
  mutex.lock();
  if( delayLineIndex >= 0 && delayLineIndex < (int) delayLines.size() )
  {
    switch( parameterIndex )
    {
    case DELAY_TIME:         
      delayLines[delayLineIndex]->setDelayTime(RAPT::rsMax(0.01,newValue));          
      break;
    case GAIN_FACTOR:        
      delayLines[delayLineIndex]->setGlobalGainFactor(newValue);   
      break;
    case FEEDBACK_FACTOR:    
      delayLines[delayLineIndex]->setFeedbackFactor(newValue);     
      break;
      /*
    case FEEDFORWARD_FACTOR: 
      delayLines[delayLineIndex]->setFeedforwardFactor(newValue);  
      break;
    case BLEND_FACTOR:       
      delayLines[delayLineIndex]->setBlendFactor(newValue);        
      break;
      */
    case PAN:                
      delayLines[delayLineIndex]->setPan(newValue);                
      break;
    case PING_PONG:                
      delayLines[delayLineIndex]->setPingPongMode(newValue >= 0.5);                
      break;
    case MUTE:                
      delayLines[delayLineIndex]->setMute(newValue >= 0.5);                
      break;
      /*
    case DELAY_MOD_CYCLE:    
      delayLines[delayLineIndex]->setDelayModulationCycleLength(newValue); 
      break;
    case DELAY_MOD_DEPTH:    
      delayLines[delayLineIndex]->setDelayModulationDepth(newValue); 
      break;
    case DELAY_MOD_PHASE_L:  
      delayLines[delayLineIndex]->setDelayModulationPhaseLeft(newValue); 
      break;
    case DELAY_MOD_PHASE_R:  
      delayLines[delayLineIndex]->setDelayModulationPhaseRight(newValue); 
      break;
    case AMP_MOD_CYCLE:    
      delayLines[delayLineIndex]->setAmplitudeModulationCycleLength(newValue); 
      break;
    case AMP_MOD_DEPTH:    
      delayLines[delayLineIndex]->setAmplitudeModulationDepth(newValue); 
      break;
    case AMP_MOD_PHASE_L:  
      delayLines[delayLineIndex]->setAmplitudeModulationPhaseLeft(newValue); 
      break;
    case AMP_MOD_PHASE_R:  
      delayLines[delayLineIndex]->setAmplitudeModulationPhaseRight(newValue); 
      break;
      */
    }
    mutex.unlock();
    return true;
  }
  else
  {
    DEBUG_BREAK;
    mutex.unlock();
    return false;
  }
}

double EchoLab::getDelayLineParameterThreadSafe(int parameterIndex, int delayLineIndex)
{
  mutex.lock();
  if( delayLineIndex < 0 || delayLineIndex >= (int) delayLines.size() )
  {
    DEBUG_BREAK;
    mutex.unlock();
    return 0.0;
  }
  double result = 0.0; 
  switch( parameterIndex )
  {
  case DELAY_TIME:         
    result = delayLines[delayLineIndex]->getDelayTime();          
    break;
  case GAIN_FACTOR:        
    result = delayLines[delayLineIndex]->getGlobalGainFactor();   
    break;
  case FEEDBACK_FACTOR:    
    result = delayLines[delayLineIndex]->getFeedbackFactor();     
    break;
    /*
  case FEEDFORWARD_FACTOR: 
    result = delayLines[delayLineIndex]->getFeedforwardFactor();  
    break;
  case BLEND_FACTOR:       
    result = delayLines[delayLineIndex]->getBlendFactor();        
    break;
    */
  case PAN:                
    result = delayLines[delayLineIndex]->getPan();                
    break;
  case PING_PONG:                
    result = (double) delayLines[delayLineIndex]->isInPingPongMode();                
    break;
  case MUTE:                
    result = (double) delayLines[delayLineIndex]->isMuted();                
    break;
    /*
  case DELAY_MOD_CYCLE:    
    result = delayLines[delayLineIndex]->getDelayModulationCycleLength(); 
    break;
  case DELAY_MOD_DEPTH:    
    result = delayLines[delayLineIndex]->getDelayModulationDepth(); 
    break;
  case DELAY_MOD_PHASE_L:  
    result = delayLines[delayLineIndex]->getDelayModulationPhaseLeft(); 
    break;
  case DELAY_MOD_PHASE_R:  
    result = delayLines[delayLineIndex]->getDelayModulationPhaseRight(); 
    break;
  case AMP_MOD_CYCLE:    
    result = delayLines[delayLineIndex]->getAmplitudeModulationCycleLength(); 
    break;
  case AMP_MOD_DEPTH:    
    result = delayLines[delayLineIndex]->getAmplitudeModulationDepth(); 
    break;
  case AMP_MOD_PHASE_L:  
    result = delayLines[delayLineIndex]->getAmplitudeModulationPhaseLeft(); 
    break;
  case AMP_MOD_PHASE_R:  
    result = delayLines[delayLineIndex]->getAmplitudeModulationPhaseRight(); 
    break;
    */
  }
  mutex.unlock();
  return result;
}

/*

Ideas:
-Instead of allowing only feedback of a delayline to itself, allow crossfeedback between all 
 delaylines. Use a (sparse) feedback matrix of coeffs and filters.

*/


