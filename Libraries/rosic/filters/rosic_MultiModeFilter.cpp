#include "rosic_MultiModeFilter.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiModeFilter::MultiModeFilter()
{
  parameters        = new MultiModeFilterParameters;
  isMaster          = true;

  freqWithKeyAndVel = 1000.0;      
  freqInstantaneous = 1000.0;
  currentKey        = 64.0;
  currentVel        = 64.0;
  //glideMode         = false;
  //glideTime         = 0.1;

  resetBuffers();

  //mode              = MultiModeFilterParameters::BYPASS


  //twoStageBiquad.setMode(FourPoleFilter::LOWPASS_RBJ);
  //twoStageBiquad.setMode(FourPoleFilter::MORPH_LP_RES_HP);
  //setMode(MultiModeFilterParameters::BYPASS);

  //calculateCoefficients();

  //markPresetAsClean();
}

MultiModeFilter::~MultiModeFilter()
{
  if( isMaster && parameters != NULL )
    delete parameters;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void MultiModeFilter::setSampleRate(double newSampleRate)
{
  ladderFilter.setSampleRate(newSampleRate);
  twoStageBiquad.setSampleRate(newSampleRate);
}

void MultiModeFilter::setMode(int newMode)
{
  if( newMode >= MultiModeFilterParameters::BYPASS 
    && newMode < MultiModeFilterParameters::NUM_FILTER_MODES )
  {
    parameters->mode = newMode;
    if( parameters->mode == MultiModeFilterParameters::BYPASS )
      parameters->filterClass = MultiModeFilterParameters::NO_FILTER;
    else if( parameters->mode == MultiModeFilterParameters::MOOGISH_LOWPASS )
      parameters->filterClass = MultiModeFilterParameters::LADDER_FILTER;
    else
    {
      parameters->filterClass = MultiModeFilterParameters::TWO_STAGE_BIQUAD;
      switch( parameters->mode )
      {
      case MultiModeFilterParameters::LOWPASS_6:  
        twoStageBiquad.setMode(FourPoleFilterParameters::LOWPASS_6);  
        break;
      case MultiModeFilterParameters::LOWPASS_RBJ:  
        twoStageBiquad.setMode(FourPoleFilterParameters::LOWPASS_12);  
        break;
      case MultiModeFilterParameters::HIGHPASS_6: 
        twoStageBiquad.setMode(FourPoleFilterParameters::HIGHPASS_6); 
        break;
      case MultiModeFilterParameters::HIGHPASS_RBJ: 
        twoStageBiquad.setMode(FourPoleFilterParameters::HIGHPASS_12); 
        break;
      case MultiModeFilterParameters::BANDPASS_RBJ: 
        twoStageBiquad.setMode(FourPoleFilterParameters::BANDPASS_RBJ); 
        break;
      case MultiModeFilterParameters::BANDREJECT_RBJ: 
        twoStageBiquad.setMode(FourPoleFilterParameters::BANDREJECT_RBJ); 
        break;
      case MultiModeFilterParameters::PEAK_OR_DIP_RBJ:
        twoStageBiquad.setMode(FourPoleFilterParameters::PEAK_OR_DIP_RBJ ); 
        break;
      case MultiModeFilterParameters::LOW_SHELV_1ST:
        twoStageBiquad.setMode(FourPoleFilterParameters::LOW_SHELV_1ST ); 
        break;
      case MultiModeFilterParameters::LOW_SHELV_RBJ:
        twoStageBiquad.setMode(FourPoleFilterParameters::LOW_SHELV_RBJ ); 
        break;
      case MultiModeFilterParameters::HIGH_SHELV_1ST:
        twoStageBiquad.setMode(FourPoleFilterParameters::HIGH_SHELV_1ST ); 
        break;
      case MultiModeFilterParameters::HIGH_SHELV_RBJ:
        twoStageBiquad.setMode(FourPoleFilterParameters::HIGH_SHELV_RBJ ); 
        break;
      case MultiModeFilterParameters::ALLPASS_1ST:
        twoStageBiquad.setMode(FourPoleFilterParameters::ALLPASS_1ST ); 
        break;
      case MultiModeFilterParameters::ALLPASS_RBJ:
        twoStageBiquad.setMode(FourPoleFilterParameters::ALLPASS_RBJ ); 
        break;


      case MultiModeFilterParameters::MORPH_LP_BP_HP: 
        twoStageBiquad.setMode(FourPoleFilterParameters::MORPH_LP_BP_HP); 
        break;
      case MultiModeFilterParameters::MORPH_LP_PK_HP: 
        twoStageBiquad.setMode(FourPoleFilterParameters::MORPH_LP_PK_HP); 
        break;

      }
    }
    calculateCoefficients();
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->setMode(newMode);
    //markPresetAsDirty();
  }
  else
    DEBUG_BREAK; // invalid filter mode index
}

void MultiModeFilter::setFrequencyNominal(double newFrequency)
{
  parameters->freqNominal = newFrequency;
  updateFreqWithKeyAndVel();

  if( isMaster == true )
  {
    setFrequencyInstantaneous(newFrequency, true); 
    for(unsigned int s = 0; s < slaves.size(); s++)
      slaves[s]->updateFreqWithKeyAndVel();
  }
  //markPresetAsDirty();

  /*
  if( updateCoefficients == true )
  {
    calculateCoefficients();
    for(int s=0; s<slaves.size(); s++)
      slaves[s]->calculateCoefficients();
  }
  */
}

void MultiModeFilter::setFrequencyByKey(double newFrequencyByKey)
{
  parameters->freqByKey = newFrequencyByKey;
  updateFreqWithKeyAndVel();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->updateFreqWithKeyAndVel();
  //markPresetAsDirty();
}

void MultiModeFilter::setFrequencyByVel(double newFrequencyByVel)
{
  parameters->freqByVel = newFrequencyByVel;
  updateFreqWithKeyAndVel();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->updateFreqWithKeyAndVel();
  //markPresetAsDirty();
}

void MultiModeFilter::setKey(double newKey)
{ 
  currentKey = newKey;
  updateFreqWithKeyAndVel();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->updateFreqWithKeyAndVel();
}

void MultiModeFilter::setKeyAndVel(double newKey, double newVel)
{ 
  currentKey = newKey;
  currentVel = newVel;
  updateFreqWithKeyAndVel();
  for(unsigned int s = 0; s < slaves.size(); s++)
    slaves[s]->updateFreqWithKeyAndVel();
}

void MultiModeFilter::setResonance(double newResonance, bool updateCoefficients)
{
  parameters->resonance = newResonance;
  ladderFilter.setResonance(0.01*parameters->resonance, updateCoefficients); 
  /*
  twoStageBiquad.setResonance(0.01*parameters->resonance);  
  if( updateCoefficients == true )
    twoStageBiquad.updateFilterCoefficients();
  */
  //markPresetAsDirty();
}

void MultiModeFilter::setQ(double newQ, bool updateCoefficients)
{
  twoStageBiquad.setQ(newQ);  
  if( updateCoefficients == true )
    twoStageBiquad.updateFilterCoefficients();
  //markPresetAsDirty();
}

void MultiModeFilter::setGain(double newGain)
{
  twoStageBiquad.setGain(newGain);
  //markPresetAsDirty();
}

void MultiModeFilter::setDrive(double newDrive)
{
  ladderFilter.setDrive(newDrive);
  //twoStageBiquad.setDrive(newDrive);
  //markPresetAsDirty();
}

void MultiModeFilter::setOrder(int newOrder)
{
  parameters->order = newOrder;
  //twoStageBiquad.setOrder(newOrder);
  ladderFilter.setOutputStage(newOrder);
  //markPresetAsDirty();
}

void MultiModeFilter::useTwoStages(bool shouldUseTwoStages)
{
  twoStageBiquad.useTwoStages(shouldUseTwoStages);
  //markPresetAsDirty();
}

void MultiModeFilter::setAllpassFreq(double newAllpassFreq)
{
  ladderFilter.setAllpassFreq(newAllpassFreq);
  //markPresetAsDirty();
}

void MultiModeFilter::setMakeUp(double newMakeUp)
{
  ladderFilter.setMakeUp(newMakeUp, true);
  //markPresetAsDirty();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

Complex MultiModeFilter::getTransferFunctionAt(Complex z)
{
  ///< \todo switch or ifs
  return ladderFilter.getTransferFunctionAt(z, true, true, true, ladderFilter.getOutputStage());
}

double MultiModeFilter::getMagnitudeAt(double frequency)
{
  switch( parameters->filterClass )
  {
  case MultiModeFilterParameters::LADDER_FILTER: 
    return ladderFilter.getMagnitudeAt(frequency, true, true, true, ladderFilter.getOutputStage());
  case MultiModeFilterParameters::TWO_STAGE_BIQUAD:
    return twoStageBiquad.getMagnitudeAt(frequency);
  default:
    return 1.0;
  }
  /*
  if( parameters->mode == BYPASS )
    return 1.0;
  else if( parameters->mode == MOOGISH_LOWPASS )
    return ladderFilter.getMagnitudeAt(frequency, true, true, true, ladderFilter.getOutputStage());
  else
    return twoStageBiquad.getMagnitudeAt(frequency);
  */
}

void MultiModeFilter::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
                                           bool inDecibels, bool accumulate)
{
  if( parameters->filterClass == MultiModeFilterParameters::NO_FILTER )
  {
    // fill magnitude-array with bypass-rspeonse:
    int k;
    if( accumulate == false && inDecibels == true )
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] = 0.0;
    }
    else if( accumulate == false && inDecibels == false )
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] = 1.0;
    }
    return;
  }
  /*
  else if( parameters->mode == MOOGISH_LOWPASS )
    ladderFilter.getMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, accumulate);
  else
    twoStageBiquad.getMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, accumulate);
  */
  switch( parameters->filterClass )
  {
  case MultiModeFilterParameters::LADDER_FILTER:    
    ladderFilter.getMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, accumulate);
    break;
  case MultiModeFilterParameters::TWO_STAGE_BIQUAD:
    twoStageBiquad.getMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, accumulate);
    break;
  }
}

int MultiModeFilter::getMode()
{
  return parameters->mode;
}

double MultiModeFilter::getFrequencyNominal()
{ 
  return parameters->freqNominal;
}

double MultiModeFilter::getFrequencyByKey()
{
  return parameters->freqByKey;
}

double MultiModeFilter::getFrequencyByVel()
{
  return parameters->freqByVel;
}

double MultiModeFilter::getFrequencyWithKeyAndVel()
{
  return freqWithKeyAndVel;
}

double MultiModeFilter::getResonance()
{
  return parameters->resonance;
}

double MultiModeFilter::getQ()
{
  return twoStageBiquad.getQ();;
}

double MultiModeFilter::getGain()
{
  return twoStageBiquad.getGain();
}

double MultiModeFilter::getDrive()
{
  return ladderFilter.getDrive();
}

int MultiModeFilter::getOrder()
{
  return parameters->order;
}

bool MultiModeFilter::usesTwoStages()
{
  return twoStageBiquad.usesTwoStages();
}

double MultiModeFilter::getMorph()
{ 
  return twoStageBiquad.getMorph();
}

double MultiModeFilter::getAllpassFreq()
{
  return ladderFilter.getAllpassFreq();
}

double MultiModeFilter::getMakeUp()
{
  return ladderFilter.getMakeUp();
}

bool MultiModeFilter::currentModeSupportsQ()
{
  if(  parameters->mode == MultiModeFilterParameters::LOWPASS_RBJ
    || parameters->mode == MultiModeFilterParameters::HIGHPASS_RBJ 
    || parameters->mode == MultiModeFilterParameters::BANDPASS_RBJ
    || parameters->mode == MultiModeFilterParameters::BANDREJECT_RBJ  
    || parameters->mode == MultiModeFilterParameters::ALLPASS_RBJ
    || parameters->mode == MultiModeFilterParameters::PEAK_OR_DIP_RBJ
    || parameters->mode == MultiModeFilterParameters::LOW_SHELV_RBJ 
    || parameters->mode == MultiModeFilterParameters::HIGH_SHELV_RBJ   
    || parameters->mode == MultiModeFilterParameters::MORPH_LP_PK_HP )
  {
    return true;
  }
  else 
    return false;
}

bool MultiModeFilter::currentModeSupportsGain()
{
  if(  parameters->mode == MultiModeFilterParameters::PEAK_OR_DIP_RBJ
    || parameters->mode == MultiModeFilterParameters::LOW_SHELV_1ST 
    || parameters->mode == MultiModeFilterParameters::LOW_SHELV_RBJ 
    || parameters->mode == MultiModeFilterParameters::HIGH_SHELV_1ST  
    || parameters->mode == MultiModeFilterParameters::HIGH_SHELV_RBJ   )
  {
    return true;
  }
  else 
    return false;
}

bool MultiModeFilter::currentModeSupportsTwoStages()
{
  if( parameters->filterClass == MultiModeFilterParameters::TWO_STAGE_BIQUAD )
    return true;
  else
    return false;
  /*
  if(  parameters->mode == MultiModeFilterParameters::LOWPASS_6
    || parameters->mode == MultiModeFilterParameters::LOWPASS_RBJ
    || parameters->mode == MultiModeFilterParameters::HIGHPASS_6
    || parameters->mode == MultiModeFilterParameters::HIGHPASS_RBJ 
    || parameters->mode == MultiModeFilterParameters::BANDPASS_RBJ
    || parameters->mode == MultiModeFilterParameters::BANDREJECT_RBJ  
    || parameters->mode == MultiModeFilterParameters::ALLPASS_1ST 
    || parameters->mode == MultiModeFilterParameters::ALLPASS_RBJ 
    || parameters->mode == MultiModeFilterParameters::MORPH_LP_PK_HP )
  {
    return true;
  }
  else 
    return false;
  */
}

//-------------------------------------------------------------------------------------------------
// master/slave config:

void MultiModeFilter::addSlave(MultiModeFilter* newSlave)
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

  // add the embedded filters in the newSlave as slaves to the respective embedded filter here in 
  // this instance:
  ladderFilter.addSlave( &(newSlave->ladderFilter) );
  twoStageBiquad.addSlave( &(newSlave->twoStageBiquad) );


  // set the isMaster-flag of the new slave to false: 
  newSlave->isMaster = false;

  // this flag will prevent the destructor of the slave from trying to delete the parameter-set 
  // which is now shared - only masters delete their parameter-set on destruction
}

//-------------------------------------------------------------------------------------------------
// others:

void MultiModeFilter::updateFreqWithKeyAndVel()
{
  freqWithKeyAndVel  = parameters->freqNominal;
  freqWithKeyAndVel *= pow(2.0, (0.01*parameters->freqByKey/12.0) * (currentKey-64.0) );
  freqWithKeyAndVel *= pow(2.0, (0.01*parameters->freqByVel/63.0) * (currentVel-64.0) );
}

void MultiModeFilter::resetBuffers()
{
  ladderFilter.reset();
  twoStageBiquad.reset();
}