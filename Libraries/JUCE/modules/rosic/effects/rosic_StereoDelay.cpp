#include "rosic_StereoDelay.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

StereoDelay::StereoDelay()
{
  inL2L                = 1.0;
  inL2R                = 0.0;
  inR2L                = 0.0;
  inR2R                = 1.0;
  fbL2L                = 0.0;
  fbL2R                = 0.0;
  fbR2L                = 0.0;
  fbR2R                = 0.0;
  outL2L               = 1.0;
  outL2R               = 0.0;
  outR2L               = 0.0;
  outR2R               = 1.0;
  dryGain              = sqrt(0.5);
  wetGain              = sqrt(0.5);
  sampleRate           = 44100.0;
  bpm                  = 120.0;
  delayInBeatsL        = 1.0;
  delayInBeatsR        = 1.0;
  delayScaleL          = 1.0;
  delayScaleR          = 1.0;
  lowpassCutoffL       = 20000.0;
  highpassCutoffL      = 20.0;
  lowpassCutoffR       = 20000.0;
  highpassCutoffR      = 20.0;
  cutoffScale          = 1.0;
  wetDelayInBeatsL     = 0.0;
  wetDelayInBeatsR     = 0.0;
  dryWet               = 50.0;

  updateDelayTimes();
  updateFilterFrequencies();
}

StereoDelay::~StereoDelay()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void StereoDelay::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;

    delayLineL.setSampleRate(sampleRate);
    delayLineR.setSampleRate(sampleRate);
    diffusorL.setSampleRate(sampleRate);
    diffusorR.setSampleRate(sampleRate);
    filterL.setSampleRate(sampleRate);
    filterR.setSampleRate(sampleRate);
    wetDelayLineL.setSampleRate(sampleRate);
    wetDelayLineR.setSampleRate(sampleRate);
  }
  else
    DEBUG_BREAK; // invalid sample-rate (must be > 0.0)
}

void StereoDelay::setBeatsPerMinute(double newBpm)
{
  if( newBpm > 0.0 )
  {
    bpm = newBpm;
    updateDelayTimes();
  }
  else
  {
    //DEBUG_BREAK; // invalid tempo (must be > 0.0)
  }
}

void StereoDelay::setDelayInBeats(double newDelay, int delayline)
{
  if( newDelay > 0.0 )
  {
    if( delayline == LEFT )
      delayInBeatsL = newDelay;
    else if( delayline == RIGHT )
      delayInBeatsR = newDelay;
    updateDelayTimes();
  }
  else
    DEBUG_BREAK; // invalid delay-time (must be > 0.0)
}

void StereoDelay::setDelayScale(double newFactor, int delayline)
{
  if( newFactor > 0.0 )
  {
    if( delayline == LEFT )
      delayScaleL = newFactor;
    else if( delayline == RIGHT )
      delayScaleR = newFactor;
    updateDelayTimes();
  }
  else
    DEBUG_BREAK; // invalid factor (must be > 0.0)
}

void StereoDelay::setInjection(double newPercentage, int sourceChannel, int delayline)
{
  if( sourceChannel == LEFT && delayline == LEFT )
    inL2L = 0.01*newPercentage;
  if( sourceChannel == LEFT && delayline == RIGHT )
    inL2R = 0.01*newPercentage;
  if( sourceChannel == RIGHT && delayline == LEFT )
    inR2L = 0.01*newPercentage;
  if( sourceChannel == RIGHT && delayline == RIGHT )
    inR2R = 0.01*newPercentage;
}

void StereoDelay::setDiffusorTimeInMilliseconds(double newDiffusorTime, int delayline)
{
  if( delayline == LEFT )
    diffusorL.setDelayInMilliseconds(newDiffusorTime);
  else if( delayline == RIGHT )
    diffusorR.setDelayInMilliseconds(newDiffusorTime);  
  updateDelayTimes();
}

void StereoDelay::setDiffusorAmount(double newDiffusorAmount, int delayline)
{
  if( delayline == LEFT )
    diffusorL.setDiffusionAmount(newDiffusorAmount);
  else if( delayline == RIGHT )
    diffusorR.setDiffusionAmount(newDiffusorAmount);
}

void StereoDelay::setLowpassCutoff(double newCutoff, int delayline)
{
  if( newCutoff >= 20.0 && newCutoff <= 20000.0 )
  {
    if( delayline == LEFT )
      lowpassCutoffL = newCutoff;
    else if( delayline == RIGHT )
      lowpassCutoffR = newCutoff;
    updateFilterFrequencies();
  }
}

void StereoDelay::setHighpassCutoff(double newCutoff, int delayline)
{
  if( newCutoff >= 20.0 && newCutoff <= 20000.0 )
  {
    if( delayline == LEFT )
      highpassCutoffL = newCutoff;
    else if( delayline == RIGHT )
      highpassCutoffR = newCutoff;
    updateFilterFrequencies();
  }
}

void StereoDelay::setCutoffScale(double newCutoffScale)
{
  if( newCutoffScale > 0.0 )
  {
    cutoffScale = newCutoffScale;
    updateFilterFrequencies();
  }
  else
    DEBUG_BREAK; // invalid factor (must be > 0.0)
}

void StereoDelay::setFeedback(double newFeedback, int source, int target)
{
  if( source == LEFT && target == LEFT )
    fbL2L = 0.01*newFeedback;
  if( source == LEFT && target == RIGHT )
    fbL2R = 0.01*newFeedback;
  if( source == RIGHT && target == LEFT )
    fbR2L = 0.01*newFeedback;
  if( source == RIGHT && target == RIGHT )
    fbR2R = 0.01*newFeedback;
}

void StereoDelay::setOutputMix(double newPercentage, int delayline, int targetChannel)
{
  if( delayline == LEFT && targetChannel == LEFT )
    outL2L = 0.01*newPercentage;
  if( delayline == LEFT && targetChannel == RIGHT )
    outL2R = 0.01*newPercentage;
  if( delayline == RIGHT && targetChannel == LEFT )
    outR2L = 0.01*newPercentage;
  if( delayline == RIGHT && targetChannel == RIGHT )
    outR2R = 0.01*newPercentage;
}

void StereoDelay::setWetDelayInBeats(double newDelay, int channel)
{
  if( newDelay >= 0.0 )
  {
    if( channel == LEFT )
      wetDelayInBeatsL = newDelay;
    else if( channel == RIGHT )
      wetDelayInBeatsR = newDelay;
    updateDelayTimes();
  }
  else
    DEBUG_BREAK; // invalid delay-time (must be > 0.0)
}

void StereoDelay::setDryWet(double newDryWet)
{
  if( newDryWet >= 0.0 && newDryWet <= 100.0 )
  {
    dryWet = newDryWet;
    equalPowerGainFactors(dryWet, &dryGain, &wetGain, 0.0, 100.0);
    //wet = 0.01 * newDryWet;
    //dry = 1.0 - wet;
  }
  else
    DEBUG_BREAK; // invalid value (must be >= 0.0 and <= 100.0)
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double StereoDelay::getDelayInBeats(int delayline)
{
  if( delayline == LEFT )
    return delayInBeatsL;
  else if( delayline == RIGHT )
    return delayInBeatsR;
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getDelayScale(int delayline)
{
  if( delayline == LEFT )
    return delayScaleL;
  else if( delayline == RIGHT )
    return delayScaleR;
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getInjection(int sourceChannel, int delayline)
{
  if( sourceChannel == LEFT && delayline == LEFT )
    return 100.0*inL2L;
  if( sourceChannel == LEFT && delayline == RIGHT )
    return 100.0*inL2R;
  if( sourceChannel == RIGHT && delayline == LEFT )
    return 100.0*inR2L;
  if( sourceChannel == RIGHT && delayline == RIGHT )
    return 100.0*inR2R;
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getDiffusorTimeInMilliseconds(int delayline)
{
  if( delayline == LEFT )
    return diffusorL.getDelayInMilliseconds();
  else if( delayline == RIGHT )
    return diffusorR.getDelayInMilliseconds();
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getDiffusorAmount(int delayline)
{
  if( delayline == LEFT )
    return diffusorL.getDiffusionAmount();
  else if( delayline == RIGHT )
    return diffusorR.getDiffusionAmount();
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getLowpassCutoff(int delayline)
{
  if( delayline == LEFT )
    return filterL.getLowpassCutoff();
  else if( delayline == RIGHT )
    return filterR.getLowpassCutoff();
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getHighpassCutoff(int delayline)
{
  if( delayline == LEFT )
    return filterL.getHighpassCutoff();
  else if( delayline == RIGHT )
    return filterR.getHighpassCutoff();
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getFeedback(int source, int target)
{
  if( source == LEFT && target == LEFT )
    return 100.0*fbL2L;
  else if( source == LEFT && target == RIGHT )
    return 100.0*fbL2R;
  else if( source == RIGHT && target == LEFT )
    return 100.0*fbR2L;
  else if( source == RIGHT && target == RIGHT )
    return 100.0*fbR2R;
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getOutputMix(int delayline, int targetChannel)
{
  if( delayline == LEFT && targetChannel == LEFT )
    return 100.0*outL2L;
  if( delayline == LEFT && targetChannel == RIGHT )
    return 100.0*outL2R;
  if( delayline == RIGHT && targetChannel == LEFT )
    return 100.0*outR2L;
  if( delayline == RIGHT && targetChannel == RIGHT )
    return 100.0*outR2R;
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getWetDelayInBeats(int channel)
{
  if( channel == LEFT )
    return wetDelayInBeatsL;
  else if( channel == RIGHT )
    return wetDelayInBeatsR;
  else
  {
    DEBUG_BREAK; // invalid delayline index
    return 0.0;
  }
}

double StereoDelay::getCutoffScale()
{
  return cutoffScale;
}

double StereoDelay::getDryWet()
{
  return dryWet;
}

//-------------------------------------------------------------------------------------------------
// others:

void StereoDelay::reset()
{
  //yL = 0.0;
  //yR = 0.0;
  delayLineL.clearDelayBuffer();
  delayLineR.clearDelayBuffer();
  diffusorL.clearDelayBuffer();  
  diffusorR.clearDelayBuffer();
  filterL.reset();
  filterR.reset();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void StereoDelay::updateDelayTimes()
{
  // calculate the desired total delay-time for the left and right delayline which includes the 
  // delay introduced by the allpass diffusors:
  double dL = delayScaleL*beatsToSeconds(delayInBeatsL, bpm);
  double dR = delayScaleR*beatsToSeconds(delayInBeatsR, bpm);

  // subtract the delay from the allpass diffusors to obtain the desired delay for the delay-lines:
  dL -= diffusorL.getDelayInSeconds();
  dR -= diffusorR.getDelayInSeconds();

  // setup the delaylines:
  delayLineL.setDelayInSeconds(dL);
  delayLineR.setDelayInSeconds(dR);

  // calculate and set the desired delay-time for the left and right wet signal delaylines:
  dL = beatsToSeconds(wetDelayInBeatsL, bpm);
  dR = beatsToSeconds(wetDelayInBeatsR, bpm);
  wetDelayLineL.setDelayInSeconds(dL);
  wetDelayLineR.setDelayInSeconds(dR);
}

void StereoDelay::updateFilterFrequencies()
{
  filterL.setLowpassCutoff(cutoffScale*lowpassCutoffL);
  filterL.setHighpassCutoff(cutoffScale*highpassCutoffL);
  filterR.setLowpassCutoff(cutoffScale*lowpassCutoffR);
  filterR.setHighpassCutoff(cutoffScale*highpassCutoffR);
}