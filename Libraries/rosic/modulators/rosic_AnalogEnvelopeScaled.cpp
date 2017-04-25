#include "rosic_AnalogEnvelopeScaled.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AnalogEnvelopeScaled::AnalogEnvelopeScaled()
{
  mu             = 0.001;
  shape          = 1.0;
  sampleCounter  = 0;
  
  holdSamples    = 0;

  attackSamples  = 0;   
  attackAccuMin  = 0.0;
  attackAccuMax  = 1.0;
  attackAccu     = 1.0;
  attackCoeff    = 0.0;

  decaySamples   = 0;   
  decayAccuMin   = 0.0;
  decayAccuMax   = 1.0;
  decayAccu      = 1.0;
  decayCoeff     = 0.0;

  releaseSamples = 0;   
  releaseAccuMin = 0.0;
  releaseAccuMax = 1.0;
  releaseAccu    = 1.0;
  releaseCoeff   = 0.0;
  

  setAttack( 10.0);
  setDecay(  50.0);
  setRelease(20.0);

  // init envelope as if it has run to its end:
  sampleCounter = attackSamples+holdSamples+decaySamples+releaseSamples;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AnalogEnvelopeScaled::setAttack(double newAttackTime)
{
  if( newAttackTime > 0.0 )
  {
    attackTime    = newAttackTime;
    attackSamples = roundToInt(0.001*attackTime*sampleRate*timeScale); 
    attackAccuMin = mu;
    attackAccuMax = pow(attackAccuMin, (1.0/(attackSamples+1)));
    attackCoeff   = (peakScale*peakLevel-startLevel) / (attackAccuMax-attackAccuMin);  // scaleFactor
    attackAccu    = attackAccuMax;
  }
  else // newAttackTime <= 0
  {
    attackTime    = 0.0;
    attackSamples = 0;
    attackAccuMin = 0.0;
    attackAccuMax = 1.0;
    attackCoeff   = 0.0;
    attackAccu    = 1.0;
  }
}

void AnalogEnvelopeScaled::setHold(double newHoldTime)
{
  if( newHoldTime >= 0 ) 
  {
    holdTime    = newHoldTime;
    holdSamples = roundToInt(0.001*holdTime*sampleRate*timeScale); 
  }
}

void AnalogEnvelopeScaled::setDecay(double newDecayTime)
{
  if( newDecayTime > 0.0 )
  {
    decayTime    = newDecayTime;
    decaySamples = roundToInt(0.001*decayTime*sampleRate*timeScale); 
    decayAccuMin = mu;
    decayAccuMax = pow(decayAccuMin, (1.0/(decaySamples+1)));
    decayCoeff   = (sustainLevel-peakScale*peakLevel) / (decayAccuMax-decayAccuMin);  // scaleFactor
    decayAccu    = decayAccuMax;
  }
  else // newDecayTime <= 0
  {
    decayTime    = 0.0;
    decaySamples = 0;
    decayAccuMin = 0.0;
    decayAccuMax = 1.0;
    decayCoeff   = 0.0;
    decayAccu    = 1.0;
  }
}

void AnalogEnvelopeScaled::setRelease(double newReleaseTime)
{
  if( newReleaseTime > 0.0 )
  {
    releaseTime    = newReleaseTime;
    releaseSamples = roundToInt(0.001*releaseTime*sampleRate*timeScale); 
    releaseAccuMin = mu;
    releaseAccuMax = pow(releaseAccuMin, (1.0/(releaseSamples+1)));
    releaseCoeff   = (endLevel-sustainLevel) / (releaseAccuMax-releaseAccuMin);  // scaleFactor
    releaseAccu    = releaseAccuMax;
  }
  else // newReleaseTime <= 0
  {
    releaseTime    = 0.0;
    releaseSamples = 0;
    releaseAccuMin = 0.0;
    releaseAccuMax = 1.0;
    releaseCoeff   = 0.0;
    releaseAccu    = 1.0;
  }
}

void AnalogEnvelopeScaled::setTimeScale(double newTimeScale)
{
  if( newTimeScale > 0 )
    timeScale = newTimeScale;

  increment  = 1000.0*timeScale/sampleRate;

  // re-calculate coefficients for the 3 filters:
  setAttack (attackTime);
  setHold   (holdTime);
  setDecay  (decayTime);
  setRelease(releaseTime);
}

void AnalogEnvelopeScaled::setShape(double newShape)
{
  if( newShape > 0.0 )
  {
    shape = newShape;
    mu    = pow(0.001, shape);
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void AnalogEnvelopeScaled::reset()
{
  time          = 0.0;
  sampleCounter = 0;
}

void AnalogEnvelopeScaled::noteOn(bool startFromCurrentLevel, int newKey, int newVel)
{
  // \todo: if( !startFromCurrentLevel ) ...
  // \todo: calculate key and velocity scale factors for duration and peak-value...

  if( startFromCurrentLevel == true )
  {
    //leftLevel   = previousOutput;
    attackCoeff = (peakScale*peakLevel-previousOutput) / (attackAccuMax-attackAccuMin); 
  }
  else
  {
    //leftLevel   = startLevel;
    attackCoeff = (peakScale*peakLevel-startLevel) / (attackAccuMax-attackAccuMin);  
  }


  // reset time for the new note:
  sampleCounter = 0;
  time          = 0.0;   
  noteIsOn      = true;
  outputIsZero  = false;
}

void AnalogEnvelopeScaled::noteOff()
{
  noteIsOn = false; 

  //leftLevel    = previousOutput;
  releaseCoeff = (endLevel-previousOutput) / (releaseAccuMax-releaseAccuMin); 


  // advance time to the beginnig of the release phase:
  time = (attackTime + holdTime + decayTime + increment);
  sampleCounter = attackSamples+holdSamples+decaySamples;
}

bool AnalogEnvelopeScaled::endIsReached()
{
  if( sampleCounter >= attackSamples+holdSamples+decaySamples+releaseSamples )
    return true;
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

