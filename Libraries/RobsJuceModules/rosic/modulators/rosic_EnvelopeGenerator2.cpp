//#include "rosic_EnvelopeGenerator2.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

EnvelopeGenerator2::EnvelopeGenerator2()
{
  sampleRate     = 44100.0;
  sampleRateRec  = 1.0/sampleRate;
  startLevel     = 0.0;
  attackTime     = 0.0;
  peakLevel      = 1.0;
  holdTime       = 0.0;
  decayTime      = 0.1;
  sustainLevel   = 0.5;
  releaseTime    = 0.01;
  endLevel       = 0.0;     
  time           = 0.0;
  timeScale      = 1.0;
  increment      = timeScale/sampleRate;
  tauScale       = 1.0;
  peakScale      = 1.0;
  noteIsOn       = false;
  outputIsZero   = true;

  previousOutput = 0.0;

  // call some internal functions to trigge the coefficient calculations:
  setAttack(attackTime);
  setDecay(decayTime);
  setRelease(releaseTime);
}

EnvelopeGenerator2::~EnvelopeGenerator2()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void EnvelopeGenerator2::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  // adjust time increment:
  increment = timeScale/sampleRate;

  //re-calculate coefficients for the 3 filters:
  setAttack (attackTime);
  setDecay  (decayTime);
  setRelease(releaseTime);
}

void EnvelopeGenerator2::setAttack(double newAttackTime)
{
  if( newAttackTime > 0.0 )
  {
    attackTime  = newAttackTime;
    double tau  = (sampleRate*attackTime) * tauScale/timeScale;
    coeffAttack = 1.0 - exp( -1.0 / tau ); 
  }
  else // newAttackTime <= 0
  {
    attackTime  = 0.0;
    coeffAttack = 1.0; 
  }
  calculateAccumulatedTimes();
}

void EnvelopeGenerator2::setHold(double newHoldTime)
{
  if( newHoldTime >= 0 ) 
    holdTime = newHoldTime;
  calculateAccumulatedTimes();
}

void EnvelopeGenerator2::setDecay(double newDecayTime)
{
  if( newDecayTime > 0.0 )
  {
    decayTime  = newDecayTime;
    double tau = (sampleRate*decayTime) * tauScale/timeScale;
    coeffDecay = 1.0 - exp( -1.0 / tau  );
  }
  else // newDecayTime <= 0
  {
    decayTime  = 0.0;
    coeffDecay = 1.0; 
  }
  calculateAccumulatedTimes();
}

void EnvelopeGenerator2::setRelease(double newReleaseTime)
{
  if( newReleaseTime > 0.0 )
  {
    releaseTime  = newReleaseTime;
    double tau   = (sampleRate*releaseTime) * tauScale/timeScale;
    coeffRelease = 1.0 - exp( -1.0 / tau  );
  }
  else // newReleaseTime <= 0
  {
    releaseTime  = 0.0;
    coeffRelease = 1.0; 
  }
  calculateAccumulatedTimes();
}

void EnvelopeGenerator2::setTimeScale(double newTimeScale)
{
  if( newTimeScale > 0 )
    timeScale = newTimeScale;

  increment  = timeScale/sampleRate;

  //re-calculate coefficients for the 3 filters:
  setAttack (attackTime);
  setDecay  (decayTime);
  setRelease(releaseTime);
}

void EnvelopeGenerator2::setTauScale(double newTauScale)
{
  if( newTauScale > 0 )
    tauScale = newTauScale;

  setAttack(attackTime);
  setDecay(decayTime);
  setRelease(releaseTime);
}

void EnvelopeGenerator2::setPeakScale(double newPeakScale)
{
  if( newPeakScale > 0 )
    peakScale = newPeakScale;
}

//-------------------------------------------------------------------------------------------------
// others:

void EnvelopeGenerator2::reset()
{
  time = 0.0;
}

void EnvelopeGenerator2::trigger(bool startFromCurrentValue)
{
  if( !startFromCurrentValue )
    previousOutput = startLevel;  // may lead to clicks

  // reset time for the new note:
  time         = 0.0;   
  noteIsOn     = true;
  outputIsZero = false;
}

void EnvelopeGenerator2::noteOff()
{
  noteIsOn = false; 

  // advance time to the beginnig of the release phase:
  time = (attackTime + holdTime + decayTime + increment);
}

bool EnvelopeGenerator2::endIsReached()
{
  if( noteIsOn == false && previousOutput < 0.0001 )
    return true;
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void EnvelopeGenerator2::calculateAccumulatedTimes()
{
  attPlusHld               = attackTime + holdTime;
  attPlusHldPlusDec        = attPlusHld + decayTime;
  attPlusHldPlusDecPlusRel = attPlusHldPlusDec + releaseTime;
}