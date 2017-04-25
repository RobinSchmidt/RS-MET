#include "rosic_Modulator.h"
using namespace rosic;

//-----------------------------------------------------------------------------
// construction/destruction:

Modulator::Modulator()
{
  //TEST:
  //fadeInOutSlewRateLimiter.setReleaseTime(5.0);

  scaleFactor      = 1.0;
  offset           = 0.0;

  currentKey       = -1;
  currentVel       = -1;
  applyMode        = ADDITIVE;

  timeScaleNominal = 1.0;
  //timeScale        = 1.0;       
  timeScaleByKey   = 0.0;  
  timeScaleByVel   = 0.0; 

  depthNominal     = 1.0;
  depth            = 1.0;        
  depthByKey       = 0.0;   
  depthByVel       = 0.0;   

  slewUp           = 0.0;
  slewUpByKey      = 0.0;
  slewUpByVel      = 0.0;

  slewDown         = 0.0;
  slewDownByKey    = 0.0;
  slewDownByVel    = 0.0;

  fadeInTime       = 0.0;
  fadeInTimeByKey  = 0.0;
  fadeInTimeByVel  = 0.0;
  fadeInIsOn       = false;

  fadeOutTime      = 0.0;
  fadeOutTimeByKey = 0.0;
  fadeOutTimeByVel = 0.0;
  fadeOutIsOn      = false;

  noRelease        = false;
}

Modulator::~Modulator()
{

}

void Modulator::copyDataFrom(const rosic::Modulator &source)
{
  breakpointModulator.copyDataFrom(source.breakpointModulator);
  sampleModulator.copyDataFrom(source.sampleModulator);

  upDownSlewRateLimiter    = source.upDownSlewRateLimiter;
  fadeInOutSlewRateLimiter = source.fadeInOutSlewRateLimiter;

  applyMode                = source.applyMode;
  currentKey               = source.currentKey;
  currentVel               = source.currentVel;
  offset                   = source.offset;
  depth                    = source.depth;
  depthByKey               = source.depthByKey;
  depthByVel               = source.depthByVel;
  depthNominal             = source.depthNominal;
  fadeInIsOn               = source.fadeInIsOn;
  fadeInTime               = source.fadeInTime;
  fadeInTimeByKey          = source.fadeInTimeByKey;
  fadeInTimeByVel          = source.fadeInTimeByVel;
  fadeOutIsOn              = source.fadeOutIsOn;
  fadeOutTime              = source.fadeOutTime;
  fadeOutTimeByKey         = source.fadeOutTimeByKey;
  fadeOutTimeByVel         = source.fadeOutTimeByVel;
  modulationSource         = source.modulationSource;
  noRelease                = source.noRelease;
  scaleFactor              = source.scaleFactor;
  slewDown                 = source.slewDown;
  slewDownByKey            = source.slewDownByKey;
  slewDownByVel            = source.slewDownByVel;
  slewUp                   = source.slewUp;
  slewUpByKey              = source.slewUpByKey;
  slewUpByVel              = source.slewUpByVel;
  timeScaleByKey           = source.timeScaleByKey;
  timeScaleByVel           = source.timeScaleByVel;
  timeScaleNominal         = source.timeScaleNominal;
}

//-----------------------------------------------------------------------------
// parameter settings:

void Modulator::setSampleRate(double newSampleRate)
{
  breakpointModulator.setSampleRate(newSampleRate);
  upDownSlewRateLimiter.setSampleRate(newSampleRate);
  fadeInOutSlewRateLimiter.setSampleRate(newSampleRate);
  sampleModulator.setSampleRate(newSampleRate);
  //riseRamp.setSampleRate(newSampleRate);
}

void Modulator::setModulationSource(int newModulationSource)
{
  if( newModulationSource >= BREAKPOINT_MODULATOR &&
      newModulationSource <= SAMPLE_MODULATOR               )
  {
    modulationSource = newModulationSource;
  }
}

void Modulator::setScaleFactor(double newScaleFactor)
{
  scaleFactor = newScaleFactor;
}

void Modulator::setOffset(double newOffset)
{
  offset = newOffset;
}

void Modulator::setBeatsPerMinute(double newBpm)
{
  breakpointModulator.setBeatsPerMinute(newBpm);
  sampleModulator.setBeatsPerMinute(newBpm);
}

void Modulator::setTimeScale(double newTimeScale)
{ 
  if( newTimeScale >= 0.0001 )
    timeScaleNominal = newTimeScale;

  double timeScale = timeScaleNominal * pow(2.0, 0.01*timeScaleByKey*(currentKey-60.0)/12.0)
                                      * pow(2.0, 0.01*timeScaleByVel*(currentVel-64.0)/63.0);
  breakpointModulator.setTimeScale(timeScale);
  sampleModulator.setNumCyclesPerTimeUnit(1.0/timeScale);
}

void Modulator::setTimeScaleByKey(double newTimeScaleByKey)
{ 
  timeScaleByKey = newTimeScaleByKey;

  double timeScale = timeScaleNominal * pow(2.0, 0.01*timeScaleByKey*(currentKey-60.0)/12.0)
                                      * pow(2.0, 0.01*timeScaleByVel*(currentVel-64.0)/63.0);
  breakpointModulator.setTimeScale(timeScale);
  sampleModulator.setNumCyclesPerTimeUnit(1.0/timeScale);
}

void Modulator::setTimeScaleByVel(double newTimeScaleByVel)
{ 
  timeScaleByVel = newTimeScaleByVel;

  double timeScale = timeScaleNominal * pow(2.0, 0.01*timeScaleByKey*(currentKey-60.0)/12.0)
                                      * pow(2.0, 0.01*timeScaleByVel*(currentVel-64.0)/63.0);
  breakpointModulator.setTimeScale(timeScale);
  sampleModulator.setNumCyclesPerTimeUnit(1.0/timeScale);
}

void Modulator::setDepth(double newDepth)
{ 
  depthNominal = newDepth;

  depth =   depthNominal 
    * pow(2.0, 0.01*depth*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*depth*(currentVel-64.0)/63.0);
}

void Modulator::setDepthByKey(double newDepthByKey)
{ 
  depthByKey = newDepthByKey;

  depth =   depthNominal 
    * pow(2.0, 0.01*depth*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*depth*(currentVel-64.0)/63.0);
}

void Modulator::setDepthByVel(double newDepthByVel)
{ 
  depthByVel = newDepthByVel;

  depth =   depthNominal 
    * pow(2.0, 0.01*depth*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*depth*(currentVel-64.0)/63.0);
}

void Modulator::setUpwardSlewRateLimit(double newLimit) 
{
  if( newLimit >= 0.0 )
    slewUp = newLimit;

  double tmp =   slewUp
    * pow(2.0, 0.01*slewUpByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewUpByVel*(currentVel-64.0)/63.0);

  upDownSlewRateLimiter.setAttackTime(tmp); 
}

void Modulator::setUpwardSlewRateLimitByKey(double newLimitByKey) 
{
  slewUpByKey = newLimitByKey;

  double tmp =   slewUp
    * pow(2.0, 0.01*slewUpByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewUpByVel*(currentVel-64.0)/63.0);

  upDownSlewRateLimiter.setAttackTime(tmp); 
}

void Modulator::setUpwardSlewRateLimitByVel(double newLimitByVel) 
{
  slewUpByVel = newLimitByVel;

  double tmp =   slewUp
    * pow(2.0, 0.01*slewUpByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewUpByVel*(currentVel-64.0)/63.0);

  upDownSlewRateLimiter.setAttackTime(tmp); 
}

void Modulator::setDownwardSlewRateLimit(double newLimit)
{ 
  if( newLimit >= 0.0 )
    slewDown = newLimit;

  double tmp =   slewDown
    * pow(2.0, 0.01*slewDownByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewDownByVel*(currentVel-64.0)/63.0);

  upDownSlewRateLimiter.setReleaseTime(tmp); 
}

void Modulator::setDownwardSlewRateLimitByKey(double newLimitByKey) 
{
  slewDownByKey = newLimitByKey;

  double tmp =   slewDown
    * pow(2.0, 0.01*slewDownByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewDownByVel*(currentVel-64.0)/63.0);

  upDownSlewRateLimiter.setReleaseTime(tmp); 
}

void Modulator::setDownwardSlewRateLimitByVel(double newLimitByVel) 
{
  slewDownByVel = newLimitByVel;

  double tmp =   slewDown
    * pow(2.0, 0.01*slewDownByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewDownByVel*(currentVel-64.0)/63.0);

  upDownSlewRateLimiter.setReleaseTime(tmp); 
}

void Modulator::setFadeInTime(double newFadeInTime)
{ 
  if( newFadeInTime >= 0.0 )
    fadeInTime = newFadeInTime;

  double tmp =   fadeInTime
    * pow(2.0, 0.01*fadeInTimeByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeInTimeByVel*(currentVel-64.0)/63.0);

  fadeInOutSlewRateLimiter.setAttackTime(tmp); 
}

void Modulator::setFadeInTimeByKey(double newFadeInTimeByKey) 
{
  fadeInTimeByKey = newFadeInTimeByKey;

  double tmp =   fadeInTime
    * pow(2.0, 0.01*fadeInTimeByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeInTimeByVel*(currentVel-64.0)/63.0);

  fadeInOutSlewRateLimiter.setAttackTime(tmp); 
}

void Modulator::setFadeInTimeByVel(double newFadeInTimeByVel) 
{
  fadeInTimeByVel = newFadeInTimeByVel;

  double tmp =   fadeInTime
    * pow(2.0, 0.01*fadeInTimeByVel*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeInTimeByVel*(currentVel-64.0)/63.0);

  fadeInOutSlewRateLimiter.setAttackTime(tmp); 
}

void Modulator::switchFadeInOnOff(bool fadeInShouldBeOn)
{
  fadeInIsOn = fadeInShouldBeOn;
}

void Modulator::setFadeOutTime(double newFadeOutTime)
{ 
  if( newFadeOutTime >= 0.0 )
    fadeOutTime = newFadeOutTime;

  double tmp =   fadeOutTime
    * pow(2.0, 0.01*fadeOutTimeByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeOutTimeByVel*(currentVel-64.0)/63.0);

  fadeInOutSlewRateLimiter.setReleaseTime(tmp); 
}

void Modulator::setFadeOutTimeByKey(double newFadeOutTimeByKey) 
{
  fadeOutTimeByKey = newFadeOutTimeByKey;

  double tmp =   fadeOutTime
    * pow(2.0, 0.01*fadeOutTimeByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeOutTimeByVel*(currentVel-64.0)/63.0);

  fadeInOutSlewRateLimiter.setReleaseTime(tmp); 
}

void Modulator::setFadeOutTimeByVel(double newFadeOutTimeByVel) 
{
  fadeOutTimeByVel = newFadeOutTimeByVel;

  double tmp =   fadeOutTime
    * pow(2.0, 0.01*fadeOutTimeByVel*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeOutTimeByVel*(currentVel-64.0)/63.0);

  fadeInOutSlewRateLimiter.setReleaseTime(tmp); 
}

void Modulator::switchFadeOutOnOff(bool fadeOutShouldBeOn)
{
  fadeOutIsOn = fadeOutShouldBeOn;
}

void Modulator::setApplyMode(int newApplyMode)
{
  applyMode = newApplyMode;
}

void Modulator::setLoopMode(int newLoopMode)
{
  breakpointModulator.setLoopMode(newLoopMode != BreakpointModulator::NO_LOOP);
  sampleModulator.setLoopMode(newLoopMode);
}

void Modulator::setLoopStart(int newLoopStart)
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    breakpointModulator.setLoopStartIndex(newLoopStart);
  else if( modulationSource == SAMPLE_MODULATOR )
    sampleModulator.setLoopStartSample(newLoopStart);
}

void Modulator::setLoopEnd(int newLoopEnd)
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    breakpointModulator.setLoopEndIndex(newLoopEnd);
  else if( modulationSource == SAMPLE_MODULATOR )
    sampleModulator.setLoopEndSample(newLoopEnd);
}

void Modulator::setNumCyclesInLoop(int newNumCycles)
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    breakpointModulator.setNumCyclesInLoop(newNumCycles);
  else if( modulationSource == SAMPLE_MODULATOR )
    sampleModulator.setNumCyclesInLoop(newNumCycles);
}

void Modulator::setSyncMode(bool shouldBeSynced)
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    breakpointModulator.setSyncMode(shouldBeSynced);
  else if( modulationSource == SAMPLE_MODULATOR )
    sampleModulator.setSyncMode(shouldBeSynced);
}

void Modulator::switchReleaseOnOff(bool shouldBeOn)
{
  noRelease = !shouldBeOn;
}

//-----------------------------------------------------------------------------
// inquiry:

double Modulator::getSampleRate()
{
  return sampleModulator.getSampleRate();
  // it doesn't matter from which of the embedded objects we return the sample-rate as they are
  // all in sync.
}

int Modulator::getModulationSource() const
{
  return modulationSource;
}

double Modulator::getScaleFactor() const
{
  return scaleFactor;
}

double Modulator::getOffset() const
{
  return offset;
}

double Modulator::getTimeScale() const
{
  return timeScaleNominal;
}

double Modulator::getTimeScaleByKey() const
{
  return timeScaleByKey;
}

double Modulator::getTimeScaleByVel() const
{
  return timeScaleByVel;
}


double Modulator::getDepth() const
{
  return depthNominal;
}


double Modulator::getDepthByKey() const
{ 
  return depthByKey;
}

double Modulator::getDepthByVel() const
{
  return depthByVel;
}

double Modulator::getUpwardSlewRateLimit() const
{ 
  return slewUp; 
}

double Modulator::getUpwardSlewRateLimitByKey() const
{ 
  return slewUpByKey; 
}

double Modulator::getUpwardSlewRateLimitByVel() const
{ 
  return slewUpByVel; 
}

double Modulator::getDownwardSlewRateLimit() const
{ 
  return slewDown; 
}

double Modulator::getDownwardSlewRateLimitByKey() const
{ 
  return slewDownByKey; 
}

double Modulator::getDownwardSlewRateLimitByVel() const
{ 
  return slewDownByVel; 
}

double Modulator::getFadeInTime() const
{ 
  return fadeInTime;
}

double Modulator::getFadeInTimeByKey() const
{ 
  return fadeInTimeByKey;
}

double Modulator::getFadeInTimeByVel() const
{ 
  return fadeInTimeByVel;
}

bool Modulator::isFadeInOn()
{
  return fadeInIsOn;
}

double Modulator::getFadeOutTime() const
{ 
  return fadeOutTime;
}

double Modulator::getFadeOutTimeByKey() const
{ 
  return fadeOutTimeByKey;
}

double Modulator::getFadeOutTimeByVel() const
{ 
  return fadeOutTimeByVel;
}

bool Modulator::isFadeOutOn()
{
  return fadeOutIsOn;
}

int Modulator::getApplyMode()
{
  return applyMode;
}

int Modulator::getLoopMode()
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    return breakpointModulator.getLoopMode();
  else if( modulationSource == SAMPLE_MODULATOR )
    return sampleModulator.getLoopMode();
  else 
    return 0;
}

int Modulator::getLoopStart()
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    return breakpointModulator.getLoopStartIndex();
  else if( modulationSource == SAMPLE_MODULATOR )
    return sampleModulator.getLoopStartSample();
  else 
    return 0;
}

int Modulator::getLoopEnd()
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    return breakpointModulator.getLoopEndIndex();
  else if( modulationSource == SAMPLE_MODULATOR )
    return sampleModulator.getLoopEndSample();
  else 
    return 0;
}

int Modulator::getNumCyclesInLoop()
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    return breakpointModulator.getNumCyclesInLoop();
  else if( modulationSource == SAMPLE_MODULATOR )
    return sampleModulator.getNumCyclesInLoop();
  else 
    return 0;
}

bool Modulator::isInSyncMode()
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    return breakpointModulator.isInSyncMode();
  else if( modulationSource == SAMPLE_MODULATOR )
    return sampleModulator.isInSyncMode();
  else 
    return false;
}

bool Modulator::isReleaseOn()
{
  return !noRelease;
}

bool Modulator::endIsReached()
{
  if( modulationSource == BREAKPOINT_MODULATOR )
    return breakpointModulator.endIsReached;
  else if( modulationSource == SAMPLE_MODULATOR )
    return sampleModulator.endIsReached;
  else 
    return false;
}


//-----------------------------------------------------------------------------
// event handling :

void Modulator::noteOn(bool startFromCurrentLevel, int newKey, int newVelocity)
{
  // keep track of the note-key and velocity:
  currentKey = newKey;
  currentVel = newVelocity;

  // set the time scale-variable for the inherited BreakpointModulator object 
  // according to key and velocity:
  double timeScale = timeScaleNominal * pow(2.0, 0.01*timeScaleByKey*(currentKey-60.0)/12.0)
                                      * pow(2.0, 0.01*timeScaleByVel*(currentVel-64.0)/63.0);

  if( modulationSource == BREAKPOINT_MODULATOR )
    breakpointModulator.setTimeScale(timeScale);
  else
    sampleModulator.setNumCyclesPerTimeUnit(1.0/timeScale);

  depth = depthNominal * pow(2.0, 0.01*depth*(currentKey-60.0)/12.0)
                       * pow(2.0, 0.01*depth*(currentVel-64.0)/63.0);

  // set up the embedded SlewRateLimiters according to the key and velocity:
  double tmp;
  tmp =   slewUp
    * pow(2.0, 0.01*slewUpByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewUpByVel*(currentVel-64.0)/63.0);
  upDownSlewRateLimiter.setAttackTime(tmp);
  tmp =   slewDown
    * pow(2.0, 0.01*slewDownByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*slewDownByVel*(currentVel-64.0)/63.0);
  upDownSlewRateLimiter.setReleaseTime(tmp);
  tmp =   fadeInTime
    * pow(2.0, 0.01*fadeInTimeByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeInTimeByVel*(currentVel-64.0)/63.0);
  fadeInOutSlewRateLimiter.setAttackTime(tmp);
  tmp =   fadeOutTime
    * pow(2.0, 0.01*fadeOutTimeByKey*(currentKey-60.0)/12.0)
    * pow(2.0, 0.01*fadeOutTimeByVel*(currentVel-64.0)/63.0);
  fadeInOutSlewRateLimiter.setReleaseTime(tmp);

  // let the inherited BreakpointModulator object do its stuff which happens on
  // note-on:
  breakpointModulator.noteOn(startFromCurrentLevel);

  // retrigger the embedded SlewRateLimiter and ExponentialRamp, if we do not 
  // start from the current level:
  if( startFromCurrentLevel == false )
  {
    upDownSlewRateLimiter.reset();
    fadeInOutSlewRateLimiter.reset();
  }
}

void Modulator::noteOff(bool startFromCurrentLevel)
{
  if( noRelease == false )
    breakpointModulator.noteOff(startFromCurrentLevel);

  currentKey = -1;
  currentVel = -1;
}



