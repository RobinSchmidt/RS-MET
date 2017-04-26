#include "rosic_LevelDetector.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

LevelDetector::LevelDetector()
{
  //mode = MEAN_ABS;   
  //mode = MEAN_SQUARE;   
  //mode = ROOT_MEAN_SQUARE;   
  mode = INSTANTANEOUS_ENVELOPE;   
}

LevelDetector::~LevelDetector()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void LevelDetector::setSampleRate(double newSampleRate)
{
  preSmoother.setSampleRate(newSampleRate);
  attackReleaseFollower.setSampleRate(newSampleRate);
}

void LevelDetector::setAttackTime(double newAttackTime)
{
  attackReleaseFollower.setAttackTime(newAttackTime);
  double tau = rmin(attackReleaseFollower.getAttackTime(), attackReleaseFollower.getReleaseTime());
  preSmoother.setTimeConstant(tau);
}

void LevelDetector::setReleaseTime(double newReleaseTime)
{
  attackReleaseFollower.setReleaseTime(newReleaseTime);
  double tau = rmin(attackReleaseFollower.getAttackTime(), attackReleaseFollower.getReleaseTime());
  preSmoother.setTimeConstant(tau);
}

void LevelDetector::setMode(int Mode)
{
  if( Mode >= MEAN_ABS && Mode < NUM_MODES )
    mode = Mode;
}

//-------------------------------------------------------------------------------------------------
// others:

void LevelDetector::reset()
{
  instEnvDetector.reset();
  preSmoother.reset();
  attackReleaseFollower.reset();
}