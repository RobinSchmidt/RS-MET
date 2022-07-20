//#include "rosic_TrackMeter.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

TrackMeter::TrackMeter()
{
  maxLeftLevel      = 0.0;
  maxRightLevel     = 0.0;
  maxMidLevel       = 0.0;
  maxSideLevel      = 0.0;
  sumOfSquaresLeft  = 0.0;
  sumOfSquaresRight = 0.0;
  sumOfProducts     = 0.0;
  sampleCounter     = 0;

  leftLevelExtractor.setMode(EnvelopeFollower::MEAN_SQUARE);
  rightLevelExtractor.setMode(EnvelopeFollower::MEAN_SQUARE);
  midLevelExtractor.setMode(EnvelopeFollower::MEAN_SQUARE);
  sideLevelExtractor.setMode(EnvelopeFollower::MEAN_SQUARE);

  productLevelExtractor.setMode(rsOnePoleFilterDD::LOWPASS_IIT);
  productLevelExtractor.setCutoff(0.25); ///< \todo setTimeConstant instead
  meanSquareExtractorLeft.setMode(rsOnePoleFilterDD::LOWPASS_IIT);
  meanSquareExtractorLeft.setCutoff(0.25);
  meanSquareExtractorRight.setMode(rsOnePoleFilterDD::LOWPASS_IIT);
  meanSquareExtractorRight.setCutoff(0.25); 

  //setAttackTimeInMilliseconds(10.0);
  //setReleaseTimeInMilliseconds(300.0);

  setAttackTimeInMilliseconds( 65.144172285487784);  // VU time-constant
  setReleaseTimeInMilliseconds(65.144172285487784);
}

TrackMeter::~TrackMeter()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void TrackMeter::setSampleRate(double newSampleRate)
{
  leftLevelExtractor.setSampleRate(newSampleRate);
  rightLevelExtractor.setSampleRate(newSampleRate);
  midLevelExtractor.setSampleRate(newSampleRate);
  sideLevelExtractor.setSampleRate(newSampleRate);
  productLevelExtractor.setSampleRate(newSampleRate);
  meanSquareExtractorLeft.setSampleRate(newSampleRate);
  meanSquareExtractorRight.setSampleRate(newSampleRate);
}

void TrackMeter::setAttackTimeInMilliseconds(double newAttackTime)
{
  leftLevelExtractor.setAttackTime(newAttackTime);
  rightLevelExtractor.setAttackTime(newAttackTime);
  midLevelExtractor.setAttackTime(newAttackTime);
  sideLevelExtractor.setAttackTime(newAttackTime);
}

void TrackMeter::setReleaseTimeInMilliseconds(double newReleaseTime)
{
  leftLevelExtractor.setReleaseTime(newReleaseTime);
  rightLevelExtractor.setReleaseTime(newReleaseTime);
  midLevelExtractor.setReleaseTime(newReleaseTime);
  sideLevelExtractor.setReleaseTime(newReleaseTime);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

SignalMeasures TrackMeter::getCurrentMeasurement(bool reset)
{
  if(sampleCounter < 1)
    return currentMeasures; // if this gets called multiple times during one sample

  // assign the measured levels:
  currentMeasures.leftLevel  = (float) RAPT::rsMax( -200.0, RAPT::rsAmpToDb(sqrt(maxLeftLevel))  );
  currentMeasures.rightLevel = (float) RAPT::rsMax( -200.0, RAPT::rsAmpToDb(sqrt(maxRightLevel)) );
  currentMeasures.midLevel   = (float) RAPT::rsMax( -200.0, RAPT::rsAmpToDb(sqrt(maxMidLevel))   );
  currentMeasures.sideLevel  = (float) RAPT::rsMax( -200.0, RAPT::rsAmpToDb(sqrt(maxSideLevel))  );

  // calculate and assign the cross-correlation:
  double factor          = 1.0 / (double)sampleCounter;
  double meanSquareLeft  = factor * sumOfSquaresLeft;
  double meanSquareRight = factor * sumOfSquaresRight;
  double meanProduct     = factor * sumOfProducts;
  double normalizer      = sqrt(meanSquareLeft * meanSquareRight);
  if(normalizer < RAPT::rsDbToAmp(-120.0))
    currentMeasures.crossCorrelation = 0.f;
  else
    currentMeasures.crossCorrelation = (float) (meanProduct / normalizer);

  // reset the accumulating variables:
  if( reset )
  {
    sampleCounter     = 0;
    maxLeftLevel      = 0.0;
    maxRightLevel     = 0.0;
    maxMidLevel       = 0.0;
    maxSideLevel      = 0.0;
    sumOfSquaresLeft  = 0.0;
    sumOfSquaresRight = 0.0;
    sumOfProducts     = 0.0;
  }

  return currentMeasures;
}

