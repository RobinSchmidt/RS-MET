//#include "rosic_PitchDetector.h"
//using namespace rosic;

// Construction/Destruction:

PitchDetector::PitchDetector() : formantRemover(30)
{
  // initialize parameters:
  sampleRate        = 44100.0;
  sampleRateRec     = 1.0/sampleRate;
  minFundamental    = 20.0;
  maxFundamental    = 10000.0;
  minPeriod         = 1.0 / maxFundamental;
  maxPeriod         = 1.0 / minFundamental;
  y0 = y1 = y2 = y3 = 0.0;
  fracOld           = 0.0;
  periodEstimate    = 0.001;
  frequencyEstimate = 1.0 / periodEstimate;
  cycleCorrelation  = 0.0;
  sampleCounter     = 0;

  formantRemover.setOrder(30);

  dcBlocker.setMode(OnePoleFilter::HIGHPASS);
  dcBlocker.setCutoff(20.0);
  dcBlocker.setSampleRate(sampleRate);

  lowpass.setMode(FourPoleFilterParameters::LOWPASS_6);
  lowpass.useTwoStages(true);
  lowpass.setFrequency(50.0);
  lowpass.setSampleRate(sampleRate);

  envFollower.setMode(EnvelopeFollower::MEAN_ABS);
  envFollower.setAttackTime(0.0);
  envFollower.setReleaseTime(20.0);
  envFollower.setSampleRate(sampleRate);
}

// Setup:

void PitchDetector::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate    = newSampleRate;
    sampleRateRec = 1.0/sampleRate;
    dcBlocker.setSampleRate(sampleRate);
    lowpass.setSampleRate(sampleRate);
    envFollower.setSampleRate(sampleRate);
  }
  else
    DEBUG_BREAK; // invalid sample-rate
}

void PitchDetector::setMinFundamental(double newMinFundamental)
{
  if( newMinFundamental >= 10.0 
    && newMinFundamental <= 5000.0 
    && newMinFundamental < maxFundamental)
  {
    minFundamental = newMinFundamental;
    maxPeriod      = 1.0 / minFundamental;
  }
}

void PitchDetector::setMaxFundamental(double newMaxFundamental)
{
  if( newMaxFundamental >= 100.0 
    && newMaxFundamental <= 20000.0 
    && newMaxFundamental > minFundamental)
  {
    maxFundamental = newMaxFundamental;
    minPeriod      = 1.0 / maxFundamental;
  }
}

