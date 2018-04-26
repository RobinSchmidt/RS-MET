using namespace RSLib;

// Construction/Destruction:

rsZeroCrossingPitchDetector::rsZeroCrossingPitchDetector() : formantRemover(30)
{
  // initialize parameters:
  sampleRate        = 44100.0;
  sampleRateRec     = 1.0/sampleRate;
  minFundamental    = 20.0;
  maxFundamental    = 10000.0;
  minPeriod         = 1.0 / maxFundamental;
  maxPeriod         = 1.0 / minFundamental;

  formantRemover.setOrder(30);

  dcBlocker.setMode(rsOnePoleFilter::HIGHPASS);
  dcBlocker.setCutoff(20.0);
  dcBlocker.setSampleRate(sampleRate);

  lowpass.setMode(rsFourPoleFilterParameters::LOWPASS_6);
  lowpass.useTwoStages(true);
  lowpass.setFrequency(50.0);
  lowpass.setSampleRate(sampleRate);

  envFollower.setMode(rsEnvelopeFollower::MEAN_ABS);
  envFollower.setAttackTime(0.0);
  envFollower.setReleaseTime(20.0);
  envFollower.setSampleRate(sampleRate);

  reset();
  /*
  y0 = y1 = y2 = y3 = 0.0;
  fracOld           = 0.0;
  periodEstimate    = 0.001;
  frequencyEstimate = 1.0 / periodEstimate;
  cycleCorrelation  = 0.0;
  sampleCounter     = 0;
  */
}

// Setup:

void rsZeroCrossingPitchDetector::setSampleRate(double newSampleRate)
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
    RS_DEBUG_BREAK; // invalid sample-rate
}

void rsZeroCrossingPitchDetector::setMinFundamental(double newMinFundamental)
{
  if( newMinFundamental >= 10.0 
    && newMinFundamental <= 5000.0 
    && newMinFundamental < maxFundamental)
  {
    minFundamental = newMinFundamental;
    maxPeriod      = 1.0 / minFundamental;
  }
}

void rsZeroCrossingPitchDetector::setMaxFundamental(double newMaxFundamental)
{
  if( newMaxFundamental >= 100.0 
    && newMaxFundamental <= 20000.0 
    && newMaxFundamental > minFundamental)
  {
    maxFundamental = newMaxFundamental;
    minPeriod      = 1.0 / maxFundamental;
  }
}

void rsZeroCrossingPitchDetector::reset(double initialEstimate)
{
  y0 = y1 = y2 = y3 = 0.0;
  fracOld           = 0.0;
  frequencyEstimate = initialEstimate;
  periodEstimate    = 1.0 / frequencyEstimate;
  cycleCorrelation  = 0.0;
  sampleCounter     = 0;

  formantRemover.reset();
  dcBlocker.reset();
  lowpass.reset();
  envFollower.reset();
}