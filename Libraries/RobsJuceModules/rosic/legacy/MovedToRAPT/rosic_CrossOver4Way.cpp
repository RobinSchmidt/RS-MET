// construction/destruction:

rsCrossOver4Way::rsCrossOver4Way() 
: lowBranchCompensationAllpass(4)
, highBranchCompensationAllpass(4)
{
 
}

rsCrossOver4Way::~rsCrossOver4Way()
{

}

// setup:

void rsCrossOver4Way::setSampleRate(double newSampleRate)
{
  stage1.setSampleRate(newSampleRate);
  for(int s=0; s<2; s++)
    stage2[s].setSampleRate(newSampleRate);
  setupCompensationAllpasses();
}

void rsCrossOver4Way::setBandActive(bool shouldBeActive, int treeLevel, int indexInLevel)
{
  if( treeLevel == 0 )
    stage1.setActive(shouldBeActive);
  else if( treeLevel == 1 )
    stage2[indexInLevel].setActive(shouldBeActive);
}

void rsCrossOver4Way::setCrossoverFrequency(double newCrossoverFrequency, int treeLevel, 
  int indexInLevel)
{
  if( treeLevel == 0 )
    stage1.setCrossoverFrequency(newCrossoverFrequency);
  else if( treeLevel == 1 )
    stage2[indexInLevel].setCrossoverFrequency(newCrossoverFrequency);
  setupCompensationAllpasses();
}

void rsCrossOver4Way::setSlope(int newSlope, int treeLevel, int indexInLevel)
{
  if( treeLevel == 0 )
    stage1.setSlope(newSlope);
  else if( treeLevel == 1 )
    stage2[indexInLevel].setSlope(newSlope);
  setupCompensationAllpasses();
}

void rsCrossOver4Way::setMonoMode(bool shouldBeMono)
{
  stage1.setMono(shouldBeMono);
  for(int s=0; s<2; s++)
    stage2[s].setMono(shouldBeMono);
}

// inquiry:

bool rsCrossOver4Way::isBandActive(int treeLevel, int indexInLevel) const
{
  if( treeLevel == 0 )
    return stage1.isActive();
  else if( treeLevel == 1 )
    return stage2[indexInLevel].isActive();
  else 
    return false;
}

double rsCrossOver4Way::getCrossoverFrequency(int treeLevel, int indexInLevel) const
{
  if( treeLevel == 0 )
    return stage1.getCrossoverFrequency();
  else if( treeLevel == 1 )
    return stage2[indexInLevel].getCrossoverFrequency();
  else 
    return 0.0;
}

int rsCrossOver4Way::getSlope(int treeLevel, int indexInLevel) const
{
  if( treeLevel == 0 )
    return stage1.getSlope();
  else if( treeLevel == 1 )
    return stage2[indexInLevel].getSlope();
  else 
    return 0;
}

void rsCrossOver4Way::getMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
  int outputChannel, bool inDecibels)
{
  RAPT::rsArrayTools::fillWithValue(magnitudes, numBins, -100.0);

  if( !stage2[0].isActive() && !stage2[1].isActive() ) 
  {
    // 2 bands:
    if( outputChannel == 0 )
      stage1.getLowpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, false);
    else if( outputChannel == 1 )
      stage1.getHighpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, false);
  }
  else if( stage2[0].isActive() && !stage2[1].isActive() ) 
  {
    // 3 bands, lower band split further:
    if( outputChannel == 0 )
    {
      stage1.getLowpassMagnitudeResponse(   frequencies, magnitudes, numBins, inDecibels, false);
      stage2[0].getLowpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
    else if( outputChannel == 1 )
    {
      stage1.getLowpassMagnitudeResponse(    frequencies, magnitudes, numBins, inDecibels, false);
      stage2[0].getHighpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
    else if( outputChannel == 2 )
      stage1.getHighpassMagnitudeResponse(   frequencies, magnitudes, numBins, inDecibels, false);
  }
  else if( !stage2[0].isActive() && stage2[1].isActive() ) 
  {
    // 3 bands, upper band split further
    if( outputChannel == 0 )
      stage1.getLowpassMagnitudeResponse(   frequencies, magnitudes, numBins, inDecibels, false);
    else if( outputChannel == 1 )
    {
      stage1.getHighpassMagnitudeResponse(  frequencies, magnitudes, numBins, inDecibels, false);
      stage2[1].getLowpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
    else if( outputChannel == 2 )
    {
      stage1.getHighpassMagnitudeResponse(   frequencies, magnitudes, numBins, inDecibels, false);
      stage2[1].getHighpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
  }
  else
  {
    // 4 bands:
    if( outputChannel == 0 )
    {
      stage1.getLowpassMagnitudeResponse(   frequencies, magnitudes, numBins, inDecibels, false);
      stage2[0].getLowpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
    else if( outputChannel == 1 )
    {
      stage1.getLowpassMagnitudeResponse(    frequencies, magnitudes, numBins, inDecibels, false);
      stage2[0].getHighpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
    else if( outputChannel == 2 )
    {
      stage1.getHighpassMagnitudeResponse(  frequencies, magnitudes, numBins, inDecibels, false);
      stage2[1].getLowpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
    else if( outputChannel == 3 )
    {
      stage1.getHighpassMagnitudeResponse(   frequencies, magnitudes, numBins, inDecibels, false);
      stage2[1].getHighpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, true);
    }
  }

  rsFilterAnalyzerD::clampValuesAboveNyquist(frequencies, magnitudes, numBins, 
    stage1.getSampleRate(), -100.0);
  RAPT::rsArrayTools::clipBuffer(magnitudes, numBins, -150.0, 10.0);
}

// audio-processing:

void rsCrossOver4Way::processBuffer(float **inOutBuffer, int length)
{
  int c, n;
  double sampleFrame[8];

  for(n=0; n<length; n++)
  {
    sampleFrame[0] = inOutBuffer[0][n];
    sampleFrame[1] = inOutBuffer[1][n];

    processSampleFrame(sampleFrame);

    for(c=0; c<8; c++)
      inOutBuffer[c][n] = (float) sampleFrame[c];
  }
}

// others:

void rsCrossOver4Way::resetBuffers()
{
  stage1.resetBuffers();
  for(int s=0; s<2; s++)
    stage2[s].resetBuffers();
  lowBranchCompensationAllpass.reset();
  highBranchCompensationAllpass.reset();
}

void rsCrossOver4Way::setupCompensationAllpasses()
{
  lowBranchCompensationAllpass.copySettingsFrom(&stage2[1].crossoverL.sumAllpass);   // we could also copy from the lowpass
  highBranchCompensationAllpass.copySettingsFrom(&stage2[0].crossoverL.sumAllpass);
  highBranchCompensationAllpass.turnIntoAllpass();

  // Interestingly, the allpass resulting from adding a lowpass- and highpass Butterworth-squared 
  // response is itself only a non-squared Butterworth-allpass response. Half of the zeros end up 
  // inside the unit circle canceling with half of the poles.
}
