// construction/destruction:

template<class TSig, class TPar>
rsCrossOver4Way<TSig, TPar>::rsCrossOver4Way() 
: lowBranchCompensationAllpass(4)
, highBranchCompensationAllpass(4)
{
 
}

// setup:

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  stage1.setSampleRate(newSampleRate);
  for(int s=0; s<2; s++)
    stage2[s].setSampleRate(newSampleRate);
  setupCompensationAllpasses();
}

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::setBandActive(bool shouldBeActive, int treeLevel, int indexInLevel)
{
  if( treeLevel == 0 )
    stage1.setActive(shouldBeActive);
  else if( treeLevel == 1 )
    stage2[indexInLevel].setActive(shouldBeActive);
}

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::setCrossoverFrequency(TPar newCrossoverFrequency, int treeLevel, 
  int indexInLevel)
{
  if( treeLevel == 0 )
    stage1.setCrossoverFrequency(newCrossoverFrequency);
  else if( treeLevel == 1 )
    stage2[indexInLevel].setCrossoverFrequency(newCrossoverFrequency);
  setupCompensationAllpasses();
}

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::setSlope(int newSlope, int treeLevel, int indexInLevel)
{
  if( treeLevel == 0 )
    stage1.setSlope(newSlope);
  else if( treeLevel == 1 )
    stage2[indexInLevel].setSlope(newSlope);
  setupCompensationAllpasses();
}

// inquiry:

template<class TSig, class TPar>
bool rsCrossOver4Way<TSig, TPar>::isBandActive(int treeLevel, int indexInLevel) const
{
  if( treeLevel == 0 )
    return stage1.isActive();
  else if( treeLevel == 1 )
    return stage2[indexInLevel].isActive();
  else 
    return false;
}

template<class TSig, class TPar>
TPar rsCrossOver4Way<TSig, TPar>::getCrossoverFrequency(int treeLevel, int indexInLevel) const
{
  if( treeLevel == 0 )
    return stage1.getCrossoverFrequency();
  else if( treeLevel == 1 )
    return stage2[indexInLevel].getCrossoverFrequency();
  else 
    return 0.0;
}

template<class TSig, class TPar>
int rsCrossOver4Way<TSig, TPar>::getSlope(int treeLevel, int indexInLevel) const
{
  if( treeLevel == 0 )
    return stage1.getSlope();
  else if( treeLevel == 1 )
    return stage2[indexInLevel].getSlope();
  else 
    return 0;
}

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::getMagnitudeResponse(TPar* frequencies, TPar* magnitudes, 
  int numBins, int outputChannel, bool inDecibels)
{
  rsArrayTools::fillWithValue(magnitudes, numBins, TPar(-100));

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

  rsFilterAnalyzer<TPar>::clampValuesAboveNyquist(frequencies, magnitudes, numBins, 
    stage1.getSampleRate(), -100.0);
  rsArrayTools::clip(magnitudes, numBins, TPar(-150), TPar(10));
}

// audio-processing:

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::processBuffer(TSig** inOutBuffer, int length)
{
  /*
  int c, n;
  double sampleFrame[8];
  for(n = 0; n < length; n++)
  {
    sampleFrame[0] = inOutBuffer[0][n];
    sampleFrame[1] = inOutBuffer[1][n];

    processSampleFrame(sampleFrame);

    for(c = 0; c < 8; c++)
      inOutBuffer[c][n] = (float) sampleFrame[c];
  }
  */
}

// others:

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::resetBuffers()
{
  stage1.resetBuffers();
  for(int s = 0; s < 2; s++)
    stage2[s].resetBuffers();
  lowBranchCompensationAllpass.reset();
  highBranchCompensationAllpass.reset();
}

template<class TSig, class TPar>
void rsCrossOver4Way<TSig, TPar>::setupCompensationAllpasses()
{
  lowBranchCompensationAllpass.copySettingsFrom(&stage2[1].sumAllpass);   // we could also copy from the lowpass
  highBranchCompensationAllpass.copySettingsFrom(&stage2[0].sumAllpass);
  highBranchCompensationAllpass.turnIntoAllpass();

  // Interestingly, the allpass resulting from adding a lowpass- and highpass Butterworth-squared 
  // response is itself only a non-squared Butterworth-allpass response. Half of the zeros end up 
  // inside the unit circle canceling with half of the poles.
}
