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
  rsWarning("rsCrossOver4Way::processBuffer() needs tests");

  // The implementation was commented out. I have uncommented it again to fix a compiler warning 
  // but I'm not sure if the implementation below is good or not. We need to set up a unit test for
  // it and then re-implement/uncomment and test it. 

  TSig sampleFrame[8];
  rsArrayTools::fillWithZeros(sampleFrame, 8);  // Not sure, if that's strictly needed.
  for(int n = 0; n < length; n++)
  {
    // Fetch:
    sampleFrame[0] = inOutBuffer[0][n];
    sampleFrame[1] = inOutBuffer[1][n];

    // Process:
    processSampleFrame(sampleFrame);

    // Store:
    for(int c = 0; c < 8; c++)
      inOutBuffer[c][n] = sampleFrame[c];
  }

  // ToDo:
  //
  // - Document why sampleFrame is a vector of 8. I think it may be because we have a stereo input 
  //   (= 2 channels) and with at most 4 frequency bands, we could get up to 2*4 = 8 output 
  //   channels, so processSampleFrame() expects in inOut vector of 8. Verify that!
  // 
  // - Figure out if we should perhaps initialize the sampleFrame vector to all zeros. If so, do 
  //   it and document why it's needed...done.. I think, if we have less than 4 bands, the call to 
  //   processSampleFrame() will not touch the upper components of the in/out vector, so without
  //   the zeroing, we may return garbage in the upper bands. Maybe in many contexts that doesn't 
  //   matter because the caller knows that we are using less bands and will also ignore the upper
  //   components. But there are also conceivable contexts where this may not be the case, so let's
  //   better play it safe.
  //
  // - Maybe unroll the loop over c manually as:
  // 
  //     inOutBuffer[0][n] = sampleFrame[0];
  //     inOutBuffer[1][n] = sampleFrame[1];
  //     ...
  //     inOutBuffer[7][n] = sampleFrame[7];
  //
  //   but I guess the compiler should be able to do that as well. Maybe inspect the generated 
  //   assembly and/or set up a benchmark with both versions.
  //
  // - Why do we even need this Fetch/Process/Store business anyway? Can't we just directly call
  //   processSampleFrame(&inOutBuffer[0][n])? Set up a unit test an try it!
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
