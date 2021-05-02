//#include "rosic_FeedbackDelayNetwork16.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FeedbackDelayNetwork16::FeedbackDelayNetwork16()
{
  numDelayLines       = 16;
  impulseAmplitude    = 0.f;
  sampleRate          = 44100.0;
  midReverbTime       = 3.0;
  lowReverbTimeScale  = 1.0;
  highReverbTimeScale = 0.3;
  lowCrossoverFreq    = 250.0;
  highCrossoverFreq   = 4000.0;
  delayOrdering       = ASCENDING;
  referenceDelayTime  = 1000*1000/sampleRate;   // -> 22.67..ms -> 1000 samples

  // initialize the pointers to the beginnings of the delaylines:
  for(int d=0; d < maxNumDelayLines; d++)
  {
    tapIns[d]  = 0;
    tapOuts[d] = 0;
  }

  //assignRelativeDelayTimesAlgorithmically(DISTANCE_DECAY, 0.4, 0.8);
  //assignRelativeDelayTimesAlgorithmically(SIMPLE_RATIO_MODES, 0.0, 0.0);
  assignRelativeDelayTimesAlgorithmically(GEOMETRIC_MEANS, 1.0, 2.4);
  setSampleRate(sampleRate);

  // initialize the damping- and correction-filters:
  //setLowReverbTimeScale(2.0);
  setLowReverbTimeScale(1.0);
  setLowCrossoverFreq(250.0);
  setMidReverbTime(3.0);
  setHighReverbTimeScale(0.3);
  //setHighReverbTimeScale(1.0);
  setHighCrossoverFreq(2000.0);

  // setup the FDN-parameters:
  setInjectionVector(IN_ALL_ONES);
  setFeedbackMatrix(HADAMARD);
  setOutputVector(OUT_04);
  allpassMode = false;

  // setup the output filters:
  setWetLpfCutoff(20000.0);
  setWetHpfCutoff(20.0);

  // setup the output parameters:
  wetPinking       = true;
  stereoSwapSwitch = false;
  dryVol           = 0.0;
  wetVol           = 1.0;

  // reset the content of the delaylines samples to zero and trigger a reset
  // in the embedded modules:
  reset();
}

FeedbackDelayNetwork16::~FeedbackDelayNetwork16()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void FeedbackDelayNetwork16::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  // adjust the delay-times in samples:
  adjustDelayTimes();

  // this, in turn, requires the injection-vector to be re-calculated:
  setInjectionVector(injectionVectorIndex);

  // the feedback-filters must be informed about the new sample-rate, too:
  for(int k=0; k<numDelayLines; k++)
    dampingFilters[k].setSampleRate(sampleRate);

  // update the correction-filters and the filters for the wet signal:
  correctionFilterL.setSampleRate(sampleRate);
  correctionFilterR.setSampleRate(sampleRate);
  wetFilterL.setSampleRate(sampleRate);
  wetFilterR.setSampleRate(sampleRate);
  updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork16::setLowReverbTimeScale(double newLowReverbTimeScale)
{
  if( newLowReverbTimeScale > 0.0 )
    lowReverbTimeScale = newLowReverbTimeScale;
  updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork16::setLowCrossoverFreq(double newLowCrossoverFreq)
{
  if( newLowCrossoverFreq >= 20.0 && newLowCrossoverFreq <= 20000.0 )
    lowCrossoverFreq = newLowCrossoverFreq;
  updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork16::setMidReverbTime(double newMidReverbTime)
{
  if( newMidReverbTime > 0.0 )
    midReverbTime = newMidReverbTime;
  updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork16::setHighReverbTimeScale(double newHighReverbTimeScale)
{
  if( newHighReverbTimeScale > 0.0 )
    highReverbTimeScale = newHighReverbTimeScale;
  updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork16::setHighCrossoverFreq(double newHighCrossoverFreq)
{
  if( newHighCrossoverFreq >= 20.0 && newHighCrossoverFreq <= 20000.0 )
    highCrossoverFreq = newHighCrossoverFreq;
  updateDampingAndCorrectionFilters();
}

/*
void FeedbackDelayNetwork16::setDensity(int newDensity)
{
  double tmp = linToLin(newDensity, 0.0, 100.0, 0.6, 0.85);
  setDelayDistributionForm(tmp);
}
*/

void FeedbackDelayNetwork16::setAllpassMode(bool shouldBeAllpass)
{
  allpassMode = shouldBeAllpass;
}

void FeedbackDelayNetwork16::setDryVolume(double newDryVolume)
{
  dryVol = RAPT::rsDbToAmp(newDryVolume);
}

void FeedbackDelayNetwork16::setWetVolume(double newWetVolume)
{
  wetVol = RAPT::rsDbToAmp(newWetVolume);
}

void FeedbackDelayNetwork16::setWetLpfCutoff(double newWetLpfCutoff)
{
  wetFilterL.setLowpassCutoff(newWetLpfCutoff);
  wetFilterR.setLowpassCutoff(newWetLpfCutoff);
}

void FeedbackDelayNetwork16::setWetHpfCutoff(double newWetHpfCutoff)
{
  wetFilterL.setHighpassCutoff(newWetHpfCutoff);
  wetFilterR.setHighpassCutoff(newWetHpfCutoff);
}

void FeedbackDelayNetwork16::setStereoSwapSwitch(bool newStereoSwapSwitch)
{
  stereoSwapSwitch = newStereoSwapSwitch;
}

void FeedbackDelayNetwork16::setWetPinkingSwitch(bool newWetPinkingSwitch)
{
  wetPinking = newWetPinkingSwitch;
}

/*
void FeedbackDelayNetwork16::setMinDelayTime(double newMinDelayTime)
{
  // should be at least 1 millisecond shorter than the maximum:
  if( newMinDelayTime >= 0.1 && newMinDelayTime <= 2000.0 && newMinDelayTime <  maxDelayTime-1.0 )
    minDelayTime = newMinDelayTime;
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::setMaxDelayTime(double newMaxDelayTime)
{
  // should be at least 1 millisecond longer than the minimum:
  if( newMaxDelayTime >= 0.1 && newMaxDelayTime <= 2000.0  && newMaxDelayTime >  minDelayTime+1.0 )
    maxDelayTime = newMaxDelayTime;
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::setDelayDistributionIndex(int newDelayDistributionIndex)
{
  if( newDelayDistributionIndex >= 0 &&
    newDelayDistributionIndex <  NUM_DELAY_DISTRIBUTIONS )
  {
    delayDistributionIndex = newDelayDistributionIndex;
  }
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::setDelayDistributionForm(double newDelayDistributionForm)
{
  if( newDelayDistributionForm >= 0.1 && newDelayDistributionForm <= 0.9 )
    delayDistributionForm = newDelayDistributionForm;
  adjustDelayTimes();
}
*/

void FeedbackDelayNetwork16::setReferenceDelayTime(double newReferenceDelayTime)
{
  if( newReferenceDelayTime >= 1.0 && newReferenceDelayTime <= 100.0 )
    referenceDelayTime = newReferenceDelayTime;
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::setRelativeDelayTime(int delayLineIndex, double newRelativeDelayTime)
{
  if( delayLineIndex >= 0 && delayLineIndex < maxNumDelayLines && newRelativeDelayTime > 0.0)
    relativeDelayTimes[delayLineIndex] = newRelativeDelayTime;
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::setAllRelativeDelayTimes(double *newRelativeDelayTimes)
{
  for(int d=0; d<maxNumDelayLines; d++)
  {
    if( newRelativeDelayTimes[d] > 0.0 )
      relativeDelayTimes[d] = newRelativeDelayTimes[d];
  }
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::assignRelativeDelayTimesAlgorithmically(int distributionIndex,
                                                                     double parameter1,
                                                                     double parameter2)
{
  int i; // for indexing the delay-time
  switch( distributionIndex )
  {
  case LINEAR:
    {
      relativeDelayTimes[0]               = 1.0;
      relativeDelayTimes[numDelayLines-1] = parameter1;  // maximum relative delaytime
      double offset                       = (parameter1-1.0) / (double) (numDelayLines-1);
      for(i=1; i<=numDelayLines-2; i++)
        relativeDelayTimes[i] = relativeDelayTimes[0] + i*offset;
    }
    break;
  case DISTANCE_DECAY:
    {
      double delta         = parameter1;    // first distance
      double factor        = parameter2;    // decay factor for successive distances
      double accu          = delta;
      double relativeDelay = 1.0;
      for(i=0; i<numDelayLines; i++)
      {
        relativeDelayTimes[i]  = relativeDelay;
        relativeDelay         += accu;
        accu                  *= factor;
      }
    }
    break;
  case SIMPLE_RATIO_MODES:
    {
      relativeDelayTimes[ 0] = 1.0;
      relativeDelayTimes[ 1] =  8.0 / 7.0;  // 1.142857...
      relativeDelayTimes[ 2] =  7.0 / 6.0;  // 1.1666...
      relativeDelayTimes[ 3] =  6.0 / 5.0;  // 1.2
      relativeDelayTimes[ 4] =  5.0 / 4.0;  // 1.25
      relativeDelayTimes[ 5] =  9.0 / 7.0;  // 1.285714...
      relativeDelayTimes[ 6] =  4.0 / 3.0;  // 1.3333...
      relativeDelayTimes[ 7] =  7.0 / 5.0;  // 1.4
      relativeDelayTimes[ 8] = 10.0 / 7.0;  // 1.42857...
      relativeDelayTimes[ 9] =  3.0 / 2.0;  // 1.5
      relativeDelayTimes[10] = 11.0 / 7.0;  // 1.571428...
      relativeDelayTimes[11] =  8.0 / 5.0;  // 1.6
      relativeDelayTimes[12] =  5.0 / 3.0;  // 1.6666...
      relativeDelayTimes[13] =  7.0 / 4.0;  // 1.75
      relativeDelayTimes[14] =  9.0 / 5.0;  // 1.8
      relativeDelayTimes[15] = 11.0 / 6.0;  // 1.8333...
    }
    break;

  case GEOMETRIC_MEANS:
    {
      relativeDelayTimes[ 0] = parameter1;
      double dMax            = parameter2;

      relativeDelayTimes[8]  = sqrt(dMax                   * relativeDelayTimes[0]  );
      relativeDelayTimes[4]  = sqrt(relativeDelayTimes[ 8] * relativeDelayTimes[0]  );
      relativeDelayTimes[12] = sqrt(dMax                   * relativeDelayTimes[8]  );
      relativeDelayTimes[2]  = sqrt(relativeDelayTimes[4]  * relativeDelayTimes[0]  );
      relativeDelayTimes[6]  = sqrt(relativeDelayTimes[8]  * relativeDelayTimes[4]  );
      relativeDelayTimes[10] = sqrt(relativeDelayTimes[12] * relativeDelayTimes[8]  );
      relativeDelayTimes[14] = sqrt(dMax                   * relativeDelayTimes[12] );
      relativeDelayTimes[1]  = sqrt(relativeDelayTimes[2]  * relativeDelayTimes[0]  );
      relativeDelayTimes[3]  = sqrt(relativeDelayTimes[4]  * relativeDelayTimes[2]  );
      relativeDelayTimes[5]  = sqrt(relativeDelayTimes[6]  * relativeDelayTimes[4]  );
      relativeDelayTimes[7]  = sqrt(relativeDelayTimes[8]  * relativeDelayTimes[6]  );
      relativeDelayTimes[9]  = sqrt(relativeDelayTimes[10] * relativeDelayTimes[8]  );
      relativeDelayTimes[11] = sqrt(relativeDelayTimes[12] * relativeDelayTimes[10] );
      relativeDelayTimes[13] = sqrt(relativeDelayTimes[14] * relativeDelayTimes[12] );
      relativeDelayTimes[15] = sqrt(dMax                   * relativeDelayTimes[14] );

      // hmm..isn't this equivalent to using a pow-function:
      // delayTimes[i] = pow(dMax, i/(numDelayLines-1)) or something?
    }
    break;





    /*
  case EXPONENTIAL:
    {
      double ratio = maxDelayTime/minDelayTime;

      for(i=1; i <= numDelayLines-2; i++)
        desiredDelays[i] = desiredDelays[0] *
        pow(ratio, (double) i / (double) (numDelayLines-1) );
    }
    break;
    */

   /*
  case PRIME_ALGO_1:
    {
      relativeDelayTimes[ 1] = (13.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 2] = (17.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 3] = (19.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 4] = (23.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 5] = (29.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 6] = (31.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 7] = (37.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 8] = (41.0/11.0) * desiredDelays[0];
      relativeDelayTimes[ 9] = (43.0/11.0) * desiredDelays[0];
      relativeDelayTimes[10] = (47.0/11.0) * desiredDelays[0];
      relativeDelayTimes[11] = (53.0/11.0) * desiredDelays[0];
      relativeDelayTimes[12] = (59.0/11.0) * desiredDelays[0];
      relativeDelayTimes[13] = (61.0/11.0) * desiredDelays[0];
      relativeDelayTimes[14] = (67.0/11.0) * desiredDelays[0];
      relativeDelayTimes[15] = (71.0/11.0) * desiredDelays[0];
    }
    break;
  case PRIME_ALGO_2:
    {
      desiredDelays[ 1] = ( 13.0/11.0) * desiredDelays[0];
      desiredDelays[ 2] = ( 19.0/11.0) * desiredDelays[0];
      desiredDelays[ 3] = ( 23.0/11.0) * desiredDelays[0];
      desiredDelays[ 4] = ( 31.0/11.0) * desiredDelays[0];
      desiredDelays[ 5] = ( 37.0/11.0) * desiredDelays[0];
      desiredDelays[ 6] = ( 43.0/11.0) * desiredDelays[0];
      desiredDelays[ 7] = ( 47.0/11.0) * desiredDelays[0];
      desiredDelays[ 8] = ( 59.0/11.0) * desiredDelays[0];
      desiredDelays[ 9] = ( 61.0/11.0) * desiredDelays[0];
      desiredDelays[10] = ( 71.0/11.0) * desiredDelays[0];
      desiredDelays[11] = ( 73.0/11.0) * desiredDelays[0];
      desiredDelays[12] = ( 83.0/11.0) * desiredDelays[0];
      desiredDelays[13] = ( 89.0/11.0) * desiredDelays[0];
      desiredDelays[14] = (101.0/11.0) * desiredDelays[0];
      desiredDelays[15] = (103.0/11.0) * desiredDelays[0];
    }
    break;
    */
   

  }
  adjustDelayTimes();

  // todo: compute delay-times according to room acoustic considerations - each delaline can the 
  // thought of modeling standing wvaes in a particular direction - maybe use the dimensions of
  // a rectangular room, the 3 different plane diagonals and the room diagonal - this gives
  // 7 out of 16 delay-times
}


void FeedbackDelayNetwork16::setDelayOrdering(int newDelayOrdering)
{
  if( newDelayOrdering >= 0 &&
    newDelayOrdering <  NUM_DELAY_ORDERINGS )
  {
    delayOrdering = newDelayOrdering;
  }
  adjustDelayTimes();
}

void FeedbackDelayNetwork16::setInjectionVector(int newInjectionVectorIndex)
{
  if( newInjectionVectorIndex >= 0 &&
    newInjectionVectorIndex <  NUM_INJECTION_VECTORS )
  {
    injectionVectorIndex = newInjectionVectorIndex;
  }

  int d;
  switch( newInjectionVectorIndex )
  {
  case IN_ALL_ONES:
    {
      for(d=0; d<maxNumDelayLines; d++)
      {
        injectionVectorL[d] = 1.0;
        injectionVectorR[d] = 1.0;
      }
    }
    break;

    /*
    case IN_ALTERNATING_PLUSMINUS_01:
    {
    injectionVectorL[0] = +1.0;
    injectionVectorR[0] = +1.0;
    injectionVectorL[1] = -1.0;
    injectionVectorR[1] = -1.0;
    injectionVectorL[2] = +1.0;
    injectionVectorR[2] = +1.0;
    injectionVectorL[3] = -1.0;
    injectionVectorR[3] = -1.0;
    injectionVectorL[4] = +1.0;
    injectionVectorR[4] = +1.0;
    injectionVectorL[5] = -1.0;
    injectionVectorR[5] = -1.0;
    injectionVectorL[6] = +1.0;
    injectionVectorR[6] = +1.0;
    injectionVectorL[7] = -1.0;
    injectionVectorR[7] = -1.0;
    }
    break;
    */
  } // end of switch



  /*
  // normalize the input gains with respect to the delayline-lengths in order to ensure that the
  // height of the modes of the individual combs start at the same initial level:
  double normalizer;
  double refDelayInSamplesDbl = (0.001*referenceDelayTime*sampleRate);
  int    refDelayInSamples    = PrimeNumbers::findClosestPrime((int) refDelayInSamplesDbl);
  for(d=0; d<numDelayLines; d++)
  {
    normalizer = (double) delaysInSamples[d] / (double) refDelayInSamples;
    injectionVectorL[d] *= normalizer;
    injectionVectorR[d] *= normalizer;
  }
  */



  /*
  double normalizer;
  for(d=1; d<numDelayLines; d++)
  {
    normalizer = (double) delaysInSamples[0] / (double) delaysInSamples[d];
    injectionVectorL[d] *= normalizer;
    injectionVectorR[d] *= normalizer;
  }
  */




}

void FeedbackDelayNetwork16::setFeedbackMatrix(int newFeedbackMatrixIndex)
{
  if( newFeedbackMatrixIndex >= 0 && newFeedbackMatrixIndex <  NUM_FEEDBACK_MATRICES )
    feedbackMatrixIndex = newFeedbackMatrixIndex;
}

void FeedbackDelayNetwork16::setOutputVector(int newOutputVectorIndex)
{
  // For good stereo effect, it seems desirable that the output vectors should be orthogonal and
  // dense (i guess)

  if( newOutputVectorIndex >= 0 && newOutputVectorIndex <  NUM_OUTPUT_VECTORS )
    outputVectorIndex = newOutputVectorIndex;

  double c;
  int    d;

  // first, we set all output weights to zero:
  for(d=0; d<maxNumDelayLines; d++)
  {
    outputVectorL[d] = 0.0;
    outputVectorR[d] = 0.0;
  }

  switch( newOutputVectorIndex )
  {

  case OUT_ALL_ONES:
    {
      outputVectorL[ 0] = +1.0;
      outputVectorL[ 1] = +1.0;
      outputVectorL[ 2] = +1.0;
      outputVectorL[ 3] = +1.0;
      outputVectorL[ 4] = +1.0;
      outputVectorL[ 5] = +1.0;
      outputVectorL[ 6] = +1.0;
      outputVectorL[ 7] = +1.0;
      outputVectorL[ 8] = +1.0;
      outputVectorL[ 9] = +1.0;
      outputVectorL[10] = +1.0;
      outputVectorL[11] = +1.0;
      outputVectorL[12] = +1.0;
      outputVectorL[13] = +1.0;
      outputVectorL[14] = +1.0;
      outputVectorL[15] = +1.0;

      outputVectorR[ 0] = +1.0;
      outputVectorR[ 1] = +1.0;
      outputVectorR[ 2] = +1.0;
      outputVectorR[ 3] = +1.0;
      outputVectorR[ 4] = +1.0;
      outputVectorR[ 5] = +1.0;
      outputVectorR[ 6] = +1.0;
      outputVectorR[ 7] = +1.0;
      outputVectorR[ 8] = +1.0;
      outputVectorR[ 9] = +1.0;
      outputVectorR[10] = +1.0;
      outputVectorR[11] = +1.0;
      outputVectorR[12] = +1.0;
      outputVectorR[13] = +1.0;
      outputVectorR[14] = +1.0;
      outputVectorR[15] = +1.0;
    }
    break;

  case OUT_01:
    {
      outputVectorL[ 0] = +1.0;
      outputVectorL[ 1] =  0.0;
      outputVectorL[ 2] = +1.0;
      outputVectorL[ 3] =  0.0;
      outputVectorL[ 4] = +1.0;
      outputVectorL[ 5] =  0.0;
      outputVectorL[ 6] = +1.0;
      outputVectorL[ 7] =  0.0;
      outputVectorL[ 8] = +1.0;
      outputVectorL[ 9] =  0.0;
      outputVectorL[10] = +1.0;
      outputVectorL[11] =  0.0;
      outputVectorL[12] = +1.0;
      outputVectorL[13] =  0.0;
      outputVectorL[14] = +1.0;
      outputVectorL[15] =  0.0;

      outputVectorR[ 0] =  0.0;
      outputVectorR[ 1] = +1.0;
      outputVectorR[ 2] =  0.0;
      outputVectorR[ 3] = +1.0;
      outputVectorR[ 4] =  0.0;
      outputVectorR[ 5] = +1.0;
      outputVectorR[ 6] =  0.0;
      outputVectorR[ 7] = +1.0;
      outputVectorR[ 8] =  0.0;
      outputVectorR[ 9] = +1.0;
      outputVectorR[10] =  0.0;
      outputVectorR[11] = +1.0;
      outputVectorR[12] =  0.0;
      outputVectorR[13] = +1.0;
      outputVectorR[14] =  0.0;
      outputVectorR[15] = +1.0;
    }
    break;

  case OUT_02:
    {
      c = (0.001*referenceDelayTime*sampleRate); // reference delaytime in samples

      outputVectorL[ 0] = c / (double)delaysInSamples[ 0];
      outputVectorL[ 1] =  0.0;
      outputVectorL[ 2] = c / (double)delaysInSamples[ 2];
      outputVectorL[ 3] =  0.0;
      outputVectorL[ 4] = c / (double)delaysInSamples[ 4];
      outputVectorL[ 5] =  0.0;
      outputVectorL[ 6] = c / (double)delaysInSamples[ 6];
      outputVectorL[ 7] =  0.0;
      outputVectorL[ 8] = c / (double)delaysInSamples[ 8];
      outputVectorL[ 9] =  0.0;
      outputVectorL[10] = c / (double)delaysInSamples[10];
      outputVectorL[11] =  0.0;
      outputVectorL[12] = c / (double)delaysInSamples[12];
      outputVectorL[13] =  0.0;
      outputVectorL[14] = c / (double)delaysInSamples[14];
      outputVectorL[15] =  0.0;

      outputVectorR[ 0] =  0.0;
      outputVectorR[ 1] = c / (double)delaysInSamples[ 1];;
      outputVectorR[ 2] =  0.0;
      outputVectorR[ 3] = c / (double)delaysInSamples[ 3];;
      outputVectorR[ 4] =  0.0;
      outputVectorR[ 5] = c / (double)delaysInSamples[ 5];;
      outputVectorR[ 6] =  0.0;
      outputVectorR[ 7] = c / (double)delaysInSamples[ 7];;
      outputVectorR[ 8] =  0.0;
      outputVectorR[ 9] = c / (double)delaysInSamples[ 9];;
      outputVectorR[10] =  0.0;
      outputVectorR[11] = c / (double)delaysInSamples[11];;
      outputVectorR[12] =  0.0;
      outputVectorR[13] = c / (double)delaysInSamples[13];;
      outputVectorR[14] =  0.0;
      outputVectorR[15] = c / (double)delaysInSamples[15];;
    }
    break;

  case OUT_03:
    {
      outputVectorL[ 0] = +1.0;
      outputVectorL[ 1] = -0.5;
      outputVectorL[ 2] = +1.0;
      outputVectorL[ 3] = +0.5;
      outputVectorL[ 4] = +1.0;
      outputVectorL[ 5] = -0.5;
      outputVectorL[ 6] = +1.0;
      outputVectorL[ 7] = +0.5;
      outputVectorL[ 8] = +1.0;
      outputVectorL[ 9] = -0.5;
      outputVectorL[10] = +1.0;
      outputVectorL[11] = +0.5;
      outputVectorL[12] = +1.0;
      outputVectorL[13] = -0.5;
      outputVectorL[14] = +1.0;
      outputVectorL[15] = +0.5;

      outputVectorR[ 0] = +0.5;
      outputVectorR[ 1] = +1.0;
      outputVectorR[ 2] = -0.5;
      outputVectorR[ 3] = +1.0;
      outputVectorR[ 4] = +0.5;
      outputVectorR[ 5] = +1.0;
      outputVectorR[ 6] = -0.5;
      outputVectorR[ 7] = +1.0;
      outputVectorR[ 8] = +0.5;
      outputVectorR[ 9] = +1.0;
      outputVectorR[10] = -0.5;
      outputVectorR[11] = +1.0;
      outputVectorR[12] = +0.5;
      outputVectorR[13] = +1.0;
      outputVectorR[14] = -0.5;
      outputVectorR[15] = +1.0;
    }
    break;



  case OUT_04:
    {
      outputVectorL[ 0] = +1.0;
      outputVectorL[ 1] = +0.5;
      outputVectorL[ 2] = -1.0;
      outputVectorL[ 3] = +0.5;
      outputVectorL[ 4] = +1.0;
      outputVectorL[ 5] = +0.5;
      outputVectorL[ 6] = -1.0;
      outputVectorL[ 7] = +0.5;
      outputVectorL[ 8] = +1.0;
      outputVectorL[ 9] = +0.5;
      outputVectorL[10] = -1.0;
      outputVectorL[11] = +0.5;
      outputVectorL[12] = +1.0;
      outputVectorL[13] = +0.5;
      outputVectorL[14] = -1.0;
      outputVectorL[15] = +0.5;

      outputVectorR[ 0] = +0.5;
      outputVectorR[ 1] = -1.0;
      outputVectorR[ 2] = +0.5;
      outputVectorR[ 3] = +1.0;
      outputVectorR[ 4] = +0.5;
      outputVectorR[ 5] = -1.0;
      outputVectorR[ 6] = +0.5;
      outputVectorR[ 7] = +1.0;
      outputVectorR[ 8] = +0.5;
      outputVectorR[ 9] = -1.0;
      outputVectorR[10] = +0.5;
      outputVectorR[11] = +1.0;
      outputVectorR[12] = +0.5;
      outputVectorR[13] = -1.0;
      outputVectorR[14] = +0.5;
      outputVectorR[15] = +1.0;
    }
    break;



    /*
    case OUT_02:
    {
    outputVectorL[0] = +1.0;
    outputVectorL[1] = +1.0;
    outputVectorL[2] = +1.0;
    outputVectorL[3] = +1.0;
    outputVectorL[4] = +1.0;
    outputVectorL[5] = +1.0;
    outputVectorL[6] = +1.0;
    outputVectorL[7] = +1.0;

    outputVectorR[0] = -1.0;
    outputVectorR[1] = -1.0;
    outputVectorR[2] = -1.0;
    outputVectorR[3] = -1.0;
    outputVectorR[4] = +1.0;
    outputVectorR[5] = +1.0;
    outputVectorR[6] = +1.0;
    outputVectorR[7] = +1.0;
    }
    break;

    case OUT_03:
    {
    outputVectorL[0] = +1.0;
    outputVectorL[1] = +0.5;
    outputVectorL[2] = +1.0;
    outputVectorL[3] = +0.5;
    outputVectorL[4] = +1.0;
    outputVectorL[5] = +0.5;
    outputVectorL[6] = +1.0;
    outputVectorL[7] = +0.5;

    outputVectorR[0] = -0.5;
    outputVectorR[1] = -1.0;
    outputVectorR[2] = -0.5;
    outputVectorR[3] = -1.0;
    outputVectorR[4] = -0.5;
    outputVectorR[5] = -1.0;
    outputVectorR[6] = -0.5;
    outputVectorR[7] = -1.0;
    }
    break;

    case OUT_04:
    {
    outputVectorL[0] = +1.0;
    outputVectorL[1] =  0.0;
    outputVectorL[2] = -1.0;
    outputVectorL[3] =  0.0;
    outputVectorL[4] = +1.0;
    outputVectorL[5] =  0.0;
    outputVectorL[6] = -1.0;
    outputVectorL[7] =  0.0;

    outputVectorR[0] =  0.0;
    outputVectorR[1] = -1.0;
    outputVectorR[2] =  0.0;
    outputVectorR[3] = +1.0;
    outputVectorR[4] =  0.0;
    outputVectorR[5] = -1.0;
    outputVectorR[6] =  0.0;
    outputVectorR[7] = +1.0;
    }
    break;

    case OUT_05:
    {
    outputVectorL[0] = +1.0;
    outputVectorL[1] =  0.0;
    outputVectorL[2] = -1.0;
    outputVectorL[3] =  0.0;
    outputVectorL[4] = +1.0;
    outputVectorL[5] =  0.0;
    outputVectorL[6] = -1.0;
    outputVectorL[7] =  0.0;

    outputVectorR[0] =  0.0;
    outputVectorR[1] = +1.0;
    outputVectorR[2] =  0.0;
    outputVectorR[3] = -1.0;
    outputVectorR[4] =  0.0;
    outputVectorR[5] = +1.0;
    outputVectorR[6] =  0.0;
    outputVectorR[7] = -1.0;
    }
    break;
    */

  } // end of switch

  /*
  // the input gains are normalized with respect to the delayline-lengths in
  // order to ensure that the height of the modes of the individual combs
  // start at the same initial level. to compensate for the undesired side
  // effect on the impulse-response we apply reciprocal factors in the
  // output-vector:
  //double c = 1.0/sqrt((double)numDelayLines);
  double c = 0.0;
  for(d=0; d<numDelayLines; d++)
  c += outputVectorL[d]*outputVectorL[d];
  c = 1/sqrt(c);

  double normalizer;
  double maxDelayInSamplesDbl = (0.001*maxDelayTime*sampleRate);
  int    maxDelayInSamples
  = PrimeNumbers::findClosestPrime((int) maxDelayInSamplesDbl);
  for(d=0; d<numDelayLines; d++)
  {
  normalizer = c * (double) maxDelayInSamples / (double) delaysInSamples[d];

  outputVectorL[d] *= normalizer;
  outputVectorR[d] *= normalizer;
  }
  */








  double normalizer;
  for(d=1; d<numDelayLines; d++)
  {
    normalizer = (double) delaysInSamples[0] / (double) delaysInSamples[d];
    outputVectorL[d] *= normalizer;
    outputVectorR[d] *= normalizer;
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

void FeedbackDelayNetwork16::getImpulseResponse(double* /*bufferL*/, double* /*bufferR*/,
                                                int /*length*/)
{

}


//-------------------------------------------------------------------------------------------------
// others:

void FeedbackDelayNetwork16::adjustReadPointer(int index)
{
  // calculate the delay in samples and put the read-pointer that number of
  // samples behind the write-pointer:
  tapOuts[index] = tapIns[index] - delaysInSamples[index];

  // now we must take care, that the read-pointer is not below zero, that is,
  // we need to do a forward-wraparound:
  while( tapOuts[index] < 0 )
    tapOuts[index] += maxDelayInSamples;
}

void FeedbackDelayNetwork16::adjustDelayTimes()
{
  int tmpDelays[maxNumDelayLines];
  bool usePrimes = true;
  double factor  = 0.001*referenceDelayTime*sampleRate;
  int i;
  for(i=0; i<numDelayLines; i++)
  {
    tmpDelays[i] = roundToInt(factor*relativeDelayTimes[i]);
    if( usePrimes == true )
      tmpDelays[i] = PrimeNumbers::findClosestPrime(tmpDelays[i]);
  }

  // tmpDelay now contains the delay-times in ascending order. now re-order the delay-times to the
  // desired ordering/permutation:
  switch( delayOrdering )
  {
  case ASCENDING:
    {
      for(i=0; i<numDelayLines; i++)
        delaysInSamples[i] = tmpDelays[i];
    }
    break;
  case DESCENDING:
    {
      for(i=0; i<numDelayLines; i++)
        delaysInSamples[i] = tmpDelays[numDelayLines-i-1];
    }
    break;
  case ALTERNATING:
    {
      delaysInSamples[ 0] = tmpDelays[ 0];
      delaysInSamples[ 1] = tmpDelays[15];
      delaysInSamples[ 2] = tmpDelays[ 1];
      delaysInSamples[ 3] = tmpDelays[14];
      delaysInSamples[ 4] = tmpDelays[ 2];
      delaysInSamples[ 5] = tmpDelays[13];
      delaysInSamples[ 6] = tmpDelays[ 3];
      delaysInSamples[ 7] = tmpDelays[12];
      delaysInSamples[ 8] = tmpDelays[ 4];
      delaysInSamples[ 9] = tmpDelays[11];
      delaysInSamples[10] = tmpDelays[ 5];
      delaysInSamples[11] = tmpDelays[10];
      delaysInSamples[12] = tmpDelays[ 6];
      delaysInSamples[13] = tmpDelays[ 9];
      delaysInSamples[14] = tmpDelays[ 7];
      delaysInSamples[15] = tmpDelays[ 8];
    }
    break;
  } // end of switch( delayOrdering )

  // the new delay-times have been calculated (in samples), we now need to update the relative
  // positions of the read-pointers
  for(i=0; i<numDelayLines; i++)
    adjustReadPointer(i);

  // we also need to re-calculate the injection- and output-vectors, as they take into account the
  // delay-times:
  setInjectionVector(injectionVectorIndex);
  setOutputVector(outputVectorIndex);

  // changing the lengths of the delay-lines also affects the overall reverberation time, so we
  // need to update the damping-filters to compensate this:
  updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork16::updateDampingAndCorrectionFilters()
{
  double T_l, T_lb, T_m, T_hb, T_h;    // desired reverb-times at 5 frequencies
  double g_l_abs, g_m_abs, g_h_abs;    // absolute gain-factors for the 3 bands
  double g_lb_abs, g_hb_abs;           // absolute gains at which the bandwidths are measured
                                       // (using the terminology from the Orfanidis-papaer)
  double g_l_rel, g_h_rel;             // gains for low and high band relative to the mid-band
  double g_lb_rel, g_hb_rel;           // relative gains at the bandwidth measurement frequencies


  T_l  = lowReverbTimeScale*midReverbTime;   // desired reverb-time at low frequencies
  T_m  = midReverbTime;                      // desired reverb-time at mid frequencies
  T_h  = highReverbTimeScale*midReverbTime;  // desired reverb-time at high frequencies
  T_lb = sqrt(T_l*T_m);                      // desired reverb-time at the low crossover frequency
  T_hb = sqrt(T_h*T_m);                      // desired reverb-time at the high crossover frequency

  for(int d=0; d<numDelayLines; d++)
  {
    // calculate the required absolute gain-factors at 5 frequencies:
    g_l_abs  = pow(10.0, -3.0*delaysInSamples[d]/(T_l* sampleRate) );
    g_lb_abs = pow(10.0, -3.0*delaysInSamples[d]/(T_lb*sampleRate) );
    g_m_abs  = pow(10.0, -3.0*delaysInSamples[d]/(T_m* sampleRate) );
    g_hb_abs = pow(10.0, -3.0*delaysInSamples[d]/(T_hb*sampleRate) );
    g_h_abs  = pow(10.0, -3.0*delaysInSamples[d]/(T_h* sampleRate) );

    // the desired absolute gains will be approximately realized by a low- and high-shelving filter
    // and a global gain, where the global gain is equal to the desired mid-frequency gain and the
    // shelves are set up to account for the relative gain-deviation:
    g_l_rel  = g_l_abs  / g_m_abs;
    g_lb_rel = g_lb_abs / g_m_abs;
    g_h_rel  = g_h_abs  / g_m_abs;
    g_hb_rel = g_hb_abs / g_m_abs;

    // set up the damping-filter for delayline d:
    dampingFilters[d].setGlobalGainFactor(g_m_abs);
    dampingFilters[d].setLowCrossoverFreq(lowCrossoverFreq);
    dampingFilters[d].setLowCrossoverGainFactor(g_lb_rel);
    dampingFilters[d].setLowGainFactor(g_l_rel);
    dampingFilters[d].setHighCrossoverFreq(highCrossoverFreq);
    dampingFilters[d].setHighCrossoverGainFactor(g_hb_rel);
    dampingFilters[d].setHighGainFactor(g_h_rel);

    // store the mid-gain factor in a member-variable because we need to know that value in the
    // allpass-mode:
    feedbackGains[d] = g_m_abs;
  }

  // set up the correction filters which decouple the overall frequency response from the frequency
  // dependent decay times:
  correctionFilterL.setGlobalGainFactor(1.0/sqrt(T_m));
  correctionFilterR.setGlobalGainFactor(1.0/sqrt(T_m));
  correctionFilterL.setLowGainFactor(1.0/sqrt(T_l/T_m));
  correctionFilterR.setLowGainFactor(1.0/sqrt(T_l/T_m));
  correctionFilterL.setLowCrossoverGainFactor(1.0/sqrt(T_lb/T_m));
  correctionFilterR.setLowCrossoverGainFactor(1.0/sqrt(T_lb/T_m));
  correctionFilterL.setHighGainFactor(1.0/sqrt(T_h/T_m));
  correctionFilterR.setHighGainFactor(1.0/sqrt(T_h/T_m));
  correctionFilterL.setHighCrossoverGainFactor(1.0/sqrt(T_hb/T_m));
  correctionFilterR.setHighCrossoverGainFactor(1.0/sqrt(T_hb/T_m));
  correctionFilterL.setLowCrossoverFreq(lowCrossoverFreq);
  correctionFilterR.setLowCrossoverFreq(lowCrossoverFreq);
  correctionFilterL.setHighCrossoverFreq(highCrossoverFreq);
  correctionFilterR.setHighCrossoverFreq(highCrossoverFreq);
}

void FeedbackDelayNetwork16::reset()
{
  int d,n;   // index for the delayline and sample-number
  for(d=0; d<maxNumDelayLines; d++)
  {
    // reset the feedback damping filters:
    dampingFilters[d].reset();

    // reset the output-slots:
    delayLineOuts[d] = 0.0;

    // set content of the delaylines to zero:
    for(n=0; n<maxDelayInSamples; n++)
      delayLines[d][n] = 0.0;
  }

  // reset the correction and output-filters:
  correctionFilterL.reset();
  correctionFilterR.reset();
  wetFilterL.reset();
  wetFilterR.reset();
}


/*
Idea for producing relative delayline lengths:
start with an irrational number r, like, for example: r = (1+sqrt(5))/2
(1) take r as 1st relative delay-time
(2) square it: r = r^2  ...or maybe r = r^2 + c for c some additive constant (like 0.2)
(3) while r > 2: r = r/2
    -use r as next delay time
(4) back to (2)

can use different seed numbers (can be user parameter)
-maybe, when all times are computed, the should be sorted ascending (but may be optional)
*/