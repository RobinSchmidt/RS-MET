#include "FeedbackDelayNetwork8.h"

//----------------------------------------------------------------------------
//construction/destruction:

FeedbackDelayNetwork8::FeedbackDelayNetwork8()
{
 numDelayLines = 8;
	setSampleRate(44100.0);      // sampleRate = 44100 Hz by default

 // initialize the pointers to the beginnings of the delaylines:
 for(int d=0; d < maxNumDelayLines; d++)
 {
  tapIns[d]  = 0;
  tapOuts[d] = 0;
 }

 // initialize the delay-times:
 minDelayTime           = 30.0;        // 30 ms
 maxDelayTime           = 45.0;        // 45 ms
 delayDistributionIndex = EXPONENTIAL;
 delayOrdering          = ASCENDING;
 delayDistributionForm  = 0.5;
 adjustDelayTimes();

 // initialize the damping- and correction-filters:
 setLowReverbTimeScale(2.0);
 setLowCrossoverFreq(250.0);
 setMidReverbTime(4.0);
 setHighReverbTimeScale(0.25);
 setHighCrossoverFreq(4000.0);

 // setup the FDN-parameters:
 setInjectionVector(IN_ALL_ONES);
 setFeedbackMatrix(HADAMARD);
 setOutputVector(OUT_04);
 allpassMode      = false;

 // setup the output filters:
 setWetLpfCutoff(20000.0);
 setWetHpfCutoff(20.0);

 // setup the output parameters:
 wetPinking       = false;
 stereoSwapSwitch = false;
 dryVol           = 0.0;
 wetVol           = 1.0;

 // reset the content of the delaylines samples to zero and trigger a reset
 // in the embedded modules:
 resetBuffers();         
}

FeedbackDelayNetwork8::~FeedbackDelayNetwork8()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void FeedbackDelayNetwork8::setSampleRate(double newSampleRate)
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

void FeedbackDelayNetwork8::setLowReverbTimeScale(double newLowReverbTimeScale)
{
 if( newLowReverbTimeScale > 0.0 )
  lowReverbTimeScale = newLowReverbTimeScale;
 updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork8::setLowCrossoverFreq(double newLowCrossoverFreq)
{
 if( newLowCrossoverFreq >= 20.0 && newLowCrossoverFreq <= 20000.0 )
  lowCrossoverFreq = newLowCrossoverFreq;
 updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork8::setMidReverbTime(double newMidReverbTime)
{
 if( newMidReverbTime > 0.0 )
  midReverbTime = newMidReverbTime;
 updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork8::setHighReverbTimeScale(double newHighReverbTimeScale)
{
 if( newHighReverbTimeScale > 0.0 )
  highReverbTimeScale = newHighReverbTimeScale;
 updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork8::setHighCrossoverFreq(double newHighCrossoverFreq)
{
 if( newHighCrossoverFreq >= 20.0 && newHighCrossoverFreq <= 20000.0 )
  highCrossoverFreq = newHighCrossoverFreq;
 updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork8::setDensity(int newDensity)
{

}

void FeedbackDelayNetwork8::setAllpassMode(bool shouldBeAllpass)
{
 allpassMode = shouldBeAllpass;
}

void FeedbackDelayNetwork8::setDryVolume(double newDryVolume)
{
 dryVol = dB2amp(newDryVolume);
}

void FeedbackDelayNetwork8::setWetVolume(double newWetVolume)
{
 wetVol = dB2amp(newWetVolume);
}

void FeedbackDelayNetwork8::setWetLpfCutoff(double newWetLpfCutoff)
{
 wetFilterL.setLpfCutoff(newWetLpfCutoff);
 wetFilterR.setLpfCutoff(newWetLpfCutoff);
}

void FeedbackDelayNetwork8::setWetHpfCutoff(double newWetHpfCutoff)
{
 wetFilterL.setHpfCutoff(newWetHpfCutoff);
 wetFilterR.setHpfCutoff(newWetHpfCutoff);
}

void FeedbackDelayNetwork8::setStereoSwapSwitch(bool newStereoSwapSwitch)
{
 stereoSwapSwitch = newStereoSwapSwitch;
}

void FeedbackDelayNetwork8::setWetPinkingSwitch(bool newWetPinkingSwitch)
{
 wetPinking = newWetPinkingSwitch;
}

void FeedbackDelayNetwork8::setMinDelayTime(double newMinDelayTime)
{
 if( newMinDelayTime >= 0.1 && 
     newMinDelayTime <= 2000.0  && 
     newMinDelayTime <  maxDelayTime-1.0 ) // should be at least 1 millisecond
                                           // shorter than the maximum
  minDelayTime = newMinDelayTime;

 adjustDelayTimes();
}

void FeedbackDelayNetwork8::setMaxDelayTime(double newMaxDelayTime)
{
 if( newMaxDelayTime >= 0.1 && 
     newMaxDelayTime <= 2000.0  && 
     newMaxDelayTime >  minDelayTime+1.0 ) // should be at least 1 millisecond
                                           // longer than the minimum
  maxDelayTime = newMaxDelayTime;

 adjustDelayTimes();
}

void FeedbackDelayNetwork8::setDelayDistributionIndex(int newDelayDistributionIndex)
{
 if( newDelayDistributionIndex >= 0 &&
     newDelayDistributionIndex <  NUM_DELAY_DISTRIBUTIONS )
 {
  delayDistributionIndex = newDelayDistributionIndex;
 }
 adjustDelayTimes();
}

void FeedbackDelayNetwork8::setDelayDistributionForm(double newDelayDistributionForm)
{
 if( newDelayDistributionForm >= 0.1 &&
     newDelayDistributionForm <= 0.9 )
 {
  delayDistributionForm = newDelayDistributionForm;
 }
 adjustDelayTimes();
}

void FeedbackDelayNetwork8::setDelayOrdering(int newDelayOrdering)
{
 if( newDelayOrdering >= 0 &&
     newDelayOrdering <  NUM_DELAY_ORDERINGS )
 {
  delayOrdering = newDelayOrdering;
 }
 adjustDelayTimes();
}

void FeedbackDelayNetwork8::setInjectionVector(int newInjectionVectorIndex)
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
 } // end of switch

 // normalize the input gains with respect to the delayline-lengths in order
 // to ensure that the height of the modes of the individual combs start at 
 // the same initial level:
 double normalizer;
 double maxDelayInSamplesDbl = (0.001*maxDelayTime*sampleRate);
 int    maxDelayInSamples    
         = primeNumbers.findClosestLowerPrime((int) maxDelayInSamplesDbl);
 for(d=0; d<numDelayLines; d++)
 {
  normalizer = (double) delaysInSamples[d] / (double) maxDelayInSamples; 

  injectionVectorL[d] *= normalizer;
  injectionVectorR[d] *= normalizer;
 }
}

void FeedbackDelayNetwork8::setFeedbackMatrix(int newFeedbackMatrixIndex)
{
 if( newFeedbackMatrixIndex >= 0 && 
     newFeedbackMatrixIndex <  NUM_FEEDBACK_MATRICES )
 {
  feedbackMatrixIndex = newFeedbackMatrixIndex;
 }
}

void FeedbackDelayNetwork8::setOutputVector(int newOutputVectorIndex)
{
 if( newOutputVectorIndex >= 0 && 
     newOutputVectorIndex <  NUM_OUTPUT_VECTORS )
 {
  outputVectorIndex = newOutputVectorIndex;
 }

 int    d;
 switch( newOutputVectorIndex )
 {
 case OUT_01:
  {
   for(d=0; d<maxNumDelayLines; d++)
   {
   outputVectorL[0] = +1.0;
   outputVectorL[1] = +1.0;
   outputVectorL[2] = +1.0;
   outputVectorL[3] = +1.0;
   outputVectorL[4] = +1.0;
   outputVectorL[5] = +1.0;
   outputVectorL[6] = +1.0;
   outputVectorL[7] = +1.0;

   outputVectorR[0] = +1.0;
   outputVectorR[1] = -1.0;
   outputVectorR[2] = +1.0;
   outputVectorR[3] = -1.0;
   outputVectorR[4] = +1.0;
   outputVectorR[5] = -1.0;
   outputVectorR[6] = +1.0;
   outputVectorR[7] = -1.0;
   }
  }
  break;

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



 } // end of switch

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
         = primeNumbers.findClosestLowerPrime((int) maxDelayInSamplesDbl);
 for(d=0; d<numDelayLines; d++)
 {
  normalizer = c * (double) maxDelayInSamples / (double) delaysInSamples[d]; 

  outputVectorL[d] *= normalizer;
  outputVectorR[d] *= normalizer;
 }
}

//----------------------------------------------------------------------------
//others:

void FeedbackDelayNetwork8::adjustReadPointer(int index)
{
 // calculate the delay in samples and put the read-pointer that number of
 // samples behind the write-pointer:
 tapOuts[index] = tapIns[index] - delaysInSamples[index];

 // now we must take care, that the read-pointer is not below zero, that is, 
 // we need to do a forward-wraparound:
 while( tapOuts[index] < 0 )
  tapOuts[index] += maxDelayInSamples;
}

void FeedbackDelayNetwork8::adjustDelayTimes()
{
 // allocate an array to hold the desired delay-times in samples which are not
 // yet 'rounded' to the next lower prime integer:
 double desiredDelays[maxNumDelayLines];

 // calculate the shortest and the longest delayline-length in samples from 
 // the specified minimum and maximum delay-time in milliseconds:
 desiredDelays[0]               = (0.001*minDelayTime*sampleRate);
 desiredDelays[numDelayLines-1] = (0.001*maxDelayTime*sampleRate);

 // calculate the desired intermediate delayline-lengths:
 int i; // for indexing the delay-time
 switch( delayDistributionIndex )
 {
 case LINEAR:
  {
   double offset = (desiredDelays[numDelayLines-1]-desiredDelays[0]) / 
                   (double) (numDelayLines-1);

   for(i=1; i <= numDelayLines-2; i++)
    desiredDelays[i] = desiredDelays[0] + i * offset;
  }
  break;
 case EXPONENTIAL:
  {
   double ratio = maxDelayTime/minDelayTime;

   for(i=1; i <= numDelayLines-2; i++)
    desiredDelays[i] = desiredDelays[0] * 
                       pow(ratio, (double) i / (double) (numDelayLines-1) );
  }
  break;
 case DIVISION_ALGO_1:
  {
   // this algorithm assigns the delay-times with some kind of 
   // recursive subdivision-algorithm

   double c = delayDistributionForm;
   double d_mid, delta;

   delta = desiredDelays[numDelayLines-1] - desiredDelays[0];
   d_mid = desiredDelays[0] + c*delta;

   delta = d_mid - desiredDelays[0];
   desiredDelays[2] = desiredDelays[0] + c*delta;

   delta = desiredDelays[2] - desiredDelays[0];
   desiredDelays[1] = desiredDelays[0] + c*delta;

   delta = d_mid - desiredDelays[2];
   desiredDelays[3] = desiredDelays[2] + c*delta;

   delta = desiredDelays[7] - d_mid;
   desiredDelays[5] = d_mid + c*delta;

   delta = desiredDelays[5] - d_mid;
   desiredDelays[4] = d_mid + c*delta;

   delta = desiredDelays[7] - desiredDelays[5];
   desiredDelays[6] = desiredDelays[5] + c*delta;
  }
  break;

 case PRIME_ALGO_1:
  {
   desiredDelays[1] = (17.0/11.0) * desiredDelays[0];
   desiredDelays[2] = (23.0/11.0) * desiredDelays[0];
   desiredDelays[3] = (31.0/11.0) * desiredDelays[0];
   desiredDelays[4] = (41.0/11.0) * desiredDelays[0];
   desiredDelays[5] = (47.0/11.0) * desiredDelays[0];
   desiredDelays[6] = (59.0/11.0) * desiredDelays[0];
   desiredDelays[7] = (67.0/11.0) * desiredDelays[0];
  }
  break;

 case PRIME_ALGO_2:
  {
   desiredDelays[1] = (13.0/11.0) * desiredDelays[0];
   desiredDelays[2] = (19.0/11.0) * desiredDelays[0];
   desiredDelays[3] = (23.0/11.0) * desiredDelays[0];
   desiredDelays[4] = (31.0/11.0) * desiredDelays[0];
   desiredDelays[5] = (37.0/11.0) * desiredDelays[0];
   desiredDelays[6] = (43.0/11.0) * desiredDelays[0];
   desiredDelays[7] = (47.0/11.0) * desiredDelays[0];
  }
  break;

 } // end of switch( delayDistributionIndex ) 

 // the desired delay-times have been calculated as double-values in ascending
 // order. now we determine the actual integer delay-times as the closest 
 // lower primes below these values, thereby we also re-order the delay-times
 // to the desired ordering:
 switch( delayOrdering )
 {
 case ASCENDING:
  {
   for(i=0; i<numDelayLines; i++)
    delaysInSamples[i] 
     = primeNumbers.findClosestLowerPrime((int) desiredDelays[i]);
  }
  break;
 case DESCENDING:
  {
   for(i=0; i<numDelayLines; i++)
    delaysInSamples[i] 
     = primeNumbers.findClosestLowerPrime
       (
        (int) desiredDelays[numDelayLines-i-1] 
       );
  }
  break;
 case ALTERNATING:
  {
   delaysInSamples[0] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[0]);
   delaysInSamples[1] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[7]);
   delaysInSamples[2] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[1]);
   delaysInSamples[3] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[6]);
   delaysInSamples[4] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[2]);
   delaysInSamples[5] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[5]);
   delaysInSamples[6] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[3]);
   delaysInSamples[7] 
    = primeNumbers.findClosestLowerPrime((int) desiredDelays[4]);
  }
  break;
 } // end of switch( delayOrdering ) 

 // the new delay-times have been calculated (in samples), we now need to
 // update the relative positions of the read-pointers
 for(i=0; i<numDelayLines; i++)
  adjustReadPointer(i);

 // we also need to re-calculate the injection- and output-vectors, as they
 // take into account the delay-times:
 setInjectionVector(injectionVectorIndex);
 setOutputVector(outputVectorIndex);

 // changing the lengths of the delay-lines also affects the overall 
 // reverberation time, so we need to update the damping-filters to compensate
 // this:
 updateDampingAndCorrectionFilters();
}

void FeedbackDelayNetwork8::updateDampingAndCorrectionFilters()
{
 double T_l, T_lb, T_m, T_hb, T_h;             
  // the desired reverb-times at 5 frequencies

 double g_l_abs, g_m_abs, g_h_abs;  
  // absolute gain-factors for the 3 bands

 double g_lb_abs, g_hb_abs;     
  // absolute gains at which the bandwidths are measured (using the 
  // terminology from the Orfanidis-papaer)

 double g_l_rel, g_h_rel;       
  // gains for low and high band relative to the mid-band

 double g_lb_rel, g_hb_rel;     
  // relative gains at the bandwidth measurement frequencies

 // assign our local reverb-time T-variables for 5 frequencies:
 T_l  = lowReverbTimeScale*midReverbTime;
 T_m  = midReverbTime;
 T_h  = highReverbTimeScale*midReverbTime;
 T_lb = sqrt(T_l*T_m);  // desired reverb-time at the low crossover frequency
 T_hb = sqrt(T_h*T_m);  // desired reverb-time at the high crossover frequency

 for(int d=0; d<numDelayLines; d++)
 {
  // calculate the required absolute gain-factors at 5 frequencies:
  g_l_abs  = pow(10.0, -3.0*delaysInSamples[d]/(T_l* sampleRate) );
  g_lb_abs = pow(10.0, -3.0*delaysInSamples[d]/(T_lb*sampleRate) );
  g_m_abs  = pow(10.0, -3.0*delaysInSamples[d]/(T_m* sampleRate) );
  g_hb_abs = pow(10.0, -3.0*delaysInSamples[d]/(T_hb*sampleRate) );
  g_h_abs  = pow(10.0, -3.0*delaysInSamples[d]/(T_h* sampleRate) );

  // the desired absolute gains will be approximately realized by a low- and 
  // high-shelving filter and a global gain, where the global gain is equal 
  // to the desired mid-frequency gain and the shelves are set up to account 
  // for the relative gain-deviation:
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

  // store the mid-gain factor in a member-variable because we need to know
  // that value in the allpass-mode:
  feedbackGains[d] = g_m_abs;
 }

 //compensationGain = 1.0/sqrt(T_m);
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

void FeedbackDelayNetwork8::resetBuffers()
{
 int d,n;   // index for the delayline and sample-number
 for(d=0; d<maxNumDelayLines; d++)
 {
  // reset the feedback damping filters:
  dampingFilters[d].resetBuffers();

  // reset the output-slots:
  delayLineOuts[d] = 0.0;

  // set content of the delaylines to zero:
  for(n=0; n<maxDelayInSamples; n++)
   delayLines[d][n] = 0.0;
 }

 // reset the correction and output-filters:
 correctionFilterL.resetBuffers();
 correctionFilterR.resetBuffers();
 wetFilterL.resetBuffers();
 wetFilterR.resetBuffers();
}
