//#include "rosic_Reverb.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

rsReverb::rsReverb(int delayMemoryInSamplesToAllocate)
{
  maxDelayInSamples = delayMemoryInSamplesToAllocate / numDelayLines;
  for(int d=0; d<numDelayLines; d++)
    delayLines[d] = new double[maxDelayInSamples];

  sampleRate          = 44100.0;
  referenceDelayTime  = 1000*1000/sampleRate;   // -> 22.67..ms -> 1000 samples
  dryVol              = 0.0;
  wetVol              = 0.25;
  midReverbTime       = 3.0;
  lowCrossoverFreq    = 250.0;
  lowReverbTimeScale  = 1.0;
  highCrossoverFreq   = 4000.0;
  highReverbTimeScale = 0.3f;  
  wetPinking          = true;
  stereoSwapSwitch    = false;

  tapIn = 0;
  for(int d=0; d<numDelayLines; d++)
  {
    tapOuts[d]       = 0;
    delayLineOuts[d] = 0.0;
  }

  setPreDelay(0.0);
  assignRelativeDelayTimes(); 
  setupOutputVector();
  setWetLowpassCutoff(20000.0);
  setWetHighpassCutoff(20.0);
  reset();         
}

rsReverb::~rsReverb()
{
  for(int d=0; d<numDelayLines; d++)
    delete[] delayLines[d];
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void rsReverb::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;

    adjustDelayTimes();

    for(int k=0; k<numDelayLines; k++)
      dampingFilters[k].setSampleRate(sampleRate);

    preDelayLineL.setSampleRate(sampleRate);
    preDelayLineR.setSampleRate(sampleRate);
    correctionFilterL.setSampleRate(sampleRate);
    correctionFilterR.setSampleRate(sampleRate);
    wetFilterL.setSampleRate(sampleRate);
    wetFilterR.setSampleRate(sampleRate);
    updateDampingAndCorrectionFilters();
  }
}

void rsReverb::setLowReverbTimeScale(double newLowReverbTimeScale)
{
  if( newLowReverbTimeScale > 0.0 && newLowReverbTimeScale != lowReverbTimeScale )
  {
    lowReverbTimeScale = newLowReverbTimeScale;
    updateDampingAndCorrectionFilters();
  }
}

void rsReverb::setLowCrossoverFreq(double newLowCrossoverFreq)
{
  if(  newLowCrossoverFreq >= 20.0 && newLowCrossoverFreq <= 20000.0 
    && newLowCrossoverFreq != lowCrossoverFreq )
  {
    lowCrossoverFreq = newLowCrossoverFreq;
    updateDampingAndCorrectionFilters();
  }
}

void rsReverb::setMidReverbTime(double newMidReverbTime)
{
  if( newMidReverbTime > 0.0 && newMidReverbTime != midReverbTime )
  {
    midReverbTime = newMidReverbTime;
    updateDampingAndCorrectionFilters();
  }
}

void rsReverb::setHighReverbTimeScale(double newHighReverbTimeScale)
{
  if( newHighReverbTimeScale > 0.0 && newHighReverbTimeScale != highReverbTimeScale )
  {
    highReverbTimeScale = newHighReverbTimeScale;
    updateDampingAndCorrectionFilters();
  }
}

void rsReverb::setHighCrossoverFreq(double newHighCrossoverFreq)
{
  if(  newHighCrossoverFreq >= 20.0 && newHighCrossoverFreq <= 20000.0 
    && newHighCrossoverFreq != highCrossoverFreq )
  {
    highCrossoverFreq = newHighCrossoverFreq;
    updateDampingAndCorrectionFilters();
  }
}

void rsReverb::setWetLowpassCutoff(double newCutoff)
{
  wetFilterL.setLowpassCutoff(newCutoff);
  wetFilterR.setLowpassCutoff(newCutoff);
}

void rsReverb::setWetHighpassCutoff(double newCutoff)
{
  wetFilterL.setHighpassCutoff(newCutoff);
  wetFilterR.setHighpassCutoff(newCutoff);
}

void rsReverb::setStereoSwapSwitch(bool newStereoSwapSwitch)
{
  stereoSwapSwitch = newStereoSwapSwitch;
}

void rsReverb::setWetPinkingSwitch(bool newWetPinkingSwitch)
{
  wetPinking = newWetPinkingSwitch;
}

void rsReverb::setReferenceDelayTime(double newReferenceDelayTime)
{
  if(  newReferenceDelayTime >= 1.0 && newReferenceDelayTime <= 100.0 
    && newReferenceDelayTime != referenceDelayTime )
  {
    referenceDelayTime = newReferenceDelayTime;
    adjustDelayTimes();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void rsReverb::assignRelativeDelayTimes()
{
  relativeDelayTimes[ 0] = 1.0;    
  double dMax            = 2.4;
  relativeDelayTimes[8]  = sqrt(dMax                   * relativeDelayTimes[ 0] ); 
  relativeDelayTimes[4]  = sqrt(relativeDelayTimes[ 8] * relativeDelayTimes[ 0] );
  relativeDelayTimes[12] = sqrt(dMax                   * relativeDelayTimes[ 8] );
  relativeDelayTimes[2]  = sqrt(relativeDelayTimes[ 4] * relativeDelayTimes[ 0] );
  relativeDelayTimes[6]  = sqrt(relativeDelayTimes[ 8] * relativeDelayTimes[ 4] );
  relativeDelayTimes[10] = sqrt(relativeDelayTimes[12] * relativeDelayTimes[ 8] );
  relativeDelayTimes[14] = sqrt(dMax                   * relativeDelayTimes[12] );
  relativeDelayTimes[1]  = sqrt(relativeDelayTimes[ 2] * relativeDelayTimes[ 0] );
  relativeDelayTimes[3]  = sqrt(relativeDelayTimes[ 4] * relativeDelayTimes[ 2] );
  relativeDelayTimes[5]  = sqrt(relativeDelayTimes[ 6] * relativeDelayTimes[ 4] );
  relativeDelayTimes[7]  = sqrt(relativeDelayTimes[ 8] * relativeDelayTimes[ 6] );
  relativeDelayTimes[9]  = sqrt(relativeDelayTimes[10] * relativeDelayTimes[ 8] );
  relativeDelayTimes[11] = sqrt(relativeDelayTimes[12] * relativeDelayTimes[10] );
  relativeDelayTimes[13] = sqrt(relativeDelayTimes[14] * relativeDelayTimes[12] );
  relativeDelayTimes[15] = sqrt(dMax                   * relativeDelayTimes[14] );
  adjustDelayTimes();
}

void rsReverb::setupOutputVector()
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

  double normalizer;
  for(int d=1; d<numDelayLines; d++)
  {
    normalizer = ((double) delaysInSamples[0] / (double) delaysInSamples[d]); 
    outputVectorL[d] *= normalizer;
    outputVectorR[d] *= normalizer;
  }
}

void rsReverb::adjustReadPointer(int index)
{
  tapOuts[index] = tapIn - delaysInSamples[index];
  if( tapOuts[index] < 0 )
    tapOuts[index] += maxDelayInSamples;
}

void rsReverb::adjustDelayTimes()
{
  int tmpDelays[numDelayLines];
  bool usePrimes = true;
  double factor  = 0.001*referenceDelayTime*sampleRate;
  int i;
  for(i=0; i<numDelayLines; i++)
  {
    tmpDelays[i] = roundToInt(factor*relativeDelayTimes[i]);
    if( usePrimes == true )
      tmpDelays[i] = PrimeNumbers::findClosestPrime(tmpDelays[i]);
    if( tmpDelays[i] > maxDelayInSamples )
      tmpDelays[i] = maxDelayInSamples;
    delaysInSamples[i] = tmpDelays[i];
  }

  // update dependent variables:
  for(i=0; i<numDelayLines; i++)
    adjustReadPointer(i);
  setupOutputVector();
  updateDampingAndCorrectionFilters();
  reset();
}

void rsReverb::updateDampingAndCorrectionFilters()
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

void rsReverb::reset()
{
  correctionFilterL.reset();
  correctionFilterR.reset();
  wetFilterL.reset();
  wetFilterR.reset();
  for(int d=0; d<numDelayLines; d++)
  {
    dampingFilters[d].reset();
    delayLineOuts[d] = 0.0;
    for(int n=0; n<maxDelayInSamples; n++)
      delayLines[d][n] = 0.0;
  }
}


/*

Ideas:
-Use allpass diffusors on left and right output. 
-Each diffusor is a chain of multi-sample delay allpass filters. That means, instead of z^-1, we 
 have z^-M delays.
-Each such allpass filter should have the same decay time, i.e. the feedback coeff should be scaled
 according to the number of delay samples M.
-The sum of all delays should be the same for the left and right chain. This may be achieved by 
 using twin primes like so:
   L: 5, 13, 17, 31, 41, 61 
   R: 7, 11, 19, 29, 43, 59
 The sum of bot rows is the same, namely 168. We alternatingly put the larger of the twin to the L 
 or R row. The number of primes in each row must be even (here it is 6). We may also prepend 
 1,2,3 to each row to make the impulse response even denser.


*/