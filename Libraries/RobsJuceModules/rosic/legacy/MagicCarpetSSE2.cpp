#include "MagicCarpetSSE2.h"
//#include "fft4g.c"

//----------------------------------------------------------------------------
// construction/destruction:

MagicCarpetSSE2::MagicCarpetSSE2()
{
 // init member variables:
 sampleRate           = 44100.0;    
 numVoices            = 16;
 //tableLength          = MAXTABLELENGTH;
 //fundamental    = 110.0; 
 //fundamentalRec = 1.0 / fundamental;
 x                    = 0.0;
 y                    = 0.0;
 randomSeed           = 17;
 masterAmplitude      = 1.0;
 midSideMix           = 0.5;
 mostRecentNote       = -1;
 mostRecentNoteVel    = 0;
 mostRecentNoteDetune = 0;

 // initialize the table contets to all zeros:
 zeroAllTables();

 x1L = 0.0;
 x1R = 0.0;

 outFilterIsOn = true;


 int i;

 // init all the tables with zeros:
 for(i=0; i<MAXTABLELENGTH+4; i++)
 {
  topLeftTable[i]     = 0.0;
  topRightTable[i]    = 0.0;
  bottomLeftTable[i]  = 0.0;
  bottomRightTable[i] = 0.0;
  xModTable[i]        = 0.0;
  yModTable[i]        = 0.0;
 }

 // pass pointers to the beginnings of the wavetables to all the voices:
 for(i=0; i<MAXVOICES; i++)
 {
  magicCarpetVoiceArray[i].topLeftSource.setTableAdress(topLeftTable);
  magicCarpetVoiceArray[i].topRightSource.setTableAdress(topRightTable);
  magicCarpetVoiceArray[i].bottomLeftSource.setTableAdress(bottomLeftTable);
  magicCarpetVoiceArray[i].bottomRightSource.setTableAdress(bottomRightTable);
 }
 xModulator.setTableAdress(xModTable);
 yModulator.setTableAdress(yModTable);


 fourierTransformer.setBlockSize(262144);


 // for performance testing:
 //for(i=0; i<1; i++)
  //magicCarpetVoiceArray[i].noteOn(64, 64);

 // test:
 /*
 for(i=0; i<tableSize; i++)
  table0[i] = sin( (2*PI*i) / 128 );
 */

}

MagicCarpetSSE2::~MagicCarpetSSE2()
{
 
}

//----------------------------------------------------------------------------
// parameter settings:

void MagicCarpetSSE2::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0 )
  sampleRate = newSampleRate;

 for(long i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].setSampleRate(sampleRate);

 // update the LFOs:


 // update the output-filters:
 outputFilter.setSampleRate(sampleRate);
 outputEqualizer.setSampleRate(sampleRate);
}

void MagicCarpetSSE2::setNumVoices(int newNumVoices)
{
 /*
 if( newNumVoices>0 && newNumVoices<=MAXVOICES )
  numVoices = newNumVoices;

 //reset all the voices:
 for(long i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].noteOn(0,0);
 */
}

void MagicCarpetSSE2::setWaveTable(int    destinationSlot, 
                               double *newWaveTable, 
                               int    newTableLength)
{
 // make sure, that the passed table is not too long:
 if( newTableLength > MAXTABLELENGTH || newTableLength < 2)
  return;

 int i;
 // remove the dc-component from the newTable:
 double dc = 0.0;
 for(i=0; i<newTableLength; i++)
  dc += newWaveTable[i];
 dc /= (double) newTableLength;
 for(i=0; i<newTableLength; i++)
  newWaveTable[i] -= dc;


 double accu = 0.0; // we want store the integrated waveform, so wee need an 
                    // accumulator for the integrator filter

 switch( destinationSlot )
 {
 case TOP_LEFT:
  {
   // copy the passed table into our member-array: 
   for(i=0; i<newTableLength; i++)
   {
    accu += newWaveTable[i];
    topLeftTable[i] = accu;
   }
   // and fill the rest with zeros, when the new waveform is shorter than our
   // array:
   if( newTableLength < MAXTABLELENGTH )
    for(i=newTableLength+1; i<MAXTABLELENGTH; i++)
     topLeftTable[i] = 0.0;

   topLeftTable[newTableLength]   = topLeftTable[0]; // for interpolation
   topLeftTable[newTableLength+1] = topLeftTable[1];
   topLeftTable[newTableLength+2] = topLeftTable[2];
   topLeftTable[newTableLength+3] = topLeftTable[3];
   topLeftTable[newTableLength+4] = topLeftTable[4];
   topLeftTable[newTableLength+5] = topLeftTable[5];

   // tell all the voices the new tableLength:
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].topLeftSource.setTableLength(newTableLength);
  }
  break;
  // and the same for the other slots:
 case TOP_RIGHT:
  {
   for(i=0; i<newTableLength; i++)
   {
    accu += newWaveTable[i];
    topRightTable[i] = accu;
   }
   if( newTableLength < MAXTABLELENGTH )
    for(i=newTableLength+1; i<MAXTABLELENGTH; i++)
     topRightTable[i] = 0.0;
   topRightTable[newTableLength]   = topRightTable[0]; // for interpolation
   topRightTable[newTableLength+1] = topRightTable[1];
   topRightTable[newTableLength+2] = topRightTable[2];
   topRightTable[newTableLength+3] = topRightTable[3];
   topRightTable[newTableLength+4] = topRightTable[4];
   topRightTable[newTableLength+5] = topRightTable[5];
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].topRightSource.setTableLength(newTableLength);
  }
  break;
 case BOTTOM_LEFT:
  {
   for(i=0; i<newTableLength; i++)
   {
    accu += newWaveTable[i];
    bottomLeftTable[i] = accu;
   }
   if( newTableLength < MAXTABLELENGTH )
    for(i=newTableLength+1; i<MAXTABLELENGTH; i++)
     bottomLeftTable[i] = 0.0;
   bottomLeftTable[newTableLength]   = bottomLeftTable[0]; // for interpolation
   bottomLeftTable[newTableLength+1] = bottomLeftTable[1];
   bottomLeftTable[newTableLength+2] = bottomLeftTable[2];
   bottomLeftTable[newTableLength+3] = bottomLeftTable[3];
   bottomLeftTable[newTableLength+4] = bottomLeftTable[4];
   bottomLeftTable[newTableLength+5] = bottomLeftTable[5];
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].bottomLeftSource.setTableLength(newTableLength);
  }
  break;
 case BOTTOM_RIGHT:
  {
   for(i=0; i<newTableLength; i++)
   {
    accu += newWaveTable[i];
    bottomRightTable[i] = accu;
   }
   if( newTableLength < MAXTABLELENGTH )
    for(i=newTableLength+1; i<MAXTABLELENGTH; i++)
     bottomRightTable[i] = 0.0;
   bottomRightTable[newTableLength]   = bottomRightTable[0]; // for interpolation
   bottomRightTable[newTableLength+1] = bottomRightTable[1];
   bottomRightTable[newTableLength+2] = bottomRightTable[2];
   bottomRightTable[newTableLength+3] = bottomRightTable[3];
   bottomRightTable[newTableLength+4] = bottomRightTable[4];
   bottomRightTable[newTableLength+5] = bottomRightTable[5];
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].bottomRightSource.setTableLength(newTableLength);
  }
  break;

 case X_MOD:
  {
   // copy the passed table into our member array:
   for(i=0; i<newTableLength; i++)
    xModTable[i] = newWaveTable[i];

   // fill the rest with zeros:
   if( newTableLength < MAXTABLELENGTH )
    for(i=newTableLength+1; i<MAXTABLELENGTH; i++)
     xModTable[i] = 0.0;

   // repeat the first few samples at the end for the interpolator:
   xModTable[newTableLength]   = xModTable[0]; // for interpolation
   xModTable[newTableLength+1] = xModTable[1];
   xModTable[newTableLength+2] = xModTable[2];
   xModTable[newTableLength+3] = xModTable[3];
   xModTable[newTableLength+4] = xModTable[4];
   xModTable[newTableLength+5] = xModTable[5];

   // tell the modulator the new table-length:
   xModulator.setTableLength(newTableLength);
  }
  break;

 case Y_MOD:
  {
   // copy the passed table into our member array:
   for(i=0; i<newTableLength; i++)
    yModTable[i] = newWaveTable[i];

   // fill the rest with zeros:
   if( newTableLength < MAXTABLELENGTH )
    for(i=newTableLength+1; i<MAXTABLELENGTH; i++)
     yModTable[i] = 0.0;

   // repeat the first few samples at the end for the interpolator:
   yModTable[newTableLength]   = yModTable[0]; // for interpolation
   yModTable[newTableLength+1] = yModTable[1];
   yModTable[newTableLength+2] = yModTable[2];
   yModTable[newTableLength+3] = yModTable[3];
   yModTable[newTableLength+4] = yModTable[4];
   yModTable[newTableLength+5] = yModTable[5];

   // tell the modulator the new table-length:
   yModulator.setTableLength(newTableLength);
  }
  break;

 } //  end of  switch( destinationSlot )
}

void MagicCarpetSSE2::setSpectrum(int    destinationSlot, 
                              double *newSpectrum, 
                              int    newSpectrumLength)
{
 // copy the spectrum into aour temporary array and fill with zeros:
 int i;
 int numBins = 131072;
 double* magnitudeSpectrum = new double[numBins];
 for(i=0; i<newSpectrumLength; i++)
  magnitudeSpectrum[i] = 0.5*(newSpectrum[i]+1.0);
 for(i=newSpectrumLength; i<numBins; i++)
  magnitudeSpectrum[i] = 0.0;

 // generate a pseudo-random phase-spectrum:
 srand(randomSeed);
 double* phaseSpectrum = new double[numBins];
 for(i=0; i<numBins; i++)
  phaseSpectrum[i] = 2.0*PI*((double)rand())/32767.0;

 // from the passed magnitude spectrum ("newSpectrum") and our pseudo-random 
 // phase-spectrum, we now render the time-domain signal via an iFFT:
 fourierTransformer.getSigFromMagAndPhs(magnitudeSpectrum, phaseSpectrum, tmpTable);

 // find the maximum for normalizing:
 double max = 0.0;
 for(i=0; i<2*numBins; i++)
 {
  if( fabs(tmpTable[i]) > max )
   max = fabs(tmpTable[i]);
 }

 // normalize:
 double normalizer = 1.0/max;
 for(i=0; i<2*numBins; i++)
  tmpTable[i] *= normalizer;

 // the temporary table will now be stored in its dedicated slot:
 setWaveTable(destinationSlot, tmpTable, 2*numBins);

 // for debug, copy the phase into the tmp-array::
 for(i=0; i<numBins; i++)
  tmpTable[i] = phaseSpectrum[i];

 // free dynamically allocated memory:
 delete[] magnitudeSpectrum;
 delete[] phaseSpectrum;
}

void MagicCarpetSSE2::setSlotLoopMode(int  whichSlot, 
                                  bool newLoopSwitch)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setLoopMode(newLoopSwitch);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setLoopMode(newLoopSwitch);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setLoopMode(newLoopSwitch);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setLoopMode(newLoopSwitch);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotSingleCycleMode(int  whichSlot, 
                                         bool newSingleCycleSwitch)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setSingleCycleMode(newSingleCycleSwitch);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setSingleCycleMode(newSingleCycleSwitch);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setSingleCycleMode(newSingleCycleSwitch);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setSingleCycleMode(newSingleCycleSwitch);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotStereoMode(int  whichSlot, 
                                    bool newStereoSwitch)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setStereoMode(newStereoSwitch);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setStereoMode(newStereoSwitch);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setStereoMode(newStereoSwitch);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setStereoMode(newStereoSwitch);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotMute(int  whichSlot, 
                              bool newMuteSwitch)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setMute(newMuteSwitch);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setMute(newMuteSwitch);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setMute(newMuteSwitch);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setMute(newMuteSwitch);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotDetune(int whichSlot, double newDetune)
{
 double newDetuneFactor = pitchOffsetToFreqFactor(newDetune);
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setDetuneFactor(newDetuneFactor);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setDetuneFactor(newDetuneFactor);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setDetuneFactor(newDetuneFactor);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setDetuneFactor(newDetuneFactor);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotTuneByKey(int whichSlot, double newTuneByKey)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftTuneByKey = newTuneByKey;
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightTuneByKey = newTuneByKey;
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftTuneByKey = newTuneByKey;
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightTuneByKey = newTuneByKey;
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotLevel(int whichSlot, double newLevel)
{
 double newAmplitude = dB2amp(newLevel);
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setAmplitude(newAmplitude);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setAmplitude(newAmplitude);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setAmplitude(newAmplitude);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setAmplitude(newAmplitude);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotLevelByKey(int whichSlot, double newLevelByKey)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftLevelByKey = newLevelByKey;
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightLevelByKey = newLevelByKey;
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftLevelByKey = newLevelByKey;
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightLevelByKey = newLevelByKey;
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotLevelByVel(int whichSlot, double newLevelByVel)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftLevelByVel = newLevelByVel;
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightLevelByVel = newLevelByVel;
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftLevelByVel = newLevelByVel;
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightLevelByVel = newLevelByVel;
   break;
  } // end switch
 } // end for
}



void MagicCarpetSSE2::setSlotHpf(int whichSlot, double newHpfCutoff)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setHpfCutoff(newHpfCutoff);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setHpfCutoff(newHpfCutoff);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setHpfCutoff(newHpfCutoff);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setHpfCutoff(newHpfCutoff);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotLpf(int whichSlot, double newLpfCutoff)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setLpfCutoff(newLpfCutoff);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setLpfCutoff(newLpfCutoff);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setLpfCutoff(newLpfCutoff);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setLpfCutoff(newLpfCutoff);
   break;
  } // end switch
 } // end for
}

void MagicCarpetSSE2::setSlotRootKey(int whichSlot, double newRootKey)
{
 double newFundamental = pitchToFreq(newRootKey);
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setFundamentalFreq(newFundamental);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setFundamentalFreq(newFundamental);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setFundamentalFreq(newFundamental);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setFundamentalFreq(newFundamental);
   break;
  } // end switch
 } // end for
}




void MagicCarpetSSE2::setXModNumerator(int newXModNumerator)
{
 xModulator.setPeriodNumerator(newXModNumerator);
}

void MagicCarpetSSE2::setXModDenominator(int newXModDenominator)
{
 xModulator.setPeriodDenominator(newXModDenominator);
}

void MagicCarpetSSE2::setXModStartPhase(double newXModStartPhase)
{
 xModulator.setStartPhase(newXModStartPhase);
}

void MagicCarpetSSE2::setXModAttack(double newXModAttack)
{
 xModulator.setAttack(newXModAttack);
}

void MagicCarpetSSE2::setXModRelease(double newXModRelease)
{
 xModulator.setRelease(newXModRelease);
}

void MagicCarpetSSE2::setXModRise(double newXModRise)
{
 xModulator.setRiseTime(newXModRise);
}

void MagicCarpetSSE2::setXModAmount(double newXModAmount)
{
 xModulator.setAmount(newXModAmount);
}

void MagicCarpetSSE2::setXModBpm(double newXModBpm)
{
 xModulator.setBpm(newXModBpm);
}

void MagicCarpetSSE2::setXModLoopMode(bool newXModLoopModeSwitch)
{
 xModulator.setLoopMode(newXModLoopModeSwitch);
}

void MagicCarpetSSE2::setYModNumerator(int newYModNumerator)
{
 yModulator.setPeriodNumerator(newYModNumerator);
}

void MagicCarpetSSE2::setYModDenominator(int newYModDenominator)
{
 yModulator.setPeriodDenominator(newYModDenominator);
}

void MagicCarpetSSE2::setYModStartPhase(double newYModStartPhase)
{
 yModulator.setStartPhase(newYModStartPhase);
}

void MagicCarpetSSE2::setYModAttack(double newYModAttack)
{
 yModulator.setAttack(newYModAttack);
}

void MagicCarpetSSE2::setYModRelease(double newYModRelease)
{
 yModulator.setRelease(newYModRelease);
}

void MagicCarpetSSE2::setYModRise(double newYModRise)
{
 yModulator.setRiseTime(newYModRise);
}

void MagicCarpetSSE2::setYModAmount(double newYModAmount)
{
 yModulator.setAmount(newYModAmount);
}

void MagicCarpetSSE2::setYModBpm(double newYModBpm)
{
 yModulator.setBpm(newYModBpm);
}

void MagicCarpetSSE2::setYModLoopMode(bool newYModLoopModeSwitch)
{
 yModulator.setLoopMode(newYModLoopModeSwitch);
}

void MagicCarpetSSE2::setFltMode(int newFltMode)
{ 
 for(int i=0; i<MAXVOICES; i++)
 {
  magicCarpetVoiceArray[i].filter.setMode(newFltMode);
  if( newFltMode == MagicCarpetFilterSSE2::BYPASS )
   magicCarpetVoiceArray[i].filterIsOn = false;
  else
   magicCarpetVoiceArray[i].filterIsOn = true;
 }
}

void MagicCarpetSSE2::setFltFreq(double newFltFreq)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltFreq = newFltFreq;
}

void MagicCarpetSSE2::setFltFreqByKey(double newFltFreqByKey)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltFreqByKey = newFltFreqByKey;
}

void MagicCarpetSSE2::setFltFreqByVel(double newFltFreqByVel)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltFreqByVel = newFltFreqByVel;
}

void MagicCarpetSSE2::setFltReso(double newFltReso)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].filter.setQ(newFltReso);
}

void MagicCarpetSSE2::setFltGain(double newFltGain)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].filter.setGain(newFltGain);
  // gain is scaled by 0.5 because the filter uses 2 stages, which should be
  // compensated......not anymore
}

void MagicCarpetSSE2::setFltEnvAttack(double newFltEnvAttack)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setAttack(0.001*newFltEnvAttack);
}

void MagicCarpetSSE2::setFltEnvPeak(double newFltEnvPeak)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setPeak(newFltEnvPeak);
}

void MagicCarpetSSE2::setFltEnvPeakByVel(double newFltEnvPeakByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnvPeakByVel = newFltEnvPeakByVel;
}

void MagicCarpetSSE2::setFltEnvDecay(double newFltEnvDecay)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setDecay(0.001*newFltEnvDecay);
}



void MagicCarpetSSE2::setFltEnvRelease(double newFltEnvRelease)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setRelease(0.001*newFltEnvRelease);
}
void MagicCarpetSSE2::setFltEnvEnd(double newFltEnvEnd)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setEnd(newFltEnvEnd);
}

void MagicCarpetSSE2::setFltEnvDurByKey(double newFltEnvDurByKey)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnvDurByKey = newFltEnvDurByKey;
}

void MagicCarpetSSE2::setFltEnvDurByVel(double newFltEnvDurByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnvDurByVel = newFltEnvDurByVel;
}

void MagicCarpetSSE2::setFltEnvSlope(double newFltEnvSlope)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setTauScale(newFltEnvSlope);
}


void MagicCarpetSSE2::setMasterVolume(double newMasterVolume)
{ 
 masterAmplitude = dB2amp(newMasterVolume);
}

void MagicCarpetSSE2::setMidSideMix(double newMidSideMix)
{ 
 midSideMix = newMidSideMix;
}

void MagicCarpetSSE2::setX(double newX)
{
 x = newX;
}

void MagicCarpetSSE2::setY(double newY)
{
 y = newY;
}

void MagicCarpetSSE2::setAmpEnvAttack(double newAmpEnvAttack)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setAttack(0.001*newAmpEnvAttack);
}

void MagicCarpetSSE2::setAmpEnvPeak(double newAmpEnvPeak)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setPeak(newAmpEnvPeak);
}

void MagicCarpetSSE2::setAmpEnvPeakByVel(double newAmpEnvPeakByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnvPeakByVel = newAmpEnvPeakByVel;
}

void MagicCarpetSSE2::setAmpEnvDecay(double newAmpEnvDecay)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setDecay(0.001*newAmpEnvDecay);
}

void MagicCarpetSSE2::setAmpEnvSustain(double newAmpEnvSustain)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setSustain(newAmpEnvSustain);
}

void MagicCarpetSSE2::setAmpEnvRelease(double newAmpEnvRelease)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setRelease(0.001*newAmpEnvRelease);
}

void MagicCarpetSSE2::setAmpEnvDurByKey(double newAmpEnvDurByKey)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnvDurByKey = newAmpEnvDurByKey;
}

void MagicCarpetSSE2::setAmpEnvDurByVel(double newAmpEnvDurByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnvDurByVel = newAmpEnvDurByVel;
}

void MagicCarpetSSE2::setAmpEnvSlope(double newAmpEnvSlope)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setTauScale(newAmpEnvSlope);
}

void MagicCarpetSSE2::setOutHpfFreq(double newOutHpfFreq)
{
 outputFilter.setHpfCutoff(newOutHpfFreq);
}
void MagicCarpetSSE2::setOutLpfFreq(double newOutLpfFreq)
{
 outputFilter.setLpfCutoff(newOutLpfFreq);
}
void MagicCarpetSSE2::setOutEqFreq(double newOutEqFreq, int whichBand)
{
 outputEqualizer.setEqFreq(newOutEqFreq, whichBand);
}
void MagicCarpetSSE2::setOutEqQ(double newOutEqQ, int whichBand)
{
 outputEqualizer.setEqQ(newOutEqQ, whichBand);
}
void MagicCarpetSSE2::setOutEqGain(double newOutEqGain, int whichBand)
{
 outputEqualizer.setEqGain(newOutEqGain, whichBand);
}


/*
void MagicCarpetSSE2::setStereoWidth(flt64 newStereoWidth)
{ 
 
}

// the amplitude-envelope settings:
void MagicCarpetSSE2::setAmpEnv1Attack(double newEnv1Attack)
{
 for(int i=0; i<numVoices; i++)
  magicCarpetVoiceArray[i].setAmpEnvAttack(newAmpEnvAttack);
}

//...

*/

//----------------------------------------------------------------------------
// event processing:

void MagicCarpetSSE2::noteOn(long NoteNumber, long Velocity)
{
 mostRecentNote       = NoteNumber;
 mostRecentNoteVel    = Velocity;
 //mostRecentNoteDetune = Detune;

 //magicCarpetVoiceArray[NoteNumber].noteOn(mostRecentNote, mostRecentNoteVel);


 long i = 0;
 
 if(mostRecentNoteVel==0)   // a note-off event occured
 {
  // loop through the voices to find the one which is playing the note
  // for which the note-off was was sent:
  for(i=0; i<numVoices; i++)
  {
   if(magicCarpetVoiceArray[i].currentNote == mostRecentNote)
    magicCarpetVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
     // a note-on with velocity 0 is sent to the voice
  }
 }

 else                         // a note-on event occured
 {
  // look at all the voices, to see if ALL of them have reached their ends -
  // only in this case, the modulator-phases are reset:
  bool aVoiceIsPlaying  = false;
  for(i=0; i<numVoices; i++)
  {
   if( !magicCarpetVoiceArray[i].ampEnv.endIsReached() )
    aVoiceIsPlaying = true;
  }
  if( aVoiceIsPlaying == false )
  {
   xModulator.reset();
   yModulator.reset();
  }

  // loop through the voices to find a free one:
  int  oldestNoteAge    = 0;
  int  oldestVoice      = 0;       // voice with the oldest note
  for(i=0; i<numVoices; i++)
  {
   // check, if voice i is currently releasing the note, which is coming in,
   // if so, use that voice again:
   // check if the voice i is free:
   if(magicCarpetVoiceArray[i].currentNote == mostRecentNote)
   {
    // voice is used for the same note (which is currently releasing) again:
    magicCarpetVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
    return; //jump out of the function
   }
   // check if the voice i is free:
   else if(magicCarpetVoiceArray[i].ampEnv.endIsReached())
   {
    // voice i is free and can be used for the new note:
    magicCarpetVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
    return; //jump out of the function
   }
   // keep track of the oldest note for the case, when no voice is
   // free for the new note (last note priority assignment):
   else if(magicCarpetVoiceArray[i].currentNoteAge > oldestNoteAge)
   {
    oldestNoteAge = magicCarpetVoiceArray[i].currentNoteAge;
    oldestVoice   = i;
   }
  } // end of for-loop


  // no free voice has been found - set the voice with the oldest note
  // to the new note:
  magicCarpetVoiceArray[oldestVoice].noteOn(mostRecentNote, mostRecentNoteVel);
 } // end of else

}


//----------------------------------------------------------------------------
// internal functions:

void MagicCarpetSSE2::zeroAllTables()
{
 int i;
 for(i=0; i<MAXTABLELENGTH+1; i++)
 {
  topLeftTable[i] = 0.0;
  topRightTable[i] = 0.0;
  bottomLeftTable[i] = 0.0;
  bottomRightTable[i] = 0.0;
 }
}