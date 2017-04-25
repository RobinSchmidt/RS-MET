#include "MagicCarpet.h"
//#include "fft4g.c"

//----------------------------------------------------------------------------
// construction/destruction:

MagicCarpet::MagicCarpet()
{
 // init member variables:
 sampleRate           = 44100.0;    
 numVoices            = 16;
 //tableLength          = MAXTABLELENGTH;
 //fundamental    = 110.0; 
 //fundamentalRec = 1.0 / fundamental;
 x                       = 0.0;
 y                       = 0.0;
 randomSeed              = 17;
 masterAmplitude         = 1.0;
 masterAmplitudeByVoices = 0.0;
 midSideMix              = 0.5;
 mostRecentNote          = -1;
 mostRecentNoteVel       = 0;
 mostRecentNoteDetune    = 0;

 // initialize the table contets to all zeros:
 //zeroAllTables();

 // init all the table-pointers with NULL-pointers:
 topLeftTableL     = NULL;
 topLeftTableR     = NULL;
 topRightTableL    = NULL;
 topRightTableR    = NULL;
 bottomLeftTableL  = NULL;
 bottomLeftTableR  = NULL;
 bottomRightTableL = NULL;
 bottomRightTableR = NULL;
 xModTable         = NULL;
 yModTable         = NULL;

 outFilterIsOn   = true;
 numActiveVoices = 0;

 xSmoother.setMode(OnePoleFilter::LOWPASS);
 xSmoother.setSampleRate(sampleRate);
 xSmoother.setCutoff(50.0);
 xSmoother.resetBuffers();

 ySmoother.setMode(OnePoleFilter::LOWPASS);
 ySmoother.setSampleRate(sampleRate);
 ySmoother.setCutoff(50.0);
 ySmoother.resetBuffers();

 fourierTransformer.setBlockSize(262144);
}

MagicCarpet::~MagicCarpet()
{
 
}

//----------------------------------------------------------------------------
// parameter settings:

void MagicCarpet::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0 )
  sampleRate = newSampleRate;

 for(long i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].setSampleRate(sampleRate);

 // update the LFOs:

 // update the output-filters:
 outputFilter.setSampleRate(sampleRate);
 outputEqualizer.setSampleRate(sampleRate);

 // update the smoothers:
 xSmoother.setSampleRate(sampleRate);
 ySmoother.setSampleRate(sampleRate);
}

void MagicCarpet::setNumVoices(int newNumVoices)
{
 if( newNumVoices>0 && newNumVoices<=MAXVOICES )
  numVoices = newNumVoices;

 // send a note-off to all the voices:
 for(long i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].noteOn(0,0);
}

int MagicCarpet::getNumActiveVoices()
{
 return numActiveVoices;
}

void MagicCarpet::invalidateTablePointers(int whichSlot)
{
 for(int i=0; i<MAXVOICES; i++)
 {
  switch( whichSlot )
  {
  case(TOP_LEFT):
   magicCarpetVoiceArray[i].topLeftSource.setTableAddresses(NULL, NULL);
   break;
  case(TOP_RIGHT):
   magicCarpetVoiceArray[i].topRightSource.setTableAddresses(NULL, NULL);
   break;
  case(BOTTOM_LEFT):
   magicCarpetVoiceArray[i].bottomLeftSource.setTableAddresses(NULL, NULL);
   break;
  case(BOTTOM_RIGHT):
   magicCarpetVoiceArray[i].bottomRightSource.setTableAddresses(NULL, NULL);
   break;
  case(X_MOD):
   xModulator.setTableAdress(NULL);
   break;
  case(Y_MOD):
   yModulator.setTableAdress(NULL);
   break;
  } // end switch
 } // end for
}

bool MagicCarpet::setWaveTables(int    destinationSlot, 
                                double *newWaveTableL, 
                                double *newWaveTableR,
                                int    newTableLength)
{
 int i;

 // remove the dc-component from the newTable (but not for the modulators):
 if( destinationSlot < X_MOD )
 {
  double dcL = 0.0;
  double dcR = 0.0;
  for(i=0; i<newTableLength; i++)
  {
   dcL += newWaveTableL[i];
   dcR += newWaveTableR[i];
  }
  dcL /= (double) newTableLength;
  dcR /= (double) newTableLength;
  for(i=0; i<newTableLength; i++)
  {
   newWaveTableL[i] -= dcL;
   newWaveTableR[i] -= dcR;
  } 
 }

 double accuL = 0.0; // we want store the integrated waveform, so wee need an 
 double accuR = 0.0; // accumulator for the integrator filter

 switch( destinationSlot )
 {
 case TOP_LEFT:
  {
   // tell the oscillators, that their table-addresses are invalid now:
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].topLeftSource.setTableAddresses(NULL, NULL);

   // free the old memory chunk and invalidate the old-pointers:
   if( topLeftTableL != NULL )
   {
    delete[] topLeftTableL;
    topLeftTableL = NULL;
   }
   if( topLeftTableR != NULL )
   {
    delete[] topLeftTableR;
    topLeftTableR = NULL;
   }

   // try to allocate new memory for the new table:
   topLeftTableL = new double[newTableLength+6];
   topLeftTableR = new double[newTableLength+6];

   // report failure if the memory-allocation failed; in the case, that only 
   // one table was allocated, delete the other one also:
   if( topLeftTableL == NULL || topLeftTableR == NULL )
   {
    if( topLeftTableL == NULL && topLeftTableR != NULL )
    {
     delete[] topLeftTableR;
     topLeftTableR = NULL;
    }
    if( topLeftTableR == NULL && topLeftTableL != NULL )
    {
     delete[] topLeftTableL;
     topLeftTableL = NULL;
    }
    return false;
   }

   // write the integrated table into our member-array: 
   for(i=0; i<newTableLength; i++)
   {
    accuL += newWaveTableL[i];
    accuR += newWaveTableR[i];
    topLeftTableL[i] = accuL;
    topLeftTableR[i] = accuR;
   }

   // repeat the first 6 samples at the end for the interpolator:
   for(i=newTableLength; i<=newTableLength+5; i++)
   {
    topLeftTableL[i] = topLeftTableL[i-newTableLength]; 
    topLeftTableR[i] = topLeftTableR[i-newTableLength]; 
   }

   // tell all the voices the new table-adresses and tableLength:
   for(i=0; i<MAXVOICES; i++)
   {
    magicCarpetVoiceArray[i].topLeftSource.setTableLength(newTableLength);
    magicCarpetVoiceArray[i].topLeftSource.setTableAddresses(topLeftTableL, 
                                                             topLeftTableR);
   }

  }
  break;

  // and the same for the other slots:
 case TOP_RIGHT:
  {
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].topRightSource.setTableAddresses(NULL, NULL);
   if( topRightTableL != NULL )
   {
    delete[] topRightTableL;
    topRightTableL = NULL;
   }
   if( topRightTableR != NULL )
   {
    delete[] topRightTableR;
    topRightTableR = NULL;
   }
   topRightTableL = new double[newTableLength+6];
   topRightTableR = new double[newTableLength+6];
   if( topRightTableL == NULL || topRightTableR == NULL )
   {
    if( topRightTableL == NULL && topRightTableR != NULL )
    {
     delete[] topRightTableR;
     topRightTableR = NULL;
    }
    if( topRightTableR == NULL && topRightTableL != NULL )
    {
     delete[] topRightTableL;
     topRightTableL = NULL;
    }
    return false;
   }
   for(i=0; i<newTableLength; i++)
   {
    accuL += newWaveTableL[i];
    accuR += newWaveTableR[i];
    topRightTableL[i] = accuL;
    topRightTableR[i] = accuR;
   }
   for(i=newTableLength; i<=newTableLength+5; i++)
   {
    topRightTableL[i] = topRightTableL[i-newTableLength]; 
    topRightTableR[i] = topRightTableR[i-newTableLength]; 
   }
   for(i=0; i<MAXVOICES; i++)
   {
    magicCarpetVoiceArray[i].topRightSource.setTableLength(newTableLength);
    magicCarpetVoiceArray[i].topRightSource.setTableAddresses(topRightTableL, 
                                                             topRightTableR);
   }
  }
  break;

 case BOTTOM_LEFT:
  {
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].bottomLeftSource.setTableAddresses(NULL, NULL);
   if( bottomLeftTableL != NULL )
   {
    delete[] bottomLeftTableL;
    bottomLeftTableL = NULL;
   }
   if( bottomLeftTableR != NULL )
   {
    delete[] bottomLeftTableR;
    bottomLeftTableR = NULL;
   }
   bottomLeftTableL = new double[newTableLength+6];
   bottomLeftTableR = new double[newTableLength+6];
   if( bottomLeftTableL == NULL || bottomLeftTableR == NULL )
   {
    if( bottomLeftTableL == NULL && bottomLeftTableR != NULL )
    {
     delete[] bottomLeftTableR;
     bottomLeftTableR = NULL;
    }
    if( bottomLeftTableR == NULL && bottomLeftTableL != NULL )
    {
     delete[] bottomLeftTableL;
     bottomLeftTableL = NULL;
    }
    return false;
   }
   for(i=0; i<newTableLength; i++)
   {
    accuL += newWaveTableL[i];
    accuR += newWaveTableR[i];
    bottomLeftTableL[i] = accuL;
    bottomLeftTableR[i] = accuR;
   }
   for(i=newTableLength; i<=newTableLength+5; i++)
   {
    bottomLeftTableL[i] = bottomLeftTableL[i-newTableLength]; 
    bottomLeftTableR[i] = bottomLeftTableR[i-newTableLength]; 
   }
   for(i=0; i<MAXVOICES; i++)
   {
    magicCarpetVoiceArray[i].bottomLeftSource.setTableLength(newTableLength);
    magicCarpetVoiceArray[i].bottomLeftSource.setTableAddresses(bottomLeftTableL, 
                                                             bottomLeftTableR);
   }
  }
  break;

 case BOTTOM_RIGHT:
  {
   for(i=0; i<MAXVOICES; i++)
    magicCarpetVoiceArray[i].bottomRightSource.setTableAddresses(NULL, NULL);
   if( bottomRightTableL != NULL )
   {
    delete[] bottomRightTableL;
    bottomRightTableL = NULL;
   }
   if( bottomRightTableR != NULL )
   {
    delete[] bottomRightTableR;
    bottomRightTableR = NULL;
   }
   bottomRightTableL = new double[newTableLength+6];
   bottomRightTableR = new double[newTableLength+6];
   if( bottomRightTableL == NULL || bottomRightTableR == NULL )
   {
    if( bottomRightTableL == NULL && bottomRightTableR != NULL )
    {
     delete[] bottomRightTableR;
     bottomRightTableR = NULL;
    }
    if( bottomRightTableR == NULL && bottomRightTableL != NULL )
    {
     delete[] bottomRightTableL;
     bottomRightTableL = NULL;
    }
    return false;
   }
   for(i=0; i<newTableLength; i++)
   {
    accuL += newWaveTableL[i];
    accuR += newWaveTableR[i];
    bottomRightTableL[i] = accuL;
    bottomRightTableR[i] = accuR;
   }
   for(i=newTableLength; i<=newTableLength+5; i++)
   {
    bottomRightTableL[i] = bottomRightTableL[i-newTableLength]; 
    bottomRightTableR[i] = bottomRightTableR[i-newTableLength]; 
   }
   for(i=0; i<MAXVOICES; i++)
   {
    magicCarpetVoiceArray[i].bottomRightSource.setTableLength(newTableLength);
    magicCarpetVoiceArray[i].bottomRightSource.setTableAddresses(bottomRightTableL, 
                                                             bottomRightTableR);
   }
  }
  break;

 case X_MOD:
  {
   // tell the modulator, that its table-addresses is invalid now:
   xModulator.setTableAdress(NULL);

   // free the old memory chunk and invalidate the old-pointers:
   if( xModTable != NULL )
   {
    delete[] xModTable;
    xModTable = NULL;
   }

   // try to allocate new memory for the new table:
   xModTable = new double[newTableLength+6];

   // report failure if the memory-allocation failed:
   if( xModTable == NULL )
    return false;

   // write the table into our member-array:
   for(i=0; i<newTableLength; i++)
    xModTable[i] = newWaveTableL[i];

   // repeat the first 6 samples at the end for the interpolator:
   for(i=newTableLength; i<=newTableLength+5; i++)
    xModTable[i] = xModTable[i-newTableLength]; 

   // tell all the modulator the new table-address and tableLength:
   xModulator.setTableLength(newTableLength);
   xModulator.setTableAdress(xModTable);
  }
  break;

 // the same for the y-modulator:
 case Y_MOD:
  {
   yModulator.setTableAdress(NULL);
   if( yModTable != NULL )
   {
    delete[] yModTable;
    yModTable = NULL;
   }
   yModTable = new double[newTableLength+6];
   if( yModTable == NULL )
    return false;
   for(i=0; i<newTableLength; i++)
    yModTable[i] = newWaveTableL[i];
   for(i=newTableLength; i<=newTableLength+5; i++)
    yModTable[i] = yModTable[i-newTableLength]; 
   yModulator.setTableLength(newTableLength);
   yModulator.setTableAdress(yModTable);
  }
  break;

 } //  end of  switch( destinationSlot )

 return true;
}

bool MagicCarpet::setSpectrum(int    destinationSlot, 
                              double *newSpectrum, 
                              int    newSpectrumLength)
{
 // copy the spectrum into aour temporary array and fill with zeros:
 int i;
 int numBins        = 131072; // numBins maybe chosen equal to or larger than
                              // newSpectrumLength - larger numBins results in 
                              // zero-padding in the freq-domain - however, it 
                              // is generally not advisable to use such a 
                              // hardcoded constant - get rid of this...
 int newTableLength = 2*numBins;

 // allocate some temporary buffers to work with:
 double* tmpTableL         = new double[newTableLength];
 double* tmpTableR         = new double[newTableLength];
 double* magnitudeSpectrum = new double[numBins];
 double* phaseSpectrum     = new double[numBins];

 // check if all temporary buffers were successfully allocated, otherwise free
 // the occupied memory (there may be the case that some buffers could be 
 // allocated and others not) and report failure:
 if( tmpTableL == NULL || 
     tmpTableR == NULL || 
     magnitudeSpectrum == NULL || 
     phaseSpectrum == NULL )
 {
  if( tmpTableL != NULL )
   delete[] tmpTableL;
  if( tmpTableR != NULL )
   delete[] tmpTableR;
  if( magnitudeSpectrum != NULL )
   delete[] magnitudeSpectrum;
  if( phaseSpectrum != NULL )
   delete[] phaseSpectrum;

  return false;
 }

 for(i=0; i<newSpectrumLength; i++)
  magnitudeSpectrum[i] = 0.5*(newSpectrum[i]+1.0);
 for(i=newSpectrumLength; i<numBins; i++)
  magnitudeSpectrum[i] = 0.0;

 // generate a pseudo-random phase-spectrum:
 srand(randomSeed);

 for(i=0; i<numBins; i++)
  phaseSpectrum[i] = 2.0*PI*((double)rand())/32767.0;

 // from the passed magnitude spectrum ("newSpectrum") and our pseudo-random 
 // phase-spectrum, we now render the time-domain signal via an iFFT:
 fourierTransformer.getSigFromMagAndPhs(magnitudeSpectrum, phaseSpectrum, 
                                        tmpTableL);

 // find the maximum for normalizing:
 double max = 0.0;
 for(i=0; i<2*numBins; i++)
 {
  if( fabs(tmpTableL[i]) > max )
   max = fabs(tmpTableL[i]);
 }

 // normalize:
 double normalizer = 1.0/max;
 for(i=0; i<2*numBins; i++)
 {
  tmpTableL[i] *= normalizer;
 }

 // obtain the table for the right channel by means of a circular shift of
 // the left table by half the table-length:
 for(i=0; i<newTableLength; i++)
  tmpTableR[i] = tmpTableL[(i+newTableLength/2)%newTableLength];

 // the temporary table will now be stored in its dedicated slot:
 bool success = setWaveTables(destinationSlot, 
                              tmpTableL, 
                              tmpTableR, 
                              newTableLength);

 // free dynamically allocated memory:
 delete[] tmpTableL;
 delete[] tmpTableR;
 delete[] magnitudeSpectrum;
 delete[] phaseSpectrum;

 return success;
}

void MagicCarpet::setSlotLoopMode(int  whichSlot, 
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

void MagicCarpet::setSlotSingleCycleMode(int  whichSlot, 
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

void MagicCarpet::setSlotStereoMode(int  whichSlot, 
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

void MagicCarpet::setSlotMute(int  whichSlot, 
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

void MagicCarpet::setSlotDetune(int whichSlot, double newDetune)
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

void MagicCarpet::setSlotTuneByKey(int whichSlot, double newTuneByKey)
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

void MagicCarpet::setSlotLevel(int whichSlot, double newLevel)
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

void MagicCarpet::setSlotLevelByKey(int whichSlot, double newLevelByKey)
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

void MagicCarpet::setSlotLevelByVel(int whichSlot, double newLevelByVel)
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



void MagicCarpet::setSlotHpf(int whichSlot, double newHpfCutoff)
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

void MagicCarpet::setSlotLpf(int whichSlot, double newLpfCutoff)
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

void MagicCarpet::setSlotRootKey(int whichSlot, double newRootKey)
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




void MagicCarpet::setXModNumerator(int newXModNumerator)
{
 xModulator.setPeriodNumerator(newXModNumerator);
}

void MagicCarpet::setXModDenominator(int newXModDenominator)
{
 xModulator.setPeriodDenominator(newXModDenominator);
}

void MagicCarpet::setXModStartPhase(double newXModStartPhase)
{
 xModulator.setStartPhase(newXModStartPhase);
}

void MagicCarpet::setXModAttack(double newXModAttack)
{
 xModulator.setAttack(newXModAttack);
}

void MagicCarpet::setXModRelease(double newXModRelease)
{
 xModulator.setRelease(newXModRelease);
}

void MagicCarpet::setXModRise(double newXModRise)
{
 xModulator.setRiseTime(newXModRise);
}

void MagicCarpet::setXModAmount(double newXModAmount)
{
 xModulator.setAmount(newXModAmount);
}

void MagicCarpet::setXModBpm(double newXModBpm)
{
 xModulator.setBpm(newXModBpm);
}

void MagicCarpet::setXModLoopMode(bool newXModLoopModeSwitch)
{
 xModulator.setLoopMode(newXModLoopModeSwitch);
}

void MagicCarpet::setYModNumerator(int newYModNumerator)
{
 yModulator.setPeriodNumerator(newYModNumerator);
}

void MagicCarpet::setYModDenominator(int newYModDenominator)
{
 yModulator.setPeriodDenominator(newYModDenominator);
}

void MagicCarpet::setYModStartPhase(double newYModStartPhase)
{
 yModulator.setStartPhase(newYModStartPhase);
}

void MagicCarpet::setYModAttack(double newYModAttack)
{
 yModulator.setAttack(newYModAttack);
}

void MagicCarpet::setYModRelease(double newYModRelease)
{
 yModulator.setRelease(newYModRelease);
}

void MagicCarpet::setYModRise(double newYModRise)
{
 yModulator.setRiseTime(newYModRise);
}

void MagicCarpet::setYModAmount(double newYModAmount)
{
 yModulator.setAmount(newYModAmount);
}

void MagicCarpet::setYModBpm(double newYModBpm)
{
 yModulator.setBpm(newYModBpm);
}

void MagicCarpet::setYModLoopMode(bool newYModLoopModeSwitch)
{
 yModulator.setLoopMode(newYModLoopModeSwitch);
}

void MagicCarpet::setFltMode(int newFltMode)
{ 
 for(int i=0; i<MAXVOICES; i++)
 {
  magicCarpetVoiceArray[i].filter.setMode(newFltMode);
  if( newFltMode == MagicCarpetFilter::BYPASS )
   magicCarpetVoiceArray[i].filterIsOn = false;
  else
   magicCarpetVoiceArray[i].filterIsOn = true;
 }
}

void MagicCarpet::setFltTwoStages(bool newTwoStagesSwitch)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].filter.setTwoStages(newTwoStagesSwitch);
}

void MagicCarpet::setFltFreq(double newFltFreq)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltFreq = newFltFreq;
}

void MagicCarpet::setFltFreqByKey(double newFltFreqByKey)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltFreqByKey = newFltFreqByKey;
}

void MagicCarpet::setFltFreqByVel(double newFltFreqByVel)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltFreqByVel = newFltFreqByVel;
}

void MagicCarpet::setFltReso(double newFltReso)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].filter.setQ(newFltReso);
}

void MagicCarpet::setFltGain(double newFltGain)
{ 
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].filter.setGain(newFltGain);
  // gain is scaled by 0.5 because the filter uses 2 stages, which should be
  // compensated......not anymore
}

void MagicCarpet::setFltEnvAttack(double newFltEnvAttack)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setAttack(0.001*newFltEnvAttack);
}

void MagicCarpet::setFltEnvPeak(double newFltEnvPeak)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setPeak(newFltEnvPeak);
}

void MagicCarpet::setFltEnvPeakByVel(double newFltEnvPeakByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnvPeakByVel = newFltEnvPeakByVel;
}

void MagicCarpet::setFltEnvDecay(double newFltEnvDecay)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setDecay(0.001*newFltEnvDecay);
}



void MagicCarpet::setFltEnvRelease(double newFltEnvRelease)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setRelease(0.001*newFltEnvRelease);
}
void MagicCarpet::setFltEnvEnd(double newFltEnvEnd)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setEnd(newFltEnvEnd);
}

void MagicCarpet::setFltEnvDurByKey(double newFltEnvDurByKey)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnvDurByKey = newFltEnvDurByKey;
}

void MagicCarpet::setFltEnvDurByVel(double newFltEnvDurByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnvDurByVel = newFltEnvDurByVel;
}

void MagicCarpet::setFltEnvSlope(double newFltEnvSlope)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].fltEnv.setTauScale(newFltEnvSlope);
}

void MagicCarpet::setMasterVolume(double newMasterVolume)
{ 
 masterAmplitude = dB2amp(newMasterVolume);
}

void MagicCarpet::setMasterVolumeByVoices(double newMasterVolumeByVoices)
{ 
 masterAmplitudeByVoices = newMasterVolumeByVoices;
}

void MagicCarpet::setMidSideMix(double newMidSideMix)
{ 
 midSideMix = newMidSideMix;
}

void MagicCarpet::setX(double newX)
{
 x = newX;
}

void MagicCarpet::setY(double newY)
{
 y = newY;
}

void MagicCarpet::setAmpEnvAttack(double newAmpEnvAttack)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setAttack(0.001*newAmpEnvAttack);
}

void MagicCarpet::setAmpEnvPeak(double newAmpEnvPeak)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setPeak(newAmpEnvPeak);
}

void MagicCarpet::setAmpEnvPeakByVel(double newAmpEnvPeakByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnvPeakByVel = newAmpEnvPeakByVel;
}

void MagicCarpet::setAmpEnvDecay(double newAmpEnvDecay)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setDecay(0.001*newAmpEnvDecay);
}

void MagicCarpet::setAmpEnvSustain(double newAmpEnvSustain)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setSustain(newAmpEnvSustain);
}

void MagicCarpet::setAmpEnvRelease(double newAmpEnvRelease)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setRelease(0.001*newAmpEnvRelease);
}

void MagicCarpet::setAmpEnvDurByKey(double newAmpEnvDurByKey)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnvDurByKey = newAmpEnvDurByKey;
}

void MagicCarpet::setAmpEnvDurByVel(double newAmpEnvDurByVel)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnvDurByVel = newAmpEnvDurByVel;
}

void MagicCarpet::setAmpEnvSlope(double newAmpEnvSlope)
{
 for(int i=0; i<MAXVOICES; i++)
  magicCarpetVoiceArray[i].ampEnv.setTauScale(newAmpEnvSlope);
}

void MagicCarpet::setOutHpfFreq(double newOutHpfFreq)
{
 outputFilter.setHpfCutoff(newOutHpfFreq);
}
void MagicCarpet::setOutLpfFreq(double newOutLpfFreq)
{
 outputFilter.setLpfCutoff(newOutLpfFreq);
}
void MagicCarpet::setOutEqFreq(double newOutEqFreq, int whichBand)
{
 outputEqualizer.setEqFreq(newOutEqFreq, whichBand);
}
void MagicCarpet::setOutEqQ(double newOutEqQ, int whichBand)
{
 outputEqualizer.setEqQ(newOutEqQ, whichBand);
}
void MagicCarpet::setOutEqGain(double newOutEqGain, int whichBand)
{
 outputEqualizer.setEqGain(newOutEqGain, whichBand);
}


/*
void MagicCarpet::setStereoWidth(flt64 newStereoWidth)
{ 
 
}

// the amplitude-envelope settings:
void MagicCarpet::setAmpEnv1Attack(double newEnv1Attack)
{
 for(int i=0; i<numVoices; i++)
  magicCarpetVoiceArray[i].setAmpEnvAttack(newAmpEnvAttack);
}

//...

*/

//----------------------------------------------------------------------------
// event processing:

void MagicCarpet::noteOn(long NoteNumber, long Velocity)
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

void MagicCarpet::setPitchBend(double transpositionInSemitones)
{
 double transpositionFactor 
        = pitchOffsetToFreqFactor(transpositionInSemitones);

 for(int i=0; i<MAXVOICES; i++)
 {
  magicCarpetVoiceArray[i].topLeftSource.
   setTranspositionFactor(transpositionFactor);
  magicCarpetVoiceArray[i].topRightSource.
   setTranspositionFactor(transpositionFactor);
  magicCarpetVoiceArray[i].bottomLeftSource.
   setTranspositionFactor(transpositionFactor);
  magicCarpetVoiceArray[i].bottomRightSource.
   setTranspositionFactor(transpositionFactor);
 }
}

//----------------------------------------------------------------------------
// internal functions:
/*
void MagicCarpet::zeroAllTables()
{
 int i;
 for(i=0; i<MAXTABLELENGTH+1; i++)
 {
  topLeftTableL[i] = 0.0;
  topLeftTableR[i] = 0.0;
  topRightTableL[i] = 0.0;
  topRightTableR[i] = 0.0;
  bottomLeftTableL[i] = 0.0;
  bottomLeftTableR[i] = 0.0;
  bottomRightTableL[i] = 0.0;
  bottomRightTableR[i] = 0.0;
 }
}
*/