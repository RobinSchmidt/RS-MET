#include "DrumSynth.h"

//----------------------------------------------------------------------------
// construction/destruction:

DrumSynth::DrumSynth()
{
 // init member variables:
 sampleRate              = 44100.0;    
 numVoices               = 1;
 numActiveVoices         = 0;
 randomSeed              = 17;
 masterAmplitude         = 1.0;
 mostRecentNote          = -1;
 mostRecentNoteVel       = 0;
 mostRecentNoteDetune    = 0;

 // init all the table-pointers with NULL-pointers:
 osc1TableL        = NULL;
 osc1TableR        = NULL;
 osc2TableL        = NULL;
 osc2TableR        = NULL;

 fourierTransformer.setBlockSize(262144);
}

DrumSynth::~DrumSynth()
{
 
}

//----------------------------------------------------------------------------
// parameter settings:

void DrumSynth::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0 )
  sampleRate = newSampleRate;

 for(long i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].setSampleRate(sampleRate);
}

void DrumSynth::setNumVoices(int newNumVoices)
{
 if( newNumVoices>0 && newNumVoices<=maxNumVoices )
  numVoices = newNumVoices;

 // send a note-off to all the voices:
 for(long i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].noteOn(0,0);
}

int DrumSynth::getNumActiveVoices()
{
 return numActiveVoices;
}

void DrumSynth::invalidateTablePointers(int whichSlot)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setTableAddresses(NULL, NULL);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setTableAddresses(NULL, NULL);
   break;
  } // end switch
 } // end for
}

bool DrumSynth::setWaveTables(int    destinationSlot, 
                              double *newWaveTableL, 
                              double *newWaveTableR,
                              int    newTableLength)
{
 int i;

 // remove the dc-component from the newTable:
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

 double accuL = 0.0; // we want store the integrated waveform, so wee need an 
 double accuR = 0.0; // accumulator for the integrator filter

 switch( destinationSlot )
 {
 case OSC1:
  {
   // tell the oscillators, that their table-addresses are invalid now:
   for(i=0; i<maxNumVoices; i++)
    drumSynthVoiceArray[i].osc1.setTableAddresses(NULL, NULL);

   // free the old memory chunk and invalidate the old-pointers:
   if( osc1TableL != NULL )
   {
    delete[] osc1TableL;
    osc1TableL = NULL;
   }
   if( osc1TableR != NULL )
   {
    delete[] osc1TableR;
    osc1TableR = NULL;
   }

   // try to allocate new memory for the new table:
   osc1TableL = new double[newTableLength+6];
   osc1TableR = new double[newTableLength+6];

   // report failure if the memory-allocation failed; in the case, that only 
   // one table was allocated, delete the other one also:
   if( osc1TableL == NULL || osc1TableR == NULL )
   {
    if( osc1TableL == NULL && osc1TableR != NULL )
    {
     delete[] osc1TableR;
     osc1TableR = NULL;
    }
    if( osc1TableR == NULL && osc1TableL != NULL )
    {
     delete[] osc1TableL;
     osc1TableL = NULL;
    }
    return false;
   }

   // write the integrated table into our member-array: 
   for(i=0; i<newTableLength; i++)
   {
    accuL += newWaveTableL[i];
    accuR += newWaveTableR[i];
    osc1TableL[i] = accuL;
    osc1TableR[i] = accuR;
   }

   // repeat the first 6 samples at the end for the interpolator:
   for(i=newTableLength; i<=newTableLength+5; i++)
   {
    osc1TableL[i] = osc1TableL[i-newTableLength]; 
    osc1TableR[i] = osc1TableR[i-newTableLength]; 
   }

   // tell all the voices the new table-adresses and tableLength:
   for(i=0; i<maxNumVoices; i++)
   {
    drumSynthVoiceArray[i].osc1.setTableLength(newTableLength);
    drumSynthVoiceArray[i].osc1.setTableAddresses(osc1TableL, 
                                                  osc1TableR);
   }

  }
  break;

  // and the same for the other slot:
 case OSC2:
  {
   for(i=0; i<maxNumVoices; i++)
    drumSynthVoiceArray[i].osc2.setTableAddresses(NULL, NULL);
   if( osc2TableL != NULL )
   {
    delete[] osc2TableL;
    osc2TableL = NULL;
   }
   if( osc2TableR != NULL )
   {
    delete[] osc2TableR;
    osc2TableR = NULL;
   }
   osc2TableL = new double[newTableLength+6];
   osc2TableR = new double[newTableLength+6];
   if( osc2TableL == NULL || osc2TableR == NULL )
   {
    if( osc2TableL == NULL && osc2TableR != NULL )
    {
     delete[] osc2TableR;
     osc2TableR = NULL;
    }
    if( osc2TableR == NULL && osc2TableL != NULL )
    {
     delete[] osc2TableL;
     osc2TableL = NULL;
    }
    return false;
   }
   for(i=0; i<newTableLength; i++)
   {
    accuL += newWaveTableL[i];
    accuR += newWaveTableR[i];
    osc2TableL[i] = accuL;
    osc2TableR[i] = accuR;
   }
   for(i=newTableLength; i<=newTableLength+5; i++)
   {
    osc2TableL[i] = osc2TableL[i-newTableLength]; 
    osc2TableR[i] = osc2TableR[i-newTableLength]; 
   }
   for(i=0; i<maxNumVoices; i++)
   {
    drumSynthVoiceArray[i].osc2.setTableLength(newTableLength);
    drumSynthVoiceArray[i].osc2.setTableAddresses(osc2TableL, 
                                                  osc2TableR);
   }
  }
  break;
 } //  end of  switch( destinationSlot )

 return true;
}

bool DrumSynth::setSpectrum(int    destinationSlot, 
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

void DrumSynth::setSlotLoopMode(int  whichSlot, 
                                bool newLoopSwitch)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setLoopMode(newLoopSwitch);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setLoopMode(newLoopSwitch);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotSingleCycleMode(int  whichSlot, 
                                         bool newSingleCycleSwitch)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setSingleCycleMode(newSingleCycleSwitch);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setSingleCycleMode(newSingleCycleSwitch);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotStereoMode(int  whichSlot, 
                                    bool newStereoSwitch)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setStereoMode(newStereoSwitch);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setStereoMode(newStereoSwitch);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotMute(int  whichSlot, 
                              bool newMuteSwitch)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setMute(newMuteSwitch);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setMute(newMuteSwitch);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotDetune(int whichSlot, double newDetune)
{
 double newDetuneFactor = pitchOffsetToFreqFactor(newDetune);
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setDetuneFactor(newDetuneFactor);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setDetuneFactor(newDetuneFactor);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotTuneByKey(int whichSlot, double newTuneByKey)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1TuneByKey = newTuneByKey;
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2TuneByKey = newTuneByKey;
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotLevel(int whichSlot, double newLevel)
{
 double newAmplitude = dB2amp(newLevel);
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setAmplitude(newAmplitude);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setAmplitude(newAmplitude);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotLevelByVel(int whichSlot, double newLevelByVel)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1LevelByVel = newLevelByVel;
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2LevelByVel = newLevelByVel;
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotHpf(int whichSlot, double newHpfCutoff)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setHpfCutoff(newHpfCutoff);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setHpfCutoff(newHpfCutoff);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotLpf(int whichSlot, double newLpfCutoff)
{
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setLpfCutoff(newLpfCutoff);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setLpfCutoff(newLpfCutoff);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setSlotRootKey(int whichSlot, double newRootKey)
{
 double newFundamental = pitchToFreq(newRootKey);
 for(int i=0; i<maxNumVoices; i++)
 {
  switch( whichSlot )
  {
  case(OSC1):
   drumSynthVoiceArray[i].osc1.setFundamentalFreq(newFundamental);
   break;
  case(OSC2):
   drumSynthVoiceArray[i].osc2.setFundamentalFreq(newFundamental);
   break;
  } // end switch
 } // end for
}

void DrumSynth::setFltMode(int newMode)
{
 /*
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].filter.setMode(newMode);
  */
}

void DrumSynth::setFltTwoStages(bool newTwoStagesSwitch)
{
 /*
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].filter.setTwoStages(newTwoStagesSwitch);
  */
}

void DrumSynth::setFltFreq(double newFltFreq)
{ 
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].fltFreq = newFltFreq;
}

void DrumSynth::setFltFreqByKey(double newFltFreqByKey)
{ 
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].fltFreqByKey = newFltFreqByKey;
}

void DrumSynth::setFltFreqByVel(double newFltFreqByVel)
{ 
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].fltFreqByVel = newFltFreqByVel;
}

void DrumSynth::setFltReso(double newFltReso)
{ 
 /*
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].filter.setQ(newFltReso);
 */
}

void DrumSynth::setFltGain(double newFltGain)
{
 /*
 for(int i=0; i<maxNumVoices; i++)
  drumSynthVoiceArray[i].filter.setGain(newFltGain);
  */
  // gain is scaled by 0.5 because the filter uses 2 stages, which should be
  // compensated......not anymore
}

void DrumSynth::setMasterVolume(double newMasterVolume)
{ 

}

//----------------------------------------------------------------------------
// event processing:

void DrumSynth::noteOn(long NoteNumber, long Velocity)
{
 mostRecentNote       = NoteNumber;
 mostRecentNoteVel    = Velocity;
 //mostRecentNoteDetune = Detune;

 //drumSynthVoiceArray[NoteNumber].noteOn(mostRecentNote, mostRecentNoteVel);


 long i = 0;
 
 if(mostRecentNoteVel==0)   // a note-off event occured
 {
  // loop through the voices to find the one which is playing the note
  // for which the note-off was was sent:
  for(i=0; i<numVoices; i++)
  {
   if(drumSynthVoiceArray[i].currentNote == mostRecentNote)
    drumSynthVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
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
   if( !drumSynthVoiceArray[i].ampEnv.endIsReached )
    aVoiceIsPlaying = true;
  }

  // loop through the voices to find a free one:
  int  oldestNoteAge    = 0;
  int  oldestVoice      = 0;       // voice with the oldest note
  for(i=0; i<numVoices; i++)
  {
   // check, if voice i is currently releasing the note, which is coming in,
   // if so, use that voice again:
   // check if the voice i is free:
   if(drumSynthVoiceArray[i].currentNote == mostRecentNote)
   {
    // voice is used for the same note (which is currently releasing) again:
    drumSynthVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
    return; //jump out of the function
   }
   // check if the voice i is free:
   else if(drumSynthVoiceArray[i].ampEnv.endIsReached)
   {
    // voice i is free and can be used for the new note:
    drumSynthVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
    return; //jump out of the function
   }
   // keep track of the oldest note for the case, when no voice is
   // free for the new note (last note priority assignment):
   else if(drumSynthVoiceArray[i].currentNoteAge > oldestNoteAge)
   {
    oldestNoteAge = drumSynthVoiceArray[i].currentNoteAge;
    oldestVoice   = i;
   }
  } // end of for-loop


  // no free voice has been found - set the voice with the oldest note
  // to the new note:
  drumSynthVoiceArray[oldestVoice].noteOn(mostRecentNote, mostRecentNoteVel);
 } // end of else

}

void DrumSynth::setPitchBend(double transpositionInSemitones)
{
 double transpositionFactor 
        = pitchOffsetToFreqFactor(transpositionInSemitones);

 for(int i=0; i<maxNumVoices; i++)
 {
  drumSynthVoiceArray[i].osc1.
   setTranspositionFactor(transpositionFactor);
  drumSynthVoiceArray[i].osc2.
   setTranspositionFactor(transpositionFactor);
 }
}

//----------------------------------------------------------------------------
// internal functions:
