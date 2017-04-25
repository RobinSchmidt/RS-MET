#include "MagicCarpetOscillator.h"

//----------------------------------------------------------------------------
// construction/destruction:

MagicCarpetOscillator::MagicCarpetOscillator()
{
 // init the pointers to the table with a NULL-pointer:
 tableL = NULL;
 tableR = NULL;    

 double phaseIncrement;

 // init parameters:
 x1_L                  = 0.0;
 x1_R                  = 0.0;
 amplitude             = 1.0;
 ampScaler             = 1.0;
 sampleRate            = 44100.0;
 sampleRateRec         = 1.0/sampleRate;
 tableLength           = MAXTABLELENGTH ;
 bitMask               = tableLength-1;
 //tableLengthDiv2    = tableLength/2 ;
 //tableLengthDbl     = (double) tableLength;
 fundamentalFreq       = 110.0; 
 fundamentalFreqRec    = 1.0/fundamentalFreq;
 detuneFactor          = 1.0;
 transpositionFactor   = 1.0;
 freq                  = 1000.0;
 finalFreq             = transpositionFactor*detuneFactor*freq;
 loop                  = true;
 singleCycleMode       = false;
 stereo                = true;
 mute                  = false;
 draftMode             = false;
 phaseIncrement        = freq  * fundamentalFreqRec;
 incrementInt          = (int) phaseIncrement;
 incrementFrac         = phaseIncrement - (double) incrementInt;
 finalAmplitude        = amplitude/phaseIncrement;
 //positionL          = 0.0;
 positionInt           = 0;
 positionFrac          = 0.0;
 //positionR          = (double) tableLength/2;

 // init embedded modules:
 filter.setSampleRate(sampleRate); 
 //filterR.setSampleRate(sampleRate);
}

MagicCarpetOscillator::~MagicCarpetOscillator()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void MagicCarpetOscillator::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.1 )
  sampleRate = newSampleRate;

 // trigger a re-calculation of the phase-increment:
 setFrequency(freq);

 // tell the embedded modules the new sample-rate:
 filter.setSampleRate(sampleRate);
 //filterR.setSampleRate(sampleRate);
}

void MagicCarpetOscillator::setAmpScaler(double newAmpScaler)
{
 ampScaler      = newAmpScaler;
 finalAmplitude = ampScaler*amplitude/((double)incrementInt+incrementFrac);
}

void MagicCarpetOscillator::setFundamentalFreq(double newFundamentalFreq)
{
 if( newFundamentalFreq > 0.01 )
  fundamentalFreq = newFundamentalFreq;
 fundamentalFreqRec = 1.0/fundamentalFreq;

 // trigger a re-calculation of the phase-increment:
 setFrequency(freq);
}

void MagicCarpetOscillator::setDetuneFactor(double newDetuneFactor)
{
 if( newDetuneFactor > 0.00001 )
  detuneFactor = newDetuneFactor;

 // trigger a re-calculation of the phase-increment:
 setFrequency(freq);
}

void MagicCarpetOscillator::setTranspositionFactor(double newTranspositionFactor)
{
 if( newTranspositionFactor > 0.00001 )
  transpositionFactor = newTranspositionFactor;

 // trigger a re-calculation of the phase-increment:
 setFrequency(freq);
}

void MagicCarpetOscillator::setLoopMode(bool newLoopMode)
{
 loop = newLoopMode;
 reset();
}

void MagicCarpetOscillator::setSingleCycleMode(bool newSingleCycleMode)
{
 singleCycleMode = newSingleCycleMode;
 reset();
 setFrequency(freq);
}

void MagicCarpetOscillator::setStereoMode(bool newStereoMode)
{
 stereo = newStereoMode;
 reset();
}

void MagicCarpetOscillator::setMute(bool shouldBeMuted)
{
 mute = shouldBeMuted;
}

void MagicCarpetOscillator::setTableAddresses(double* newTableAddressL, 
                                              double* newTableAddressR)
{
 tableL = newTableAddressL;
 if( newTableAddressR == NULL )
  tableR = tableL;
 else
  tableR = newTableAddressR;
}

void MagicCarpetOscillator::setTableLength(int newTableLength)
{
 if( newTableLength <= MAXTABLELENGTH )
  tableLength = newTableLength;
 bitMask = tableLength-1;

 // reset the phase-pointers:
 reset();

 // trigger a re-calculation of the phase-increment:
 setFrequency(freq);
}

void MagicCarpetOscillator::setOffsetBetweenLeftAndRight(int newOffset)
{

}

void MagicCarpetOscillator::setAmplitude(double newAmplitude)
{
 amplitude      = newAmplitude;
 finalAmplitude = ampScaler*amplitude/((double)incrementInt+incrementFrac);
}

void MagicCarpetOscillator::setHpfCutoff(double newHpfCutoff)
{
 filter.setHpfCutoff(newHpfCutoff);
 //filterR.setHpfCutoff(newHpfCutoff);
}

void MagicCarpetOscillator::setLpfCutoff(double newLpfCutoff)
{
 filter.setLpfCutoff(newLpfCutoff);
 //filterR.setLpfCutoff(newLpfCutoff);
}

//----------------------------------------------------------------------------
// others:

/*
void MagicCarpetOscillator::reset()
{
 positionL = 0.0;
 
 filterL.resetBuffers();
 filterR.resetBuffers();
}
*/

void MagicCarpetOscillator::reset(double offset)
{
 double positionL = offset;

 // do wraparound:
 while( positionL >= (double) tableLength )
  positionL -= (double) tableLength;

 // do forward wraparound, if necessary:
 while( positionL < 0.0 )
  positionL += (double) tableLength;

 positionInt  = (int) positionL;
 positionFrac = positionL - (double) positionInt;

 x1_L = 0.0;
 x1_R = 0.0;

 filter.resetBuffers();
 //filterR.resetBuffers();
}

void MagicCarpetOscillator::decrementPhase()
{
 double positionL  =  (double) positionInt  + positionFrac;
 positionL        -= ((double) incrementInt + incrementFrac);

 // do forward wraparound, if necessary:
 while( positionL < 0.0 )
  positionL += (double) tableLength;

 positionInt  = (int) positionL;
 positionFrac = positionL - (double) positionInt;
}
