#include "Oscillator.h"

#include "BlepInitialization.h"
 // this initializes our static blep member-variable

//construction/destruction:
Oscillator::Oscillator()
{
 //init member variables:
 tableLengthDbl       = (double) WaveTable::tableLength;  // typecasted versions
 tableLengthDiv2Dbl   = 0.5 * tableLengthDbl;
 tableLengthMinus1Dbl = tableLengthDbl - 1.0;
	tableLengthRec       = 1.0/tableLengthDbl;         //reciprocal
 sampleRate           = 44100.0; 
 freq                 = 440.0;
	basicIncrement       = (tableLengthDbl*freq)/sampleRate;
 increment1           = basicIncrement;
	increment2           = basicIncrement;
	pulseWidth           = 0.5;
	pulseFactor1         = 1/(pulseWidth*2);
 pulseFactor2         = 1/((1-pulseWidth)*2);
 phaseIndex           = 0.0;
 startIndex           = 0.0;

 //tableLength          = WaveTable::tableLength;



 rec16384             = 1.0/16384.0; //reciprocal of 16384

 //blepTimeScaler       = 1.0/((double)blepOversampling);
 blepOversampling = 128.0;
 blepLengthDbl    = (double) blepLength;
 blepIncrement    = blepOversampling/2.0;
 blepPosition     = 8192.0;
 blepSizeScaler   = 0.0;
 blepStartLevel   = 0.0;
 previousOutput   = 0.0;

 waveTable = NULL;
 waveTable = new WaveTable(); // for test purposes only

	//somewhat redundant:
	setSampleRate(44100.0);  // sampleRate = 44100 Hz by default
	setFreq      (440.0);    // frequency = 440 Hz by default
	setStartPhase(0.0);      // sartPhase = 0 by default
	setWaveForm  (3);        // fills the table with the default waveform (sine)
	resetPhase   ();
}

Oscillator::~Oscillator()
{
 if(waveTable)
  delete waveTable; // for test purposes only
}

//----------------------------------------------------------------------------
//parameter settings:

void Oscillator::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;
	sampleRateRec  = 1.0 / sampleRate;

 //calculate the increment as if there were no pulseWidth-parameter:
 basicIncrement = tableLengthDbl*freq*sampleRateRec; 

	//the increment for the 1st and 2nd half-wave are derived via the pulse-factor:
	increment1     = basicIncrement*pulseFactor1;
	increment2     = basicIncrement*pulseFactor2;
}

void Oscillator::setWaveForm(int newWaveForm)
{
	waveForm = newWaveForm;
 if( waveTable != NULL )
  waveTable->setWaveform(waveForm-2);
}

void Oscillator::setWaveTable(WaveTable* newWaveTable)
{
 waveTable = newWaveTable;
}

void Oscillator::setStartPhase(double StartPhase)
{
 if( (StartPhase>=0) && (StartPhase<=360) )
  startIndex = (StartPhase/360.0)*tableLengthDbl;
}

//----------------------------------------------------------------------------
//event processing:

void Oscillator::resetPhase()
{
 phaseIndex = startIndex;
}

void Oscillator::setPhase(double PhaseIndex)
{
 phaseIndex = startIndex+PhaseIndex;
}


void Oscillator::triggerSync(double timeStamp)
{
 static intA    tableNumber;
 static doubleA tableValueHere;

 phaseIndex = startIndex + timeStamp*basicIncrement;

 // wraparound if necessary (actually this should occur only if the 
 // initialPhaseAdvance has an unreasonably large value):
 while ( phaseIndex>=tableLengthMinus1Dbl )
  phaseIndex = phaseIndex - tableLengthMinus1Dbl;

	// forward-warparound for cases when phaseIndex is negative (can occur due
 // to sync):
 while ( phaseIndex<0.0 )
  phaseIndex = phaseIndex + tableLengthMinus1Dbl;

 // select the table and calculate the offset:
 tableNumber  = ((int)EXPOFDBL(basicIncrement)); 
 tableNumber += 1; // generates frequencies up to nyquist/2 only
 //tableNumber += 3; // test
 if( tableNumber<=0 )
  tableNumber = 0;
 else if ( tableNumber>11 )
  tableNumber = 11;

 // calculate, what the output would be at this index in order to determine 
 // the scaling factor for the blep:
 tableValueHere = waveTable->getValueLinear(phaseIndex, tableNumber);

 // determine the scaling factor for the minBlep:
 blepSizeScaler = tableValueHere - previousOutput;
 blepStartLevel = previousOutput;
 blepPosition   = blepIncrement*timeStamp;
}














