//#include "MagicCarpetModulator.h"

//----------------------------------------------------------------------------
// construction/destruction:

MagicCarpetModulator::MagicCarpetModulator()
{
 // init the pointers to the table with a NULL-pointer:
 table      = NULL;

 // init parameters:
 sampleRate         = 44100.0;
 sampleRateRec      = 1.0/sampleRate;
 tableLength        = MAXTABLELENGTH ;
 tableLengthDiv2    = tableLength/2 ;
 tableLengthDbl     = (double) tableLength;
 loop               = true;
 isOff              = false;
 phaseIncrement     = 1.0;
 amount             = 1.0;
 position           = 0.0;
 startPhase         = 0.0;
 bpm                = 140.0;
 periodNumerator    = 1;
 periodDenominator  = 4;

 calculateIncrement();

 // init embedded modules:
 slewRateLimiter.setSampleRate(sampleRate);


 riseRamp.setSampleRate(sampleRate);
 riseRamp.setStart(0.0);
 riseRamp.setEnd(1.0);
 riseRamp.setTime(0.0);
}

MagicCarpetModulator::~MagicCarpetModulator()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void MagicCarpetModulator::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.1 )
  sampleRate = newSampleRate;

 calculateIncrement();

 // tell the embedded modules the new sample-rate:
 slewRateLimiter.setSampleRate(sampleRate);
 riseRamp.setSampleRate(sampleRate);
}

void MagicCarpetModulator::setPeriodNumerator(int newPeriodNumerator)
{
 if( newPeriodNumerator >= 1 )
  periodNumerator = newPeriodNumerator;
 calculateIncrement();
}

void MagicCarpetModulator::setPeriodDenominator(int newPeriodDenominator)
{
 if( newPeriodDenominator >= 1 )
  periodDenominator = newPeriodDenominator;
 calculateIncrement();
}

void MagicCarpetModulator::setBpm(double newBpm)
{
 if( newBpm >= 1.0 )
  bpm = newBpm;
 calculateIncrement();
}

void MagicCarpetModulator::setLoopMode(bool newLoopMode)
{
 loop = newLoopMode;
 reset();
}

void MagicCarpetModulator::setStartPhase(double newStartPhase)
{
 if( newStartPhase >= 0.0 && newStartPhase <= 360.0 )
  startPhase = newStartPhase;
}

void MagicCarpetModulator::setAttack(double newAttack)
{
 slewRateLimiter.setAttackTime(0.001*newAttack);
}

void MagicCarpetModulator::setRelease(double newRelease)
{
 slewRateLimiter.setReleaseTime(0.001*newRelease);
}

void MagicCarpetModulator::setRiseTime(double newRiseTime)
{
 riseRamp.setTime(0.001*newRiseTime);
}

void MagicCarpetModulator::setAmount(double newAmount)
{
 amount = newAmount;
 riseRamp.setEnd(amount);
 if( abs(amount) < 0.01 )
  isOff = true;
 else
  isOff = false;
}

void MagicCarpetModulator::setTableAdress(double* newTableAdress)
{
 table = newTableAdress;
}

void MagicCarpetModulator::setTableLength(int newTableLength)
{
 if( tableLength <= MAXTABLELENGTH )
  tableLength = newTableLength;
 tableLengthDbl = (double) tableLength;
 tableLengthDiv2    = tableLength/2 ;

 // reset the phase-pointers:
 reset();
 calculateIncrement();
}


//----------------------------------------------------------------------------
// others:

void MagicCarpetModulator::calculateIncrement()
{
 phaseIncrement = ((double)periodDenominator*tableLengthDbl*bpm) / 
                  ((double)periodNumerator*sampleRate*240.0);
}

void MagicCarpetModulator::reset()
{
 position = (startPhase/360.0)*tableLengthDbl;

 slewRateLimiter.reset();
 riseRamp.trigger();
}

