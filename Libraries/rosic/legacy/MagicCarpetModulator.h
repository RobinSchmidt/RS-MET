#ifndef MagicCarpetModulator_h
#define MagicCarpetModulator_h

#include "stdlib.h"          // for the NULL macro
#include "MagicCarpetDefinitions.h"
#include "MoreMath.h"
#include "ExponentialRamp.h"
#include "SlewRateLimiter.h"

/**

This class implements the hybrid envelope/LFO-modulator as used in the 
MagicCarpet synthesizer.

*/

class MagicCarpetModulator
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 MagicCarpetModulator();
 ~MagicCarpetModulator();

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< sets the sample-rate

 void setPeriodNumerator(int newPeriodNumerator);
 /**< Sets the numerator of the note-length which the modulator needs to 
      complete one cycle through the table, in units of whole notes 
      (= 4 beats) */

 void setPeriodDenominator(int newPeriodDenominator);
 /**< Sets the denominator of the note-length which the modulator needs to 
      complete one cycle through the table, in units of whole notes 
      (= 4 beats) */

 void setBpm(double newBpm);
 /**< Sets the song-tempo as BPM-value, this is needed as reference to 
      calculate the actual speed according to the target period measured in 
      whole notes. */

 void setLoopMode(bool newLoopMode);
 /**< Set this to "true" when the table should be played as loop. */

 bool isInLoopMode() {return loop;}

 void setStartPhase(double newStartPhase);
 /**< Sets the position in the table to which the modulator will be reset, 
      when it is triggered. The assumption here is, that the modulator has a 
      single cycle waveform loaded and the unit expected here is degrees. */

 void setAttack(double newAttack);
 /**< Sets the attack time of a slew-rate limiter which is in series to the 
      actual table-lookup oscillator. */

 void setRelease(double newRelease);
 /**< Sets the release time of a slew-rate limiter which is in series to the 
      actual table-lookup oscillator. */

 void setRiseTime(double newRiseTime);
 /**< Sets the time which it needs for the modulator to build up to its full 
      modulation depth (given by "amount") */

 void setAmount(double newAmount);
 /**< This sets the amount of modulation (in percent) */

 void setTableAdress(double* newTableAdress);
 /**< Sets the adress of the (integrated) wave-table. Make very sure that
      you call this function BEFORE the first call to getSample() - otherwise
      the getSample()-function will have to deal with an invalid pointer. */

 void setTableLength(int newTableLength);
 ///< Sets the size of the table in samples.

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE double getSample();
 /**< Calculates one output-sample at a time. */

 //---------------------------------------------------------------------------
 // others:

 void reset();
  // resets the table-pointer to its start-position 

 //===========================================================================

protected:

 double  tableLengthDbl; // table-length as double value
 double  position;       // current sample position in the table for left
                         // channel

 double  phaseIncrement; // phase increment

 double  amount;

 int     tableLength;         // size of one table (actually it is one sample
                              // longer for interpolation)
 int     tableLengthDiv2;     // half of the table-length

 bool    loop;            // indicates, if sample should be looped
 bool    isOff;           // indicates, that the modulator is off

 // embedded audio-modules:
 ExponentialRamp riseRamp;
 SlewRateLimiter slewRateLimiter;

 double* table;               // pointer to the wave-table

 double  startPhase;

 double  bpm;
 int     periodNumerator;
 int     periodDenominator;

 double  sampleRate;
 double  sampleRateRec;

 void    calculateIncrement();
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE double MagicCarpetModulator::getSample()
{
 static double truncatedPosition;  
  // truncated position-pointers (integer parts as float)

 static int intPosition;  
  // type-casted integer part of the position-pointers

 static double frac;          
  // fractional parts of the position-pointers

 static double out; 
  // for the actual output-sample

 // catch some special conditions:
 if( isOff == true ||
     (loop == false && position >= tableLengthDbl) ||
     table == NULL )
 {
  return 0.0;
 }

 // wraparound the position-pointer if necesarry:
 while( position >= tableLengthDbl )
  position -= tableLengthDbl;


 // truncate position-pointer and convert to integer:
 truncatedPosition = floor(position);
 intPosition       = (int) truncatedPosition;  
 frac              = position - truncatedPosition;  // OPTIMIZE HERE!

 // calculate ouput by means of linear interpolation:
 out = (1.0-frac)*table[intPosition] + frac*table[intPosition+1]; // OPTIMIZE HERE!

 // increment position-pointer:
 position += phaseIncrement;

 //out *= amount;

 // apply the slew-rate-limiter:
 out = slewRateLimiter.getSample(out);

 // apply the rise-ramp (this includes the amount-factor):
 out *= riseRamp.getSample();

 // return the output-sample:
 return out;
}

#endif // MagicCarpetModulator_h
