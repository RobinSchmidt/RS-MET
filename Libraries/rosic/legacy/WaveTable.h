#ifndef WaveTable_h
#define WaveTable_h

#include "FourierTransformer.h"
#include "Interpolator.h"
#include "Tools.h"
#include "stdlib.h"             // for rand()
#include "White2PinkFilter.h"
//#include "Oscillator.h"

/**

This is a class for generating and storing a single-cycle-waveform in a 
lookup-table and retrieving values form it at arbitrary positions by means of
interpolation.

*/

class WaveTable  
{

 // Oscillator and SuperOscillator classes need access to certain protected 
 // member-variables (namely the tableLength and related quantities), so we 
 // declare them as friend-classes:
 friend class Oscillator;
 friend class SuperOscillator;

public:
 
 //---------------------------------------------------------------------------
 // construction/destruction:

	WaveTable();          ///< Default constructor.
	~WaveTable();         ///< Destructor.

 //---------------------------------------------------------------------------
 // parmeter-settings:

 void setWaveform(int newWaveform);
 /**< Selects a waveform from the set of built-in wavforms. The object 
      generates the prototype-waveform by some algorithmic rules and renders 
      various bandlimited version of it via FFT/iFFT. */

 void setWaveform(double* newWaveform, int lengthInSamples);
 /**< Overloaded function to set the waveform form outside this class. This 
      function expects a poniter to the prototype-waveform to be handed over 
      along with the length of this waveform. It copies the values into the 
      internal buffers and renders various bandlimited version via FFT/iFFT.
      still ToDo: Interpolation for the case that lengthInSamples does not 
      match the length of the internal table-length. */

 //---------------------------------------------------------------------------
 //audio processing:

 INLINE double getValueLinear(double phaseIndex, int tableIndex);
 /**< Returns the value at position "phaseIndex" of table "tableIndex" with 
      linear interpolation */

 INLINE double getValueHermite(double phaseIndex, int tableIndex);
 /**< Returns the value at position "phaseIndex" of table "tableIndex" with 
      hermite interpolation */

protected:

 // data:
 //static const intA tableLength = 2048;   
 static const intA tableLength = 8192;   
  /*   Length of the lookup-table. The actual length of the allocated
       memory is 4 samples longer, to store additional samples for the
       interpolator (which are the same values as at the beginning of the
       buffer) */

 // some related quantities:
 static const intA tableLengthDiv2   = tableLength/2;
 static const intA tableLengthMinus1 = tableLength-1;
 doubleA tableLengthDbl;
 doubleA tableLengthDiv2Dbl;
 doubleA tableLengthMinus1Dbl;
	doubleA tableLengthRec;        
 doubleA rec16384;   // reciprocal of 16384
  // for some reason, the doubles cannot be declared as static const and cannot 
  // be initialized here therefore - so we have to do this in the contructor 
  // instead


 static const intA numTables = 12;
  /*   The Oscillator class uses a one table-per octave multisampling to avoid
       aliasing. With a table-size of 8192 and a sample-sample rate of
       44100, the 12th table will have a fundamental frequency (the
       frequency where the increment is 1) of 11025 which is good for the
       highest frequency. */

 int    waveform;   // index of the currently chosen native waveform
 double sampleRate; // the sampleRate

	doubleA prototypeTable[tableLength];
  /* this is the prototype-table with full bandwidth.
	  one additional sample (same as prototypeTable[0]) for linear interpolation
   without need for table wraparound at the last sample (-> saves one
   if-statement each audio-cycle) ...and a three further addtional samples
   for more elaborate interpolations like cubic (not implemented yet, also:
   the fillWith...()-functions don't support these samples yet). */

 doubleA tableSet[numTables][tableLength+4]; 
  /* The multisample for anti-aliased waveform generation. The 4 additional 
     values are equal to the first 4 values in the table for easier 
     interpolation. The first index is for the table-number - index 0 accesses
     the first version which has full bandwidth, index 1 accesses the second 
     version which is bandlimited to Nyquist/2, 2->Nyquist/4, 3->Nyquist/8, 
     etc. */


 // functions to fill table with the built-in waveforms (these functions are
 // called from setWaveform(int newWaveform):
 void fillWithSine();
 void fillWithTriangle();
 void fillWithPulse();
 void fillWithSawUp();
 void fillWithSawDown();
 void fillWithPeakUp();
 void fillWithPeakDown();
 void fillWithMoogSaw();
 void fillWithSqrMinusSin();
 void fillWithTriangleModulation();
 void fillWithSquareSweep1();
 void fillWithSineSweepUp(double startFreq, double endFreq);
 void fillWithSineSweepDown(double startFreq, double endFreq);

 void fillWithWhiteSineOcts(); // octave-spaced sine-waves 
                               // (amplitude=const) up to nyquist

 void fillWithPinkSineOcts();  // octave-spaced sine-waves 
                               // (amplitude = 1/sqrt(n)) up to nyquist

 void fillWithBrownSineOcts(); // octave-spaced sine-waves 
                               // (amplitude = 1/n) up to nyquist

 void fillWithSawOcts();

 // noise cycles of different color:
 void fillWithNoiseCycle(int noiseColor); // 0:white, 1:pink, 2:brown


 void fillWithSineOcts(int NumOcts);    // fills the table with 
                                        // octave-spaced sine waves

 void fillWithSineOctsAlt(int NumOcts); // fills the table with
                                        // octave-spaced sine waves
                                        // and alternating signs

 void fillWithTriangleOcts   (int NumOcts);  
 void fillWithTriangleOctsAlt(int NumOcts); 
 void fillWithPulseOcts      (int NumOcts);  
 void fillWithPulseOctsAlt   (int NumOcts); 
 void fillWithSawUpOcts      (int NumOcts);  
 void fillWithSawUpOctsAlt   (int NumOcts); 
 void fillWithSawDownOcts    (int NumOcts);  
 void fillWithSawDownOctsAlt (int NumOcts);

 void fillWithH1_4_7();   // wave with 1st, 4th, 7th, ... harmonic
 void fillWithH2_5_8();   // wave with 2nd, 5th, 8th, ... harmonic
 void fillWithH3_6_9();   // wave with 3rd, 6th, 9th, ... harmonic


	// internal helper methods:

 void initPrototypeTable(); 
  // fills the "prototypeTable"-variable with all zeros

 void initTableSet();  
  // fills the "tableSet"-variable with all zeros

	void removeDC();      
  // removes dc-component from the waveform in the prototype-table

	void normalize();     
  // normalizes the amplitude of the prototype-table to 1.0

 void reverseTime();   
  // time-reverses the prototype-table

 void generateMultiSample();  
  // generates a multisample from the prototype table, where each of the
  // successive tables contains one half of the spectrum of the previous one

 // embedded objects:
 White2PinkFilter   white2pink;
 FourierTransformer fourierTransformer;
 Interpolator       interpolator;
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE double WaveTable::getValueLinear(double phaseIndex, int tableIndex)
{
	doubleA floorOfIndex, fracIndex;//, out;
 intA    intIndex;

 // ensure, that the table index is in the valid range:
 if( tableIndex<=0 )
  tableIndex = 0;
 else if ( tableIndex>numTables )
  tableIndex = 11;

 // calculate integer and fractional part of the phaseIndex:
 floorOfIndex  = floor(phaseIndex);         // integer part as double
 fracIndex     = phaseIndex - floorOfIndex; // fractional part  
 intIndex      = (int) floorOfIndex;        // integer part a int

  // lookup value in the table with linear interpolation and return it:
 return interpolator.getSampleLinear(fracIndex, 
                                     &(tableSet[tableIndex][intIndex]) );
}

INLINE double WaveTable::getValueHermite(double phaseIndex, int tableIndex)
{
	static doubleA floorOfIndex, fracIndex, out;
 static intA    intIndex;

 // ensure, that the table index is in the valid range:
 if( tableIndex<=0 )
  tableIndex = 0;
 else if ( tableIndex>numTables )
  tableIndex = 11;

 // calculate integer and fractional part of the phaseIndex:
 floorOfIndex  = floor(phaseIndex);         // integer part as double
 fracIndex     = phaseIndex - floorOfIndex; // fractional part  
 intIndex      = (int) floorOfIndex;        // integer part a int

  // lookup value in the table with linear interpolation and return it:
 return interpolator.getSampleHermite4p3o(fracIndex, 
                                          &(tableSet[tableIndex][intIndex]) );
}

#endif // WaveTable_h
