#ifndef rosic_SuperOscillator_h
#define rosic_SuperOscillator_h

//// rosic-indcludes:
//#include "rosic_Oscillator.h"
////#include "../filters/rosic_OnePoleFilter.h"

namespace rosic
{

 /**

 This class is derived from the Osc class. As opposed to the
 Oscillator-class, the SuperOscillator-class can produce several voices at
 once (just like the well known SuperSaw waveform - but this "Super"-feature
 is available for all waveforms). Therefore it provides as additional
 parameters "numVoices", "detune" and "freqSpacing".

 */

 class SuperOscillator : public Oscillator
 {

 public:

  //---------------------------------------------------------------------------
	 // construction destruction:

	 SuperOscillator();  ///< Constructor.
	 ~SuperOscillator();  ///< Destructor.

  //---------------------------------------------------------------------------
	 // new parameters:

	 void setNumVoices(int newNumVoices);
  /**< Set the number of voices to be generated. The maximum is 127. */

	 void setFreqSpacing(int newFreqSpacing);
  /**< Switch between linear and exponential frequency spacing between the
       individual voices. Linear spacing means, that the voices are seperated
       by equal frequency differences, exponential spacing means, that the
       voices are separated by the same musical interval. */

  void setDetuneRatio(double newDetuneRatio);
  /**< Sets a ratio factor between adjacent frequencies in unison-mode */

  void setDetunePhaseSpread(double newDetunePhaseSpread);
  /**< Sets the spreading of the start-phases of the unison-voices.*/

  //---------------------------------------------------------------------------
	 // overriden parameter settings:

  void setSampleRate(double newSampleRate);
  ///< Overrides the setSampleRate() method from the Oscillator base-class.

  void setStartPhase(double newStartPhase);
  ///< Overrides the setStartPhase() method from the Oscillator base-class.

  /*
	 void modulateIncrement (sample IncOffset);   //for feedback-FM
  */

  // parameters that are supposed to be updated at sample-rate
  // (inlined access-functions):

  INLINE void setFreq(double newFreq);
  ///< Overrides the setFreq() method from the Oscillator base-class.

  INLINE void setPulseWidth(double newPulseWidth);
  ///< Overrides the setPulseWidth() method from the Oscillator base-class.

  INLINE void setDetune(double newDetune);
  /**< Sets the detune factor for the detuning between the voices. */

  //__forceinline void calcIncrements();
  // calculates the phase-increments for first and second half-period
  // according to freq and pulseWidth


  //---------------------------------------------------------------------------
  // audio processing:

  INLINE double getSampleDraft();
  ///< Overrides the getSample() method from the Oscillator base-class.

  INLINE double getSample();
  /**< Overrides the getSampleAntiAliased() method from the Oscillator
       base-class. */


  //---------------------------------------------------------------------------
  // others:
  void resetPhase();

  static const intA maxVoices = 17; // maximum number of simultaneous voices

  // the phaseindices, phaseIndices[0] has to be accessed from outside for
  // sync (such that the synced osc not just resets it's phase but sets it
  // to the value phaseIndices[0]:
	 doubleA phaseIndices[maxVoices];

  //---------------------------------------------------------------------------
	 // info for an outlying class:

	 bool wraparoundOccurred;
  /**< This flag is true only once per cycle - it is set, everytime the
       phaseIndex is wrapped around. It can be accessed from outside and is
       meant to be used for syncing other oscillators. */

 protected:

  // embedded audio-modules:
  //OnePoleFilter levelDetector; // for detecting th current output-level - used
  //                             // to scale output in unison-mode

	 intA       numVoices;  // current number of simultaneous voices

	 doubleA detune;  // detuning between the center frequency and the most
                  // "left" frequency in semitones

	 doubleA detuneFactor;  // the same as factor for direct multiplication
                        // (=PITCHOFFSET2FREQFACTOR(detune))

  int freqSpacing;

  doubleA detuneRatio;
  doubleA detuneDamp;
  doubleA phaseSpread;

  bool freqsAreDirty; // indicates, if the frequencies of the individual voices
                      // must be recalculated even when the center frequency is
                      // unchanged (this happens due to new detune-settings)

	 doubleA increments1[maxVoices],
          increments2[maxVoices]; // increments for 1st and 2nd halfwave
	                                 // for all oscillator voices

  doubleA incShuffle[maxVoices];  // holds the random values for the increments

	 doubleA ampScale;  // compensates for the gain in amplitude when adding
                     // the voices ( = 1 / sqrt(numVoices) )

  intA sampleCount;

 };

 //-----------------------------------------------------------------------------
 //from here: definitions of the functions to be inlined, i.e. all functions
 //which are supposed to be called at audio-rate (they can't be put into
 //the .cpp file):

 INLINE void SuperOscillator::setFreq(double newFreq)
 {
  static int i;
  doubleA freqFactor, freqFactorRec, accuFactor, accuFactorRec;
  doubleA maxFreqDeviation, freqOffset, incOffset, maxIncOffset;
  doubleA c;

  // save CPU, when the freq didn't change:
  if( newFreq == freq && !freqsAreDirty)
   return;

  freq           = newFreq;
  basicIncrement = tableLengthDbl*freq*sampleRateRec;
  increments1[0] = basicIncrement*pulseFactor1;
  increments2[0] = basicIncrement*pulseFactor2;

  switch(freqSpacing)
  {
  case 0: // frequencies will be spaced exponentially
   {
    // calculate ratio between to adjacent frequenies:
    freqFactor    = pow( detuneFactor , 1.0/((numVoices-1)/2) );

    // sample freqFactor    = pow( numVoices+1 , 1/PITCHOFFSET2FREQFACTOR(12) ); //test
    freqFactorRec = 1.0/freqFactor;

    // init accumulating freq-factors:
    accuFactor    = freqFactor;
    accuFactorRec = freqFactorRec;

    for(i=1; i<(numVoices-1); i+=2 )   //numVoices has to be odd! , i is always odd
    {
     // the odd indices will rise the frequency:
     increments1[i]   = accuFactor * increments1[0];
     increments2[i]   = accuFactor * increments2[0];

     // and the even indices will lower the frequency:
     increments1[i+1] = accuFactorRec * increments1[0];
     increments2[i+1] = accuFactorRec * increments2[0];

     // accumulate the frequency-factors:
     accuFactor    *= freqFactor;
     accuFactorRec *= freqFactorRec;
    } //end of "for(long i=1; i<((numVoices-1)/2); i+=2 )"
   } break;
  case 1:
   {
    // calculate offset bewtween center and lowest frequency:
    maxFreqDeviation = freq - freq/detuneFactor;

    // calculate offset between 2 adjacent frequencies:
    freqOffset = maxFreqDeviation / ((numVoices-1)/2);

    // from this, calculate the offset bewtween 2 adjacent increments:
    incOffset  = tableLengthDbl*freqOffset*sampleRateRec;

    // the first voice above the center frequency (has to look one position
    // back and add the incOffset to the increment which it finds there):
    increments1[1]   = increments1[0] + incOffset;
    increments2[1]   = increments2[0] + incOffset;

    // all other voices look 2 positions back and add or subtract the incOffset:
    for(i=2; i<(numVoices); i++ )   //numVoices has to be odd!
    {
     if( i%2 == 1) //i is odd
     {
      //the odd indices will rise the frequency:
      increments1[i]   = increments1[i-2] + incOffset;
      increments2[i]   = increments2[i-2] + incOffset;
     }
     else
     {
      //and the even indices will lower the frequency:
      increments1[i]   = increments1[i-2] - incOffset;
      increments2[i]   = increments2[i-2] - incOffset;
     }
    } //end of "for(long i=1; i<(numVoices); i++ )"
   } break;
  case 2: // self-similar frequency spacing (variation 1)
   {
    // setup the constant:
    c = detuneRatio;

    // calculate offset bewtween center and lowest frequency:
    maxFreqDeviation = freq - freq/detuneFactor;

    // from this, calculate the maximum increment-offset:
    maxIncOffset  = tableLengthDbl*maxFreqDeviation*sampleRateRec;

    // the first voice above the center frequency (has to look one position
    // back and add the incOffset to the increment which it finds there):
    increments1[1]   = increments1[0] + maxIncOffset;
    increments2[1]   = increments2[0] + maxIncOffset;

    // and so on:
    increments1[2] = increments1[0] - maxIncOffset;
    increments2[2] = increments2[0] - maxIncOffset;

    increments1[3] = increments1[0] + c*(increments1[1]-increments1[0]);
    increments2[3] = increments2[0] + c*(increments2[1]-increments2[0]);

    increments1[4] = increments1[0] + c*(increments1[2]-increments1[0]);
    increments2[4] = increments2[0] + c*(increments2[2]-increments2[0]);

    increments1[5] = increments1[0] + c*(increments1[3]-increments1[0]);
    increments2[5] = increments2[0] + c*(increments2[3]-increments2[0]);

    increments1[6] = increments1[0] + c*(increments1[4]-increments1[0]);
    increments2[6] = increments2[0] + c*(increments2[4]-increments2[0]);

    increments1[7] = increments1[3] + c*(increments1[1]-increments1[3]);
    increments2[7] = increments2[3] + c*(increments2[1]-increments2[3]);

    increments1[8] = increments1[4] + c*(increments1[2]-increments1[4]);
    increments2[8] = increments2[4] + c*(increments2[2]-increments2[4]);

    increments1[9] = increments1[0] + c*(increments1[5]-increments1[0]);
    increments2[9] = increments2[0] + c*(increments2[5]-increments2[0]);

    increments1[10] = increments1[0] + c*(increments1[6]-increments1[0]);
    increments2[10] = increments2[0] + c*(increments2[6]-increments2[0]);

    increments1[11] = increments1[5] + c*(increments1[3]-increments1[5]);
    increments2[11] = increments2[5] + c*(increments2[3]-increments2[5]);

    increments1[12] = increments1[6] + c*(increments1[4]-increments1[6]);
    increments2[12] = increments2[6] + c*(increments2[4]-increments2[6]);

    increments1[13] = increments1[3] + c*(increments1[7]-increments1[3]);
    increments2[13] = increments2[3] + c*(increments2[7]-increments2[3]);

    increments1[14] = increments1[4] + c*(increments1[8]-increments1[4]);
    increments2[14] = increments2[4] + c*(increments2[8]-increments2[4]);

    increments1[15] = increments1[7] + c*(increments1[1]-increments1[7]);
    increments2[15] = increments2[7] + c*(increments2[1]-increments2[7]);

    increments1[16] = increments1[8] + c*(increments1[2]-increments1[8]);
    increments2[16] = increments2[8] + c*(increments2[2]-increments2[8]);
   } break;
  } // end of swicth(freqSpacing)

  freqsAreDirty = false; // the frequencies are all correct now
 }

 INLINE void SuperOscillator::setPulseWidth(double newPulseWidth)
 {
  if( (newPulseWidth>0.0) && (newPulseWidth<1.0) )
   pulseWidth = newPulseWidth;

  //calculate factors for the increments for 1st and 2nd half-wave
  pulseFactor1   = 1.0/(pulseWidth*2.0);
  pulseFactor2   = 1.0/((1-pulseWidth)*2.0);

  freqsAreDirty = true;

  //setFreq has to be called for the new setting to take effect!
 }

 INLINE void SuperOscillator::setDetune(double newDetune)
 {
  if( newDetune >= 0.0 )
   detune = newDetune;

  detuneFactor = PITCHOFFSET2FREQFACTOR(detune);

  freqsAreDirty = true;

  // re-calculation of the increments is not performed here - you will have
  // to call setFreq for the new detuning setting to become active
 }

 INLINE double SuperOscillator::getSampleDraft()
 {
	 return 0.0;
 }

 INLINE double SuperOscillator::getSample()
 {
	 static doubleA increment, outputSample, tmpNoise; //, level;
  static intA i; // j;             //counter for the loop through the voices
  static intA tableNumber;
  //static intA offset;

  if( waveTable == NULL )
   return 0.0;

	 if(waveForm>=3)  // oscillator is really playing a waveform from the
                   // lookup-table for lower values of waveForm it outputs
                   // noise (1,2) or nothing (0)
	 {
		 outputSample = 0.0;

	  // check if in first or second half-wave and choose
   // the appropriate increment:
	  if( phaseIndices[0] < tableLengthDiv2Dbl ) //first half period
	 	 increment = increments1[0];
	  else
	 	 increment = increments2[0];

   // from this increment, decide which table is to be used:
   tableNumber  = ((int)EXPOFDBL(increment)); // & 0x000B;
   tableNumber += 1; // generates frequencies up to nyquist/2 only
   if( tableNumber<=0 )
    tableNumber = 0;
   else if ( tableNumber>11 )
    tableNumber = 11;

   // calculate new phase-index:
   phaseIndices[0] = phaseIndices[0] + increment;

	  // init wraparound flag:
	  wraparoundOccurred = false;

   // wraparound if necessary:
   while ( phaseIndices[0]>=tableLengthDbl )
	  {
    phaseIndices[0] = phaseIndices[0] - tableLengthDbl;
	 	 wraparoundOccurred = true; //set flag
	  }

		 // forward-warparound for cases when increment is negative
   // (can occur due to frequency-modulation):
   while ( phaseIndices[0]<0.0 )
	  {
    phaseIndices[0] = phaseIndices[0] + tableLengthDbl;
	  }

   // O.K. phaseIndices[0] is assured to be valid now - read out the table:
   outputSample = waveTable->getValueLinear(phaseIndices[0], tableNumber);

		 for(i=1; i<numVoices; i++) // loop through the voices 2 - numVoices, the
                              // first voice has to be outside the loop,
                              // because it changes the wraparound-flag
																													 // while the other voices do not
	  {
	   // check if in first or second half-wave and choose the
    // appropriate increment:
	   if( phaseIndices[i] < tableLengthDiv2Dbl ) //first half-period
	  	 increment = increments1[i];
	   else
	  	 increment = increments2[i];

    //from this increment, decide which table is to be used:
    tableNumber = ((int)EXPOFDBL(increment)); // & 0x000B;
    tableNumber += 1; // generates frequencies up to nyquist/2 only
    if( tableNumber<=0 )
     tableNumber = 0;
    else if ( tableNumber>11 )
     tableNumber = 11;

    // calculate new phase-index:
    phaseIndices[i] = phaseIndices[i] + increment;

    //wraparound if necessary:
    while ( phaseIndices[i]>=tableLengthDbl )
	   {
     phaseIndices[i] = phaseIndices[i] - tableLengthDbl;
	   }

	 	 // forward-warparound for cases when increment is negative
    // (can occur due to frequency-modulation):
    while ( phaseIndices[i]<0.0 )
	   {
     phaseIndices[i] = phaseIndices[i] + tableLengthDbl;
	   }

    // O.K. phaseIndices[i] is assured to be valid now - read out the table and
    // add it to the output:
    outputSample += waveTable->getValueLinear(phaseIndices[i], tableNumber);

		 } // end of "for(long i=0; i<numVoices; i++)"

   sampleCount++;

   outputSample *= ampScale; // scales the outut by 1/sqrt(numVoices)

		 return outputSample;
	 } // end of "if(waveForm>=3)"


	 else if(waveForm==1) //white noise
		 return ( (double) (rand()-16384) * rec16384 ); // rand() produces integer
                                                  // numbers between
	                                                 // 0 and 32768

	 else if(waveForm==2) //pink noise
  {
   tmpNoise = ( (double) (rand()-16384) * rec16384 ); // white noise
                                                      // between -1 and 1
   tmpNoise = white2pink.getSample(tmpNoise);
		 return tmpNoise;
  }

	 else return 0.0;
 }

} // end namespace rosic

#endif // rosic_SuperOscillator_h
