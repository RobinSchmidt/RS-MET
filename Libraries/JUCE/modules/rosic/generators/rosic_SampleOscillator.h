#ifndef rosic_SampleOscillator_h
#define rosic_SampleOscillator_h

// standard-library includes:
#include <stdlib.h>          // for the NULL macro

// rosic-indcludes:
//#include "../basics/rosic_Definitions.h"
#include "../basics/rosic_Interpolator.h"
#include "../filters/rosic_LowpassHighpassStereo.h"
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

 /**

 This class implements the spectral-cluster-oscillator as used in the 
 MagicCarpet synthesizer.

 */

 class SampleOscillator
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

  SampleOscillator();
  ~SampleOscillator();

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< sets the sample-rate

  void setAmpScaler(double newAmpScaler);
  /**< Set's a scale-factor for the amplitude. */

  void setFundamentalFreq(double newFundamentalFreq);
  /**< Set's the fundamental frequency of the table. */

  void setDetuneFactor(double newDetuneFactor);
  /**< Set's the detune-factor for this oscillator. */

  void setTranspositionFactor(double newTranspositionFactor);
  /**< Set's the transposition-factor for this oscillator which comes from the 
       pitch-wheel. */

  void setLoopMode(bool newLoopMode);
  /**< Set this to "true" when the table should be played as loop. */

  bool isInLoopMode() {return loop;}

  void setSingleCycleMode(bool newSingleCycleMode);
  /**< Set this to "true" when the table represents a single cycle of a 
       waveform. */

  // TODO: implement MultiCycleMode -> the loop represents an arbitrary
  // (integer) number of cycles, function setCyclesPerLoop or something
  // will be required -> can be consolidated in a single status variable
  // cycleMode for which a zero value indicates "off"

  void setStereoMode(bool newStereoMode);
  /**< Set this to "true" when it makes sense to offset the read-pointers for
       left and right output by tableLength/2 to create a stereo-effect. */

  void setMute(bool shouldBeMuted);
  /**< Set this to "true" when this instance should be muted. */

  void setTableAddresses(double* newTableAddressL, 
                         double* newTableAddressR = NULL);
  /**< Sets the adresses of the (integrated) wave-tables. Make very sure that
       you call this function BEFORE the first call to getSample() - otherwise
       the getSample()-function will have to deal with an invalid pointer. When
       no second pointer or a NULL-pointer is passed for the right channel, it
       will be set to the same address to which the left-channel pointer points
       to, such that both channels play the same signal (mono). */

  void setTableLength(int newTableLength);
  ///< Sets the size of the table in samples.

  void setOffsetBetweenLeftAndRight(int newOffset);
  /**< Set the offset of the phase pointers between left and right channel. 
       Set this to tableLength/2 for maximum stereo width. */

  void setAmplitude(double newAmplitude);
  ///< Sets the output amplitude of this instance.

  void setHpfCutoff(double newHpfCutoff);
  ///< Sets the cutoff-frequency of the static highpass-filter.

  void setLpfCutoff(double newLpfCutoff);
  ///< Sets the cutoff-frequency of the static lowpass-filter.


  INLINE void setFrequency(double newFrequency);
  /**< Sets the nominal frequency to be played. The actual playback-frequency 
       will additionally take into account the detune-factor as set by 
       setDetuneFactor(). */

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE void getSampleFrameStereo(double* outL, double* outR);
  /**< Calculates one stereo sample-frame at a time. */

  //---------------------------------------------------------------------------
  // others:
  //void reset();
   // resets the table-pointer to its start-position and resets the filters.

  void reset(double offset = 0.0);
   // resets the table-pointer to its start-position minus some offset (this is
   // needed to initialize the differentiators properly in the MagicCarpetVoice
   // class);

  void decrementPhase();

  //void resetFilters();
   // resets the filters only
   

 protected:

  double positionFrac;
  double incrementFrac;
  double finalAmplitude;  // final factor for the output (including reciprocal 
                          // of the phaseIncrement)

  double x1_L, x1_R;      // state variables for the differentiators

  int    tableLength;     // size of one table (actually it is one sample
                          // longer for interpolation)
  int    positionInt;
  int    incrementInt;

  bool   mute;            // indicates, if output should be quiet
  bool   loop;            // indicates, if sample should be looped
  bool   stereo;          // indicates if left and right output pointer should
                          // be spaced tableLength/2 apart
  bool   singleCycleMode; // indicates, if the table is a single cycle
  bool   draftMode;


  int     bitMask; // not used at the moment

  // embedded audio-modules:
  Interpolator          interpolator;
  LowpassHighpassStereo filter;

  double* tableL;              // pointer to the (integrated) wave-table for L
  double* tableR;              // pointer to the (integrated) wave-table for R

  double  amplitude;           // output amplitude (without the scaling with the
                               // ampScaler and the recirocal increment)
  double  ampScaler;           // addtional amplitude scaler (used for 
                               // slot-volume by velocity)

  double  detuneFactor;        // a detune-factor
  double  transpositionFactor; // transpose from pitch-wheel 
  double  freq;                // nominal frequency (without detune factor)
  double  finalFreq;           // final frequency including detuning

  double  sampleRate;
  double  sampleRateRec;

  //

 public:

  double  fundamentalFreq;     // fundamental frequency of the tables
  double  fundamentalFreqRec;  // reciprocal of the fundamental frequency
 };

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):

 INLINE void SampleOscillator::setFrequency(double newFrequency)
 {
  if( newFrequency > 0 )
   freq = newFrequency;

  // take transposition and detuning into account for the calcluation of the
  // final frequency:
  finalFreq = transpositionFactor*detuneFactor*freq;

  // calculate the new phase increment:
  double phaseIncrement;
  if( singleCycleMode == true )
   phaseIncrement = finalFreq * (double) tableLength * sampleRateRec;
  else
   phaseIncrement = finalFreq * fundamentalFreqRec; 

  incrementInt  = (int) phaseIncrement;
  incrementFrac = phaseIncrement - (double) incrementInt;

  // for click-TEST:
  //phaseIncrement = 1.0;

  // calculate te final amplitude factor - this should include the reciprocal
  // of the phase-increment to accomodate for the gain differencing-filter at
  // high frequencies (which is applied in the MagicCarpetVoice-class to the 
  // summed output of the 4 oscillators);
  finalAmplitude     = ampScaler*amplitude/phaseIncrement;
 }

 INLINE void SampleOscillator::getSampleFrameStereo(double* outL, 
                                                         double* outR)
 {
  double tmpL, tmpR;

  // catch some special conditions:
  if( (mute == true) || 
      //(loop == false && positionInt >= tableLength) ||
      (tableL == NULL) ||
      (tableR == NULL) )
  {
   *outL = 0.0;
   *outR = 0.0;
   return;
  }

  // wraparound the integer part of the position-pointer if necesarry:
  if( loop == true )
  {
   //positionInt &= bitMask;
   while( positionInt >= tableLength )
    positionInt -= tableLength;
     // replace this with while (position >= loopEnd) position -= loopLength;
     // for introducing arbitrary loops
  }
  else
  {
   // if we are not in loop-mode, we keep the position pointer fixed at the end
   // of  the table, when the end is reached:
   if( positionInt >= tableLength-3 )
   {
    positionInt  = tableLength-3;
    positionFrac = 0.0;
   }
  }

  if( stereo )
  {
   // read out the wavetable (with interpolation) for left and right
   // channel:
   tmpL = interpolator.getSampleParabolic4p2o2x(positionFrac, &(tableL[positionInt]) );
   tmpR = interpolator.getSampleParabolic4p2o2x(positionFrac, &(tableR[positionInt]) );

   // apply differentiator and update its state-veriables:
   *outL = tmpL - x1_L;
   *outR = tmpR - x1_R;
   x1_L  = tmpL;
   x1_R  = tmpR;

   // apply lowpass and highpass filtering and global volume control:
   filter.getSampleFrameStereo(outL, outR, outL, outR);
   *outL *= finalAmplitude;
   *outR *= finalAmplitude;
  }
  else // mono mode
  {
   tmpL = interpolator.getSampleParabolic4p2o2x(positionFrac, &(tableL[positionInt]) );

   *outL = tmpL - x1_L;
   x1_L  = tmpL;

   *outL *= finalAmplitude;
   *outR  = *outL;
   filter.getSampleFrameStereo(outL, outR, outL, outR);
  }

  // increment position-pointer:
  positionInt  += incrementInt;
  positionFrac += incrementFrac;
  if( positionFrac >= 1.0 )
  {
   positionFrac -= 1.0;
   positionInt  += 1;
  }
 }

} // end namespace rosic

#endif // rosic_SampleOscillator_h
