#ifndef DrumSynth_h
#define DrumSynth_h

#include "Definitions.h"
#include "FourierTransformer.h"
#include "MoreMath.h"
#include "AudioModule.h"
#include "LowpassHighpassStereo.h"
#include "DrumSynthVoice.h"
#include "OnePoleFilter.h"

/**

DrumSynth is a synthesizer based on a mix of special audio-samples. blah...

*/

enum sampleSlots
{
 OSC1 = 0,
 OSC2
};


class DrumSynth
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 DrumSynth();   ///< Constructor.
 ~DrumSynth();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Sets the sample-rate.

 void setNumVoices(int newNumVoices); 
 ///< Sets the maximum number of voices (polyphony).

 int getNumActiveVoices(); 
 ///< Returns the number of currently playnig voices (polyphony).

 //void clearWaveTables(int destinationSlot);

 void invalidateTablePointers(int whichSlot);
 /**< This function is called from the plugIn before a new table is passed to
      set the pointers temporarily to NULL. */

 bool setWaveTables(int     destinationSlot,
                    double* newWaveTableL, 
                    double* newWaveTableR,
                    int     newTableLength);
 /**< Passes a pointers to the wave-tables which should be used for left and 
      right channel. The return-value reporst, if the required memory could be
      allocated.*/

 bool setSpectrum(int     destinationSlot,
                  double* newSpectrum, 
                  int     newSpectrumLength);

 // switches for sample-playback options:
 void setSlotLoopMode       (int whichSlot, bool newLoopModeSwitch);
 void setSlotSingleCycleMode(int whichSlot, bool newSingleCycleSwitch);
 void setSlotStereoMode     (int whichSlot, bool newStereoSwitch);
 void setSlotMute           (int whichSlot, bool newMuteSwitch);

 // numeric parameters for the sample-playback:
 void setSlotDetune    (int whichSlot, double newDetune);
 void setSlotTuneByKey (int whichSlot, double newTuneByKey);
 void setSlotLevel     (int whichSlot, double newLevel);
 void setSlotLevelByKey(int whichSlot, double newLevelByKey);
 void setSlotLevelByVel(int whichSlot, double newLevelByVel);
 void setSlotHpf       (int whichSlot, double newHpfCutoff);
 void setSlotLpf       (int whichSlot, double newLpfCutoff);
 void setSlotRootKey   (int whichSlot, double newRootKey);

 // the general filter envelope:
 //void switchFilterOn    (bool   shouldBeOn); 
 void setFltMode        (int    newFltMode);
 void setFltTwoStages   (bool   newFltTwoStagesSwitch);
 void setFltFreq        (double newFltFreq);
 void setFltFreqByKey   (double newFltFreqByKey);
 void setFltFreqByVel   (double newFltFreqByVel);
 void setFltReso        (double newFltReso);
 void setFltGain        (double newFltGain);

 // the amplitude and vector-mixer settings:
 void setMasterVolume        (double newMasterVolume);  

 //---------------------------------------------------------------------------
 // event processing:

 void noteOn(long NoteNumber, long Velocity);
 void setPitchBend(double transpositionInSemitones);

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double *outL, double *outR);   
 /**< Calculates the output-samples for both channels and stores them at the
      adresses of *outL and *outR. */

//============================================================================

protected:


 static const int maxNumVoices = 32;

 // embedded audio-modules:
 FourierTransformer  fourierTransformer;
 DrumSynthVoice      drumSynthVoiceArray[maxNumVoices];
 BreakpointModulator ampEnv;

 // parameters:
 double sampleRate;
 int    numVoices;
 double masterAmplitude;

 // MIDI-data and voice-management:
 int mostRecentNote, mostRecentNoteVel, mostRecentNoteDetune;
 //....

 // allocate the arrays for the actual wavetables (the +6 is needed for
 // the interpolator):
 double* osc1TableL;
 double* osc1TableR;
 double* osc2TableL;
 double* osc2TableR;

 int    randomSeed;
 int    numActiveVoices; // number of currently playing voices
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE void DrumSynth::getSampleFrameStereo(double *outL, double *outR)
{
 double accuL, accuR;


 // initialize the output slots with zero in order to accumulate over the 
 // voices and then loop through the (active) voices:
 accuL             = 0.0;
 accuR             = 0.0;
 numActiveVoices   = 0;
 for(int i=0; i<numVoices; i++)
 {
  if(!drumSynthVoiceArray[i].ampEnv.endIsReached)
  {
   // let voice i generate its output-sample and accumulate it:
   drumSynthVoiceArray[i].getSampleFrameStereo(&accuL, &accuR);

   // increment the variable which keeps track of the active voices:
   numActiveVoices++;
  }
 }
 *outL = masterAmplitude*accuL;
 *outR = masterAmplitude*accuR;
}

#endif // DrumSynth_h
