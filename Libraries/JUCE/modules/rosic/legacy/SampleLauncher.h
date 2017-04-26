#ifndef SampleLauncher_h
#define SampleLauncher_h

#include "Definitions.h"
#include "SampleLauncherDefinitions.h"
#include "SampleLauncherVoice.h"
#include "MoreMath.h"


/**

SampleLauncher is a a DrumSampler...

*/


class SampleLauncher
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 SampleLauncher();   ///< Constructor.
 ~SampleLauncher();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Sets the sample-rate.

 void setNumVoices(int newNumVoices); 
 ///< Sets the maximum number of voices (polyphony).

 void setMasterVolume(double newMasterVolume); 
 ///< Sets the master output volume in dB.

 void clearSample(int whichSlot);

 void setSample(int    whichSlot,
                float* newSampleL, 
                float* newSampleR,
                int    newSampleLength,
                double newFileSampleRate);

 // numeric parameters for the sample-playback:
 void setSlotDetune       (int whichSlot, double newDetune);
 void setSlotLevel        (int whichSlot, double newLevel);
 void setSlotOutputChannel(int whichSlot, int    newOutputChannel);

 //---------------------------------------------------------------------------
 // event processing:

 void noteOn(long NoteNumber, long Velocity);

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameMultiChannel(float *outL_1, float *outR_1,
                                        float *outL_2, float *outR_2,
                                        float *outL_3, float *outR_3,
                                        float *outL_4, float *outR_4,
                                        float *outL_5, float *outR_5,
                                        float *outL_6, float *outR_6,
                                        float *outL_7, float *outR_7,
                                        float *outL_8, float *outR_8);   
 /**< Calculates the output-samples for the 8 stereo-pairs (16 channels) */

//============================================================================

protected:

 // embedded audio-modules:
 SampleLauncherVoice sampleLauncherVoiceArray[MAXVOICES];

 // parameters:
 double playbackSampleRate;
 int    numVoices;
 float  masterAmplitude;

 // MIDI-data and voice-management:
 int mostRecentNote, mostRecentNoteVel, mostRecentNoteDetune;

 // allocate the arrays for the actual samples, as they are:
 float theSamplesOriginal[128][MAX_SAMPLE_LENGTH][2];

 // allocate the arrays for the samples, transposed to the desired pitch:
 float theSamplesTransposed[128][MAX_SAMPLE_LENGTH][2];

   // in the above declarations, the first index refers to the key to which 
   // the sample belongs, the second index is the position in the sample and 
   // the third is the stereo-channel ( -> samples are stored L/R-interleaved
   // in memory)

 // for the individual samples file-sample rates and lengths
 int    sampleLengths[128];
 double fileSampleRates[128];

 void zeroAllTables();
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE void SampleLauncher::getSampleFrameMultiChannel(float *outL_1, float *outR_1, 
                                                       float *outL_2, float *outR_2, 
                                                       float *outL_3, float *outR_3, 
                                                       float *outL_4, float *outR_4, 
                                                       float *outL_5, float *outR_5, 
                                                       float *outL_6, float *outR_6, 
                                                       float *outL_7, float *outR_7, 
                                                       float *outL_8, float *outR_8)
{



}

#endif // SampleLauncher_h
