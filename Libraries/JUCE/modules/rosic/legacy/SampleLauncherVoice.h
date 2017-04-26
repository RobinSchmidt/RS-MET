#ifndef SampleLauncherVoice_h
#define SampleLauncherVoice_h

#include "stdlib.h"          // for the NULL macro
#include "Definitions.h"
#include "SampleLauncherDefinitions.h"

/**

This class ...

*/

class SampleLauncherVoice
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 SampleLauncherVoice();
 ~SampleLauncherVoice();

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleAddresses(float* newSampleAddressL, 
                         float* newSampleAddressR = NULL);
 /**< Sets the adresses of the sample to be played. Make very sure that
      you call this function BEFORE the first call to getSample() - otherwise
      the getSample()-function will have to deal with an invalid pointer. When
      no second pointer or a NULL-pointer is passed for the right channel, it
      will be set to the same address to which the left-channel pointer points
      to, such that both channels play the same signal (mono). */

 void setSampleLength(int newSampleLength);
 ///< Sets the size of the table in samples.

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(float* outL, float* outR);
 /**< Calculates one stereo sample-frame at a time. */

 //---------------------------------------------------------------------------
 // others:
 void reset();
  // resets the table-pointer to its start-position and resets the filters.

protected:

 int position;         // cuurent position in the sample
 int sampleLength;     // size of the sample

 float* sampleL;       // pointers to the samples
 float* sampleR;
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

INLINE void SampleLauncherVoice::getSampleFrameStereo(float* outL, 
                                                      float* outR)
{



}

#endif // SampleLauncherVoice_h
