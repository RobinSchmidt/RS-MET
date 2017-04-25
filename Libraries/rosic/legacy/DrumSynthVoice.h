#ifndef DrumSynthVoice_h
#define DrumSynthVoice_h

#include "stdlib.h"
#include "math.h"
#include "BreakpointModulator.h"
#include "Definitions.h"
#include "MoreMath.h"
#include "SampleOscillator.h"


/**

This class realizes a single voice for the DrumSynth-Synthesizer.

*/

class DrumSynthVoice  
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

	DrumSynthVoice();
	~DrumSynthVoice();

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate      (double newSampleRate); // set the sample rate
 void resetDifferentiator();

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double *outL, double *outR);   
 /**< Calculates the output-samples for both channels and adds them at the
      adresses of *outL and *outR. */

 //---------------------------------------------------------------------------
 // event processing:

 void noteOn(int NoteNumber, int Velocity);

 //---------------------------------------------------------------------------
 // public member variables (to be accessed by the outlying DrumSynth 
 // class):

 // MIDI-data:
 int  currentNote, 
      currentNoteVel, 
      currentNoteDetune;
 int  currentNoteAge;  // age of the note in samples (for last note priority
                       // voice assignment)
 bool isPlaying;


 // embedded public audio-modules (they have to be accessed by the outlying
 // DrumSynth class):
 SampleOscillator    osc1, osc2;
 BreakpointModulator ampEnv;


 // parameters:
 // nominal source amplitudes (without velocity and key-scaling):
 double osc1Amplitude;
 double osc2Amplitude;

 // key- and vel-dependencies:
 double osc1TuneByKey;
 double osc2TuneByKey;
 double osc1LevelByVel;
 double osc2LevelByVel;

 // filter settings:
 bool   filterIsOn;
 double fltFreq;
 double fltFreqByKey;
 double fltFreqByVel;
 double fltFreqScaler;

 //===========================================================================

protected:

 // parameters:
 double sampleRate;        // sample-rate
 bool   stereo;
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

//----------------------------------------------------------------------------
// audio processing:

INLINE void DrumSynthVoice::getSampleFrameStereo(double *outL, double *outR)
{
 double tmpL, tmpR;
 double ampEnvOut;
 //double fltEnvOut;

 // generate the source-signal:
 osc1.getSampleFrameStereo(&tmpL, &tmpR);

 // read and apply the amplitude-envelope:
 ampEnvOut = ampEnv.getSample();
 tmpL *= ampEnvOut;
 tmpR *= ampEnvOut;

 // add the signal to the output-slots:
 *outL += tmpL;
 *outR += tmpR;

 // increment the note age:
 currentNoteAge++;
}

#endif //  DrumSynthVoice_h