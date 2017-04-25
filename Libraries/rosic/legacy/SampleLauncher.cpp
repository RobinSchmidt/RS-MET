#include "SampleLauncher.h"

//----------------------------------------------------------------------------
// construction/destruction:

SampleLauncher::SampleLauncher()
{
 // init member variables:
 playbackSampleRate   = 44100.0f;    
 numVoices            = 16;
 masterAmplitude      = 1.0f;
 mostRecentNote       = -1;
 mostRecentNoteVel    = 0;
 mostRecentNoteDetune = 0;

 int i;
 for(i=0; i<127; i++)
 {
  fileSampleRates[i] = 44100.0f;
  sampleLengths[i]   = 0;
 }

 // initialize the table contets to all zeros:
 zeroAllTables();
}

SampleLauncher::~SampleLauncher()
{
 
}

//----------------------------------------------------------------------------
// parameter settings:

void SampleLauncher::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0 )
  playbackSampleRate = newSampleRate;

 // resample the samples according to their file sample-rates and the current 
 // playback sample-rate:
 //........
}

void SampleLauncher::setNumVoices(int newNumVoices)
{
 if( newNumVoices>0 && newNumVoices<=MAXVOICES )
  numVoices = newNumVoices;
}

void SampleLauncher::clearSample(int whichSlot)
{
 //int i;
}

void SampleLauncher::setSample(int    whichSlot, 
                               float* newSampleL, 
                               float* newSampleR,
                               int    newSampleLength,
                               double newFileSampleRate)
{
 // make sure, that the passed table is not too long:
 if( newSampleLength > MAX_SAMPLE_LENGTH || newSampleLength < 1)
  return;

 //int i;
}




void SampleLauncher::setSlotDetune(int whichSlot, double newDetune)
{

}

void SampleLauncher::setSlotLevel(int whichSlot, double newLevel)
{

}

void SampleLauncher::setMasterVolume(double newMasterVolume)
{ 
 masterAmplitude = (float) MoreMath::dB2amp(newMasterVolume);
}


//----------------------------------------------------------------------------
// event processing:

void SampleLauncher::noteOn(long NoteNumber, long Velocity)
{
 /*
 mostRecentNote       = NoteNumber;
 mostRecentNoteVel    = Velocity;
 //mostRecentNoteDetune = Detune;

 //magicCarpetVoiceArray[NoteNumber].noteOn(mostRecentNote, mostRecentNoteVel);


 long i = 0;
 
 if(mostRecentNoteVel==0)   // a note-off event occured
 {
  // loop through the voices to find the one which is playing the note
  // for which the note-off was was sent:
  for(i=0; i<numVoices; i++)
  {
   if(magicCarpetVoiceArray[i].currentNote == mostRecentNote)
    magicCarpetVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
     // a note-on with velocity 0 is sent to the voice
  }
 }

 else                         // a note-on event occured
 {
  // look at all the voices, to see if ALL of them have reached their ends -
  // only in this case, the modulator-phases are reset:
  bool aVoiceIsPlaying  = false;
  for(i=0; i<numVoices; i++)
  {
   if( !magicCarpetVoiceArray[i].ampEnv.endIsReached() )
    aVoiceIsPlaying = true;
  }
  if( aVoiceIsPlaying == false )
  {
   xModulator.reset();
   yModulator.reset();
  }

  // loop through the voices to find a free one:
  int  oldestNoteAge    = 0;
  int  oldestVoice      = 0;       // voice with the oldest note
  for(i=0; i<numVoices; i++)
  {
   // check, if voice i is currently releasing the note, which is coming in,
   // if so, use that voice again:
   // check if the voice i is free:
   if(magicCarpetVoiceArray[i].currentNote == mostRecentNote)
   {
    // voice is used for the same note (which is currently releasing) again:
    magicCarpetVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
    return; //jump out of the function
   }
   // check if the voice i is free:
   else if(magicCarpetVoiceArray[i].ampEnv.endIsReached())
   {
    // voice i is free and can be used for the new note:
    magicCarpetVoiceArray[i].noteOn(mostRecentNote, mostRecentNoteVel);
    return; //jump out of the function
   }
   // keep track of the oldest note for the case, when no voice is
   // free for the new note (last note priority assignment):
   else if(magicCarpetVoiceArray[i].currentNoteAge > oldestNoteAge)
   {
    oldestNoteAge = magicCarpetVoiceArray[i].currentNoteAge;
    oldestVoice   = i;
   }
  } // end of for-loop


  // no free voice has been found - set the voice with the oldest note
  // to the new note:
  magicCarpetVoiceArray[oldestVoice].noteOn(mostRecentNote, mostRecentNoteVel);
 } // end of else
 */
}


//----------------------------------------------------------------------------
// internal functions:

void SampleLauncher::zeroAllTables()
{
 int i;
 for(i=0; i<MAX_SAMPLE_LENGTH+1; i++)
 {

 }
}