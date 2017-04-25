#include "rosic_KeyShotVoice.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

KeyShotVoice::KeyShotVoice()
{
  ampEnv.fixEndLevelAtZero(true);
  for(int b=0; b<ampEnv.lastBreakpointIndex(); b++)
    ampEnv.setBreakpointLevel(b, 1.0);
}

KeyShotVoice::~KeyShotVoice()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings: 

void KeyShotVoice::setSampleRate(double newSampleRate)
{
  PolyphonicInstrumentVoice::setSampleRate(newSampleRate);
  samplePlayer.setSampleRate(2.0*sampleRate);
  ampEnv.setSampleRate(sampleRate);
}

void KeyShotVoice::setMasterTuneA4(double newMasterTuneA4)
{
  PolyphonicInstrumentVoice::setMasterTuneA4(newMasterTuneA4);

  // update the oscillators frequencies:
  double noteFreq = 440.0;
  if( isReleasing )
    noteFreq = getNoteFrequency(noteBeingReleased);
  else if( getCurrentNoteKey() != -1 )
    noteFreq = getNoteFrequency(getCurrentNoteKey());

  samplePlayer.setPlaybackFrequencyNominal(noteFreq);
}

void KeyShotVoice::setBeatsPerMinute(double newBeatsPerMinute)
{
  ampEnv.setBeatsPerMinute(newBeatsPerMinute);
}

//-------------------------------------------------------------------------------------------------
// event processing:

void KeyShotVoice::triggerNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::triggerNote(newKey, newVelocity, newDetune);

  // set up the oscillators frequencies and trigger:
  samplePlayer.setPlaybackFrequencyNominal(currentFrequency);
  samplePlayer.reset();

  MidiNoteEvent currentNoteEvent;
  mutex.lock();
  if( !noteList.empty() )
    currentNoteEvent = noteList.front();
  mutex.unlock();

  ampEnv.noteOn(true, newKey, newVelocity);

  // this voices has just been triggered with a note-on, so it's not silent:
  isSilent = false;
}

void KeyShotVoice::triggerRelease(int noteToBeReleased, int noteToBeReleasedVel)
{
  PolyphonicInstrumentVoice::triggerRelease(noteToBeReleased, noteToBeReleasedVel);
  ampEnv.noteOff(true);
}
