#include "rosic_WorkhorseVoice.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WorkhorseVoice::WorkhorseVoice()
{
  ampEnv.fixEndLevelAtZero(true);

  // set the filter and envelopes neutral:
  filter.setMode(MultiModeFilterParameters::BYPASS);
  int b;
  for(b=0; b<=pitchEnv.lastBreakpointIndex(); b++)
    pitchEnv.setBreakpointLevel(b, 1.0);
  for(b=0; b<ampEnv.lastBreakpointIndex(); b++)
    ampEnv.setBreakpointLevel(b, 1.0);
  for(b=0; b<=filterEnv.lastBreakpointIndex(); b++)
    filterEnv.setBreakpointLevel(b, 1.0);

  // add the embedded audio modules as child PresetRememberers:
  //addChildRememberer(&oscSection);
  //addChildRememberer(&filter);
}

WorkhorseVoice::~WorkhorseVoice()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings: 

void WorkhorseVoice::setSampleRate(double newSampleRate)
{
  PolyphonicInstrumentVoice::setSampleRate(newSampleRate);
  oscSection.setSampleRate(2.0*sampleRate);
  filter.setSampleRate(2.0*sampleRate);
  pitchEnv.setSampleRate(sampleRate);
  filterEnv.setSampleRate(sampleRate);
  ampEnv.setSampleRate(sampleRate);
}

void WorkhorseVoice::setMasterTuneA4(double newMasterTuneA4)
{
  PolyphonicInstrumentVoice::setMasterTuneA4(newMasterTuneA4);

  // update the oscillators frequencies:
  double noteFreq = 440.0;
  if( isReleasing )
    noteFreq = getNoteFrequency(noteBeingReleased);
  else if( getCurrentNoteKey() != -1 )
    noteFreq = getNoteFrequency(getCurrentNoteKey());

  oscSection.setPlaybackFrequencyNominal(noteFreq);
}

void WorkhorseVoice::setBeatsPerMinute(double newBeatsPerMinute)
{
  pitchEnv.setBeatsPerMinute(newBeatsPerMinute);
  filterEnv.setBeatsPerMinute(newBeatsPerMinute);
  ampEnv.setBeatsPerMinute(newBeatsPerMinute);

}

//-------------------------------------------------------------------------------------------------
// event processing:

void WorkhorseVoice::triggerNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::triggerNote(newKey, newVelocity, newDetune);

  // set up the oscillators frequencies and trigger:
  oscSection.setPlaybackFrequencyNominal(currentFrequency);
  oscSection.reset();

  MidiNoteEvent currentNoteEvent;
  mutex.lock();
  if( !noteList.empty() )
    currentNoteEvent = noteList.front();
  mutex.unlock();

  // tell the filter the current key and velocity, so it can set up the key and velocity 
  // scaled characteristic frequency:
  filter.setKeyAndVel((double) newKey, (double) newVelocity);

  // trigger envelopes:
  if( ampEnv.endIsReached )
  {
    pitchEnv.noteOn(false, newKey, newVelocity);
    filterEnv.noteOn(false, newKey, newVelocity);
    filter.resetBuffers();
    ampEnv.noteOn(false, newKey, newVelocity);
  }
  else
  {
    pitchEnv.noteOn(false, newKey, newVelocity);
    filterEnv.noteOn(false, newKey, newVelocity);
    ampEnv.noteOn(true, newKey, newVelocity);
  }

  // this voices has just been triggered with a note-on, so it's not silent:
  isSilent = false;
}

void WorkhorseVoice::glideToNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::glideToNote(newKey, newVelocity, newDetune);
  filter.setKeyAndVel((double) newKey, (double) newVelocity);
}

void WorkhorseVoice::triggerRelease(int noteToBeReleased, int noteToBeReleasedVel)
{
  PolyphonicInstrumentVoice::triggerRelease(noteToBeReleased, noteToBeReleasedVel);

  // trigger the releases of the envlopes:
  pitchEnv.noteOff(true);
  ampEnv.noteOff(true);
  filterEnv.noteOff(true);
}
