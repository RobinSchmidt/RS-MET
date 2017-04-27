#include "rosic_QuadrigaVoice.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

QuadrigaVoice::QuadrigaVoice()
{

}

QuadrigaVoice::~QuadrigaVoice()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void QuadrigaVoice::setSampleRate(double newSampleRate)
{
  PolyphonicInstrumentVoice::setSampleRate(newSampleRate);
  quadrigen.setSampleRate(sampleRate);
  ampEnv.setSampleRate(sampleRate);
}

void QuadrigaVoice::setMasterTuneA4(double newMasterTuneA4)
{
  PolyphonicInstrumentVoice::setMasterTuneA4(newMasterTuneA4);

  // update the oscillators frequencies:
  double noteFreq = 440.0;
  if( isReleasing )
    noteFreq = getNoteFrequency(noteBeingReleased);
  else if( getCurrentNoteKey() != -1 )
    noteFreq = getNoteFrequency(getCurrentNoteKey());
  //oscSection.setPlaybackFrequencyNominal(noteFreq);
}

void QuadrigaVoice::setBeatsPerMinute(double newBeatsPerMinute)
{
  quadrigen.setTempoInBPM(newBeatsPerMinute);
  ampEnv.setBeatsPerMinute(newBeatsPerMinute);
}

//-------------------------------------------------------------------------------------------------
// event processing:

void QuadrigaVoice::triggerNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::triggerNote(newKey, newVelocity, newDetune);

  // set up the oscillators frequencies and trigger:
  //oscSection.setPlaybackFrequencyNominal(currentFrequency);
  //oscSection.reset();
  quadrigen.reset();

  MidiNoteEvent currentNoteEvent;
  mutex.lock();
  if( !noteList.empty() )
    currentNoteEvent = noteList.front();
  mutex.unlock();

  // tell the filter the current key and velocity, so it can set up the key and velocity
  // scaled characteristic frequency:
  //oscSection.setKeyAndVel(newKey, newVelocity);

  // trigger envelopes:
  if( ampEnv.endIsReached )
  {
    quadrigen.reset();
    //quadrigen.noteOn(false, newKey, newVelocity);
    ampEnv.noteOn(false, newKey, newVelocity);
  }
  else
  {
    //quadrigen.noteOn(true, newKey, newVelocity);
    ampEnv.noteOn(true, newKey, newVelocity);
  }

  // this voices has just been triggered with a note-on, so it's not silent:
  isSilent = false;
}

void QuadrigaVoice::glideToNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::glideToNote(newKey, newVelocity, newDetune);
  //filter.setKeyAndVel((double) newKey, (double) newVelocity); // to be re-activated (?)
}

void QuadrigaVoice::triggerRelease(int noteToBeReleased, int noteToBeReleasedVel)
{
  PolyphonicInstrumentVoice::triggerRelease(noteToBeReleased, noteToBeReleasedVel);

  // trigger the releases of the envlopes:
  // quadrigen.noteOff();
  ampEnv.noteOff();
}
