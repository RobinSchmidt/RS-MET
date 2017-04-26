#include "rosic_StraightlinerVoice.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

StraightlinerVoice::StraightlinerVoice()
{
  ampEnv.fixEndLevelAtZero(true);
  cutoffSmoother.setTimeConstant(30.0);
}

StraightlinerVoice::~StraightlinerVoice()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings: 

void StraightlinerVoice::setSampleRate(double newSampleRate)
{
  PolyphonicInstrumentVoice::setSampleRate(newSampleRate);
  oscSection.setSampleRate(2.0*sampleRate);
  filter.setSampleRate(2.0*sampleRate);

  pitchEnv.setSampleRate(sampleRate);
  filterEnv.setSampleRate(sampleRate);
  ampEnv.setSampleRate(sampleRate);

  cutoffSmoother.setSampleRate(sampleRate);
}

void StraightlinerVoice::setMasterTuneA4(double newMasterTuneA4)
{
  PolyphonicInstrumentVoice::setMasterTuneA4(newMasterTuneA4);

  // update the oscillators frequencies:
  double noteFreq = 440.0;
  if( isReleasing )
    noteFreq = getNoteFrequency(noteBeingReleased);
  else if( getCurrentNoteKey() != -1 )
    noteFreq = getNoteFrequency(getCurrentNoteKey());

  oscSection.osc1.setFrequencyNominal(noteFreq);
  oscSection.osc2.setFrequencyNominal(noteFreq);
  oscSection.osc3.setFrequencyNominal(noteFreq);
  oscSection.osc4.setFrequencyNominal(noteFreq);
}

void StraightlinerVoice::setBeatsPerMinute(double newBeatsPerMinute)
{
  pitchEnv.setBeatsPerMinute(newBeatsPerMinute);
  filterEnv.setBeatsPerMinute(newBeatsPerMinute);
  ampEnv.setBeatsPerMinute(newBeatsPerMinute);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// event processing:

void StraightlinerVoice::triggerNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::triggerNote(newKey, newVelocity, newDetune);

  // set up the oscillators frequencies:
  oscSection.osc1.setFrequencyNominal(currentFrequency);
  oscSection.osc2.setFrequencyNominal(currentFrequency);
  oscSection.osc3.setFrequencyNominal(currentFrequency);
  oscSection.osc4.setFrequencyNominal(currentFrequency);


  MidiNoteEvent currentNoteEvent;
  mutex.lock();
  if( !noteList.empty() )
    currentNoteEvent = noteList.front();
  mutex.unlock();

  // tell the filter and oscillators the current key and velocity, so that they can set up their
  // key and velocity dependent parameters:
  filter.setKeyAndVel((double) newKey, (double) newVelocity);
  cutoffSmoother.setState(filter.getFrequencyWithKeyAndVel()); // avoid smoothing on note-on
  oscSection.osc1.setKeyAndVel(newKey, newVelocity);
  oscSection.osc2.setKeyAndVel(newKey, newVelocity);
  oscSection.osc3.setKeyAndVel(newKey, newVelocity);
  oscSection.osc4.setKeyAndVel(newKey, newVelocity);

  // trigger envelopes and oscillators:
  if( ampEnv.endIsReached )
  {
    pitchEnv.noteOn( false, newKey, newVelocity);
    filterEnv.noteOn(false, newKey, newVelocity);
    filter.resetBuffers();
    oscSection.osc1.reset();
    oscSection.osc2.reset();
    oscSection.osc3.reset();
    oscSection.osc4.reset();
    ampEnv.noteOn(false, newKey, newVelocity);
  }
  else
  {
    pitchEnv.noteOn(false, newKey, newVelocity);
    filterEnv.noteOn(false, newKey, newVelocity);
    ampEnv.noteOn(true, newKey, newVelocity);
  }

  // this voices has just been triggered with a note-on, so it's not silent:
  //isSilent = false;  //moved to baseclass method
}

void StraightlinerVoice::glideToNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::glideToNote(newKey, newVelocity, newDetune);
  filter.setKeyAndVel((double) newKey, (double) newVelocity);
  oscSection.osc1.setKeyAndVel(newKey, newVelocity);
  oscSection.osc2.setKeyAndVel(newKey, newVelocity);
  oscSection.osc3.setKeyAndVel(newKey, newVelocity);
  oscSection.osc4.setKeyAndVel(newKey, newVelocity);
}

void StraightlinerVoice::triggerRelease(int noteToBeReleased, int noteToBeReleasedVel)
{
  PolyphonicInstrumentVoice::triggerRelease(noteToBeReleased, noteToBeReleasedVel);

  // trigger the releases of the envlopes:
  ampEnv.noteOff(true);
  filterEnv.noteOff(true);
  pitchEnv.noteOff(true);
}
