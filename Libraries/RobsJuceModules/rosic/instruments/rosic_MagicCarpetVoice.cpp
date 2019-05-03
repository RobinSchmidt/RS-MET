//#include "rosic_MagicCarpetVoice.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MagicCarpetVoice::MagicCarpetVoice()
{
  filter.setMode(FourPoleFilterParameters::BYPASS);
}

MagicCarpetVoice::~MagicCarpetVoice()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings: 

void MagicCarpetVoice::setSampleRate(double newSampleRate)
{
  PolyphonicInstrumentVoice::setSampleRate(newSampleRate);
  oscSection.setSampleRate(sampleRate);
  filter.setSampleRate(sampleRate);
  filterEnv.setSampleRate(sampleRate);
  ampEnv.setSampleRate(sampleRate);
}

void MagicCarpetVoice::setMasterTuneA4(double newMasterTuneA4)
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

void MagicCarpetVoice::setBeatsPerMinute(double newBeatsPerMinute)
{
  oscSection.setBeatsPerMinute(newBeatsPerMinute);
  //pitchEnv.setBeatsPerMinute(newBeatsPerMinute);
  //filterEnv.setBeatsPerMinute(newBeatsPerMinute);
  //ampEnv.setBeatsPerMinute(newBeatsPerMinute);
}

//-------------------------------------------------------------------------------------------------
// event processing:

void MagicCarpetVoice::triggerNote(int newKey, int newVelocity, int newDetune)
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
  oscSection.setKeyAndVel(newKey, newVelocity);
  //filter.setKeyAndVel((double) newKey, (double) newVelocity);  
    // to be re-activated or handle these dependencies here?

  // trigger envelopes:
  if( ampEnv.endIsReached )
  {
    filter.reset();
    filterEnv.noteOn(false, newKey, newVelocity);
    ampEnv.noteOn(false, newKey, newVelocity);
  }
  else
  {
    filterEnv.noteOn(false, newKey, newVelocity);
    ampEnv.noteOn(true, newKey, newVelocity);
  }

  // this voices has just been triggered with a note-on, so it's not silent:
  isSilent = false;
}

void MagicCarpetVoice::glideToNote(int newKey, int newVelocity, int newDetune)
{
  PolyphonicInstrumentVoice::glideToNote(newKey, newVelocity, newDetune);
  //filter.setKeyAndVel((double) newKey, (double) newVelocity); // to be re-activated (?)
}

void MagicCarpetVoice::triggerRelease(int noteToBeReleased, int noteToBeReleasedVel)
{
  PolyphonicInstrumentVoice::triggerRelease(noteToBeReleased, noteToBeReleasedVel);

  // trigger the releases of the envlopes:
  ampEnv.noteOff();
  filterEnv.noteOff();
}
