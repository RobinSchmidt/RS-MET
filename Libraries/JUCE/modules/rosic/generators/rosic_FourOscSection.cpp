//#include "rosic_FourOscSection.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FourOscSection::FourOscSection()
{
  /*
  addChildRememberer(&osc1);
  addChildRememberer(&osc2);
  addChildRememberer(&osc3);
  addChildRememberer(&osc4);
  */
}

FourOscSection::~FourOscSection()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings: 

void FourOscSection::setSampleRate(double newSampleRate)
{
  osc1.setSampleRate(newSampleRate);
  osc2.setSampleRate(newSampleRate);
  osc3.setSampleRate(newSampleRate);
  osc4.setSampleRate(newSampleRate);
}

/*
void FourOscSection::setMasterTuneA4(double newMasterTuneA4)
{
  PolyphonicInstrumentVoice::setMasterTuneA4(newMasterTuneA4);

  // update the oscillators frequencies:
  double noteFreq = 440.0;
  if( isReleasing )
    noteFreq = pitchToFreq((double) noteBeingReleased, masterTuneA4);
  else if( getCurrentNoteKey() != -1 )
    noteFreq = pitchToFreq((double) getCurrentNoteKey(), masterTuneA4);

  //double noteFreq = pitchToFreq((double) getCurrentNoteKey(), masterTuneA4);
  osc1.setFrequencyNominal(noteFreq);
  osc2.setFrequencyNominal(noteFreq);
  osc3.setFrequencyNominal(noteFreq);
  osc4.setFrequencyNominal(noteFreq);
}
*/