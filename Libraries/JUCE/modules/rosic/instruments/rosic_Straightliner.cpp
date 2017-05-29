//#include "rosic_Straightliner.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

Straightliner::Straightliner()
{
  // allocate the voices and assign the pointers to the IntrumentVoice-obejcts which are inherited
  // from the PolyphonicInstrument base-class:
  voiceArray    = new StraightlinerVoice[numAllocatedVoices];
  voicePointers = new PolyphonicInstrumentVoice*[numAllocatedVoices];
  int v;
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v] = &(voiceArray[v]);

  // switch the wavetables into sum/difference (mid/side) mode
  //waveTableForOsc1.setChannelSumAndDifferenceMode(true);
  //waveTableForOsc2.setChannelSumAndDifferenceMode(true);
  //waveTableForOsc3.setChannelSumAndDifferenceMode(true);
  //waveTableForOsc4.setChannelSumAndDifferenceMode(true);

  // tell the oscillators of all the voices, where the wavetables are:
  //masterVoice[0].osc1.setWaveTableToUse(&waveTableForOsc1);
  for(v=0; v<numAllocatedVoices; v++)
  {
    voiceArray[v].oscSection.osc1.setWaveTableToUse(&waveTableForOsc1);
    voiceArray[v].oscSection.osc2.setWaveTableToUse(&waveTableForOsc2);
    voiceArray[v].oscSection.osc3.setWaveTableToUse(&waveTableForOsc3);
    voiceArray[v].oscSection.osc4.setWaveTableToUse(&waveTableForOsc4);
  }

  // configure the embedded objects of all voices other than [0] as slaves to the respective 
  // objects in voice [0]:
  for(v=1; v<numAllocatedVoices; v++)
  {
    voiceArray[0].oscSection.osc1.addSlave(      &voiceArray[v].oscSection.osc1      );
    voiceArray[0].oscSection.osc2.addSlave(      &voiceArray[v].oscSection.osc2      );
    voiceArray[0].oscSection.osc3.addSlave(      &voiceArray[v].oscSection.osc3      );
    voiceArray[0].oscSection.osc4.addSlave(      &voiceArray[v].oscSection.osc4      );
    voiceArray[0].filter.addSlave(    &voiceArray[v].filter    );
    voiceArray[0].pitchEnv.addSlave(  &voiceArray[v].pitchEnv  );
    voiceArray[0].filterEnv.addSlave( &voiceArray[v].filterEnv );
    voiceArray[0].ampEnv.addSlave(    &voiceArray[v].ampEnv    );
  }
    
  // setup the anti-alias filters:
  antiAliasFilterL.setSubDivision(2.0);
  antiAliasFilterR.setSubDivision(2.0);

  // setup the dc-blocker:
  dcBlocker.setMode(OnePoleFilterStereo::HIGHPASS);
  dcBlocker.setCutoff(10.0);

  // assign the inherited TuningTable (from PolyphonicInstrument) to all the voices:
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v]->setTuningTable(&tuningTable);

  setSampleRate(44100.0);
}

Straightliner::~Straightliner()
{
  if( voiceArray != NULL )
    delete[] voiceArray;
  if( voicePointers != NULL )
    delete[] voicePointers;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void Straightliner::setSampleRate(double newSampleRate)
{
  PolyphonicInstrument::setSampleRate(newSampleRate);
  dcBlocker.setSampleRate(newSampleRate);
}
