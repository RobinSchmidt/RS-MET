//#include "rosic_Workhorse.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Workhorse::Workhorse()
{
  // allocate the voices and assign the pointers to the IntrumentVoice-obejcts which are inherited
  // from the PolyphonicInstrument base-class:
  voiceArray    = new WorkhorseVoice[numAllocatedVoices];
  voicePointers = new PolyphonicInstrumentVoice*[numAllocatedVoices];
  int v;
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v] = &(voiceArray[v]);

  // switch the wavetables into sum/difference (mid/side) mode
  //sampleBuffer1.setChannelSumAndDifferenceMode(true);


  // tell the oscillators of all the voices, where the wavetables are:
  //masterVoice[0].osc1.setWaveTableToUse(&waveTableForOsc1);
  for(v=0; v<numAllocatedVoices; v++)
  {
    voiceArray[v].oscSection.samplePlayerTopLeft.setSampleBufferToUse(&sampleBufferTopLeft);
    voiceArray[v].oscSection.samplePlayerTopRight.setSampleBufferToUse(&sampleBufferTopRight);
    voiceArray[v].oscSection.samplePlayerBottomLeft.setSampleBufferToUse(&sampleBufferBottomLeft);
    voiceArray[v].oscSection.samplePlayerBottomRight.setSampleBufferToUse(&sampleBufferBottomRight);
  }

  // configure the embedded objects of all voices other than [0] as slaves to the respective
  // objects in voice [0]:
  for(v=1; v<numAllocatedVoices; v++)
  {
    voiceArray[0].oscSection.samplePlayerTopLeft.addSlave(&voiceArray[v].oscSection.samplePlayerTopLeft );
    voiceArray[0].oscSection.samplePlayerTopRight.addSlave(&voiceArray[v].oscSection.samplePlayerTopRight );
    voiceArray[0].oscSection.samplePlayerBottomLeft.addSlave(&voiceArray[v].oscSection.samplePlayerBottomLeft );
    voiceArray[0].oscSection.samplePlayerBottomRight.addSlave(&voiceArray[v].oscSection.samplePlayerBottomRight );
    voiceArray[0].pitchEnv.addSlave(       &voiceArray[v].pitchEnv                          );
    voiceArray[0].filter.addSlave(         &voiceArray[v].filter                            );
    voiceArray[0].filterEnv.addSlave(      &voiceArray[v].filterEnv                         );
    voiceArray[0].ampEnv.addSlave(         &voiceArray[v].ampEnv                            );
  }

  // add the 0th voice as child preset-rememberer:
  //addChildRememberer( &(voiceArray[0]) );

  // setup the anti-alias filters:
  antiAliasFilterL.setSubDivision(2.0);
  antiAliasFilterR.setSubDivision(2.0);

  // assign the inherited TuningTable (from PolyphonicInstrument) to all the voices:
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v]->setTuningTable(&tuningTable);

  setSampleRate(44100.0);
}

Workhorse::~Workhorse()
{
  if( voiceArray != NULL )
    delete[] voiceArray;
  if( voicePointers != NULL )
    delete[] voicePointers;

  // delete all the AutomatableParameter objects that have been created:
  //....nope...this is already done by the destructor of base-class AutomatableModule

}

//-------------------------------------------------------------------------------------------------
// parameter settings:


//-------------------------------------------------------------------------------------------------
// event processing:


