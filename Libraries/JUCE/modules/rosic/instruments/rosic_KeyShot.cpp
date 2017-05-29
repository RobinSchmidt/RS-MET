//#include "rosic_KeyShot.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

KeyShot::KeyShot()
{
  // allocate the voices and assign the pointers to the IntrumentVoice-obejcts which are inherited
  // from the PolyphonicInstrument base-class:
  voiceArray    = new KeyShotVoice[numAllocatedVoices];
  voicePointers = new PolyphonicInstrumentVoice*[numAllocatedVoices];
  int v;
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v] = &(voiceArray[v]);

  // tell the oscillators of all the voices, where the wavetables are:
  for(v=0; v<numAllocatedVoices; v++)
  {
    voiceArray[v].samplePlayer.setSampleBufferToUse(&sampleBuffer);
      // to be revised
  }

  // configure the embedded objects of all voices other than [0] as slaves to the respective
  // objects in voice [0]:
  for(v=1; v<numAllocatedVoices; v++)
  {
    voiceArray[0].samplePlayer.addSlave(&voiceArray[v].samplePlayer );
    voiceArray[0].ampEnv.addSlave(      &voiceArray[v].ampEnv       );
  }

  // assign the inherited TuningTable (from PolyphonicInstrument) to all the voices:
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v]->setTuningTable(&tuningTable);

  setSampleRate(44100.0);
}

KeyShot::~KeyShot()
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


