#include "rosic_MagicCarpet.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MagicCarpet::MagicCarpet()
{
  // allocate the voices and assign the pointers to the IntrumentVoice-obejcts which are inherited
  // from the PolyphonicInstrument base-class:
  voiceArray    = new MagicCarpetVoice[numAllocatedVoices];
  voicePointers = new PolyphonicInstrumentVoice*[numAllocatedVoices];
  int v;
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v] = &(voiceArray[v]);

  // tell the LFOs and oscillators of all the voices, where the wavetables are:
  for(v=0; v<numAllocatedVoices; v++)
  {
    voiceArray[v].oscSection.xLfo.setWaveTableToUse(&xLfoWaveTable);
    voiceArray[v].oscSection.yLfo.setWaveTableToUse(&yLfoWaveTable);
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
    voiceArray[0].filter.addSlave(         &voiceArray[v].filter                            );
    voiceArray[0].filterEnv.addSlave(      &voiceArray[v].filterEnv                         );
    voiceArray[0].ampEnv.addSlave(         &voiceArray[v].ampEnv                            );
  }

  // assign the inherited TuningTable (from PolyphonicInstrument) to all the voices:
  for(v=0; v<numAllocatedVoices; v++)
    voicePointers[v]->setTuningTable(&tuningTable);

  setSampleRate(44100.0);
}

MagicCarpet::~MagicCarpet()
{
  /*
  // NULL the shared memory pointers int the voices (so their detructors don't try to free memory
  // that they do not own):
  for(v=0; v<numAllocatedVoices; v++)
  {
    voiceArray[v].oscSection.xLfo.setWaveTableToUse(&xLfoWaveTable);
    voiceArray[v].oscSection.yLfo.setWaveTableToUse(&yLfoWaveTable);
    voiceArray[v].oscSection.samplePlayerTopLeft.setSampleBufferToUse(&sampleBufferTopLeft);
    voiceArray[v].oscSection.samplePlayerTopRight.setSampleBufferToUse(&sampleBufferTopRight);
    voiceArray[v].oscSection.samplePlayerBottomLeft.setSampleBufferToUse(&sampleBufferBottomLeft);
    voiceArray[v].oscSection.samplePlayerBottomRight.setSampleBufferToUse(&sampleBufferBottomRight);
  }
  */


  if( voiceArray != NULL )
    delete[] voiceArray;
  if( voicePointers != NULL )
    delete[] voicePointers;

  // delete all the AutomatableParameter objects that have been created:
  //....nope...this is already done by the destructor of base-class AutomatableModule

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void MagicCarpet::setSampleRate(double newSampleRate)
{
  PolyphonicInstrument::setSampleRate(newSampleRate);
  delayPhaser.setSampleRate(newSampleRate);
}

void MagicCarpet::setBeatsPerMinute(double newBpm)
{
  PolyphonicInstrument::setBeatsPerMinute(newBpm);
  delayPhaser.setTempoInBPM(newBpm);
}

//-------------------------------------------------------------------------------------------------
// event processing:


