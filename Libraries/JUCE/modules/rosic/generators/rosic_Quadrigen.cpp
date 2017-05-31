//#include "rosic_Quadrigen.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Quadrigen::Quadrigen() : mixMatrix(6, 3)
{
  mutex.lock();

  sampleRate  = 44100.0;
  bpm         = 120.0;
  for(int i=0; i<numGeneratorSlots; i++)
  {
    generatorAlgorithmIndices[i] = MUTE;
    generatorModules[i]          = new MuteModule();
      // don't change this - the constructor of QuadrigenAudioModule relies on them to be of type
      // MuteModule in the beginning (using static_cast) -  TODO....!!!!
  }
  reset();

  mutex.unlock();
}

Quadrigen::~Quadrigen()
{
  mutex.lock();
  for(int i=0; i<numGeneratorSlots; i++)
    delete generatorModules[i];
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Quadrigen::setSampleRate(double newSampleRate)
{
  if( newSampleRate <= 0.0 )
  {
    DEBUG_BREAK;
    return;
  }
  mutex.lock();
  sampleRate = newSampleRate;
  for(int i=0; i<numGeneratorSlots; i++)
    generatorModules[i]->setSampleRate(sampleRate);
  mutex.unlock();
}

void Quadrigen::setTempoInBPM(double newTempoInBPM)
{
  if( newTempoInBPM != bpm )
  {
    bpm = newTempoInBPM;
    mutex.lock();
    for(int i=0; i<numGeneratorSlots; i++)
      generatorModules[i]->setTempoInBPM(bpm);
    mutex.unlock();
  }
}

void Quadrigen::setGeneratorAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  if( slotIndex < 0 || slotIndex >= numGeneratorSlots )
    return;

  mutex.lock();

  // delete the old GeneratorModule:
  delete generatorModules[slotIndex];

  // create and set up the new GeneratorModule for the slot:
  generatorAlgorithmIndices[slotIndex] = newAlgorithmIndex;
  Module *newModule = NULL;
  switch( generatorAlgorithmIndices[slotIndex] )
  {
  case MUTE:                 newModule = new MuteModule;                      break;
  case OSCILLATOR_STEREO:    newModule = new OscillatorStereoModule;          break;

    //......

  default:
    {
      newModule                            = new MuteModule;
      generatorAlgorithmIndices[slotIndex] = MUTE;
    }
  }
  generatorModules[slotIndex] = newModule;
  generatorModules[slotIndex]->setSampleRate(sampleRate);
  generatorModules[slotIndex]->setTempoInBPM(bpm);
  generatorModules[slotIndex]->reset();

  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int Quadrigen::getGeneratorAlgorithmIndex(int slotIndex) const
{
  if( slotIndex < 0 || slotIndex >= numGeneratorSlots )
    return MUTE;
  else
    return generatorAlgorithmIndices[slotIndex];
}

//-------------------------------------------------------------------------------------------------
// others:

void Quadrigen::trigger()
{
  mutex.lock();
  for(int i=0; i<numGeneratorSlots; i++)
    generatorModules[i]->trigger();
  mutex.unlock();
}

void Quadrigen::reset()
{
  for(int i=0; i<5; i++)
  {
    inputsL[i]  = 0.0;
    inputsR[i]  = 0.0;
    outputsL[i] = 0.0;
    outputsR[i] = 0.0;
  }
  mutex.lock();
  for(int i=0; i<numGeneratorSlots; i++)
    generatorModules[i]->reset();
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// audio processing:

void Quadrigen::processBlock(float *inOutL, float *inOutR, int numFrames)
{
  mutex.lock();

  double tmp1L = 0.0;
  double tmp1R = 0.0;
  double tmp2L = 0.0;
  double tmp2R = 0.0;
  double tmp3L = 0.0;
  double tmp3R = 0.0;
  double tmp4L = 0.0;
  double tmp4R = 0.0;
  //double inL   = inOutL[0];
  //double inR   = inOutR[0];

  for(int n=0; n<numFrames; n++)
  {
    generatorModules[0]->processSampleFrame(&tmp1L, &tmp1R);
    generatorModules[1]->processSampleFrame(&tmp1L, &tmp1R);
    generatorModules[2]->processSampleFrame(&tmp1L, &tmp1R);
    generatorModules[3]->processSampleFrame(&tmp1L, &tmp1R);

    // todo: filtering and (possibly) voice-mixing, modulation

    inOutL[n] = (float) (tmp1L+tmp2L+tmp3L+tmp4L);
    inOutR[n] = (float) (tmp1R+tmp2R+tmp3R+tmp4R);
  }

  mutex.unlock();
}






