#include "rosof_VectorSamplePlayerAudioModule.h"
using namespace rosof;

//=================================================================================================
// class VectorSamplePlayerAudioModule

//-------------------------------------------------------------------------------------------------
// construction/destruction:

VectorSamplePlayerAudioModule::VectorSamplePlayerAudioModule(CriticalSection *newPlugInLock, 
                                                             rosic::VectorSamplePlayer *vectorSamplePlayerToWrap)                                                             
                                                             : AudioModule(newPlugInLock)
{
  jassert(vectorSamplePlayerToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedVectorSamplePlayer = vectorSamplePlayerToWrap;
  moduleName                = juce::String(T("VectorSamplePlayer"));

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/VectorSamplePlayerPresets")) );

  vectorMixerModule = new VectorMixerAudioModule(plugInLock, &wrappedVectorSamplePlayer->vectorMixer);
  vectorMixerModule->setModuleName(juce::String(T("VectorMixer")));
  addChildAudioModule(vectorMixerModule);

  samplePlayerTopLeftModule = new SamplePlayerAudioModule(plugInLock, &wrappedVectorSamplePlayer->samplePlayerTopLeft);
  samplePlayerTopLeftModule->setModuleName(juce::String(T("SamplePlayerTopLeft")));
  addChildAudioModule(samplePlayerTopLeftModule);

  samplePlayerTopRightModule = new SamplePlayerAudioModule(plugInLock, &wrappedVectorSamplePlayer->samplePlayerTopRight);
  samplePlayerTopRightModule->setModuleName(juce::String(T("SamplePlayerTopRight")));
  addChildAudioModule(samplePlayerTopRightModule);

  samplePlayerBottomLeftModule = new SamplePlayerAudioModule(plugInLock, &wrappedVectorSamplePlayer->samplePlayerBottomLeft);
  samplePlayerBottomLeftModule->setModuleName(juce::String(T("SamplePlayerBottomLeft")));
  addChildAudioModule(samplePlayerBottomLeftModule);

  samplePlayerBottomRightModule = new SamplePlayerAudioModule(plugInLock, &wrappedVectorSamplePlayer->samplePlayerBottomRight);
  samplePlayerBottomRightModule->setModuleName(juce::String(T("SamplePlayerBottomRight")));
  addChildAudioModule(samplePlayerBottomRightModule);

  xLfoModule = new LowFrequencyOscillatorAudioModule(plugInLock, &wrappedVectorSamplePlayer->xLfo);
  xLfoModule->setModuleName(juce::String(T("X-LFO")));
  addChildAudioModule(xLfoModule);

  yLfoModule = new LowFrequencyOscillatorAudioModule(plugInLock, &wrappedVectorSamplePlayer->yLfo);
  yLfoModule->setModuleName(juce::String(T("Y-LFO")));
  addChildAudioModule(yLfoModule);
}

