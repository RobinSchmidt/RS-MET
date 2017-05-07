#include "rosof_WorkhorseAudioModule.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WorkhorseAudioModule::WorkhorseAudioModule(CriticalSection *newPlugInLock, rosic::Workhorse *workhorseToWrap)
: PolyphonicInstrumentAudioModule(newPlugInLock, workhorseToWrap)
{
  jassert(workhorseToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedWorkhorse        = workhorseToWrap;
  underlyingRosicInstrument = workhorseToWrap;
  moduleName                = juce::String(T("Workhorse"));

  // initialize the current directory for preset loading and saving:
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/WorkhorsePresets")) );

  vectorMixerModule = new VectorMixerAudioModule(plugInLock, &wrappedWorkhorse->vectorMixer);
  vectorMixerModule->setModuleName(juce::String(T("VectorMixer")));
  addChildAudioModule(vectorMixerModule);

  samplePlayerTopLeftModule = new SamplePlayerAudioModule(plugInLock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerTopLeft);
  samplePlayerTopLeftModule->setModuleName(juce::String(T("SamplePlayerTopLeft")));
  addChildAudioModule(samplePlayerTopLeftModule);

  samplePlayerTopRightModule = new SamplePlayerAudioModule(plugInLock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerTopRight);
  samplePlayerTopRightModule->setModuleName(juce::String(T("SamplePlayerTopRight")));
  addChildAudioModule(samplePlayerTopRightModule);

  samplePlayerBottomLeftModule = new SamplePlayerAudioModule(plugInLock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerBottomLeft);
  samplePlayerBottomLeftModule->setModuleName(juce::String(T("SamplePlayerBottomLeft")));
  addChildAudioModule(samplePlayerBottomLeftModule);

  samplePlayerBottomRightModule = new SamplePlayerAudioModule(plugInLock, 
    &wrappedWorkhorse->voiceArray[0].oscSection.samplePlayerBottomRight);
  samplePlayerBottomRightModule->setModuleName(juce::String(T("SamplePlayerBottomRight")));
  addChildAudioModule(samplePlayerBottomRightModule);

  filterModule = new MultiModeFilterAudioModule(plugInLock, &wrappedWorkhorse->voiceArray[0].filter);
  filterModule->setModuleName(juce::String(T("Filter")));
  addChildAudioModule(filterModule);

  pitchEnvModule = new BreakpointModulatorAudioModule(plugInLock, &wrappedWorkhorse->voiceArray[0].pitchEnv);
  pitchEnvModule->setModuleName(juce::String(T("PitchEnvelope")));
  addChildAudioModule(pitchEnvModule);

  filterEnvModule = new BreakpointModulatorAudioModule(plugInLock, &wrappedWorkhorse->voiceArray[0].filterEnv);
  filterEnvModule->setModuleName(juce::String(T("FilterEnvelope")));
  addChildAudioModule(filterEnvModule);

  ampEnvModule = new BreakpointModulatorAudioModule(plugInLock, &wrappedWorkhorse->voiceArray[0].ampEnv);
  ampEnvModule->setModuleName(juce::String(T("AmpEnvelope")));
  addChildAudioModule(ampEnvModule);
}