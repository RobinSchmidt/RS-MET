#ifndef jura_AudioModuleCreation_h
#define jura_AudioModuleCreation_h


/** A do-nothing dummy AudioModule to be used as placeholder. 
todo:
-rename to DummyAudioModule
-maybe move somewhere else in the library
-maybe override create editor to create a special kind of editor that says something like
"Select Module - Editor will appear here"
-or maybe get rid and use the Bypass module from Quadrifex
*/

class JUCE_API DummyModule : public jura::AudioModule
{
public:
  DummyModule(CriticalSection *lockToUse) : AudioModule(lockToUse) 
  {
    setModuleTypeName("None");
  }
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override 
  {
    //// for debug:
    //std::vector<double> left(numSamples), right(numSamples);
    //for(int n = 0; n < numSamples; n++)
    //{
    //  left[n]  = inOutBuffer[0][n];
    //  right[n] = inOutBuffer[1][n];
    //}
    //int dummy = 0;
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DummyModule)
};

//=================================================================================================

// module creator functions (maybe move into separate file)

class JUCE_API AudioModuleCreator
{

public:

  /*
  static DummyModule*      createDummyModule(CriticalSection* lock) { return new DummyModule(lock); }
  static DebugAudioModule* createDebugModule(CriticalSection* lock) { return new DebugAudioModule(lock); }

  // Analysis:
  static PhaseScope*               createScope(CriticalSection* lock)         { return new PhaseScope(lock); }
  static MultiAnalyzerAudioModule* createMultiAnalyzer(CriticalSection* lock) { return new MultiAnalyzerAudioModule(lock); }
  static TrackMeterAudioModule*    createTrackMeter(CriticalSection* lock)    { return new TrackMeterAudioModule(lock); }
  static MidiMonitorAudioModule*   createMidiMonitor(CriticalSection* lock)   { return new MidiMonitorAudioModule(lock); }
  // hmm...maybe use std::function and lambda functions
  */



  static void populateFactoryForToolChain(AudioModuleFactory* factory, CriticalSection* lock);

};

#endif