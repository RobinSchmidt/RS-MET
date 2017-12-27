#ifndef jura_DebugModule_h
#define jura_DebugModule_h

/** A module for framework debugging puposes */

class JUCE_API DebugAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  DebugAudioModule(CriticalSection *lockToUse);
    
  virtual void createParameters();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;
  virtual void setMidiController(int controllerNumber, float controllerValue) override;

  // callback target functions:
  void setLeftValue( double newValue) { values[0] = newValue;  }
  void setRightValue(double newValue) { values[1] = newValue;  }

protected:


  static const int numValues = 2;
  double values[numValues] = { 0, 0 };

  MetaControlledParameter *leftParam, *rightParam;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(DebugAudioModule)
};

//=================================================================================================

/** Editor for the debugging AudioModule. */

#endif 