#ifndef jura_PhasorFilter_h
#define jura_PhasorFilter_h
  
/** Wraps a RAPT::PhasorFilter into an AudioModule. */

class JUCE_API PhasorFilter : public jura::AudioModule
{

public:

  PhasorFilter(CriticalSection *lockToUse);
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createStaticParameters();

  /** Creates the GUI editor (returns an object of an appropriate subclass of AudioModuleEditor) */
  AudioModuleEditor *createEditor() override;

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

protected:

  RAPT::PhasorFilter<double, double> filterCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhasorFilter)
};

#endif 