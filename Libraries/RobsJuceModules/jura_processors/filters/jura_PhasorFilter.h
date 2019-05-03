#ifndef jura_PhasorFilter_h
#define jura_PhasorFilter_h
  
typedef RAPT::rsPhasorFilter<double, double> RAPTPhasorFilter;
typedef RAPT::rsPhasorStateMapper<double> RAPTPhasorMapper;

/** Wraps a RAPT::PhasorFilter into an AudioModule. */

class JUCE_API PhasorFilter : public jura::AudioModule
{

public:

  PhasorFilter(CriticalSection *lockToUse);
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createStaticParameters();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

protected:

  RAPTPhasorFilter filterCore;
  RAPTPhasorMapper stateMapper;

  double inOld;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhasorFilter)
};

#endif 