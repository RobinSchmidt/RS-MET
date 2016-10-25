#ifndef jura_LadderFilter_h
#define jura_LadderFilter_h
  
// We want to use two RAPT::LadderFilter instances that use double precision numbers for the 
// signal as well as for the parameters/coefficients. For convenience, we make a typedef for the
// particular template instantiation:
typedef RAPT::LadderFilter<double, double> RAPTLadder;

//=================================================================================================

/** Wraps two RAPT LadderFilter instances into a single object, for stereo processing. In order to
make it more interesting than just instantiating the RAPT::LadderFilter template for a kind of 
stereo signal type (which would be also possible), we allow the two cutoff frequencies of both 
channels to be different by introducing a stereo spread parameter. */

class JUCE_API Ladder : public jura::AudioModule
{

public:

  Ladder();
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createStaticParameters();

  /** Creates the GUI editor (returns an object of an appropriate subclass of 
  AudioProcessorEditor) */
  AudioProcessorEditor *createEditor() override;

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

  // target functions for callbacks that are called on parameter changes:
  void setCutoff(double newCutoff);
  void setResonance(double newResonance);
  void setMode(int newMode);
  void setStereoSpread(double newSpreadInSemitones);
  void setMidSideMode(bool shouldBeInMidSideMode);

protected:

  // some members to store the parameter values:
  double cutoff      = 1000.0;
  double spread      = 0.0;
  double freqFactorL = 1.0;
  double freqFactorR = 1.0;
  bool   midSideMode = false;

  // embedded core DSP objects from the RAPT library:
  RAPTLadder ladderL, ladderR;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Ladder)
};

//=================================================================================================

/** This is the GUI editor class for the jura::Ladder audio processor. */

class JUCE_API LadderEditor : public AudioModuleEditor
{

public:

  LadderEditor(jura::Ladder *newLadderToEdit);
  virtual void resized();

protected:

  Ladder *ladderToEdit;

  RSlider *cutoffSlider, *resonanceSlider, *spreadSlider;
  RComboBox *modeComboBox;
  //RButton *invertButton;  

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LadderEditor)
};

#endif 