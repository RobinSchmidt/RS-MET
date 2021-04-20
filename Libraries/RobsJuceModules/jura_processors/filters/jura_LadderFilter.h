#ifndef jura_LadderFilter_h
#define jura_LadderFilter_h
  
//=================================================================================================

/** Wraps two RAPT LadderFilter instances into a single object, for stereo processing. In order to
make it more interesting than just instantiating the RAPT::LadderFilter template for a kind of 
stereo signal type (which would be also possible), we allow the two cutoff frequencies of both 
channels to be different by introducing a stereo spread parameter. 

\todo: 
-rename the class to LadderStereo
-make the GUI editor work for the basic ladder filter (not only for this stereo version)
 (maybe we have to make this a subclass of the underlying RAPT::LadderFilter)
-let the user freely adjust the mixing coefficients for the filters stages (via the DraggableNumber
 widget) - the mode-combobox then selects between various predefined settings
*/

class JUCE_API Ladder : public jura::ModulatableAudioModule
  // : public jura::AudioModuleWithMidiIn // why WithMidiIn?
{

public:

  using LDR = RAPT::rsLadderFilter<double, double>;

  Ladder(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr);
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createStaticParameters();
    // maybe rename to createParameters

  /** Creates the GUI editor (returns an object of an appropriate subclass of AudioModuleEditor) */
  AudioModuleEditor *createEditor(int type) override;

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

  // target functions for callbacks that are called on parameter changes:
  void setCutoff(double newCutoff);
  void setResonance(double newResonance);
  void setMode(int newMode);
  void setStereoSpread(double newSpreadInSemitones);
  void setMidSideMode(bool shouldBeInMidSideMode);

  // for the magnitude plot/editor later...
  //double getCutoff()    { return cutoff; }
  //double getResonance() { return 0.0;    }  // preliminary

  /** Returns the magnitude response expressed in decibels at the given frequency. */
  double getDecibelsAt(double frequency); 


protected:

  // some members to store the parameter values:
  double cutoff      = 1000.0;
  double spread      = 0.0;
  double freqFactorL = 1.0;
  double freqFactorR = 1.0;
  bool   midSideMode = false;

  // embedded core DSP objects from the rosic library:
  rosic::rsLadderFilterStereo wrappedLadder;
    // use the one with different cutoffs for the channels later (maybe)

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Ladder)
};

//=================================================================================================

class JUCE_API rsLadderPlotEditor : public ColourSchemeComponent, public ParameterObserver
{

public:

  rsLadderPlotEditor(jura::Ladder* ladderModuleToEdit);
  virtual ~rsLadderPlotEditor();

  virtual void parameterChanged(Parameter* p) override;
  virtual void resized() override;

protected:

  jura::Ladder* ladderToEdit; 

  Parameter *cutoffParam, *resoParam;

  rsFunctionPlot* freqRespPlot; // frequency response plot (maybe factor out a pointer to rsPlot)
  rsVectorPad*    vectorPad;    // for the handle, overlaid over the plot

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsLadderPlotEditor)
};

//=================================================================================================

/** This is the GUI editor class for the jura::Ladder audio processor. */

class JUCE_API LadderEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  LadderEditor(jura::Ladder *newLadderToEdit);
  virtual void resized() override;
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;

protected:

  Ladder *ladderToEdit;
  rsLadderPlotEditor* plotEditor;
  rsModulatableSlider *cutoffSlider, *resonanceSlider, *spreadSlider;
  rsAutomatableComboBox *modeComboBox;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LadderEditor)
};

#endif 