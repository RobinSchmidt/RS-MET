#ifndef jura_FuncShaper_h
#define jura_FuncShaper_h



/** This class wraps rosic::FuncShaper into a jura::AudioModule to facilitate its use as plugIn. */

class FuncShaperAudioModule : public ModulatableAudioModule, public ChangeBroadcaster
{

  friend class FuncShaperModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  FuncShaperAudioModule(CriticalSection *newPlugInLock, rosic::FuncShaper *funcShaperToWrap, 
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);

  FuncShaperAudioModule(CriticalSection *newPlugInLock, 
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr);

  void init();

  virtual ~FuncShaperAudioModule();

  AudioModuleEditor* createEditor() override;

  //---------------------------------------------------------------------------------------------
  // automation and state management:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

  /** Used internally by setStateFromXml to recall the a,b,c,d parameters along with their 
  min/max values. "name" is any of "a".."d". */
  void recallFormulaParameterFromXml(const XmlElement& xml, const juce::String& name);

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate) override
  {
    wrappedFuncShaper->setSampleRate(newSampleRate);
  }

  // formula parameter range settings, name can be any of a,b,c,d: 
  void setFormulaParameterMin(const juce::String& name, double newMin);
  void setFormulaParameterMax(const juce::String& name, double newMax);
  void setFormulaParameterAndRange(const juce::String& name, double newValue, 
    double newMin, double newMax);

  // callback targets:
  void setA(double newA) { wrappedFuncShaper->setA(newA, autoRecalcTable); }
  void setB(double newB) { wrappedFuncShaper->setB(newB, autoRecalcTable); }
  void setC(double newC) { wrappedFuncShaper->setC(newC, autoRecalcTable); }
  void setD(double newD) { wrappedFuncShaper->setD(newD, autoRecalcTable); }

  void setMinA(double newMin) { setFormulaParameterMin("a", newMin); }
  void setMinB(double newMin) { setFormulaParameterMin("b", newMin); }
  void setMinC(double newMin) { setFormulaParameterMin("c", newMin); }
  void setMinD(double newMin) { setFormulaParameterMin("d", newMin); }

  void setMaxA(double newMax) { setFormulaParameterMax("a", newMax); }
  void setMaxB(double newMax) { setFormulaParameterMax("b", newMax); }
  void setMaxC(double newMax) { setFormulaParameterMax("c", newMax); }
  void setMaxD(double newMax) { setFormulaParameterMax("d", newMax); }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedFuncShaper->getSampleFrameStereo(inOutL, inOutR, inOutL, inOutR);
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    for(int n = 0; n < numSamples; n++)
    {
      double dL = left[n];
      double dR = right[n];
      wrappedFuncShaper->getSampleFrameStereo(&dL, &dR, &dL, &dR);
      left[n] = (float)dL;
      right[n] = (float)dR;
    }
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedFuncShaper->getSampleFrameStereo(
        &inOutBuffer[0][n], &inOutBuffer[1][n],   // inputs
        &inOutBuffer[0][n], &inOutBuffer[1][n]);  // outputs
    // provide a function that we cann call like:
    // wrappedFuncShaper->processBlockStereo(inOutBuffer[0], inOutBuffer[1], numSamples);
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    wrappedFuncShaper->getSampleFrameStereo(left, right, left, right);
  }


protected:

  void createParameters();

  rosic::FuncShaper *wrappedFuncShaper;
  bool wrappedFuncShaperIsOwned = false;
  bool autoRecalcTable = true;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class FuncShaperModuleEditor : public AudioModuleEditor, public RTextEntryFieldObserver
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  FuncShaperModuleEditor(CriticalSection *newPlugInLock, 
    FuncShaperAudioModule* newFuncShaperAudioModule);

  /** Destructor. */
  virtual ~FuncShaperModuleEditor();

  //---------------------------------------------------------------------------------------------
  // callbacks:

  /** Gets called from the funcShapeAudioModule when the formula or one of the a,b,c,d parameters
  changes (in which case the plot must be updated). */
  virtual void changeListenerCallback(ChangeBroadcaster *source) override;

  /** Implements the purely virtual method inherited from RTextEntryFieldObserver. */
  virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged);

  /** Overrides paint(). */   
  virtual void paint(Graphics &g);

  /** Overrides resized(). */    
  virtual void resized();

protected:

  /** Overrides the method inherited from RPolyphonicInstrumentEditor. */
  virtual void updateWidgetsAccordingToState();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  FuncShaperAudioModule *funcShaperAudioModule;

  // some rectangles to subdivide the GUI:
  juce::Rectangle<int> formulaRectangle, inputRectangle, outputRectangle;

  // the widgets:
  RTextField *formulaLabel, *inputLabel, *outputLabel; //, *aLabel, *bLabel, *cLabel, *dLabel;

  RTextEntryField *formulaField, *aMinField, *aMaxField, *bMinField, *bMaxField, 
    *cMinField, *cMaxField, *dMinField, *dMaxField;

  RSlider *aSlider, *bSlider, *cSlider, *dSlider, *driveSlider, *dcSlider, *oversamplingSlider, 
    *dryWetSlider, *outVolumeSlider, *inHighpassSlider, *inLowpassSlider, *outHighpassSlider,
    *outLowpassSlider;
  RButton *preFilterButton, *postFilterButton;

  // the function plot which plots the characteristic line:
  CurveFamilyPlotOld* shaperPlot;
  double *xValues;
  double **yValues;

  juce_UseDebuggingNewOperator;
};

#endif 
