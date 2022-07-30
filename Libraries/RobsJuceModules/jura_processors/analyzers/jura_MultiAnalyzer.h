#ifndef jura_MultiAnalyzer_h
#define jura_MultiAnalyzer_h

////#include "../../../rosic/analysis/rosic_OscilloscopeBufferOld.h"
//#include "../../../rosic/analysis/rosic_WaveformDisplayBuffer.h"
//#include "../../../rosic/analysis/rosic_SpectrumAnalyzer.h"
//using namespace rosic;
//
//#include "../rosof_AudioModule.h"
//
//#include "../../../rojue/components/widgets/rojue_PreDefinedParameters.h"


//=======================================================================================================================================

/** This class wraps rosic::Oscilloscope into a rosof::AudioModule to facilitate its use as 
plugIn. */

class OscilloscopeAudioModule : public AudioModule
{

  friend class OscilloscopeModuleEditor;
  friend class MultiAnalyzerAudioModule;
  friend class MultiAnalyzerModuleEditor;

public:

  OscilloscopeAudioModule(CriticalSection *newPlugInLock, 
    rosic::OscilloscopeBufferOld* displayBufferToUse);

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setSampleRate(double newSampleRate) override
  { 
    waveformBuffer->setSampleRate(newSampleRate); 
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
    {
      //waveformBuffer->feedInputBuffer(&inOutBuffer[0][n], numSamples);
      // preliminary - feeds only left channel into the buffer - make stereo version

      waveformBuffer->bufferSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
    }
  }

  //virtual void processBlockStereo(float *left, float *right, int numSamples)
  //{
  //  waveformBuffer->feedInputBuffer(left, numSamples);
  //  // later: feed both signals (left and right) into a stereo version...
  //}

  juce_UseDebuggingNewOperator;

protected:

  virtual void initializeAutomatableParameters();

  rosic::OscilloscopeBufferOld* waveformBuffer;
  //rosic::SyncedWaveformDisplayBuffer* waveformDisplayBuffer;

  bool displayIsFrozen;

};

//=================================================================================================

/** This class wraps rosic::SpectrumAnalyzer into a rosof::AudioModule to facilitate its use as 
plugIn. */

class SpectrumAnalyzerAudioModule : public AudioModule
{

  friend class SpectrumAnalyzerModuleEditor;
  friend class MultiAnalyzerAudioModule;
  friend class MultiAnalyzerModuleEditor;

public:

  SpectrumAnalyzerAudioModule(CriticalSection *newPlugInLock, 
    rosic::SpectrumAnalyzer *spectrumAnalyzerToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setSampleRate(double newSampleRate) override
  { 
    wrappedSpectrumAnalyzer->setSampleRate(newSampleRate); 
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedSpectrumAnalyzer->measureSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  //{
  //  wrappedSpectrumAnalyzer->measureSampleFrameStereo(inOutL, inOutR);
  //}
  //virtual void processBlockStereo(float *left, float *right, int numSamples)
  //{
  //  for(int n = 0; n < numSamples; n++)
  //  {
  //    // this conversion sucks - fix it
  //    double l = left[n];
  //    double r = right[n];
  //    wrappedSpectrumAnalyzer->measureSampleFrameStereo(&l, &r);
  //  }
  //}

  juce_UseDebuggingNewOperator;

protected:

  virtual void initializeAutomatableParameters();

  rosic::SpectrumAnalyzer *wrappedSpectrumAnalyzer;

  bool displayIsFrozen;

};

//=================================================================================================

/** This class is a analyzer with different analysis-modes (currently oscilloscope and spectrum 
analysis) 
todo: add xy-scope, metering, spectrogram, pitch-detection, etc.
*/

class MultiAnalyzerAudioModule : public AudioModule
{

  friend class MultiAnalyzerModuleEditor;

public:

  enum modes
  {
    WAVEFORM = 0,
    SPECTRUM
    // METERS,
    // SPECTROGRAM,
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  //MultiAnalyzerAudioModule(CriticalSection *newPlugInLock, 
  //  OscilloscopeAudioModule *newOscilloscopeModule,
  //  SpectrumAnalyzerAudioModule *newSpectrumAnalyzerModule);

  MultiAnalyzerAudioModule(CriticalSection *newPlugInLock);

  virtual ~MultiAnalyzerAudioModule();

  AudioModuleEditor* createEditor(int type) override;

  //-----------------------------------------------------------------------------------------------
  // audio-setup and -processing:

  /** Selects one of the analysis modes. @see: modes */
  virtual void setMode(int newMode)
  {
    getParameterByName(juce::String("Mode"))->setValue(newMode, true, true);
  }

  /** Returns the current analysis mode. @see: modes */
  virtual int getMode() const { return mode; }

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setSampleRate(double newSampleRate) override
  {
    // preliminary - implement functions in the classes of the emebedded objects and rely on 
    // baseclass implementation to iterate over child audio modules
    spectrumAnalyzerModule->wrappedSpectrumAnalyzer->setSampleRate(newSampleRate);
    oscilloscopeModule->waveformBuffer->setSampleRate(newSampleRate);
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    switch(mode)
    {
    case WAVEFORM: 
      oscilloscopeModule->processBlock(inOutBuffer, numChannels, numSamples);     break;
    case SPECTRUM:
      spectrumAnalyzerModule->processBlock(inOutBuffer, numChannels, numSamples); break;
    }
  }

  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  //{
  //  switch( mode )
  //  {
  //  case OSCILLOSCOPE:
  //    oscilloscopeModule->getSampleFrameStereo(inOutL, inOutR);      break;
  //  case SPECTRUM_ANALYZER:
  //    spectrumAnalyzerModule->getSampleFrameStereo(inOutL, inOutR);  break;
  //  }
  //}

  //virtual void processBlockStereo(float *left, float *right, int numSamples)
  //{
  //  switch(mode)
  //  {
  //  case OSCILLOSCOPE:
  //    oscilloscopeModule->processBlockStereo(left, right, numSamples);      break;
  //  case SPECTRUM_ANALYZER:
  //    spectrumAnalyzerModule->processBlockStereo(left, right, numSamples);  break;
  //  }
  //}

protected:

  OscilloscopeAudioModule*     oscilloscopeModule;
  SpectrumAnalyzerAudioModule* spectrumAnalyzerModule;

  //rosic::SyncedWaveformDisplayBuffer* wrappedDisplayBuffer;
  rosic::OscilloscopeBufferOld* waveformBuffer;
  rosic::SpectrumAnalyzer*      wrappedSpectrumAnalyzer; // rename to spectrumBuffer

  virtual void initializeAutomatableParameters();

  int mode;

  
  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class SpectrumAnalyzerDisplay	: virtual public rsObservablePlot, 
  virtual public rsSpectrumPlot
{

public:

  SpectrumAnalyzerDisplay(const juce::String& name = juce::String("SpectrumAnalyzerDisplay"));   
  /**< Constructor. */

  virtual ~SpectrumAnalyzerDisplay(); 
  /**< Destructor. */

  virtual void useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBase = 2.0);
  /**< Overrides the inherited function from CoordinateSystem in order to adjust the 
  grid-spacing. */

  virtual void paint(Graphics &g);
  /**< Overrides the paint-function of the base-classes. */

  virtual void updateBackgroundImage();
  /**< Overrides the inherited method from the CoordinateSystem base-class. */


protected:

  virtual void plotCurveFamily(Graphics &g, juce::Image* targetImage = NULL, 
    XmlElement *targetSVG = NULL) override;
  /**< Overrides the plotCurveFamily()-function of the CurveFamilyPlot base-class. */

  virtual bool getRepresentingBins(double lowFreq, double highFreq, int k, 
    int &minBin, int &maxBin);
  /**< Returns true and stores the index of the bins which best represent a range of frequencies
  between lowFreq and highFreq in minBin and maxBin.  If there is no bin satisfying the constraint 
  of being (strictly) lower than highFreq and higher or equal to lowFreq, it will return false. If 
  there are several bins between lowFreq and highFreq which satisfy this constraint, it will 
  return the bins with the largest and smallest magnitude. */


  juce_UseDebuggingNewOperator;
};

//===============================================================================================

class OscilloscopeDisplay	: virtual public rsObservablePlot, virtual public rsDataPlot
{

public:

  /** Constructor. */
  OscilloscopeDisplay(const juce::String& name = juce::String("OscilloscopeDisplay"));   

  /** Accepts values for a family (i.e. several channels) of waveform data.  */
  virtual void setWaveformData(int newNumSamples, int newNumChannels, float** newData, 
    double* newTimeAxis);

  /** Override, to fix the currentMinX to zero when zooming in. */
  virtual void setCurrentRangeX(double newMinX, double newMaxX);

  /** Override, to trigger re-rendering. */
  virtual void setCurrentRangeY(double newMinY, double newMaxY);

  /** Override, to fix the currentMinX to zero when zooming in. */
  virtual void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);


protected:

  /** Overrides the plotCurveFamily()-function of the CurveFamilyPlot base-class. */
  virtual void plotCurveFamily(Graphics &g, juce::Image* targetImage = NULL, 
    XmlElement *targetSVG = NULL);

  /** Adjusts the grid-settings acording to the zoom-factor. */
  virtual void adjustGrid();

  double* timeAxis;
  float** peakData;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

// todo - maybe move this class to another file - might be useful for EQ with spectrum display 
// later also
class AudioModuleEditorAnimated : public AudioModuleEditor, public Timer, public RSliderListener
{

public:

  AudioModuleEditorAnimated(CriticalSection *newPlugInLock, AudioModule* newAudioModuleToEdit);

  // overrides:
  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChamged);
  virtual void timerCallback();
  virtual void setVisible(bool  shouldBeVisible);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

  /** Sets the number of frames per second. */
  virtual void setFrameRate(double newFrameRate);

  /** Freezes the display. */
  virtual void setFreeze(bool shouldBeFrozen);

  juce_UseDebuggingNewOperator;

protected:

  /** This is the method you must override in order to update the content of the editor when a new frame should be rendered. */
  virtual void updateEditorContent() = 0;

  /** Updates the time settings according to the state of frameRateSlider, freezeButton and the Component's visibility. */
  virtual void updateTimerSettings();

  RSlider *frameRateSlider;
  RButton *freezeButton; 

};

//=======================================================================================================================================

/** GUI editor for the Oscilloscope */

class OscilloscopeModuleEditor : public AudioModuleEditorAnimated, public rsPlotObserver //, public RComboBoxObserver
{

public:

  OscilloscopeModuleEditor(CriticalSection *newPlugInLock, OscilloscopeAudioModule* newOscilloscopeAudioModule);

  virtual ~OscilloscopeModuleEditor();


  virtual void coordinateSystemChanged(rsObservablePlot *coordinateSystemThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

  juce_UseDebuggingNewOperator;

protected:

  virtual void updateEditorContent();

  OscilloscopeAudioModule *oscilloscopeAudioModule;

  RButton        *midSideButton;
  RNamedComboBox *syncModeComboBox;

  OscilloscopeDisplay       *oscilloscopeDisplay;
  rsPlotZoomer *oscilloscopeZoomer;


  //// signal buffers for display:
  //double *timeAxis;
  //float  *xL, *xR;
  //float  **px;  // not used yet, later - use it like: px[0] = xL, px[1] = xR

};


//=======================================================================================================================================

/** GUI editor for the SpectrumAnalyzer */

class SpectrumAnalyzerModuleEditor : public AudioModuleEditorAnimated, public rsPlotObserver
{

public:

  SpectrumAnalyzerModuleEditor(CriticalSection *newPlugInLock, SpectrumAnalyzerAudioModule* newSpectrumAnalyzerAudioModule);

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void coordinateSystemChanged(rsObservablePlot *coordinateSystemThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

  juce_UseDebuggingNewOperator;

protected:

  virtual void updateEditorContent();

  SpectrumAnalyzerAudioModule *spectrumAnalyzerAudioModule;

  RButton        *midSideButton, *linearFrequencyAxisButton;
  RNamedComboBox *fftSizeComboBox;

  SpectrumAnalyzerDisplay   *spectrumDisplay;
  rsPlotZoomer *spectrumZoomer;

};


//=======================================================================================================================================

/** GUI editor for the MultiAnalyzer */

class MultiAnalyzerModuleEditor : virtual public AudioModuleEditor 
{

public:

  MultiAnalyzerModuleEditor(CriticalSection *newPlugInLock, MultiAnalyzerAudioModule* newMultiAnalyzerAudioModule);

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

  juce_UseDebuggingNewOperator;

protected:

  virtual void updateSubEditorVisibilitiesAndTabButtonStates();
  // Boah - what a monstrous name! Find something shorter! maybe updateVisibilities

  MultiAnalyzerAudioModule     *multiAnalyzerAudioModule;

  OscilloscopeModuleEditor     *oscilloscopeEditor;
  SpectrumAnalyzerModuleEditor *spectrumAnalyzerEditor;

  RButton *oscilloscopeButton, *spectrumAnalyzerButton;

};

#endif 
