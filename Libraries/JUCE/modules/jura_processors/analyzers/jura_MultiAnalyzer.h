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

  /**

  This class wraps rosic::Oscilloscope into a rosof::AudioModule to facilitate its use as plugIn.

  */

class OscilloscopeAudioModule : public AudioModule
{

  friend class OscilloscopeModuleEditor;
  friend class MultiAnalyzerAudioModule;
  friend class MultiAnalyzerModuleEditor;

public:

  OscilloscopeAudioModule(CriticalSection *newPlugInLock, rosic::SyncedWaveformDisplayBuffer *displayBufferToUse);
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void setSampleRate(double newSampleRate) { waveformDisplayBuffer->setSampleRate(newSampleRate); }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    waveformDisplayBuffer->feedInputBuffer(left, numSamples);
    // later: feed both signals (left and right) into a stereo version...
  }

  juce_UseDebuggingNewOperator;

protected:

  virtual void initializeAutomatableParameters();

  //rosic::OscilloscopeBufferOld *wrappedOscilloscope;

  rosic::SyncedWaveformDisplayBuffer *waveformDisplayBuffer;
  bool displayIsFrozen;

};


//=======================================================================================================================================

/**

This class wraps rosic::SpectrumAnalyzer into a rosof::AudioModule to facilitate its use as plugIn.

*/

class SpectrumAnalyzerAudioModule : public AudioModule
{

  friend class SpectrumAnalyzerModuleEditor;
  friend class MultiAnalyzerAudioModule;
  friend class MultiAnalyzerModuleEditor;

public:

  SpectrumAnalyzerAudioModule(CriticalSection *newPlugInLock, rosic::SpectrumAnalyzer *spectrumAnalyzerToWrap);
  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void setSampleRate(double newSampleRate) { wrappedSpectrumAnalyzer->setSampleRate(newSampleRate); }
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedSpectrumAnalyzer->measureSampleFrameStereo(inOutL, inOutR);
  }
  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    for(int n = 0; n < numSamples; n++)
    {
      // this conversion sucks - fix it
      double l = left[n];
      double r = right[n];
      wrappedSpectrumAnalyzer->measureSampleFrameStereo(&l, &r);
    }
  }

  juce_UseDebuggingNewOperator;

protected:

  virtual void initializeAutomatableParameters();

  rosic::SpectrumAnalyzer *wrappedSpectrumAnalyzer;

  bool displayIsFrozen;

};

//=======================================================================================================================================

/**

This class is a analyzer with different analysis-modes (currently oscilloscope and spectrum analysis)

*/

class MultiAnalyzerAudioModule : public AudioModule
{

  friend class MultiAnalyzerModuleEditor;

public:

  enum modes
  {
    OSCILLOSCOPE = 0,
    SPECTRUM_ANALYZER
  };

  //-------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  MultiAnalyzerAudioModule(CriticalSection *newPlugInLock, OscilloscopeAudioModule *newOscilloscopeModule,
    SpectrumAnalyzerAudioModule *newSpectrumAnalyzerModule);

  //-------------------------------------------------------------------------------------------------------------------------------------
  // audio-setup and -processing:

  /** Selects one of the analysis modes. @see: modes */
  virtual void setMode(int newMode)
  {
    getParameterByName(juce::String("Mode"))->setValue(newMode, true, true);
  }

  /** Returns the current analysis mode. @see: modes */
  virtual int getMode() const { return mode; }

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setSampleRate(double newSampleRate)
  {
    // preliminary - implement functions in the classes of the emebedded objects and rely on baseclass implementation to iterate over
    // child audio modules
    spectrumAnalyzerModule->wrappedSpectrumAnalyzer->setSampleRate(newSampleRate);
    oscilloscopeModule->waveformDisplayBuffer->setSampleRate(newSampleRate);
  }

  /*
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    switch( mode )
    {
    case OSCILLOSCOPE:
      oscilloscopeModule->getSampleFrameStereo(inOutL, inOutR);      break;
    case SPECTRUM_ANALYZER:
      spectrumAnalyzerModule->getSampleFrameStereo(inOutL, inOutR);  break;
    }
  }
  */

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    switch(mode)
    {
    case OSCILLOSCOPE:
      oscilloscopeModule->processBlockStereo(left, right, numSamples);      break;
    case SPECTRUM_ANALYZER:
      spectrumAnalyzerModule->processBlockStereo(left, right, numSamples);  break;
    }

    /*
    // preliminary - implement processBlockStereo in the OscilloscopeAudioModule and SpectrumAnalyzerAudioModule
    // classes later
    for(int n=0; n<numSamples; n++)
    {
      double dL = left [n];
      double dR = right[n];

      //wrappedFuncShaper->getSampleFrameStereo(&dL, &dR, &dL, &dR);
      switch( mode )
      {
      case OSCILLOSCOPE:
        oscilloscopeModule->getSampleFrameStereo(&dL, &dR);      break;
      case SPECTRUM_ANALYZER:
        spectrumAnalyzerModule->getSampleFrameStereo(&dL, &dR);  break;
      }
      //left [n] = (float) dL;
      //right[n] = (float) dR;
    }
    */
  }

protected:

  OscilloscopeAudioModule     *oscilloscopeModule;
  SpectrumAnalyzerAudioModule *spectrumAnalyzerModule;

  virtual void initializeAutomatableParameters();

  int mode;

  
  juce_UseDebuggingNewOperator;
};

//=================================================================================================

//class SpectrumAnalyzerDisplay	: virtual public MessengingCoordinateSystemOld, 
//  virtual public SpectrumDisplayOld
//{
//
//public:
//
//  SpectrumAnalyzerDisplay(const juce::String& name = juce::String(T("SpectrumAnalyzerDisplay")));   
//  /**< Constructor. */
//
//  virtual ~SpectrumAnalyzerDisplay(); 
//  /**< Destructor. */
//
//  virtual void useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBase = 2.0);
//  /**< Overrides the inherited function from CoordinateSystem in order to adjust the 
//  grid-spacing. */
//
//  virtual void paint(Graphics &g);
//  /**< Overrides the paint-function of the base-classes. */
//
//  virtual void updateBackgroundImage();
//  /**< Overrides the inherited method from the CoordinateSystem base-class. */
//
//
//protected:
//
//  virtual void plotCurveFamily(Graphics &g, Image* targetImage = NULL, XmlElement *targetSVG = NULL);
//  /**< Overrides the plotCurveFamily()-function of the CurveFamilyPlot base-class. */
//
//  virtual bool getRepresentingBins(double lowFreq, double highFreq, int k, 
//    int &minBin, int &maxBin);
//  /**< Returns true and stores the index of the bins which best represent a range of frequencies
//  between lowFreq and highFreq in minBin and maxBin.  If there is no bin satisfying the constraint 
//  of being (strictly) lower than highFreq and higher or equal to lowFreq, it will return false. If 
//  there are several bins between lowFreq and highFreq which satisfy this constraint, it will 
//  return the bins with the largest and smallest magnitude. */
//
//
//  juce_UseDebuggingNewOperator;
//};

//===============================================================================================


#endif 
