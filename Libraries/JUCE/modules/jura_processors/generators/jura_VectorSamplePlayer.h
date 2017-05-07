#ifndef jura_VectorSamplePlayer_h
#define jura_VectorSamplePlayer_h


class VectorSamplePlayerAudioModule : public AudioModule
{

  friend class VectorSamplePlayerEditor;

public:


  VectorSamplePlayerAudioModule(CriticalSection *newPlugInLock, 
    rosic::VectorSamplePlayer *vectorSamplePlayerToWrap);

  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedVectorSamplePlayer != NULL)
      wrappedVectorSamplePlayer->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedVectorSamplePlayer->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset()
  {
    if(wrappedVectorSamplePlayer != NULL)
      wrappedVectorSamplePlayer->reset();
  }

protected:

  // we maintain wrappped versions (into jura::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  SamplePlayerAudioModule *samplePlayerTopLeftModule, *samplePlayerTopRightModule,
    *samplePlayerBottomLeftModule, *samplePlayerBottomRightModule;
  VectorMixerAudioModule *vectorMixerModule;
  LowFrequencyOscillatorAudioModule *xLfoModule, *yLfoModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::VectorSamplePlayer *wrappedVectorSamplePlayer;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Subclass of SamplePlayerModuleEditor to customize the arrangement of widgets in resized. */

class VectorSamplePlayerSampleEditor : public SamplePlayerModuleEditor
{
public:
  enum layouts { TOP_LEFT, TOP_RIGHT, BOTTOM_LEFT, BOTTOM_RIGHT };
  VectorSamplePlayerSampleEditor(CriticalSection *newPlugInLock,
    SamplePlayerAudioModule* newSamplePlayerAudioModule);
  virtual void setLayout(int newLayout);
  virtual void resized();
  static const int verticalIntrusion   = 40;
  static const int horizontalIntrusion = 90;
  juce_UseDebuggingNewOperator;
protected:
  int layout, headlineWidth;
};

/** Subclass of LowFrequencyOscillatorEditor to customize the arrangement of widgets in 
resized. */
class VectorSamplePlayerLfoEditor : public LowFrequencyOscillatorEditor
{
public:
  enum layouts { LEFT, RIGHT };
  VectorSamplePlayerLfoEditor(CriticalSection *newPlugInLock, 
    LowFrequencyOscillatorAudioModule* newLowFrequencyOscillatorAudioModule);
  virtual void setLayout(int newLayout);
  virtual void paint(Graphics &g);
  virtual void resized();
protected:
  int layout, headlineWidth;
  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Editor for the VectorSamplerAudioModule. */

class VectorSamplePlayerEditor : public AudioModuleEditor
{

public:
  VectorSamplePlayerEditor(CriticalSection *newPlugInLock, 
    VectorSamplePlayerAudioModule* newVectorSamplePlayerAudioModule);

  virtual void resized() override;
  virtual void updateWidgetsAccordingToState() override;

protected:

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  VectorSamplePlayerAudioModule *vectorSamplePlayerAudioModule;

  // the sub-editors:
  VectorMixerModuleEditor *vectorMixerPad;
  VectorSamplePlayerSampleEditor *samplePlayerTopLeftEditor, *samplePlayerTopRightEditor, 
    *samplePlayerBottomLeftEditor, *samplePlayerBottomRightEditor;
  VectorSamplePlayerLfoEditor    *xLfoEditor, *yLfoEditor;

  juce_UseDebuggingNewOperator;
};

#endif 
