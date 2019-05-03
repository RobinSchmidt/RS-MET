#ifndef jura_Quadriga_h
#define jura_Quadriga_h

class QuadrigaAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class QuadrigaModuleEditor;

public:

  QuadrigaAudioModule(CriticalSection *newPlugInLock);

  virtual ~QuadrigaAudioModule();

  AudioModuleEditor* createEditor(int type) override;


  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate) override
  {
    if(wrappedQuadriga != NULL)
      wrappedQuadriga->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    if(wrappedQuadriga != NULL)
      wrappedQuadriga->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void reset() override
  {
    if(wrappedQuadriga != NULL)
      wrappedQuadriga->resetAllVoices();
  }

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  QuadrigenAudioModule *quadrigenModule;
  QuadrifexAudioModule *quadrifexModule;
  EqualizerAudioModule *equalizerModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::Quadriga *wrappedQuadriga;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class QuadrigaModuleEditor : public PolyphonicInstrumentEditor
{

public:

  QuadrigaModuleEditor(CriticalSection *newPlugInLock, QuadrigaAudioModule* newQuadrigaAudioModule);

  // callbacks:
  virtual void rButtonClicked(RButton *buttonThatWasClicked);

  /** Overrides changeListenerCallback() in order to start or stop the timers in the embedded
  module editors when the user selects another tab.*/
  virtual void changeListenerCallback(ChangeBroadcaster* objectThatHasChanged);

  virtual void resized();

protected:

  /** Overrides the method inherited from AudioModuleEditor. */
  virtual void updateWidgetsAccordingToState();

  /** Updates the visibilities of the sub-editors based on which button is down. */
  virtual void updateSubEditorVisibility();

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  QuadrigaAudioModule *quadrigaAudioModule;

  // buttons to switch between the pages:
  RButton *performanceButton, *generatorsButton, *filtersButton, *modulatorsButton, 
    *effectsButton;

  // the sub-editors:
  QuadrigenModuleEditor *quadrigenEditor;
  QuadrifexModuleEditor *quadrifexEditor;

  juce_UseDebuggingNewOperator;
};

#endif 
