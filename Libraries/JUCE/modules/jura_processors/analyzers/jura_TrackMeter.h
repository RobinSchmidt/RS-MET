#ifndef jura_TrackMeter_h
#define jura_TrackMeter_h

class TrackMeterAudioModule : public AudioModule
{

  friend class TrackMeterModuleEditor;

public:

  TrackMeterAudioModule(CriticalSection *newPlugInLock, 
    rosic::TrackMeter *trackMeterToWrap = nullptr);

  virtual ~TrackMeterAudioModule();

  AudioModuleEditor* createEditor() override;

  virtual void parameterChanged(Parameter* parameterThatHasChanged);

  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedTrackMeter != NULL)
      wrappedTrackMeter->setSampleRate(newSampleRate);
  }

  virtual void getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    wrappedTrackMeter->measureSampleFrameStereo(inOutL, inOutR);
  }

protected:

  virtual void initializeAutomatableParameters();

  rosic::TrackMeter *wrappedTrackMeter;
  bool wrappedTrackMeterIsOwned= false;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class TrackMeterModuleEditor : public AudioModuleEditor, public Timer
{

public:

  TrackMeterModuleEditor(CriticalSection *newPlugInLock, 
    TrackMeterAudioModule* newTrackMeterAudioModule);

  // callbacks:
  virtual void rButtonClicked(RButton* buttonThatWasClicked);
  virtual void paint(Graphics &g);
  virtual void resized();
  virtual void timerCallback();

protected:

  /** Draws the numerical scales next to the meters. */
  void drawMeterScales(Graphics &g);

  TrackMeterAudioModule *trackMeterModuleToEdit;

  RSlider *riseSlider, *fallSlider;
  RButton *vuButton, *ppmButton;   // maybe use RClickButton
  RTextField *leftLevelLabel, *rightLevelLabel, *midLevelLabel, *sideLevelLabel, *correlationLabel;
  MeteringDisplay *leftLevelMeter, *rightLevelMeter, *midLevelMeter, *sideLevelMeter, 
    *correlationMeter;

  double rangeMin, rangeMax;

  juce_UseDebuggingNewOperator;
};

#endif 
