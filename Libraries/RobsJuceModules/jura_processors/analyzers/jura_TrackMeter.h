#ifndef jura_TrackMeter_h
#define jura_TrackMeter_h

class TrackMeterAudioModule : public AudioModule
{

  friend class TrackMeterModuleEditor;

public:

  TrackMeterAudioModule(CriticalSection *newPlugInLock, 
    rosic::TrackMeter *trackMeterToWrap = nullptr);

  virtual ~TrackMeterAudioModule();

  AudioModuleEditor* createEditor(int type) override;

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setSampleRate(double newSampleRate) override
  {
    if(wrappedTrackMeter != NULL)
      wrappedTrackMeter->setSampleRate(newSampleRate);
  }

  /*
  virtual void getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    wrappedTrackMeter->measureSampleFrameStereo(inOutL, inOutR);
  }
  */

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    jassert(numChannels == 2); if(numChannels != 2) return;
    for(int n = 0; n < numSamples; n++)
      wrappedTrackMeter->measureSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
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

  // Inquiry:

  /** Returns true, iff the meter is shown in vertical mode. The meter switches automatically 
  between vertical and horizontal display modes based on the apsect ratio of the GUI editor. */
  bool isVertical() const { return getHeight() >= getWidth(); }

  // Callbacks:
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

  float rangeMin, rangeMax;

  juce_UseDebuggingNewOperator;
};

#endif 
