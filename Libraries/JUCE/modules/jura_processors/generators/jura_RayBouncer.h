#ifndef jura_RayBouncer_h
#define jura_RayBouncer_h
  
//=================================================================================================

/**   */

class JUCE_API RayBouncerAudioModule : public jura::AudioModule
{

public:

  RayBouncerAudioModule(CriticalSection *lockToUse);
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createParameters();

  /** Creates the GUI editor (returns an object of an appropriate subclass of AudioModuleEditor) */
  //AudioModuleEditor *createEditor() override;

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;

  // target functions for callbacks that are called on parameter changes:
  void setFrequency(double newFrequency);
  void setEllipseSize(double newSize);
  void setEllipseAspectRatio(double newRatio);
  void setEllipseAngleDegrees(double newAngle);
  void setEllipseCenterX(double newX);
  void setEllipseCenterY(double newY);
  void setStartX(double newX);
  void setStartY(double newY);
  void setLaunchAngle(double newAngle);

  void setBending( double newValue);
  void setBendX2Y( double newValue);
  void setBendY2X( double newValue);
  void setBendXX2X(double newValue);
  void setBendXX2Y(double newValue);
  void setBendXY2X(double newValue);
  void setBendXY2Y(double newValue);
  void setBendYY2X(double newValue);
  void setBendYY2Y(double newValue);



  void setAutoReset(bool shouldReset);

protected:

  /** Resets the embedded rsRayBouncer object if the flag autoReset is true. to be called from 
  parameter setters because the sound strongly depends on intial conditions, so it is advisable
  to automatically reset on parameter change during sound desgin. */
  void autoResetIfDesired();

  // embedded core DSP objects from the RAPT library:
  //RAPT::rsRayBouncerDriver<double> rayBouncer;
  rsRayBouncerDriverD rayBouncer;
  //rsRayBouncerDriverD rayBouncer;

  double frequency = 100, sampleRate = 44100;

  bool autoReset = false;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RayBouncerAudioModule)
};

//=================================================================================================

/** This is the GUI editor class for the jura::RayBouncerAudioModule. */
/*
class JUCE_API RayBouncerAudioModuleEditor : public AudioModuleEditor
{

public:

  RayBouncerAudioModuleEditor(jura::RayBouncerAudioModule *newRayBouncerToEdit);
  virtual void resized() override;

protected:

  RayBouncerAudioModule *rayBouncerToEdit;

  AutomatableSlider *frequencySlider;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RayBouncerAudioModuleEditor)
};
*/

#endif 