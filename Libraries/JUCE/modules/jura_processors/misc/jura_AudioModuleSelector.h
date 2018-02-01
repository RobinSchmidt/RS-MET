#ifndef jura_AudioModuleSelector_h
#define jura_AudioModuleSelector_h

/** A widget class for selecting a specific type of AudioModule. This can be used for dynamically
creating ("plugging in") AudioModules (like in ToolChain). */

class JUCE_API AudioModuleSelector : public RComboBox
{

public:

  /** Constructor. You must pass a pointer to an AudioModuleFactory object that will be used to 
  populate the popup menu with the AudioModule types that can be created. */
  AudioModuleSelector(AudioModuleFactory* factoryToUse);

  /** Sets the AudioModuleFactor object to be used. That determines the content of the dropdown 
  menu. */
  void setAudioModuleFactory(AudioModuleFactory* newFactory);

  /** Should be set true when the selector corresponds to the active slot. */
  void drawHighlighted(bool shouldBeHighlighted);

  virtual void paint(Graphics& g) override;

protected:

  bool highlighted = false;

  AudioModuleFactory* moduleFactory = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleSelector)
};


#endif 