#ifndef jura_AudioModuleSelector_h
#define jura_AudioModuleSelector_h

/** A widget class for selecting a specific type of AudioModule. */

class JUCE_API AudioModuleSelector : public RComboBox
{

public:

  AudioModuleSelector();

  /** Should be set true when the selector corresponds to the active slot */
  void drawHighlighted(bool shouldBeHighlighted);

  virtual void paint(Graphics& g) override;

protected:

  bool highlighted = false;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioModuleSelector)
};


#endif 