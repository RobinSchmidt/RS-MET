#ifndef PLUGINEDITOR_H_INCLUDED
#define PLUGINEDITOR_H_INCLUDED

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"

/** Editor for the simple ladder filter test plugin */

class TestPluginAudioProcessorEditor  : public AudioProcessorEditor, public Slider::Listener, 
  public ComboBox::Listener
{

public:
    TestPluginAudioProcessorEditor (TestPluginAudioProcessor&);

    virtual void paint (Graphics&) override;
    virtual void resized() override;
    virtual void sliderValueChanged (Slider* slider) override;
    virtual void comboBoxChanged(ComboBox *comboBoxThatHasChanged) override;

private:

    TestPluginAudioProcessor& processor; // reference to the "owner" AudioProcessor

    Label labelMode, labelCutoff,  labelReso;
    ComboBox boxMode;
    Slider sliderCutoff, sliderReso;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TestPluginAudioProcessorEditor)
};


#endif
