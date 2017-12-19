#include "PluginProcessor.h"
#include "PluginEditor.h"

using namespace RAPT;

TestPluginAudioProcessorEditor::TestPluginAudioProcessorEditor (TestPluginAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
  labelMode.setText("Mode", NotificationType::dontSendNotification);
  labelMode.setJustificationType(Justification::centredLeft);
  addAndMakeVisible(&labelMode);

  typedef rsLadderFilter<float, float> Ladder; // for accessing the enum of modes conveniently
  boxMode.addItem("Flat",                  Ladder::FLAT+1);    // +1 because the ID cannot be 0, when
  boxMode.addItem("Lowpass 6 dB/oct",      Ladder::LP_6+1);    // retieving the ID, we'll use the
  boxMode.addItem("Lowpass 12 dB/oct",     Ladder::LP_12+1);   // returned value -1
  boxMode.addItem("Lowpass 18 dB/oct",     Ladder::LP_18+1);
  boxMode.addItem("Lowpass 24 dB/oct",     Ladder::LP_24+1);
  boxMode.addItem("Highpass 6 dB/oct",     Ladder::HP_6+1);
  boxMode.addItem("Highpass 12 dB/oct",    Ladder::HP_12+1);
  boxMode.addItem("Highpass 18 dB/oct",    Ladder::HP_18+1);
  boxMode.addItem("Highpass 24 dB/oct",    Ladder::HP_24+1);
  boxMode.addItem("Bandpass 6/6 dB/oct",   Ladder::BP_6_6+1);
  boxMode.addItem("Bandpass 6/12 dB/oct",  Ladder::BP_6_12+1);
  boxMode.addItem("Bandpass 6/18 dB/oct",  Ladder::BP_6_18+1);
  boxMode.addItem("Bandpass 12/6 dB/oct",  Ladder::BP_12_6+1);
  boxMode.addItem("Bandpass 12/12 dB/oct", Ladder::BP_12_12+1);
  boxMode.addItem("Bandpass 18/6 dB/oct",  Ladder::BP_18_6+1);
  //// experimental modes, harvested from KVR thread:
  //boxMode.addItem("KVR_BP2",  Ladder::KVR_BP2+1);
  //boxMode.addItem("KVR_BP4",  Ladder::KVR_BP4+1);
  //boxMode.addItem("KVR_NF2",  Ladder::KVR_NF2+1);
  //boxMode.addItem("KVR_NF4",  Ladder::KVR_NF4+1);
  //boxMode.addItem("KVR_PF2",  Ladder::KVR_PF2+1);
  //boxMode.addItem("KVR_PF4",  Ladder::KVR_PF4+1);
  boxMode.addListener(this);
  boxMode.setSelectedId(Ladder::LP_24+1);
  addAndMakeVisible(&boxMode);

  labelCutoff.setText("Cutoff", NotificationType::dontSendNotification);
  labelCutoff.setJustificationType(Justification::centredLeft);
  addAndMakeVisible(&labelCutoff);

  sliderCutoff.setSliderStyle(Slider::LinearBar);
  sliderCutoff.setRange(0.0, 1.0, 0.001); // we need conversion here
  sliderCutoff.setTextBoxStyle(Slider::TextBoxRight, false, 40, 20);
  sliderCutoff.addListener(this);
  sliderCutoff.setValue(1.0);
  addAndMakeVisible(&sliderCutoff);

  labelReso.setText("Reso", NotificationType::dontSendNotification);
  labelReso.setJustificationType(Justification::centredLeft);
  addAndMakeVisible(&labelReso);

  sliderReso.setSliderStyle(Slider::LinearBar);
  sliderReso.setRange(0.0, 1.0, 0.001);
  sliderReso.setTextBoxStyle(Slider::TextBoxRight, false, 40, 20);
  sliderReso.addListener(this);
  sliderReso.setValue(0.0);
  addAndMakeVisible(&sliderReso);

  setResizeLimits(250, 100, 500, 200);
  setSize(250, 100);
}

void TestPluginAudioProcessorEditor::paint(Graphics& g)
{
  g.fillAll(Colours::grey);
}

void TestPluginAudioProcessorEditor::resized()
{
  // Lay out positions of subcomponents (widgets, etc.) in the editor:

  int margin = 10;        // margin for widgets with respect to editor's edges
  int thickness = 20;     // height ("thickness") of the widgets
  int labelWidth = 50;    // width of name-label in front of sliders
  int spacing = 5;        // spacing between widgets

  int x = margin;
  int y = margin;
  int deltaY = thickness + spacing;
  int sliderWidth = getWidth() - labelWidth - spacing - 2*margin;
  int sliderX = x + labelWidth + spacing;


  labelMode.setBounds(   x,       y, labelWidth,  thickness);
  boxMode.setBounds(     sliderX, y, sliderWidth, thickness);
  y += deltaY;
  labelCutoff.setBounds( x,       y, labelWidth,  thickness);
  sliderCutoff.setBounds(sliderX, y, sliderWidth, thickness);
  y += deltaY;
  labelReso.setBounds(   x,       y, labelWidth,  thickness);
  sliderReso.setBounds(  sliderX, y, sliderWidth, thickness);
}

void TestPluginAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
  ScopedLock(processor.getCallbackLock());

  if(slider == &sliderCutoff)
    processor.setParameter(TestPluginAudioProcessor::CUTOFF, (float)sliderCutoff.getValue());
  else if(slider == &sliderReso)
    processor.setParameter(TestPluginAudioProcessor::RESO, (float)sliderReso.getValue());

  // code below triggers a breakpoint - we probably need to add the parameters to an array of 
  // automatable parameters first - otherwise we get some index-out-of-range error:
  //if(slider == &sliderCutoff)
  //  processor.setParameterNotifyingHost(TestPluginAudioProcessor::CUTOFF, (float)sliderCutoff.getValue());
  //else if(slider == &sliderReso)
  //  processor.setParameterNotifyingHost(TestPluginAudioProcessor::RESO, (float)sliderReso.getValue());
}

void TestPluginAudioProcessorEditor::comboBoxChanged(ComboBox *box)
{
  ScopedLock(processor.getCallbackLock());

  if(box == &boxMode)
    processor.setFilterMode(boxMode.getSelectedId()-1);
}
