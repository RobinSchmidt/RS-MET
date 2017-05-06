#ifndef jura_CustomComboBoxes_h
#define jura_CustomComboBoxes_h

/** This file contains custom comboboxes for special purposes.  */

class FourPoleFilterModeComboBox : public RComboBox
{
public:
  FourPoleFilterModeComboBox(const juce::String& componentName) : RComboBox(componentName) {  }
  juce_UseDebuggingNewOperator;
protected:
  //virtual void openPopUp();
};

class WaveformComboBox : public RComboBox
{
public:
  WaveformComboBox(const juce::String& componentName);
  juce_UseDebuggingNewOperator;
protected:
  //virtual void openPopUp();
};





#endif  
