#ifndef jura_PreDefinedParameters_h
#define jura_PreDefinedParameters_h

/** This file contains a couple of frequently used parameters such that one doesn't have to
intialize all the values, default values etc. everytime one needs a parameter of the kind.

todo: get rid of these - instead make initialization functions for various kinds of parameters
such that they work with all subclasses of Parameter - see ParameterProfiles .h/cpp

*/


class ParameterPowersOfTwo : public Parameter
{
public:

  ParameterPowersOfTwo(CriticalSection *criticalSectionToUse, const juce::String& newName,
    double newMinValue     = 256.0,
    double newMaxValue     = 8192.0,
    double newInterval     = 1.0,
    double newDefaultValue = 1024.0)
    : Parameter(criticalSectionToUse, newName, newMinValue, newMaxValue, newInterval, newDefaultValue, STRING)
  {
    jassert(isPowerOfTwo((unsigned int)newMinValue));
    jassert(isPowerOfTwo((unsigned int)newMaxValue));  // min- and maxValue must be powers of two
    jassert(newMinValue <= newMaxValue);

    int currentValue = (int)newMinValue;
    while(currentValue <= newMaxValue)
    {
      addStringValue(juce::String(currentValue));
      currentValue *= 2;
    }

    interval     = 1;
    defaultValue = log2(newDefaultValue/newMinValue);
    value        = defaultValue;
  }

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterPowersOfTwo)
};

class ParameterTwoPoleFilterMode : public ModulatableParameter
{
public:
  ParameterTwoPoleFilterMode(CriticalSection *criticalSectionToUse,
    const juce::String& newName   = juce::String("Mode"),
    double newLowerLimit          = 0.0,
    double newUpperLimit          = 9.0,
    double newInterval            = 1.0,
    double newDefaultValue        = 0.0,
    Parameter::Scaling newScaling = Parameter::Scaling::STRING,
    int /*newDefaultMidiController*/ = -1,
    bool /*newSaveAndRecall*/        = true) // change parameter order, remove CriticalSection
    : ModulatableParameter(newName, newLowerLimit, newUpperLimit, newDefaultValue, 
      newScaling, newInterval)
  {
    addStringValue("Bypass");
    addStringValue("Peak/Dip");
    addStringValue("Low Shelving");
    addStringValue("High Shelving");
    addStringValue("Lowpass 6 dB/oct");
    addStringValue("Lowpass 12 dB/oct");
    addStringValue("Highpass 6 dB/oct");
    addStringValue("Highpass 12 dB/oct");
    addStringValue("Notch 2*6 dB/oct");
    addStringValue("Bandpass 2*6 dB/oct");
    addStringValue("Pole");
    addStringValue("Zero");
    addStringValue("Pole Pair");
    addStringValue("Zero Pair");

    setValue(0.0, false, false); // init to Bypass
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterTwoPoleFilterMode)
};

class ParameterFourPoleFilterMode : public ModulatableParameter
{
public:
  ParameterFourPoleFilterMode(CriticalSection *criticalSectionToUse,
    const juce::String& newName   = juce::String("Mode"),
    double newLowerLimit          = 0.0,
    double newUpperLimit          = 9.0,
    double newInterval            = 1.0,
    double newDefaultValue        = 0.0,
    Parameter::Scaling newScaling = Parameter::Scaling::STRING,
    int /*newDefaultMidiController*/ = -1,
    bool /*newSaveAndRecall*/        = true)
    : ModulatableParameter(newName, newLowerLimit, newUpperLimit, newDefaultValue, newScaling, 
      newInterval)
  {
    addStringValue("Bypass");
    addStringValue("Lowpass 6 dB/oct");
    addStringValue("Lowpass 12 dB/oct");
    addStringValue("Lowpass 18 dB/oct");
    addStringValue("Lowpass 24 dB/oct");

    addStringValue("Highpass 6 dB/oct");
    addStringValue("Highpass 12 dB/oct");
    addStringValue("Highpass 18 dB/oct");
    addStringValue("Highpass 24 dB/oct");

    addStringValue("Bandpass 6+6 dB/oct");
    addStringValue("Bandpass 12+12 dB/oct");
    addStringValue("Bandpass 6+12 dB/oct");
    addStringValue("Bandpass 12+6 dB/oct");
    addStringValue("Bandpass 6+18 dB/oct");
    addStringValue("Bandpass 18+6 dB/oct");

    addStringValue("Notch (2nd Order)");
    addStringValue("Notch (4th Order)");
    addStringValue("Two Notches");

    addStringValue("Peak/Dip (2nd Order)");
    addStringValue("Peak/Dip (4th Order)");
    addStringValue("Peak/Dip (Flat Top)");
    addStringValue("Two Peaks");

    addStringValue("Low Shelving 1st Order");
    addStringValue("Low Shelving 2nd Order");
    //addStringValue("Low Shelving 3rd order");
    addStringValue("Low Shelving 4th Order");

    addStringValue("High Shelving 1st Order");
    addStringValue("High Shelving 2nd Order");
    addStringValue("High Shelving 4th Order");

    setValue(0.0, false, false); // init to Bypass
  }

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterFourPoleFilterMode)
};

//=================================================================================================

/** For grids like in the node editor for the BreakpointModulator */ 

class ParameterGridInterval : public Parameter
{
public:

  ParameterGridInterval(const juce::String& name);

  virtual void setStringValue(const juce::String& newString, bool sendNotification,
    bool callCallbacks) override;

  virtual double getValue() const override;

  static const std::vector<String> gridIntervalStringArray;
  static const std::vector<double> gridIntervalValueArray;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterGridInterval)
};





/** A parameter that has attached key- and velocity scaling parameters. 
...this class is not usable yet. i think, it should have pointer parameters for keyParam and 
velParam and these need to be added to the AudioModule (pointers, because the AudioModule will 
delete them in its destructor) */

class ParameterWithKeyVelScaling : public Parameter
{

public:

  ParameterWithKeyVelScaling(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, Parameter::Scaling scaling = Parameter::Scaling::LINEAR, 
    double interval = 0.0);

  void setKeyScaleRange(double minValue, double maxValue, double defaultValue);
  void setVelScaleRange(double minValue, double maxValue, double defaultValue);

  void setKeyScaleCallback();
  void setVelScaleCallback();

  Parameter* getKeyScalingParameter() { return &keyParam; }
  Parameter* getVelScalingParameter() { return &velParam; }


protected:

  Parameter keyParam, velParam;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterWithKeyVelScaling)
};

#endif
