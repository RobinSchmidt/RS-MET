#ifndef jura_ColorMapLoader_h
#define jura_ColorMapLoader_h


/**   */

class JUCE_API ColorMapLoader : public StateLoadSaveWidgetSet
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ColorMapLoader(const juce::String& newName = juce::String("ColorMapLoader"));




  //-----------------------------------------------------------------------------------------------
  // setup:


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ColorMapLoader)
};

#endif  