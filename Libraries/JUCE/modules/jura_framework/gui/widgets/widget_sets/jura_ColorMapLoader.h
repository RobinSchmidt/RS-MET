#ifndef jura_ColorMapLoader_h
#define jura_ColorMapLoader_h

/** A componenent to preview a ColorMap as a rectangle showing the color gradient from lowest
to highest value. If the component is wider than high, the gradient will be shown horizontally 
(suitable for the loader widget set), otherwise vertically (suitable for a color bar next to a 
plot). */

class JUCE_API ColorMapPreviewer : public Component
{

public:

  ColorMapPreviewer(ColorMap *mapToPreview);

  virtual void paint(Graphics& g) override;

protected:

  ColorMap *colorMap;

};

//=================================================================================================

/** A widget set for loading and previewing a ColorMap.  */

class JUCE_API ColorMapLoader : public StateLoadSaveWidgetSet
{

public:

  /** Constructor. You must pass a valid pointer to a ColorMap object. This is the object that 
  will be updated, when the user loads a new colormap xml file. */
  ColorMapLoader(ColorMap *mapToUpdate);

  virtual void resized() override;

protected:

  ColorMapPreviewer previewer;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ColorMapLoader)
};

#endif  