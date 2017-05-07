#ifndef jura_PlotPreviewButton_h
#define jura_PlotPreviewButton_h

/** This class is a button that can be used to switch between various editors and preview some plot
from the editor. Like when there's a filter editor in a synth, the button could preview the 
frequency response plot. */

class PlotPreviewButton : public RButton
{

public:

  PlotPreviewButton(const juce::String& name, const CoordinateSystemOld* plotToPreview = NULL);

  virtual ~PlotPreviewButton();

  virtual void paint(Graphics &g);

protected:

  Image* plotPreviewImage;

  juce_UseDebuggingNewOperator;
};


#endif  
