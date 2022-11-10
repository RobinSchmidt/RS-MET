#ifndef jura_ColourSchemeComponent_h
#define jura_ColourSchemeComponent_h

//#include "../widgets/rojue_WidgetSet.h"
//#include "../coordinate_systems/rojue_CoordinateSystemOld.h"

/** This class serves as baseclass for GUI-components that need support of colour-schemes. It
provides functions to set up colours for various things (backgrounds, widgets, etc) and a means to 
pass this info through to child-component of class ColourSchemeComponent. */


class RWidget;
class WidgetSet;
class rsPlot;
// maybe we can get rid of them at some point

class JUCE_API ColourSchemeComponent : public Component, virtual public DescribedItem
{

public:

  enum defaultColourSchemes
  {
    BLUE,
    PURPLE,
    RED,

    NUM_DEFAULT_COLOUR_SCHEMES
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ColourSchemeComponent(const juce::String& newColourSchemeComponentName
    = juce::String("ColourSchemeComponent"));

  /** Destructor. */
  virtual ~ColourSchemeComponent();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the juce::Label in which the descriptions for the widgets will appear. */
  virtual void setDescriptionField(RTextField* newDescriptionField, 
    bool callAlsoForChildColourSchemeComponents);

  /** Selects the general appearance for the editor as one of ColourScheme::appearances. */
  virtual void setEditorAppearance(int newAppearance);

  /** Selects the general appearance for the widgets as one of ColourScheme::appearances. */
  virtual void setWidgetAppearance(int newAppearance);

  /** Selects the general appearance for the plots as one of ColourScheme::appearances. */
  virtual void setPlotAppearance(int newAppearance);

  /** Sets the central hue for the colours. */
  virtual void setCentralHue(float newHue);

  /** Sets an offset for a hue which can be used to use different hues for different parts of the 
  editor. */
  virtual void setHueOffset(int index, float newOffset);

  /** Sets the overall multiplier for the saturations of the colours. */
  virtual void setSaturationMultiplier(float newMultiplier);

  /** Sets the gamma value for the luminances of the colours. */
  virtual void setBrightnessGamma(float newGamma);

  /** Sets up the graph colours for the plots and optionally passes this through to children. */
  virtual void setGraphColours(const juce::Array<Colour> &newColours, 
    bool callAlsoForChildColourSchemeComponents = true);

  /** Copies the colourscheme settings from the passed component to this one. */
  virtual void copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom);

  /** Sets the visibility of all widgets. This is useful for adapting the visibility of widgets 
  according to the state - we can make them all invisible at once and selectively make some of them
  visible again. */
  virtual void setWidgetsVisible(bool shouldBeVisible);


  virtual void setColourScheme(const WidgetColourScheme& newColourScheme); 
  //virtual void setColourScheme(const WidgetColourScheme& newColourScheme);
    // i think, this needs to be added

  /** Sets up this colorscheme from an XmlElement. */
  virtual void setColourSchemeFromXml(const XmlElement& xmlColorScheme);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the EditorColourScheme of this ColourSchemeComponent (which defines the general 
  appearance of this component). */
  virtual const EditorColourScheme& getEditorColourScheme() const { return editorColourScheme; }

  /** Returns the WidgetColourScheme of this ColourSchemeComponent (which defines the appearance of 
  the widgets in this component). */
  virtual const WidgetColourScheme& getWidgetColourScheme() const { return widgetColourScheme; }

  /** Returns the PlotColourScheme of this ColourSchemeComponent (which defines the appearance of 
  the function-plots and visual editors in this component). */
  virtual const PlotColourScheme& getPlotColourScheme() const { return plotColourScheme; }

  /** Returns this colorscheme as XmlElement - the caller is responsible for deleting it. */
  virtual XmlElement* getColourSchemeAsXml();

  /** Returns the colour used for the outline. */
  virtual Colour getOutlineColour() const { return editorColourScheme.outline; }

  /** Returns the colour used for the top-left corner. */
  virtual Colour getTopLeftColour() const { return editorColourScheme.topLeft; }

  /** Returns the color used for regular text on the editor. */
  virtual Colour getTextColour() const { return editorColourScheme.text; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void paint(Graphics &g);
  virtual void paintOverChildren(Graphics &g);

  /** Updates the widgets according to the state of the assignedParameter (if any). */
  virtual void updateWidgetsAccordingToState();

protected:

  /** Adds a widget to our array (and optionally also adds the widget as child-component to this 
  one). Having the widget added by this method rather than directly using Component's 
  addChildComponent allows for looping through all the widgets to set up their colours, 
  description-fields etc. */
  virtual void addWidget(RWidget* widgetToAdd, bool addAsChildComponent = true, 
    bool makeVisible = true);

  /** Adds a set of widgets */
  virtual void addWidgetSet(WidgetSet* widgetSetToAdd, bool addAsChildComponent = true, 
    bool makeVisible = true);

  /** Adds a CoordinateSystem (most often used for function-plots, etc.) to our array (and 
  optionally also adds the widget as child-component to this one). Similar to addWidget(). */
  virtual void addPlot(rsPlot* plotToAdd, bool addAsChildComponent = true, 
    bool makeVisible = true);

  /** Adds a child-colourschem-component to this one. Using this method rather than Component's 
  addChildComponent allows for looping through all the sub-editors widgets to set up their colours, 
  description-fields etc. \todo: get rid of the bool flags */
  virtual void addChildColourSchemeComponent(ColourSchemeComponent* editorToAdd, 
    bool addAsChildComponent = true, bool makeVisible = true);

  /** Removes a child-editor from this one and optionally deletes the object */
  virtual void removeChildColourSchemeComponent(ColourSchemeComponent *childToRemove, 
    bool deleteObject);

  /** Removes the widget from our widgets array and optinally also removes it as child-component 
  and deletes it the object. */
  virtual void removeWidget(RWidget* widgetToRemove, bool removeAsChildComponent, 
    bool deleteObject);

  /** Passes the current colour-schemes to the embedded widgets, plots, child-editors, etc. and 
  triggers a repaint. */
  virtual void updateEmbeddedObjectsAndRepaint();

  // colour-schemes for various aspects of components:
  WidgetColourScheme widgetColourScheme;
  EditorColourScheme editorColourScheme;
  PlotColourScheme   plotColourScheme;
  // ToDo: Use pointers or referencese to allow sharing these among sub-editors. This requires to
  // elect some higher/highest level Component to be the maintainer of the actual objects which
  // will need some changes that ripple through the whole framework
  // ...but the mid-term plan is to replace this by a class "rsSkinnableComponent" soon. This 
  // should then just hold a pointer to an rsSkin object, the skin may then also countain 
  // color-schemes (among other things such as fonts, graphics, etc.). When done, this class 
  // should be made obsolete

  // we maintain components that should adhere to the colour-schemes above as arrays:
  std::vector<RWidget*>   widgets;
  std::vector<WidgetSet*> widgetSets;
  std::vector<rsPlot*>    plots;
  std::vector<ColourSchemeComponent*> childColourSchemeComponents;
  CriticalSection arrayLock; // used to access all the arrays
  // use std::vector, maybe have also arrays of non-owned widgets, etc. in order to allow client 
  // code to use non-pointer variables ...maybe use references

  bool drawWithEnclosingRectangle;

  juce_UseDebuggingNewOperator;
};

#endif  
