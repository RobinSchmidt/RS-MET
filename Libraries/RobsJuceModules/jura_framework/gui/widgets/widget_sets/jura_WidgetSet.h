#ifndef jura_WidgetSet_h
#define jura_WidgetSet_h

//#include "rojue_RWidget.h"

/** This class assembles some widgets int a set that can be treated as one entity then.

\todo: we should store the pointer to the description field here and pass it to widgets 
when they are added

Why do we need this call at all? the baseclass ColourSchemeComponent alredady does have
a widgets array. Seems like we are shadowing it here? Why? Does that make any sense? If so, 
document it, if not, try to get rid of this class

*/


class JUCE_API WidgetSet : public ColourSchemeComponent
{
  // Temporary, during trying to get rid of the WidgetSet class

  // Empty overrides for paint/OverChildren to avoid background and outline drawing:
  void paint(Graphics &g) override {}
  void paintOverChildren(Graphics &g) override {}

};

// OK - this was the old implementation, used until 2022/11/03. It's not used anymore and can 
// probably e deleted soon
class JUCE_API WidgetSetOld : public ColourSchemeComponent
{

  friend class Editor;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WidgetSetOld();

  /** Destructor. */
  virtual ~WidgetSetOld();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the colour-scheme. */
  virtual void setColourScheme(const WidgetColourScheme& newColourScheme); /* override;*/
    // should this not be an override - also available in ColourSchemeComponent?

  /** Sets up the colour-scheme from an XmlElement. */
  //void setColourSchemeFromXml(const XmlElement* widgetColours) override;
  virtual void setColourSchemeFromXml(const XmlElement* widgetColours); /* override; ?*/
  // hmm - the baseclass method takes a reference, we take a pointe here - why?
  //   virtual void setColourSchemeFromXml(const XmlElement& xmlColorScheme);

  /** Sets the juce::Label in which the descriptions for the widgets will appear. */
  void setDescriptionField(RTextField* newDescriptionField) override;

  // Empty overrides for paint/OverChildren to avoid background and outline drawing:
  void paint(Graphics &g) override {}
  void paintOverChildren(Graphics &g) override {}


protected:

  /** Adds a widget to our array (and optionally also adds the widget as child-component to this 
  one). Having the widget added by this method rather than directly using Component's 
  addChildComponent allows for looping through all the widgets to set up their colours, 
  description-fields etc. */
  void addWidget(RWidget* widgetToAdd, bool addAsChildComponent = true, 
    bool makeVisible = true) override;

  /** Analog to addWidget() but specifically for labels - they are kept in a separate array because 
  they should be treated separately. */
  //virtual void addWidget(RLabel* labelToAdd, bool addAsChildComponent = true, 
  // bool makeVisible = true);

  /** An array of all the widgets that have been added by addWidget(). */
  juce::Array<RWidget*, CriticalSection> widgets;

  /** An array of all the labels that have been added by addWidget(). */
  //juce::Array<RLabel*, CriticalSection> labels;

  juce_UseDebuggingNewOperator;
};


#endif  