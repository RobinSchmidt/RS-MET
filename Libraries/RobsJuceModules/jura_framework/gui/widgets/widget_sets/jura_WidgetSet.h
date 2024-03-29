#ifndef jura_WidgetSet_h
#define jura_WidgetSet_h

//#include "rojue_RWidget.h"

/** This class assembles some widgets int a set that can be treated as one entity then.

\todo: we should store the pointer to the description field here and pass it to widgets 
when they are added

*/

class JUCE_API WidgetSet : public ColourSchemeComponent
{

  friend class Editor;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  WidgetSet();

  /** Destructor. */
  virtual ~WidgetSet();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the colour-scheme. */
  virtual void setColourScheme(const WidgetColourScheme& newColourScheme); /* override;*/
    // should this not be an override - also available in ColourSchemeComponent?

  /** Sets up the colour-scheme from an XmlElement. */
  virtual void setColourSchemeFromXml(const XmlElement* widgetColours); /* override; ?*/
    //

  /** Sets the juce::Label in which the descriptions for the widgets will appear. */
  virtual void setDescriptionField(RTextField* newDescriptionField) override;

  // exmpty overrides ofr paint/OverChildren to avoid background and outline drawing
  virtual void paint(Graphics &g) override {}
  virtual void paintOverChildren(Graphics &g) override {}


protected:

  /** Adds a widget to our array (and optionally also adds the widget as child-component to this 
  one). Having the widget added by this method rather than directly using Component's 
  addChildComponent allows for looping through all the widgets to set up their colours, 
  description-fields etc. */
  virtual void addWidget(RWidget* widgetToAdd, bool addAsChildComponent = true, 
    bool makeVisible = true);

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