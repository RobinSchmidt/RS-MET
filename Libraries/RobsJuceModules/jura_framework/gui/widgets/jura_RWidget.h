#ifndef jura_RWidget_h
#define jura_RWidget_h

/** This class serves as base class for various GUI widgets.

\todo: maybe make a class RWidgetObserver and get rid of all the custom observer/listener classes
*/

class ColourSchemeComponent;

class JUCE_API RWidget : public DescribedComponent, public ParameterObserver/*, public AsyncUpdater*/
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RWidget(const juce::String& newDescription = juce::String(""));

  /** Destructor. */
  virtual ~RWidget();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Assigns a Parameter object for observation and manipulation. */
  virtual void assignParameter(Parameter* parameterToAssign);

  /** Sets up the colour-scheme. */
  virtual void setColourScheme(const WidgetColourScheme& newColourScheme);

  /** Sets up the colour-scheme from an XmlElement. */
  virtual void setColourSchemeFromXml(const XmlElement& xml);

  /** Adds a child-widget to this one. */
  virtual void addChildWidget(RWidget *newChild, bool addAsChildComponent = true,
    bool makeVisible = true);

  /** Removes the widget from our childWidgets array and optinally also removes it as
  child-component and deletes it the object. */
  virtual void removeChildWidget(RWidget* widgetToRemove, bool removeAsChildComponent,
    bool deleteObject);

  /** Returns a pointer the ColourSchemeComponent that contains this RWidget (if any). This works
  also recursively in the sense that when the parent is not a ColourSchemeComponent, it will
  investigate the grand-parent and possibly return it (if it is a ColourSchemeComponent) ...and so
  on. If it reaches the root of the component embedding hierarchy and still has not found a
  ColourSchemeComponent, it will return NULL. */
  virtual ColourSchemeComponent* getOutlyingColourSchemeComponent();

  /** Tell the widget that it should display itself without an outline or background - that is, it
  just draws its text on whatever background is already there. */
  virtual void setNoBackgroundAndOutline(bool shouldBeWithout);

  /** Overrides Component::enablementChanged() in order to update our alphaMultiplier and trigger a
  repaint. */
  virtual void enablementChanged() override;

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns a pointer to the assigned parameter (nullptr, if none is assigned). */
  virtual Parameter* getAssignedParameter() const { return assignedParameter; }

  /** Returns true when this widget has a parameter assigned to it, false otherwise. */
  inline bool hasAssignedParameter() { return getAssignedParameter() != nullptr; }

  virtual Colour getBackgroundColour()      const { return colourScheme.background.withMultipliedAlpha(alphaMultiplier); };
  virtual Colour getOutlineColour()         const { return colourScheme.outline.withMultipliedAlpha(alphaMultiplier); };
  virtual Colour getHandleColour()          const { return colourScheme.handle.withMultipliedAlpha(alphaMultiplier); };
  virtual Colour getTextColour()            const { return colourScheme.text.withMultipliedAlpha(alphaMultiplier); };
  virtual Colour getWeakHighlightColour()   const { return colourScheme.weakHighlight.withMultipliedAlpha(alphaMultiplier); };
  virtual Colour getStrongHighlightColour() const { return colourScheme.strongHighlight.withMultipliedAlpha(alphaMultiplier); };


  //virtual Colour getTextColour()       const { return colourScheme.text;                                            };
  virtual Colour getSpecialColour1()   const { return colourScheme.special.withMultipliedAlpha(alphaMultiplier); };
  //virtual Colour getSpecialColour2()   const { return colourScheme.specialColour2;   };

  virtual WidgetColourScheme getColourScheme() const { return colourScheme; }

  /** Returns the state of the widget as a string. Subclasses may override this in order to create
  a string from their internal variables. They should then also override setValueFromString to 
  restore their internal variables from a string by parsing it. */
  //virtual juce::String getStateAsString() const = 0;
  virtual juce::String getStateAsString() const { return juce::String(); }

  /** Returns a string representing the value for gui/display purposes using our string conversion 
  function on the value of the assigned parameter (if any, otherwise the empty string will be 
  returned). */
  virtual juce::String getValueDisplayString() const
  {
    if(assignedParameter != nullptr)
      return stringConversionFunction(assignedParameter->getValue());
    else
      return "";
  }

  /** Returns a pointer to the (bitmap) font that is used for drawing text on widgets. */
  inline const BitmapFont* getFont() { return font; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides the purely virtual parameterChanged() method inherited from ParameterObserver. */
  virtual void parameterChanged(Parameter* p) override;

  /** Overrides the virtual parameterRangeChanged() method inherited from ParameterObserver. */
  virtual void parameterRangeChanged(Parameter* p) override;

  /** Overrides the purely virtual method of the ParameterObserver base class in order to
  invalidate our pointer-member
  'assignedParameter'. */
  virtual void parameterWillBeDeleted(Parameter* p) override;

  /** Overrides the changeListenerCallback in order to receive messages which this object sends to
  itself. */
  //virtual void changeListenerCallback (void *objectThatHasChanged);

  /** Overrides handleAsyncUpdate to call updateWidgetFromAssignedParameter from there. */
  //virtual void handleAsyncUpdate();

  /** This method is called when the assigned Parameter has been changed - override it in the
  subclasses to do the actual GUI update. */
  virtual void updateWidgetFromAssignedParameter(bool sendChangeMessage = false);

  /** Sets the state of the widget from the passed string. Subclasses may override this in order to
  interpret the string in different ways (in buttons as bool, in sliders as double, etc.). */
  //virtual void setStateFromString(const juce::String &valueString, bool sendChangeMessage = true) = 0;
  virtual void setStateFromString(const juce::String& /*valueString*/,
  bool sendChangeMessage = true) {}
    // maybe make purely virtual

  /** Paints the widget. The baseclass implementations just fills a rectangle with the background
  color and draw an outline with the outline color. */
  virtual void paint(Graphics& g);

  //virtual double getValue();
  // \todo: virtual void grayOutIfDisabled(Graphics& g);

  /** Thickness for the outline of the widget (in pixels). */
  static const int outlineThickness = 2;

protected:

  /** Opens a modal field for manually entering a value and returns the value entered. */
  virtual double openModalNumberEntryField(double numberToShowInitially);

  /** The assigned Parameter object. Typically, a widget is for displaying and manipulating a 
  parameter object (but it can be used also without attaching a parameter). */
  Parameter* assignedParameter = nullptr;

  /** A pointer to the function which converts a value into a juce::String. */
  juce::String (*stringConversionFunction) (double valueToConvert);

  juce::Array<RWidget*> childWidgets;

  bool  noBackgroundAndOutline;
  float alphaMultiplier;

  static const BitmapFont* font;

private:

  // colours for various aspects of the widget:
  WidgetColourScheme colourScheme;
    // \TODO (IMPORTANT): use a pointer such that a number of widgets can share the colour-scheme
    // - if NULL we may use a global default colorscheme object

  friend class rsAutomatableWidget;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RWidget)
};

#endif
