#ifndef jura_DescribedComponent_h
#define jura_DescribedComponent_h

class RTextField;

/** This class ...  */

class JUCE_API DescribedItem 
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  DescribedItem(const juce::String& newDescription = juce::String("some item"));   

  /** Destructor. */
  virtual ~DescribedItem(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets a description for this widget - this should be short enough to fit into the Label which 
  is assigned to the descriptionField member. @see setDescriptionField() */
  virtual void setDescription(const juce::String &newDescription);

  /** Sets the juce::Label in which the description will appear. */
  virtual void setDescriptionField(RTextField* newDescriptionField);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the description for this widget. */
  virtual juce::String getDescription() const;

  /** Returns the description field for this item (or NULL if none). */
  virtual RTextField* getDescriptionField() const;

protected:

  /** A description of the parameter. */
  juce::String description;

  /** A label where the description will appear when the mouse is over the parameter. */
  RTextField* descriptionField;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This class serves as base class for various GUI-objects that may have a description which 
can be made to appear in some dedicated Label.  */

class JUCE_API DescribedComponent : virtual public Component, public DescribedItem
{

public:

  /** Constructor. */
  DescribedComponent(const juce::String& newDescription = juce::String("some item")) 
    : DescribedItem(newDescription) {}

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides the mouseEnter callback in order to show the description in the dedicated field 
  when the mouse enters the widget. */
  virtual void mouseEnter(const MouseEvent &e);

  /** Overrides the mouseExit callback in order to make the description disappear when the mouse 
  leaves the widget. */
  virtual void mouseExit(const MouseEvent &e);

  juce_UseDebuggingNewOperator;
};

#endif  