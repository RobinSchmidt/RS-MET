#ifndef jura_DescribedComponent_h
#define jura_DescribedComponent_h

class RTextField;

/** This class is used as baseclass for items that have a description text. This applies mainly to 
GUI elements. It enables tooltips and description text in a status line. 


ToDo:

- Maybe this class is now obsolete because juce::Component itself has now also a setDescription
  function. Maybe it's intention is the same so we could potentially just use that. Maybe I'm
  now just replicating functionality that is now included in JUCE itself. Figure out!

*/

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
  virtual void setItemDescription(const juce::String &newDescription);

  /** Sets the juce::Label in which the description will appear. */
  virtual void setDescriptionField(RTextField* newDescriptionField);
  // Maybe rename to setItemDescriptionField

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the description for this widget. */
  virtual juce::String getItemDescription() const;

  /** Returns the description field for this item (or NULL if none). */
  virtual RTextField* getDescriptionField() const;

protected:

  /** A description of the parameter. */
  juce::String description;  // Maybe rename to itemDescription

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
  virtual void mouseEnter(const MouseEvent &e) override;

  /** Overrides the mouseExit callback in order to make the description disappear when the mouse 
  leaves the widget. */
  virtual void mouseExit(const MouseEvent &e) override;

  /** If this is called from the message thread, it just calls repaint. If it is called from some 
  other thread, it will schedule an asynchronous call to repaint on the message thread. Sometimes, 
  repaints need to be triggered from the audio thread (for example when a plugin receives midi 
  controllers or automation and should repaint a slider), but triggering repaint from the audio 
  thread (or any other thread than the message thread) isn't allowed in juce. */
  void repaintOnMessageThread();
  // Maybe turn this into a free function taking a pointer to a juce::Component

  juce_UseDebuggingNewOperator;
};

#endif  