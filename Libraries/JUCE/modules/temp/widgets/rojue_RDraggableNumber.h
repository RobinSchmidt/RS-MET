#ifndef rojue_RDraggableNumber_h
#define rojue_RDraggableNumber_h

#include "rojue_RSlider.h"

namespace rojue
{

  class RDraggableNumber : public RSlider
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    RDraggableNumber(const juce::String& componentName);

    /** Destructor. */
    virtual ~RDraggableNumber();

    //---------------------------------------------------------------------------------------------
    // others:

    /** Overrides the way in which the slider is painted. */
    virtual void paint(Graphics& g);

    /** Overrides the behaviour on mouse-down events. */
    virtual void mouseDown (const MouseEvent& e);

    //virtual void mouseUp (const MouseEvent& e);

    /** Overrides the behaviour on mouse-drag events. */
    virtual void mouseDrag (const MouseEvent& e);  

    /** Overrides the behaviour on mouse-double click events. */
    //virtual void mouseDoubleClick (const MouseEvent& e);

    //virtual void mouseWheelMove (const MouseEvent& e, float wheelIncrementX, float wheelIncrementY);

    //===============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    double valueOnMouseDown;

    RDraggableNumber (const RDraggableNumber&);
    const RDraggableNumber& operator= (const RDraggableNumber&);

  };

}

#endif   
