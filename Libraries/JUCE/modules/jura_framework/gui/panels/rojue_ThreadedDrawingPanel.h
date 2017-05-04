#ifndef rojue_ThreadedDrawingPanel_h
#define rojue_ThreadedDrawingPanel_h

#include "rojue_ThreadedDrawingComponent.h"
#include "rojue_Panel.h"

namespace rojue
{

  /**

  This class combines the functionality of Panel with the feature of redrawing itself in a separate 
  thread inherited from ThreadedDrawingPanel

  */

  class ThreadedDrawingPanel: virtual public Panel, virtual public ThreadedDrawingComponent
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ThreadedDrawingPanel(TimeSliceThread* newThreadToUse = NULL);   

    //---------------------------------------------------------------------------------------------
    // overrides:
    virtual void resized();
    virtual void setDirty(bool shouldBeSetToDirty = true);

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    virtual void paint(Graphics &g);

  };

}

#endif  
