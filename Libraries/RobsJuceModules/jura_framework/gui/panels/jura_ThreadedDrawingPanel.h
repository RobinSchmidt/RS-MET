#ifndef jura_ThreadedDrawingPanel_h
#define jura_ThreadedDrawingPanel_h

/** This class combines the functionality of Panel with the feature of redrawing itself in a
separate thread inherited from ThreadedDrawingPanel */

class ThreadedDrawingPanel : virtual public Panel, virtual public ThreadedDrawingComponent
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ThreadedDrawingPanel(TimeSliceThread* newThreadToUse = NULL);

  //-----------------------------------------------------------------------------------------------
  // overrides:
  virtual void resized();
  virtual void setDirty(bool shouldBeSetToDirty = true);

protected:

  virtual void paint(Graphics &g);

  juce_UseDebuggingNewOperator;
};


#endif  
