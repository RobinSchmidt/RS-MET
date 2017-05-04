#ifndef rojue_ThreadedDrawingComponent_h
#define rojue_ThreadedDrawingComponent_h

#include "../../misc/rojue_MessageBoxes.h"
#include "rojue_DrawingThread.h"

namespace rojue
{

  /**

  This class is a component, intended to be used as base-class for all components that do some
  heavy work when drawing onto their client area. By default, juce calls the Component's paint 
  method in the message thread - the same thread where all the user interface interactions 
  (mouse-clicks, etc.) take place. Doing heavy graphics rendering in the same thread can make the 
  user interface unresponsive. Therefore, this class provides its own infrastructure for handling 
  redrawing: it has a member clientAreaImage (pointer to image) onto which subclasses should 
  perform their drawing operations - this image will then be shown in paint method. For subclasses 
  of ThreadedDrawingComponent, it is typically sufficient to to call setDirty() after each 
  operation which invalidates the currently shown image and to override drawComponent() to perform 
  their drawing operations on the passed juce::Image. The class ThreadedDrawingComponent together 
  with the class ThreadedDrawingComponentDrawingThread will take care of the rest. The passed image 
  will be a pointer to our member clientAreaImage and this baseclass here makes sure that 
  drawComponent() is called inside a proper CriticalSection to make it thread safe.

  */

  class ThreadedDrawingComponent	: virtual public Component, public TimeSliceClient, public Timer
  {

    friend class ThreadedDrawingPanel;

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ThreadedDrawingComponent(TimeSliceThread* newThreadToUse = NULL);   

    /** Destructor. */
    virtual ~ThreadedDrawingComponent();            

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the TimeSliceThread that will we used for redrawing. */
    virtual void setDrawingThread(TimeSliceThread* newThreadToUse);

    //---------------------------------------------------------------------------------------------
    // callbacks:

    /** Responds to a timer-callback by triggering a repaint - the ThreadedDrawingComponent will 
    itself start the timer when it has finished rendering. We use this indirect way (instead of 
    just calling repaint after rendering) because juce doesn't like repaint to be called from 
    threads other than the message thread. */
    virtual void timerCallback();

    /** Overrides the resized()-function of the component base-class. */
    virtual void resized();

    /** Implements the purely virtual method of the TimeSliceClient baseclass in order to render 
    the clientAreaImage memeber inside this function. This callback will be called from a 
    ThreadedDrawingComponentDrawingThread to which the ThreadedDrawingComponent will register 
    itself as client upon construction. This thread exists (or will be created) as a global 
    object. */
    virtual bool useTimeSlice();

    //---------------------------------------------------------------------------------------------
    // others:

    /** This function should be called by member functions that change some parameter such that a 
    re-drawing is required. */
    virtual void setDirty(bool shouldSetToDirty = true);

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** That is the function that should be overriden by your subclass to perform the actual 
    drawing operations. The baseclass implementation here will only draw a placeholder. */
    virtual void drawComponent(Image* imageToDrawOnto);

  private:

    /** Overrides the paint-function of the component base-class in order to fill the client area 
    with our member clientAreaImage. Moved to the private area because subclasses should not 
    override this - use drawComponent() instead. */
    void paint(Graphics &g);

    /** Takes care of memory (re-)allocation for the clientAreaImage and calls drawComponent() in a 
    thread-safe manner. */
    virtual void renderClientAreaImageInternal();

    /** (Re-) Allocates memeory for the client area image. The return value tells, if there was 
    actually memory re-allocated (i.e. it returns false, when re-allocation was not necessary or 
    failed). */
    virtual bool allocateClientAreaImage(int desiredWidth, int desiredHeight);

    /** The time lag that will be used by the inherited Timer object to trigger a repaint after 
    rendering has finished. */
    int repaintDelayInMilliseconds;

    /** This image will be used for drawing the client area of the panel - it will be rendered by 
    the member function renderClientAreaImage which you can override by your subclasses. This 
    rendering function can be called automatically whenever an event occurs that invalidates the 
    client area - to do that, we need the autoReRenderClientArea flag to be true (which it is by 
    default). */
    Image*          clientAreaImage; 
    bool            clientAreaImageIsDirty;
    CriticalSection clientAreaImageLock;
    TimeSliceThread *threadToUse;

  };

}

#endif  
