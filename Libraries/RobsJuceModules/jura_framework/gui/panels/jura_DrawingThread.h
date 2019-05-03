#ifndef jura_DrawingThread_h
#define jura_DrawingThread_h

 /** This class implements a thread to be used for the drawing operation of objects of class Panel 
 (or subclasses there of). It derives from TimeSliceThread and overrides the constructor in order 
 to start itself inside this constructor. Panels register themselves as TimeSliceClient to the 
 global instance 'drawingThread' - this already taken care of in the constructor of class Panel 
 (and also the de-registering in the destructor). */

class DrawingThread : public TimeSliceThread
{

public:

  static DrawingThread* getInstance();

protected:

  DrawingThread();

private:

  static DrawingThread* soleInstance;

  juce_UseDebuggingNewOperator;
};

#endif  