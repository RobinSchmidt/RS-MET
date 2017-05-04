#ifndef DualWaveformDisplay_h
#define DualWaveformDisplay_h

#include "rojue_WaveformDisplay.h"

namespace rojue
{

  /**

  This class is a component which can show two waveform displays one on top of the other for 
  displaying stereo material. It encapsulates two instances of WaveformDisplay, but can be
  used from the outside basically the same way as the WaveformDisplay class - it will 
  automatically detect whether or not the buffer is stereo and will automatically switch between 
  showing one display using the full client area and showing two displays using half of the area 
  each.

  */

  class DualWaveformDisplay : public Component, public AudioFileBufferUser
  {  

    friend class AudioClipComponent;

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    DualWaveformDisplay(AudioFileBuffer* newBuffer = NULL);  

    /** Destructor. */
    virtual ~DualWaveformDisplay();                             

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Assigns a new underlying AudioBFileuffer to this DualWaveformDisplay  */
    virtual void assignAudioFileBuffer(AudioFileBuffer* newBuffer = NULL);

    //---------------------------------------------------------------------------------------------
    // component-overrides:

    /** Overrides the setBounds()-method of the Component base-class in order to arrange the 
    wav-displays according to the size - for best viewing it is advisable to use even heights. */
    virtual void resized();

    //---------------------------------------------------------------------------------------------
    // we mimic some of the WaveformDisplay's functions and pass the respective arguments through 
    // to both of our encapsulated WaveformDisplays:

    virtual double getMaximumRangeMinX() const;
    virtual double getMaximumRangeMaxX() const;
    virtual double getCurrentRangeMinX() const;
    virtual double getCurrentRangeMaxX() const;
    virtual void setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
    virtual void setMaximumRangeX(double newMinX, double newMaxX);
    virtual void setMaximumRangeY(double newMinY, double newMaxY);
    virtual void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
    virtual void setCurrentRangeX(double newMinX, double newMaxX);
    virtual void setCurrentRangeMinX(double newMinX);
    virtual void setCurrentRangeMaxX(double newMaxX);
    virtual void setCurrentRangeY(double newMinY, double newMaxY);
    virtual void setVisibleTimeRange(double newMinTimeInSeconds, double newMaxTimeInSeconds);

    virtual void setDrawingThread(TimeSliceThread* newDrawingThread);
    virtual void updatePlotImage();

    virtual void lockUsedBufferPointer();
    virtual void unlockUsedBufferPointer();
    virtual void acquireAudioFileBufferWriteLock();
    virtual void releaseAudioFileBufferWriteLock();

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    //AudioFileBuffer *bufferToUse;
    WaveformDisplay *waveDisplayL, *waveDisplayR; 

  };

}

#endif  