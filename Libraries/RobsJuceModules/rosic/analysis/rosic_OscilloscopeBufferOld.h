#ifndef rosic_OscilloscopeBufferOld_h
#define rosic_OscilloscopeBufferOld_h

//// rosic-indcludes:
//#include "../filters/rosic_FourPoleFilter.h"

namespace rosic
{

  /**

  This is a class for temporarily storing sample values and retrieving buffers of sample-values
  for display on an oscilloscope.

  \todo use block-processing, introduce a flag to suspend processing
  \todo generate peak-data on the fly, sync only when it makes sense (don't sync for arbitrarily long buffers)

  */

  class OscilloscopeBufferOld
  {

  public:

    enum syncModes
    {
      FREE_RUNNING = 0,
      LOWPASS_ZEROS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. Requires the width of the display (in pixels) to be passed in order to 
    correctly set up the buffers for the peak-data. */
    OscilloscopeBufferOld(int displayWidthInPixels);   

    /** Destructor. */
    ~OscilloscopeBufferOld();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Switches the analysis into Mid/Side mode. */
    void setMidSideMode(bool shouldUseMidSideMode) { midSideMode = shouldUseMidSideMode; }

    /** Sets the length of the time window to be displayed (in seconds). */
    void setTimeWindowLength(double newTimeWindowLength);

    /** Sets the starting time of the time window (in seconds) - thsi will result in a constant 
    offset in the timeAxisValues array. */
    void setTimeWindowStart(double newTimeWindowStart);

    /** Sets the width of the display on which the waveform will be displayed (in pixels). This 
    information is required in order to determine the correct decimation factor for the peak-data 
    to be displayed. */
    void setDisplayWidth(int newDisplayWidth);

    /** Chooses a mode for the syncronization of the display to the signal. @see: syncModes */
    void setSyncMode(int newSyncMode);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the current mode for the syncronization. @see: syncModes */
    int getSyncMode() const { return syncMode; }

    /** Returns the number of signal channels. */
    int getNumChannels() const { return numChannels; }

    /** Returns the number of data-points in the peak buffer. */
    int getViewBufferLength() const { return viewBufferLength; }

    /** Returns a pointer to an array which holds the values for the time-axis (in seconds) which 
    correspond to the current peak-values. */
    double* getTimeAxis() const { return timeAxisValues; }

    /** Returns a pointer to the current buffers (decimated or interpolated for viewing) - the first index indicates the 
    channel, the second index indicates the bin. */
    float** getCurrentDisplayBuffers() const { return viewBuffer; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Buffers a stereo-ouput frame. */
    INLINE void bufferSampleFrameStereo(float* inL, float* inR);


    INLINE void bufferSampleFrameStereo(double* inL, double* inR);


    //---------------------------------------------------------------------------------------------
    // others:

    /** Clears all the internal buffers. */
    void clearBuffers();

    /** Causes and update of all the display-related buffers. */
    void updateDisplayBuffers();

    /** Returns a sample index which is wrapped around version of the input index. */
    INLINE int wrapAround(int index);

    //=============================================================================================

    // for debug:
    int getDecimationFactor() { return decimationFactor; }

  protected:

    /** Fills the array timeAxisValues with values according to the current settings of sampleRate
    timeWindowLength and timeWindowStart. */
    void calculateTimeAxisValues();

    /** Returns the number of samples back into the buffer, at which the display should begin to 
    read out the signal according to the syncMode.  */
    int getHindsight();

    static const int maxNumChannels = 2;
    static const int maxBufferSize  = 131072;

    // buffers and pointers:
    float*  circularBufferFlat;  
    float** circularBuffer;            
    float*  linearBufferFlat;
    float** linearBuffer;
    float*  viewBufferFlat;               // we must use min/max values when zoomed out far
    float** viewBuffer;
    float*  filteredCircularBuffer;
    double* timeAxisValues;

    // parameters:
    double  sampleRate;
    double  timeWindowStart;
    double  timeWindowLength;
    double  syncThreshold;       // threshold (in seconds) above which sync is deactivated
    int     displayWidth;
    int     decimationFactor;
    int     syncMode;
    int     numChannels;
    int     viewBufferLength;
    int     sampleCounter;
    bool    midSideMode;

    // embedded objects:
    FourPoleFilter lowpass; // for zero-crossing detection on lowpassed signal

  };

  //---------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void OscilloscopeBufferOld::bufferSampleFrameStereo(float* inL,  float* inR)
  {
    // wraparound circular buffer if necesarry:
    while( sampleCounter >= maxBufferSize )
      sampleCounter -= maxBufferSize;

    // establish the 2 input channel signals:
    float tmp1, tmp2;
    if( midSideMode == true )
    {
      tmp1 = (float) SQRT2_INV * (*inL + *inR);
      tmp2 = (float) SQRT2_INV * (*inL - *inR);
    }
    else
    {
      tmp1 = *inL;
      tmp2 = *inR;
    }

    // buffer the incoming samples:
    circularBuffer[0][sampleCounter] = tmp1;
    circularBuffer[1][sampleCounter] = tmp2;

    // store also a filtered version in case we want to sync to the zero-crossings:
    if( syncMode == LOWPASS_ZEROS )
      filteredCircularBuffer[sampleCounter] = (float) lowpass.getSample((double) tmp1);

    sampleCounter++;
  }

  INLINE void OscilloscopeBufferOld::bufferSampleFrameStereo(double* inL, double* inR)
  {
    float fL = (float) (*inL);
    float fR = (float) (*inR);
    bufferSampleFrameStereo(&fL, &fR);
  }

  INLINE int OscilloscopeBufferOld::wrapAround(int index)
  {
    while( index < 0 )
      index += maxBufferSize;
    while( index >= maxBufferSize )
      index -= maxBufferSize;
    return index;
  }

} // end namespace rosic

#endif 
