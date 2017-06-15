#ifndef rosic_WaveformDisplayBuffer_h
#define rosic_WaveformDisplayBuffer_h

//// rosic-indcludes:
//#include "../filters/rosic_FourPoleFilter.h"

namespace rosic
{

  /**

  This is a class for temporarily storing sample values and retrieving buffers of sample-values for display on a scrolling realtime
  waveform display.

  \todo introduce a flag to suspend processing - use this when the GUI is closed
  \todo make the multichannel functionality available via a wrapper that maintains an array of one-channel displaybuffers

  */

  class WaveformDisplayBuffer
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WaveformDisplayBuffer();

    /** Destructor. */
    virtual ~WaveformDisplayBuffer();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    virtual void setSampleRate(double newSampleRate);

    /** Switches the analysis into Mid/Side mode. */
    //void setMidSideMode(bool shouldUseMidSideMode) { midSideMode = shouldUseMidSideMode; }

    /** Sets the length of the time window to be displayed (in seconds). */
    void setTimeWindowLength(double newTimeWindowLength);

    /** Sets the width of the display on which the waveform will be displayed (in pixels). This
    information is required in order to determine the correct decimation factor for the peak-data
    to be displayed. */
    void setDisplayWidth(int newDisplayWidth);

    // setNumChannels

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of data-points in the peak buffer. */
    int getDisplayBufferLength() const { return displayBufferLength; }

    /** Returns a pointer to the current buffers (decimated or interpolated for viewing). */
    virtual double* getDisplayBuffer();

    /** Returns a pointer to an array which holds the values for the time-axis (in seconds) which
    correspond to the current peak-values. */
    double* getTimeAxis() const { return timeAxisValues; }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Accepts a buffer of input samples. */
    void feedInputBuffer(double* inBuffer, int length);

    void feedInputBuffer(float* inBuffer, int length);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Clears all the internal buffers. */
    void clearBuffers();

    //=====================================================================================================================================

  protected:

    /** Copies the waveform data from the input buffer into the display buffer. */
    void updateDisplayBuffer();

    /** Accepts one input sample at a time. */
    virtual void feedInputSample(double in);


    void allocateBuffers();

    void freeBuffers();

    void updateTimeVariables();

    // buffers:
    double* timeAxisValues;
    double* inputBuffer;   // used to store the incoming data
    double* displayBuffer; // used to display the data

    // parameters:
    double  sampleRate;
    double  timeWindowLength;
    double  numSamplesShown;
    int     displayWidth;
    int     displayBufferLength;

    // factor into a stereo-version that embedds two mono ones:
    //int     numChannels;
    //bool    midSideMode;

    double  xMin, xMax, xOld;
    double  dw;
    int     nMin, nMax;
    int     n;   // sample counter inside the current chunk
    int     w;   // counter in the input buffer
    double  inc; // per-sample increment for index in the display buffer
    double  wd;  // counter in the input buffer as double (used in the interpolating case)

  };


  //=======================================================================================================================================

  /**

  Adds syncronization capabilities to WaveformDisplayBuffer for usage as oscilloscope.

  */

  class SyncedWaveformDisplayBuffer : public WaveformDisplayBuffer
  {

  public:

    enum syncModes
    {
      FREE_RUNNING = 0,
      ZEROS,
      LOWPASS_ZEROS
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SyncedWaveformDisplayBuffer();

    /** Destructor. */
    virtual ~SyncedWaveformDisplayBuffer();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Chooses a mode for the syncronization of the display to the signal. @see: syncModes */
    void setSyncMode(int newSyncMode);

    /** Sets a maximum time-window length above which synchronization will be deactivated - it doesn't seem to make sense to sync
    arbitrarily long windows. */
    void setSyncThreshold(double newSyncThreshold) { syncThreshold = newSyncThreshold; }

    /** When passing "true" to this, the display-buffer will be synced in a way to make the last pixel fall on a sync-point (i.e. a
    zero-crossing). If false (the default), the first pixel will fall on a sync-point. */
    void setSyncBackward(bool shouldSyncBackward) { syncBackward = shouldSyncBackward; }

    /** Sets the sample-rate. */
    virtual void setSampleRate(double newSampleRate);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Informs whether the display-buffer is synchronized in some way - for this to be the case, sync must be activated and the length of
    the time-window must be below some threshold (in seconds). */
    bool isSyncing() { return (syncMode != FREE_RUNNING) && (timeWindowLength <= syncThreshold); }

    /** Overriden to avoid the internal updateDisplayBuffer() call in case of a synced display (in which case the updateDisplayBuffer()
    function is called from feedInputSample() instead). */
    virtual double* getDisplayBuffer();

    //=====================================================================================================================================

  protected:

    /** Overriden to incorporate the sync-functionality. */
    virtual void feedInputSample(double in);

    /** Accepts an input sample and syncs the display to another signal that is supposedly somehow related to the input-signal. */
    void feedInputSampleWithSync(double in, double syncSignal);

    /** Returns true when the next index inside the input buffer is beyond the end of the buffer. Used internally to check, if the input
    buffer is already filled completely. */
    bool isNextIndexBehindBufferEnd();

    // embedded objects:
    FourPoleFilter lowpass; // for zero-crossing detection on lowpassed signal

    // data members:
    double xOld2;
    double yOld;            // previous sample of the signal to sync to
    double syncThreshold;   // threshold (in seconds) above which sync is deactivated
    int    syncMode;
    bool   bufferFull;      // flag to indicate the input buffer is completely filled
    bool   syncBackward;

  };

} // end namespace rosic

#endif
