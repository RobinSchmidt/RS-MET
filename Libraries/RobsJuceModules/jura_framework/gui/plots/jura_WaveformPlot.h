#ifndef jura_WaveformPlot_h
#define jura_WaveformPlot_h

/** This class is intended to be used as a graphical display for a waveform.  */

class rsWaveformPlot : virtual public rsDataPlot
{

public:

  enum timeFormats
  {
    HOUR_MINUTE_SECOND = 1,
    SAMPLES
    // SECONDS
    // BEATS
  };

  rsWaveformPlot(const juce::String& name = juce::String("rsWaveformPlot"));
  /**< Constructor. */

  virtual ~rsWaveformPlot();
  /**< Destructor. */

  virtual rsPlotRange getMaximumMeaningfulRange(double relativeMarginLeft = 0.0,
    double relativeMarginRight  = 0.0, double relativeMarginTop  = 0.0,
    double relativeMarginBottom = 0.0);
  /**< Returns an Rectangle wich encloses the curve. Optionally, a margin can
  be specified in percent. */

  virtual bool setWaveform(double** newWaveformData, int newNumSampleFrames, int newNumChannels);
  /**< Sets the actual waveform-data and informs about success. */

  virtual bool setWaveform(float** newWaveformData, int newNumSampleFrames, int newNumChannels);
  /**< Sets the actual waveform-data and informs about success. */

  virtual bool setWaveform(const AudioSampleBuffer& newWaveformBuffer);
  /**< Sets the waveform-data from an AudioSampleBuffer and informs about success. */

  virtual void setSampleRate(double newSampleRate);
  /**< Updates the sample rate for this display - this is need to display the time-axis in
  seconds, when desired. */

  virtual void updateWaveformCurve();
  /**< Sets up the CoordinateSystem and redraws the waveform. */

  //virtual void paint(Graphics &g);
  /**< Overrides the paint-function of the base-classes. */

protected:

  virtual void createDecimatedData();
  /**< Creates the peak-data from the original data by means of minimum/maximum extraction and
  decimation. The peak data will be needed when one sample-unit in screen-coordinates is
  smaller than a pixel. */

  virtual void plotCurveFamily(Graphics &g, Image *targetImage = NULL,
    XmlElement *targetSVG = NULL);
  /**< Overrides rsDataPlot::plotCurveFamily in order to call the plotWaveform() function. */

  virtual void plotWaveform(Graphics &g, Image *targetImage = NULL,
    XmlElement *targetSVG = NULL);
  /**< Plots the actual waveform data on a graphics-canvas. */

  // two-dimensional array for the (multichannel) waveform-data for full 
  // detail display:
  double   sampleRate;
  int      timeFormat;
  int      numChannels;
  int      numSampleFrames;
  float*   peakData;

  juce_UseDebuggingNewOperator;
};


#endif 
