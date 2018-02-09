#ifndef jura_PlotSettings_h
#define jura_PlotSettings_h

class rsPlotSettings;

/** Baseclass for objects that must keep track of the state of an rsPlotSettings object - perhaps
in order to reapint themselves when somthing changes. */

class JUCE_API rsPlotSettingsObserver
{

public:

  /** Called when something about the appearance of the plot has changed such as the postioning of
  the axes, the grid settings, etc. */
  virtual void rsPlotAppearanceChanged(rsPlotSettings* plotSettings) = 0;

  /** Called when the currently visible range was changed, for example due to zooming or 
  scrolling. */
  virtual void rsPlotVisibleRangeChanged(rsPlotSettings* plotSettings) = 0;

  /** Called when the maximum available range has changed (i.e. then range when the plot is fully
  zoomed out). */
  virtual void rsPlotMaximumRangeChanged(rsPlotSettings* plotSettings) = 0;

};

//=================================================================================================

/** This class represents the settings for a plot, such as the ranges, the visibility of grids, 
etc. */

class JUCE_API rsPlotSettings
{

public:

  rsPlotSettings();

  //-----------------------------------------------------------------------------------------------
  // \name Enumerations

  /** Positions where the caption may appear. */
  enum captionPositions
  {
    NO_CAPTION = 0,
    TOP_CENTER,
    CENTER
  };

  /** Positions where the x or y axis may appear. */
  enum axisPositions
  {
    INVISIBLE = 0,
    ZERO,
    LEFT,
    RIGHT,
    TOP,
    BOTTOM
  };

  /** Positions where the axis annotations (such as tick-marks) may appear. */
  enum axisAnnotationPositions
  {
    NO_ANNOTATION = 0,
    LEFT_TO_AXIS,
    RIGHT_TO_AXIS,
    ABOVE_AXIS,
    BELOW_AXIS
  };


  //-----------------------------------------------------------------------------------------------
  // \name Setup


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry




  //-----------------------------------------------------------------------------------------------
  // \name Oberservation

  void registerObserver(rsPlotSettingsObserver* obs) { appendIfNotAlreadyThere(observers, obs); }

  void deregisterObserver(rsPlotSettingsObserver* obs) { removeFirstOccurrence(observers, obs); }

  void sendAppearenceChangeNotification();

  void sendVisibleRangeChangeNotification();

  void sendMaximumRangeChangeNotification();


  //-----------------------------------------------------------------------------------------------
  // \name State

  /** Creates an XmlElement from the current state and returns it. */
  virtual XmlElement* getStateAsXml(const juce::String& name = juce::String("PlotSettings")) const;

  /** Restores a state based on an XmlElement which should have been created
  with the getStateAsXml()-function. */
  virtual void setStateFromXml(const XmlElement &xml);





  //-----------------------------------------------------------------------------------------------
  // \name Data

  // make protected:

  int captionPosition;
  juce::String captionString;

  int axisPositionX;
  int axisPositionY;
  int axisLabelPositionX;
  int axisLabelPositionY;
  int axisValuesPositionX;
  int axisValuesPositionY;

  juce::String axisLabelX;
  juce::String axisLabelY;

  // actually, they are redundant - we can use values of 0 for the corresponding intervals to turn 
  // them off - but no, it's more convenient to be able to set the spacing separately from the 
  // visibility (it will be separate on a gui anyway):
  bool horizontalCoarseGridIsVisible;
  bool horizontalFineGridIsVisible;
  bool verticalCoarseGridIsVisible;
  bool verticalFineGridIsVisible;
  bool radialCoarseGridIsVisible;
  bool radialFineGridIsVisible;
  bool angularCoarseGridIsVisible;
  bool angularFineGridIsVisible;

  double horizontalCoarseGridInterval;
  double horizontalFineGridInterval;
  double verticalCoarseGridInterval;
  double verticalFineGridInterval;
  double radialCoarseGridInterval;
  double radialFineGridInterval;
  double angularCoarseGridInterval;
  double angularFineGridInterval;

  bool logScaledX;
  bool logScaledY;
  bool logScaledRadius;

  // functions for converting coordinates into strings:
  juce::String (*stringConversionForAxisX) (double valueToConvert);
  juce::String (*stringConversionForAxisY) (double valueToConvert);

  /** The currently visible range and maximum range object for the plot. */
  rsPlotRange currentRange, maximumRange;

protected:

  /** Performs a sanity check */
  //void sanityCheckGridSpacing(double& coarse, bool logScaled);



  std::vector<rsPlotSettingsObserver*> observers;

  //bool notificationsOff = false;

};

#endif