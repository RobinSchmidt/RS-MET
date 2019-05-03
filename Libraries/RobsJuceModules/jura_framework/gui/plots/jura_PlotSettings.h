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
etc. 

\todo: 
-for the range setup, we may want to check, if the new settings are actually different from the old
 settings and return early (and don't send notifications) to avoid unnecessary updates/computations
 ->check when this happens, maybe it can be avoided at a higher level
 
*/

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

  /** Sets the maximum for the currently visible range. For logarithmic x- and/or y-axis-scaling,
  make sure that the respective minimum value is greater than zero! */
  void setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY);

  /** Sets the maximum for the currently visible range. */
  void setMaximumRange(rsPlotRange newMaximumRange);

  /** Sets the maximum visible range for the y-axis. */
  void setMaximumRangeX(double newMinX, double newMaxX);

  /** Sets the maximum visible range for the y-axis. */
  void setMaximumRangeY(double newMinY, double newMaxY);

  /** Sets the minimum value for the range of x. */
  void setMaximumRangeMinX(double newMinX);

  /** Sets the maximum value for the range of x. */
  void setMaximumRangeMaxX(double newMaxX);

  /** Sets the minimum value for the range of y. */
  void setMaximumRangeMinY(double newMinY);

  /** Sets the maximum value for the range of y. */
  void setMaximumRangeMaxY(double newMaxY);

  /** Sets the currently visible range. For logarithmic x- and/or y-axis-scaling, make sure that
  the respective minimum value is greater than zero! */
  void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);

  /** Sets the currently visible range. */
  void setCurrentRange(rsPlotRange newRange);

  /** Sets the currently visible range for the y-axis. */
  void setCurrentRangeX(double newMinX, double newMaxX);

  /** Sets the currently visible range for the y-axis. */
  void setCurrentRangeY(double newMinY, double newMaxY);

  /** Sets the minimum value of x. */
  void setCurrentRangeMinX(double newMinX);

  /** Sets the maximum value of x. */
  void setCurrentRangeMaxX(double newMaxX);

  /** Sets the minimum value of y. */
  void setCurrentRangeMinY(double newMinY);

  /** Sets the maximum value of y. */
  void setCurrentRangeMaxY(double newMaxY);

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  rsPlotRange getMaximumRange() const { return maximumRange; }
  rsPlotRange getCurrentRange() const { return currentRange; }
  double getMaximumRangeMinX() const { return maximumRange.getMinX(); }
  double getMaximumRangeMaxX() const { return maximumRange.getMaxX(); }
  double getMaximumRangeMinY() const { return maximumRange.getMinY(); }
  double getMaximumRangeMaxY() const { return maximumRange.getMaxY(); }
  double getCurrentRangeMinX() const { return currentRange.getMinX(); }
  double getCurrentRangeMaxX() const { return currentRange.getMaxX(); }
  double getCurrentRangeMinY() const { return currentRange.getMinY(); }
  double getCurrentRangeMaxY() const { return currentRange.getMaxY(); }

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



protected:

  /** The currently visible range and maximum range object for the plot. */
  rsPlotRange currentRange, maximumRange;

  /** Clips the visible range to the maximum available range. If the visible range actually changes
  due to this, it will also send out a rsPlotVisibleRangeChanged notification. */
  void clipVisibleToMaximumRange();

  /** Performs a sanity check */
  //void sanityCheckGridSpacing(double& coarse, bool logScaled);

  std::vector<rsPlotSettingsObserver*> observers;

  //bool notificationsOff = false;

};

#endif