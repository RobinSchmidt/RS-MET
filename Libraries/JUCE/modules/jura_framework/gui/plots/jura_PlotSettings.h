#ifndef jura_PlotSettings_h
#define jura_PlotSettings_h

/** This class represents the settings for a plot, such as the ranges, the visibility of grids, 
etc. */

class JUCE_API rsPlotSettings
{

public:



  rsPlotSettings();


  /** Creates an XmlElement from the current state and returns it. */
  virtual XmlElement* getStateAsXml(const juce::String& name = juce::String("PlotSettings")) const;

  /** Restores a state based on an XmlElement which should have been created
  with the getStateAsXml()-function. */
  virtual void setStateFromXml(const XmlElement &xml);


  //-----------------------------------------------------------------------------------------------
  // \name Enumerations

  enum captionPositions
  {
    NO_CAPTION = 0,
    TOP_CENTER,
    CENTER
  };

  enum axisPositions
  {
    INVISIBLE = 0,
    ZERO,
    LEFT,
    RIGHT,
    TOP,
    BOTTOM
  };

  enum axisAnnotationPositions
  {
    NO_ANNOTATION = 0,
    LEFT_TO_AXIS,
    RIGHT_TO_AXIS,
    ABOVE_AXIS,
    BELOW_AXIS
  };

  //-----------------------------------------------------------------------------------------------
  // \name Data

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
  // them off:
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
  juce::String (*stringConversionForAxisX)     (double valueToConvert);
  juce::String (*stringConversionForAxisY)     (double valueToConvert);

  /** The currently visible range and maximum range object for the plot. */
  rsPlotRange currentRange, maximumRange;

};

#endif