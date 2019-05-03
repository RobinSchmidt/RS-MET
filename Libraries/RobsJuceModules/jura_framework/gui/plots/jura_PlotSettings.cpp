rsPlotSettings::rsPlotSettings()
{
  maximumRange.setRangeX(-2.2, 2.2);
  maximumRange.setRangeY(-2.2, 2.2);
  currentRange.setRangeX(-2.2, 2.2);
  currentRange.setRangeY(-2.2, 2.2);

  // move initializations into header file:

  captionPosition     =  NO_CAPTION;
  axisPositionX       =  ZERO;
  axisPositionY       =  ZERO;
  axisLabelPositionX  =  ABOVE_AXIS;
  axisLabelPositionY  =  RIGHT_TO_AXIS;
  axisValuesPositionX =  BELOW_AXIS;
  axisValuesPositionY =  LEFT_TO_AXIS;

  axisLabelX = String("x");
  axisLabelY = String("y");

  horizontalCoarseGridIsVisible =  false;
  horizontalFineGridIsVisible	  =  false;
  verticalCoarseGridIsVisible   =  false;
  verticalFineGridIsVisible	    =  false;
  radialCoarseGridIsVisible     =  false;
  radialFineGridIsVisible       =  false;
  angularCoarseGridIsVisible    =  false;
  angularFineGridIsVisible      =  false;

  horizontalCoarseGridInterval  =  1.0;
  horizontalFineGridInterval    =  0.1;
  verticalCoarseGridInterval    =  1.0;
  verticalFineGridInterval      =  0.1;
  radialCoarseGridInterval      =  1.0;
  radialFineGridInterval        =  0.1;
  angularCoarseGridInterval     =  15.0;  // 15 degrees
  angularFineGridInterval       =  5.0;   // 5 degrees

  logScaledX	    =  false;
  logScaledY	    =  false;
  logScaledRadius =  false;

  stringConversionForAxisX = &valueToString0;
  stringConversionForAxisY = &valueToString0;
}

// setup:

void rsPlotSettings::setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY)
{
  maximumRange.setRangeX(newMinX, newMaxX);
  maximumRange.setRangeY(newMinY, newMaxY);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRange(rsPlotRange newMaximumRange)
{
  maximumRange = newMaximumRange;
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRangeX(double newMinX, double newMaxX)
{
  maximumRange.setRangeX(newMinX, newMaxX);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRangeY(double newMinY, double newMaxY)
{
  maximumRange.setRangeY(newMinY, newMaxY);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRangeMinX(double newMinX)
{
  maximumRange.setMinX(newMinX);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRangeMaxX(double newMaxX)
{
  maximumRange.setMaxX(newMaxX);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRangeMinY(double newMinY)
{
  maximumRange.setMinY(newMinY);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setMaximumRangeMaxY(double newMaxY)
{
  maximumRange.setMaxY(newMaxY);
  sendMaximumRangeChangeNotification();
  clipVisibleToMaximumRange();
}

void rsPlotSettings::setCurrentRange(double newMinX, double newMaxX, 
  double newMinY, double newMaxY)
{
  currentRange.setRangeX(newMinX, newMaxX);
  currentRange.setRangeY(newMinY, newMaxY);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRange(rsPlotRange newRange)
{
  currentRange = newRange;
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRangeX(double newMinX, double newMaxX)
{
  currentRange.setRangeX(newMinX, newMaxX);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRangeY(double newMinY, double newMaxY)
{
  currentRange.setRangeY(newMinY, newMaxY);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRangeMinX(double newMinX)
{
  currentRange.setMinX(newMinX);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRangeMaxX(double newMaxX)
{
  currentRange.setMaxX(newMaxX);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRangeMinY(double newMinY)
{
  currentRange.setMinY(newMinY);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

void rsPlotSettings::setCurrentRangeMaxY(double newMaxY)
{
  currentRange.setMaxY(newMaxY);
  currentRange.clipRange(maximumRange);
  sendVisibleRangeChangeNotification();
}

// observation:

void rsPlotSettings::sendAppearenceChangeNotification()
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->rsPlotAppearanceChanged(this);
}

void rsPlotSettings::sendVisibleRangeChangeNotification()
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->rsPlotVisibleRangeChanged(this);
}

void rsPlotSettings::sendMaximumRangeChangeNotification()
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->rsPlotMaximumRangeChanged(this);
}

// state:

XmlElement* rsPlotSettings::getStateAsXml(const juce::String& name) const
{
  XmlElement* xml = new XmlElement(name); 
  // the XmlElement which stores all the releveant state-information

  xml->setAttribute("MinX", currentRange.getMinX());
  xml->setAttribute("MaxX", currentRange.getMaxX());
  xml->setAttribute("MinY", currentRange.getMinY());
  xml->setAttribute("MaxY", currentRange.getMaxY());

  xml->setAttribute("HorizontalCoarseGridIsVisible", horizontalCoarseGridIsVisible);
  xml->setAttribute("HorizontalCoarseGridInterval",  horizontalCoarseGridInterval);
  xml->setAttribute("HorizontalFineGridIsVisible",   horizontalFineGridIsVisible);
  xml->setAttribute("HorizontalFineGridInterval",    horizontalFineGridInterval); 
  xml->setAttribute("VerticalCoarseGridIsVisible",   verticalCoarseGridIsVisible);
  xml->setAttribute("VerticalCoarseGridInterval",    verticalCoarseGridInterval);
  xml->setAttribute("VerticalFineGridIsVisible",     verticalFineGridIsVisible);
  xml->setAttribute("VerticalFineGridInterval",      verticalFineGridInterval);

  return xml;
}

void rsPlotSettings::setStateFromXml(const XmlElement &xml)
{
  currentRange.setMinX( xml.getDoubleAttribute("MinX", currentRange.getMinX()) );
  currentRange.setMaxX( xml.getDoubleAttribute("MaxX", currentRange.getMaxX()) );
  currentRange.setMinY( xml.getDoubleAttribute("MinY", currentRange.getMinY()) );
  currentRange.setMaxY( xml.getDoubleAttribute("MaxY", currentRange.getMaxY()) );

  horizontalCoarseGridIsVisible = xml.getBoolAttribute(  "HorizontalCoarseGridIsVisible", false);
  horizontalCoarseGridInterval  = xml.getDoubleAttribute("HorizontalCoarseGridInterval",  1);
  horizontalFineGridIsVisible   = xml.getBoolAttribute(  "HorizontalFineGridIsVisible",   false);
  horizontalFineGridInterval    = xml.getDoubleAttribute("HorizontalFineGridInterval",    0.1);
  verticalCoarseGridIsVisible   = xml.getBoolAttribute(  "VerticalCoarseGridIsVisible",   false);
  verticalCoarseGridInterval    = xml.getDoubleAttribute("VerticalCoarseGridInterval",    1);
  verticalFineGridIsVisible     = xml.getBoolAttribute(  "VerticalFineGridIsVisible",     false);
  verticalFineGridInterval      = xml.getDoubleAttribute("VerticalFineGridInterval",      0.1);
}

void rsPlotSettings::clipVisibleToMaximumRange()
{
  if(currentRange.clipRange(maximumRange))
    sendVisibleRangeChangeNotification();
}