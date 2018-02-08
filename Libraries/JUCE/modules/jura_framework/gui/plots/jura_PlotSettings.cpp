rsPlotSettings::rsPlotSettings()
{
  // move initializations into header file:

  captionPosition               =  NO_CAPTION;
  axisPositionX                 =  ZERO;
  axisPositionY                 =  ZERO;
  axisLabelPositionX            =  ABOVE_AXIS;
  axisLabelPositionY            =  RIGHT_TO_AXIS;
  axisValuesPositionX           =  BELOW_AXIS;
  axisValuesPositionY           =  LEFT_TO_AXIS;

  axisLabelX                    =  String("x");
  axisLabelY                    =  String("y");

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

  logScaledX	                  =  false;
  logScaledY	                  =  false;
  logScaledRadius               =  false;

  stringConversionForAxisX     = &valueToString0;
  stringConversionForAxisY     = &valueToString0;
}

XmlElement* rsPlotSettings::getStateAsXml(const juce::String& name) const
{
  XmlElement* xml = new XmlElement(name); 
  // the XmlElement which stores all the releveant state-information

  /*
  xml->setAttribute(String("MinX"), currentRange.getMinX());
  xml->setAttribute(String("MaxX"), currentRange.getMaxX());
  xml->setAttribute(String("MinY"), currentRange.getMinY());
  xml->setAttribute(String("MaxY"), currentRange.getMaxY());
  */

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
  /*
  currentRange.setMinX( xml.getDoubleAttribute("MinX", getCurrentRangeMinX()) );
  currentRange.setMaxX( xml.getDoubleAttribute("MaxX", getCurrentRangeMaxX()) );
  currentRange.setMinY( xml.getDoubleAttribute("MinY", getCurrentRangeMinY()) );
  currentRange.setMaxY( xml.getDoubleAttribute("MaxY", getCurrentRangeMaxY()) );
  */

  horizontalCoarseGridIsVisible = xml.getBoolAttribute(  "HorizontalCoarseGridIsVisible", false);
  horizontalCoarseGridInterval  = xml.getDoubleAttribute("HorizontalCoarseGridInterval",  1);
  horizontalFineGridIsVisible   = xml.getBoolAttribute(  "HorizontalFineGridIsVisible",   false);
  horizontalFineGridInterval    = xml.getDoubleAttribute("HorizontalFineGridInterval",    0.1);
  verticalCoarseGridIsVisible   = xml.getBoolAttribute(  "VerticalCoarseGridIsVisible",   false);
  verticalCoarseGridInterval    = xml.getDoubleAttribute("VerticalCoarseGridInterval",    1);
  verticalFineGridIsVisible     = xml.getBoolAttribute(  "VerticalFineGridIsVisible",     false);
  verticalFineGridInterval      = xml.getDoubleAttribute("VerticalFineGridInterval",      0.1);
}