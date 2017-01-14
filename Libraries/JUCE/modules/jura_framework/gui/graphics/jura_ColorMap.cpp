ColorMap::ColorMap()
{
  int initialSize = 2048;

  colors.resize(initialSize);
  lastIndex = initialSize-1;
  setDefaultMap(grayScale);
}

void ColorMap::setFromColourGradient(const ColourGradient &g)
{
  gradient = g;
  updateArray();
}

void ColorMap::setDefaultMap(int index)
{
  setFromColourGradient(getDefaultGradient(index));
  // with the named pre-defined maps, we are not even restricted to using maps based on a 
  // juce::ColourGradient object (which uses linear interpolation RGB space, which might not be
  // ideal). We can use custom formulas or use the ColorAHSL class to create the map
}

void ColorMap::setSize(int newSize)
{
  colors.resize(newSize);
  lastIndex = newSize-1;
  updateArray();
}

void ColorMap::setFromXml(const XmlElement& xml)
{
  // not yet implemented
}

XmlElement* ColorMap::getAsXml()
{
  XmlElement* xml = new XmlElement("ColorMap"); 

  // \todo have a mechanism to determine whether this current colormap is one of the default 
  // hard-coded ones or a custom map. if it's a hard coded one, we just store the name of the map
  // as one attribute like DefaultMap="Fire" or something like that

  // we can have a defaultMapIndex member that we set to -1 when it's a custom map

  for(int i = 0; i < gradient.getNumColours(); i++)
  {
    double p = gradient.getColourPosition(i);
    Colour c = gradient.getColour(i);
    String s = c.toString();
    xml->setAttribute(s, String(p));
  }

  return xml;
}

ColourGradient ColorMap::getDefaultGradient(int index)
{
  juce::ColourGradient g;

  switch(index)
  {
  case fire:
  {
    g.addColour(0.0, Colour(  0,   0,   0));  // black
    g.addColour(0.4, Colour(255,   0,   0));  // red
    g.addColour(0.6, Colour(255, 255,   0));  // yellow
    g.addColour(1.0, Colour(255, 255, 255));  // white
  } break;
  case ice:
  {
    g.addColour(0.0, Colour(  0,   0,   0));  // black
    g.addColour(0.4, Colour(  0,   0, 255));  // blue
    g.addColour(0.6, Colour(  0, 255, 255));  // cyan
    g.addColour(1.0, Colour(255, 255, 255));  // white
  } break;
  default: // gray-scale is the default map
  {
    g.addColour(0.0, Colour(  0,   0,   0));  // black
    g.addColour(1.0, Colour(255, 255, 255));  // white
  }
  }

  return g;
}

void ColorMap::updateArray()
{
  double scaler = 1.0 / lastIndex;
  for(int i = 0; i < colors.size(); i++)
    colors[i] = gradient.getColourAtPosition(scaler*i).getARGB();
}