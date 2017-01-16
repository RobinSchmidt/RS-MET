ColorMap::ColorMap()
{
  fillDefaultMapNameArray();
  int initialSize = 2048;
  colors.resize(initialSize);
  lastIndex = initialSize-1;
  setDefaultMap(gray);
}

void ColorMap::setFromColourGradient(const ColourGradient &g)
{
  gradient = g;
  updateArray();
  defaultMapIndex = -1; // this is, in general, not a predefined map
}

void ColorMap::setDefaultMap(int index)
{
  setFromColourGradient(getDefaultGradient(index));
  defaultMapIndex = index;

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
  jassert(xml.getTagName() == "ColorMap");  // seems to be the wrong type of XML element

  String s = xml.getStringAttribute("Predefined", String::empty);
  if(s != String::empty)
    setDefaultMap(defaultMapNames.indexOf(s));
  else
  {
    ColourGradient g;
    Colour c;
    double p;
    int numColors = xml.getNumAttributes();
    for(int i = 0; i < numColors; i++)
    {
      c = Colour::fromString(xml.getAttributeName(i));
      p = xml.getAttributeValue(i).getDoubleValue();
      g.addColour(p, c);
    }
    setFromColourGradient(g);
  }
}

XmlElement* ColorMap::getAsXml()
{
  XmlElement* xml = new XmlElement("ColorMap"); 

  // \todo: allow the user to define colormaps in terms of HSL breakpoints - the xml-element 
  // should have and attribute "Type" which can be RGBA, HLSA or maybe some other color space

  if(defaultMapIndex != -1)
    xml->setAttribute("Predefined", defaultMapNames[defaultMapIndex]);
  else
  {
    for(int i = 0; i < gradient.getNumColours(); i++)
    {
      double p = gradient.getColourPosition(i);
      Colour c = gradient.getColour(i);
      String s = c.toString();
      xml->setAttribute(s, String(p));
    }
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
  case rainbow:
  {
    g.addColour(0.0, Colour(  0,   0,   0));  // black
    g.addColour(0.2, Colour(  0,   0, 255));  // blue
    g.addColour(0.4, Colour(  0, 255, 255));  // cyan
    g.addColour(0.6, Colour(  0, 255,   0));  // green
    g.addColour(0.8, Colour(255, 255,   0));  // yellow
    g.addColour(1.0, Colour(255,   0,   0));  // red
  } break;
  default: // gray-scale is the default map
  {
    g.addColour(0.0, Colour(  0,   0,   0));  // black
    g.addColour(1.0, Colour(255, 255, 255));  // white
  }
  }

  // todo: maybe define the default color maps in terms of HSL parameters. L goes from 0 to 1, 
  // H goes between two fixed values and S is a constant. For example, the fire colormap should go
  // from some red tone to a yellow, ice from purpleish blue to cyanish green..

  return g;
}

void ColorMap::fillDefaultMapNameArray()
{
  defaultMapNames.add("Gray");
  defaultMapNames.add("Fire");
  defaultMapNames.add("Ice");
  defaultMapNames.add("Rainbow");
}

void ColorMap::updateArray()
{
  double scaler = 1.0 / lastIndex;
  for(int i = 0; i < colors.size(); i++)
    colors[i] = gradient.getColourAtPosition(scaler*i).getARGB();
}