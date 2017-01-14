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
}

void ColorMap::setSize(int newSize)
{
  colors.resize(newSize);
  lastIndex = newSize-1;
  updateArray();
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
  Colour c;
  double scaler = 1.0 / lastIndex;
  for(int i = 0; i < colors.size(); i++)
  {
    c = gradient.getColourAtPosition(scaler * i);
    colors[i] = c.getARGB();
  }
}