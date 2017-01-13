ColorMap::ColorMap()
{
  colors.resize(2048);
  lastIndex = 2047;
  setDefaultMap(grayScale);
}

void ColorMap::setFromColourGradient(const ColourGradient &g)
{

}

void ColorMap::setDefaultMap(int index)
{

}

//void ColorMap::setSize(int newSize)
//{
//
//}