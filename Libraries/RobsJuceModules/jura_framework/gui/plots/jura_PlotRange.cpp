rsPlotRange::rsPlotRange(double initMinX, double initMaxX,
  double initMinY, double initMaxY)
{
  minX = initMinX;
  maxX = initMaxX;
  minY = initMinY;
  maxY = initMaxY;
}

void rsPlotRange::setMinX(double newMinX) 
{ 
  jassert(newMinX < maxX);
  if( newMinX < maxX )
    minX = newMinX;
}

void rsPlotRange::setMaxX(double newMaxX) 
{ 
  jassert(newMaxX > minX);
  if( newMaxX > minX )
    maxX = newMaxX;
}

void rsPlotRange::setRangeX(double newMinX, double newMaxX)
{
  jassert(newMaxX > newMinX);
  if( newMaxX > newMinX )
  {
    minX = newMinX;
    maxX = newMaxX;
  }
}

void rsPlotRange::setMinY(double newMinY) 
{ 
  jassert(newMinY < maxY);
  if( newMinY < maxY )
    minY = newMinY;
}

void rsPlotRange::setMaxY(double newMaxY) 
{ 
  jassert(newMaxY > minY);
  if( newMaxY > minY )
    maxY = newMaxY;
}

void rsPlotRange::setRangeY(double newMinY, double newMaxY)
{
  jassert(newMaxY > newMinY);
  if( newMaxY > newMinY )
  {
    minY = newMinY;
    maxY = newMaxY;
  }
}

bool rsPlotRange::clipRange(rsPlotRange rangeToClipTo)
{
  bool clipped = false;
  if(minX < rangeToClipTo.getMinX()) { minX = rangeToClipTo.getMinX(); clipped = true;  }
  if(maxX > rangeToClipTo.getMaxX()) { maxX = rangeToClipTo.getMaxX(); clipped = true;  }
  if(minY < rangeToClipTo.getMinY()) { minY = rangeToClipTo.getMinY(); clipped = true;  }
  if(maxY > rangeToClipTo.getMaxY()) { maxY = rangeToClipTo.getMaxY(); clipped = true;  }
  return clipped;
}
