template<class TPix>
AlphaMask<TPix>::AlphaMask()
{
  bell.setCenter(0);
  bell.setWidth(2);
  bell.setFlatTopWidth(0.5);
  bell.setPrototypeBell(&rsPositiveBellFunctions<double>::cubic);
    // make this user adjustable later

  setSize(5);
}

template<class TPix>
void AlphaMask<TPix>::setSize(double newSize)
{
  int pixelSize = (int)ceil(newSize);
  ImageResizable::setSize(pixelSize, pixelSize);
  bell.setWidth(newSize);
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::setTransitionWidth(double newWidth)
{
  bell.setFlatTopWidth(newWidth);
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::renderMask()
{
  // render circular alpha mask - alpha value depends on distance from center:
  double cx = 0.5 * width;  // x-coordinate of center
  double cy = 0.5 * height; // y-coordinate of center
  for(int y = 0; y < height; y++)
  {
    for(int x = 0; x < width; x++)
    {
      double dx = x - cx;
      double dy = y - cy;
      double distance = sqrt(dx*dx + dy*dy);
      TPix alpha = bell.getValue(distance);
      setPixelColor(x, y, alpha);
    }
  }
  // maybe the code can be generalized to an elliptic mask: divide dx by width and dy by height 
  // and use a bell with width = 1
}