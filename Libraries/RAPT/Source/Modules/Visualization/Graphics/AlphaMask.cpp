template<class TPix>
AlphaMask<TPix>::AlphaMask()
{
  flat   = 0.0;
  slope0 = 0.0;
  slope1 = 0.0;
  setSize(5);
}

template<class TPix>
void AlphaMask<TPix>::setSize(double newSize)
{
  int pixelSize = (int)ceil(newSize);
  ImageResizable<TPix>::setSize(pixelSize, pixelSize);
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::setTransitionWidth(double newWidth)
{
  flat = 1-newWidth;
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::setInnerSlope(double newSlope)
{
  slope0 = newSlope;
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::setOuterSlope(double newSlope)
{
  slope1 = newSlope;
  renderMask();
}

template<class TPix>
template<class T>
void AlphaMask<TPix>::copyShapeParametersFrom(const AlphaMask<T>& other)
{
  flat   = (TPix)other.getTransitionWidth();
  slope0 = (TPix)other.getInnerSlope();
  slope1 = (TPix)other.getOuterSlope();
  renderMask();
}

template<class TPix>
double AlphaMask<TPix>::cubicBell(double x, double steepnessAt0, double steepnessAt1)
{
  double s0 = -steepnessAt0;
  double s1 = -steepnessAt1;
  double a2 = -s1 - 2*s0 - 3;
  double a3 =  s1 +   s0 + 2;
  // todo: precompute the coeffs

  double y = 1 + s0*x + a2*x*x + a3*x*x*x;  // optimize, limit to 0..1
  return y;

  // we use a 3rd order polynomial with adjustable slopes at start- and endpoint.
  // f(x) = a0 + a1*x + a2*x^2 + a3*x^3 with f(0) = 1, f'(0) = s0, f(1) = 0, f'(1) = s1
  // solving it yields: a0 = 1, a1 = s0, a2 = -s1 - 2*s0 - 3, a3 = 2 + s0 + s1
  // maybe the user parameters slope1, slope2 are -s1, -s2
  // maybe we should figure out the condition for not having a local minimum or maximum
  // between 0..1 and perhaps restrict the parameter range
}

template<class TPix>
void AlphaMask<TPix>::renderMask()
{
  int w = this->width;
  int h = this->height;

  // render circular alpha mask - alpha value depends on distance from center:
  double cx = 0.5 * w;      // x-coordinate of center
  double cy = 0.5 * h;      // y-coordinate of center
  double sx = 1.0 / cx;     // x-scaler
  double sy = 1.0 / cy;     // y-scaler

  for(int y = 0; y < h; y++)
  {
    for(int x = 0; x < w; x++)
    {
      double dx = sx * (x - cx);
      double dy = sy * (y - cy);
      double distance = sqrt(dx*dx + dy*dy);
      TPix alpha = (TPix) abs(getAlphaForDistance(distance));
      ImageResizable<TPix>::setPixelColor(x, y, alpha);
    }
  }
  // maybe the code can be generalized to an elliptic mask: divide dx by width and dy by height
  // and use a bell with width = 1
}

template<class TPix>
double AlphaMask<TPix>::getAlphaForDistance(double d)
{
  if(d < flat)
    return 1.0;
  if(d > 1.0)
    return 0.0;

  double x;
  if(flat != 1.0)
    x = (d-flat) / (1-flat); // todo: catch div-by-zero when flat==1
  else
    x = d;

  return cubicBell(x, slope0, slope1);
}
