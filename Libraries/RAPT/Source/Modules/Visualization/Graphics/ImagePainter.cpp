template<class TPix, class TWgt, class TCor>
ImagePainter<TPix, TWgt, TCor>::ImagePainter(Image<TPix> *imageToPaintOn, AlphaMask<TWgt> *maskToUse)
{
  antiAlias = true;
  useMask = false;
  setNeighbourWeightsForSimpleDot(0, 0);

  setImageToPaintOn(imageToPaintOn);
  setAlphaMaskForDot(maskToUse);
}

// setup

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setImageToPaintOn(Image<TPix> *imageToPaintOn)
{
  image = imageToPaintOn;
  //if(image != nullptr)
  //{
  //  wi = image->getWidth();
  //  hi = image->getHeight();
  //}
  //else
  //  wi = hi = 0;
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setAlphaMaskForDot(AlphaMask<TWgt> *maskToUse)
{
  mask = maskToUse;
  //wb = mask->getWidth();
  //hb = mask->getHeight();
  // todo: handle nullptr case
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setNeighbourWeightsForSimpleDot(TWgt straight, TWgt diagonal)
{
  straightNeighbourWeight = straight;
  diagonalNeighbourWeight = diagonal;
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setAntiAlias(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setUseAlphaMask(bool shouldUseMask)
{
  useMask = shouldUseMask;
}

// painting

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::paintDot(TCor x, TCor y, TPix color)
{
  // todo: get rid of this dispatcher code - use a member function pointer instead.

  if(useMask)
  {
    if(antiAlias)
      paintDotViaMask(x, y, color);
    else
      paintDotViaMask((int) round(x), (int) round(y), color);
  }
  else
  {
    if(antiAlias)
      paintDot3x3(x, y, color, straightNeighbourWeight, diagonalNeighbourWeight);
    else
      paintDot3x3((int) round(x), (int) round(y), color, 
        straightNeighbourWeight, diagonalNeighbourWeight);
  }
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::paintDot3x3(int x, int y, TPix color, TWgt weightStraight, 
  TWgt weightDiagonal)
{
  int wi = image->getWidth();
  int hi = image->getHeight();

  if(x >= 0 && x < wi && y >= 0 && y < hi)
    blend(x, y, color);

  // apply thickness:
  if(weightStraight > 0.f && x >= 1 && x < wi-1 && y >= 1 && y < hi-1)
  {
    TPix a, ta, sa;
    a  = color;
    ta = a * (TPix)weightStraight;
    sa = a * (TPix)weightDiagonal;

    blend(x-1, y-1, sa);
    blend(x,   y-1, ta);
    blend(x+1, y-1, sa);

    blend(x-1, y,   ta);
    blend(x+1, y,   ta);

    blend(x-1, y+1, sa);
    blend(x,   y+1, ta);
    blend(x+1, y+1, sa);
  }
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::paintDot3x3(TCor x, TCor y, TPix color, TWgt weightStraight, 
  TWgt weightDiagonal)
{
  int wi = image->getWidth();
  int hi = image->getHeight();

  int xi = (int)floor(x);  // integer part of x
  int yi = (int)floor(y);  // integer part of y
  x -= xi;                 // fractional part of x
  y -= yi;                 // fractional part of y

  // compute weights for bilinear deinterpolation (maybe factor out):
  TPix a, b, c, d;
  d = TPix(x*y);
  c = TPix(y)-d;
  b = TPix(x)-d;
  a = d+TPix(1-x-y);

  // compute values to accumulate into the 4 pixels:
  a *= color;  // (xi,   yi)
  b *= color;  // (xi+1, yi)
  c *= color;  // (xi,   yi+1)
  d *= color;  // (xi+1, yi+1)

  // accumulate values into the pixels:
  if(xi >= 0 && xi < wi-1 && yi >= 0 && yi < hi-1)
  {
    blend(xi,   yi,   a);
    blend(xi+1, yi,   b);
    blend(xi,   yi+1, c);
    blend(xi+1, yi+1, d);
  }

  // apply thickness:
  if(weightStraight > 0.f && xi >= 1 && xi < wi-2 && yi >= 1 && yi < hi-2)
  {
    TPix t, s, sa, sb, sc, sd, ta, tb, tc, td;
    t = (TPix)weightStraight;      // weight for direct neighbour pixels
    s = (TPix)weightDiagonal;      // weight for diagonal neighbour pixels

    sa = s*a;
    sb = s*b;
    sc = s*c;
    sd = s*d;

    ta = t*a;
    tb = t*b;
    tc = t*c;
    td = t*d;

    blend(xi-1, yi-1, sa);
    blend(xi,   yi-1, ta+sb);
    blend(xi+1, yi-1, tb+sa);
    blend(xi+2, yi-1, sb);

    blend(xi-1, yi,   ta+sc);
    blend(xi,   yi,   sd+tb+tc);
    blend(xi+1, yi,   sc+ta+td);
    blend(xi+2, yi,   tb+sd);

    blend(xi-1, yi+1, tc+sa);
    blend(xi,   yi+1, sb+ta+td);
    blend(xi+1, yi+1, sa+tb+tc);
    blend(xi+2, yi+1, td+sb);

    blend(xi-1, yi+2, sc);
    blend(xi,   yi+2, tc+sd);
    blend(xi+1, yi+2, td+sc);
    blend(xi+2, yi+2, sd);
  }
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::paintDotViaMask(int x, int y, TPix color)
{
  int wi = image->getWidth();
  int hi = image->getHeight();
  int wm = mask->getWidth();
  int hm = mask->getHeight();

  // write coordinates in target image:
  int xs = x  - wm/2;    // start x-coordinate
  y      = y  - hm/2;    // start y coordinate
  int xe = xs + wm-1;    // end x coordinate
  int ye = y  + hm-1;    // end y coordinate

  // read coordinates in mask image:
  int mxs = 0;          // start x
  int my  = 0;          // start y

  // checks to not write beyond image bounds:
  if(xs < 0)
  {
    mxs = -xs;
    xs  = 0;
  }
  if(y < 0)
  {
    my = -y;
    y  =  0;
  }
  if(xe >= wi)
    xe = wi-1;
  if(ye >= hi)
    ye = hi-1;

  // the actual painting loop over the pixels in target image and brush image:
  int mx;
  while(y <= ye)
  {
    x  = xs;
    mx = mxs;
    while(x <= xe)
    {
      //// debug:
      //rsAssert(x  >= 0 && x  < wi);
      //rsAssert(y  >= 0 && y  < hi);
      //rsAssert(mx >= 0 && mx < wm);
      //rsAssert(my >= 0 && my < hm);

      blend(x, y, color, (*mask)(mx, my));
      x++;
      mx++;
    }
    y++;
    my++;
  }
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::paintDotViaMask(TCor x, TCor y, TPix color)
{
  //paintDotViaMask((int)round(x), (int)round(y), color);
  //return;
  //  // preliminary - calls the non anti-aliased version

  int wi = image->getWidth();
  int hi = image->getHeight();
  int wm = mask->getWidth();
  int hm = mask->getHeight();

  int  xi = (int)floor(x);  // integer part of x
  int  yi = (int)floor(y);  // integer part of y
  TCor xf = x - xi;         // fractional part of x
  TCor yf = y - yi;         // fractional part of y

  // compute pixel coverages (these are the same as the weights for bilinear deinterpolation):
  TWgt a, b, c, d;
  d = TWgt(xf*yf);
  c = TWgt(yf)-d;
  b = TWgt(xf)-d;
  a = d+TWgt(1-xf-yf);

  // write coordinates in target image:
  int xs = xi - wm/2;    // start x-coordinate in image
  int ys = yi - hm/2;    // start y coordinate in image
  int xe = xs + wm;      // end x coordinate in image
  int ye = ys + hm;      // end y coordinate in image
  //int xs = floor(x - 0.5*wm);  // start x-coordinate in image
  //int ys = floor(y - 0.5*hm);  // start y coordinate in image
  //int xe = ceil (x + 0.5*wm);  // end x coordinate in image
  //int ye = ceil (y + 0.5*hm);  // end y coordinate in image

  // read coordinates in mask image:
  int xms = 0;   // start x
  int yms = 0;   // start y

  TWgt w; // weight for current pixel

  // flags to indicate whether or not we need to draw the 4 edges of the mask:
  bool leftEdge, rightEdge, topEdge, bottomEdge;
  leftEdge   = true;
  rightEdge  = true;
  topEdge    = true;
  bottomEdge = true;


  // checks to not write beyond image bounds:
  //// old version:
  //if(xs < 0)
  //{
  //  xms = -xs;
  //  xs  = 0;
  //  leftEdge = false;
  //}
  //if(ys < 0)
  //{
  //  yms = -ys;
  //  ys  =  0;
  //  topEdge = false;
  //}
  //if(xe >= wi)
  //{
  //  xe = wi-1;
  //  rightEdge = false;
  //}
  //if(ye >= hi)
  //{
  //  ye = hi-1;
  //  bottomEdge = false;
  //}

  // new version:
  if(xs < -1)
  {
    xms = -xs;
    xs  = 0;
    leftEdge = false;
  }
  if(ys < -1)
  {
    yms = -ys;
    ys  =  0;
    topEdge = false;
  }
  if(xe > wi)
  {
    xe = wi;
    rightEdge = false;
  }
  if(ye > hi)
  {
    ye = hi;
    bottomEdge = false;
  }


  // paint edges and corners:
  int xm, ym;   // x- and y-index in mask
  if(leftEdge)
  {
    ym = yms;
    for(yi = ys+1; yi <= ye-1; yi++)
    {
      w = a * (*mask)(0, ym+1) + c * (*mask)(0, ym); 
      blend(xs, yi, color, w);
      ym++;
    }
    if(topEdge)
      blend(xs, ys, color, a * (*mask)(0, 0));  // top left corner
  }
  if(topEdge)
  {
    //xm = 0;
    xm = xms;
    for(xi = xs+1; xi <= xe-1; xi++) 
    {
      w = a * (*mask)(xm+1, 0) + b * (*mask)(xm, 0); 
      blend(xi, ys, color, w);
      xm++;
    }
    if(rightEdge)
      blend(xe, ys, color, b * (*mask)(wm-1, 0));  // top right corner
  }
  if(rightEdge)
  {
    //ym = 0;
    ym = yms;
    for(yi = ys+1; yi <= ye-1; yi++) 
    {
      w = b * (*mask)(wm-1, ym+1) + d * (*mask)(wm-1, ym); 
      blend(xe, yi, color, w);
      ym++;
    }
  }
  if(bottomEdge)
  {
    //xm = 0;
    xm = xms;
    for(xi = xs+1; xi <= xe-1; xi++) 
    {
      w = c * (*mask)(xm+1, hm-1) + d * (*mask)(xm, hm-1); 
      blend(xi, ye, color, w);
      xm++;
    }
    if(leftEdge)
      blend(xs, ye, color, c * (*mask)(0,    hm-1));  // bottom left corner
    if(rightEdge)
      blend(xe, ye, color, d * (*mask)(wm-1, hm-1));  // bottom right corner
  }

  // there's still something wrong with the edge code - there's a 1 pixel wide zone border in the 
  // image onto which nothing is drawn


  // paint interior rectangle:
  ym = yms;
  for(yi = ys+1; yi <= ye-1; yi++)
  {
    xm = xms; // why not +1?
    for(xi = xs+1; xi <= xe-1; xi++)
    {
      w = a * (*mask)(xm+1, ym+1) + b * (*mask)(xm, ym+1)
        + c * (*mask)(xm+1, ym)   + d * (*mask)(xm, ym);   // check this carefully!
      blend(xi, yi, color, w); 
      xm++;
    }
    ym++;
  }
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::drawDottedLine(TCor x1, TCor y1, TCor x2, TCor y2, TPix color, 
  TCor density)
{
  TCor dx = x2-x1;
  TCor dy = y2-y1;
  TCor pixelDistance = sqrt(dx*dx + dy*dy);
  int  numDots = rsMax(1, (int)floor(density*pixelDistance));
  TPix scaledColor = (TPix) (color / (TPix)numDots); // maybe make this scaling optional
  TCor scaler = (TCor)(1.0 / numDots);
  TCor k;
  for(int i = 1; i <= numDots; i++)
  {
    k = scaler * i;  // == i / numDots
    paintDot(x1 + k*dx, y1 + k*dy, scaledColor);
  }
}

